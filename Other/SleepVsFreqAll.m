% Define directories for different experimental conditions
baseDirs = {'/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/300Lux', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/1000LuxWk1', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/1000LuxWk4'};

saveDir = '/Users/noahmuscat/Desktop';

% Define frequency bands and their corresponding names
bandNames = {'delta', 'theta', 'spindle', 'lowbeta', 'highbeta'};
startStopFreqs = {[0.5, 4], [5, 10], [11, 19], [20, 30], [30, 40]};

% Possible sleep states: Wake (1), NREM (3), REM (5)
sleepStatesUnique = [1, 3, 5];
sleepStateNames = {'WAKE', 'NREM', 'REM'};

channelsToAnalyze = 80;

% Initialize plot
figure;

% Loop over each experimental condition
for b = 1:length(baseDirs)
    baseFolder = baseDirs{b};
    [~, condition] = fileparts(baseFolder); % Extract condition name
    
    % Initialize storage for pooled data
    allSpectrogramData = [];
    allSleepStates = [];
    allFreqs = [];

    % Get subfolders for each condition
    subFolders = dir(baseFolder);
    subFolders = subFolders([subFolders.isdir]);  % Keep only directories
    subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % Remove '.' and '..'

    % Process each subfolder, which contains the .mat files
    for k = 1:length(subFolders)
        currentSubFolder = fullfile(baseFolder, subFolders(k).name);
        
        % Load frequency data
        freqFiles = dir(fullfile(currentSubFolder, '*.eegstates.mat'));
        for m = 1:length(freqFiles)
            data = load(fullfile(freqFiles(m).folder, freqFiles(m).name));

            % Get the channel numbers
            channels = zeros(1, length(data.StateInfo.fspec));
            for i = 1:length(channels)
                channels(i) = data.StateInfo.fspec{i}.info.Ch;
            end

            % Logical index for the channels to analyze
            channelIndicesToAnalyze = find(ismember(channels, channelsToAnalyze));

            if ~isempty(channelIndicesToAnalyze)
                % Use the first match assuming there's only one channel of interest
                selectedChannelIndex = channelIndicesToAnalyze(1);
                specData = data.StateInfo.fspec{selectedChannelIndex}.spec;
                freqs = data.StateInfo.fspec{selectedChannelIndex}.fo;
                allSpectrogramData = [allSpectrogramData; specData];
                allFreqs = freqs; % Assume frequencies are consistent
            end
        end
        
        % Load sleep state data
        sleepFiles = dir(fullfile(currentSubFolder, '*.SleepState.states.mat'));
        for m = 1:length(sleepFiles)
            % Load the .mat file and extract sleep states
            fullFilePath = fullfile(sleepFiles(m).folder, sleepFiles(m).name);
            data = load(fullFilePath);
            sleepStates = data.SleepState.idx.states;
            allSleepStates = [allSleepStates; sleepStates];
        end
    end

    % Only keep valid sleep states
    validIndices = ismember(allSleepStates, sleepStatesUnique);
    allSpectrogramData = allSpectrogramData(validIndices, :);
    allSleepStates = allSleepStates(validIndices);

    % Initialize storage for percent power values for each sleep state
    percentPowers = zeros(length(sleepStatesUnique), length(bandNames));

    % Compute percent power for each sleep state
    for s = 1:length(sleepStatesUnique)
        state = sleepStatesUnique(s);
        stateIndices = (allSleepStates == state);
        
        % Extract spectral data for current sleep state
        stateSpectrogram = allSpectrogramData(stateIndices, :);

        % Calculate total power across the specified bands for normalization
        totalPowerSelectedBands = 0;
        for bandIdx = 1:length(startStopFreqs)
            bandFreqs = startStopFreqs{bandIdx};
            % Find the indices of frequencies that fall within the current band
            freqIndices = (allFreqs >= bandFreqs(1)) & (allFreqs <= bandFreqs(2));
            
            if any(freqIndices)
                % Compute the power for the current band
                bandPower = mean(stateSpectrogram(:, freqIndices), 'all');
                percentPowers(s, bandIdx) = bandPower;
                totalPowerSelectedBands = totalPowerSelectedBands + bandPower;
            else
                % If no frequencies fall within the current band, set power to zero
                percentPowers(s, bandIdx) = 0;
            end
        end

        % Normalize the power across the selected bands
        if totalPowerSelectedBands > 0
            percentPowers(s, :) = (percentPowers(s, :) / totalPowerSelectedBands) * 100;
        else
            percentPowers(s, :) = 0; % Handle case where no power is observed
        end
    end

    % Plot the data for the current condition
    subplot(length(baseDirs), 1, b);
    bar(percentPowers, 'stacked');
    title(['Condition: ', condition]);
    legend(bandNames, 'Location', 'BestOutside');
    xlabel('Sleep State');
    xticks(1:3);
    xticklabels(sleepStateNames);
    ylabel('% Power');
    ylim([0, 100]);
end

formattedTitle = 'SleepVsFreqAll';
fullPath = fullfile(saveDir, [formattedTitle, '.png']);
saveas(gcf, fullPath);