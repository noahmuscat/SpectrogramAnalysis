% Define directories for different experimental conditions
baseDirs = {'/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/300Lux', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/1000LuxWk1', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/1000LuxWk4'};

saveDir = '/Users/noahmuscat/Desktop';

allConditions = {'300Lux', '1000LuxWk1', '1000LuxWk4'};

% Define frequency bands and their corresponding names
bandNames = {'delta', 'theta', 'spindle', 'lowbeta', 'highbeta'};
startStopFreqs = {[0.5, 4], [5, 10], [11, 19], [20, 30], [30, 40]};

% Possible sleep states: Wake (1), NREM (3), REM (5)
sleepStatesUnique = [1, 3, 5];
sleepStateNames = {'W', 'NR', 'R'};

channelsToAnalyze = 80;

% Store raw power data for each band, condition, and sleep state
bandConditionPowers = zeros(length(bandNames), length(baseDirs), length(sleepStatesUnique));

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

    % Initialize raw power storage for each band and sleep state
    for s = 1:length(sleepStatesUnique)
        state = sleepStatesUnique(s);
        stateIndices = (allSleepStates == state);
        
        % Extract spectral data for current sleep state
        stateSpectrogram = allSpectrogramData(stateIndices, :);

        % Compute raw power for each frequency band
        for bandIdx = 1:length(startStopFreqs)
            bandFreqs = startStopFreqs{bandIdx};
            % Logical index to find columns that fall within the current frequency band
            freqIndices = (allFreqs >= bandFreqs(1)) & (allFreqs <= bandFreqs(2));
            
            if any(freqIndices)
                % Sum over time (rows) and frequencies (columns) within the band
                bandPower = mean(stateSpectrogram(:, freqIndices), 'all');
                bandConditionPowers(bandIdx, b, s) = bandPower;
            else
                % If no frequencies fall within the band, set power to zero
                bandConditionPowers(bandIdx, b, s) = 0;
            end
        end
    end
end

%% Generate figures for each frequency band
for bandIdx = 1:length(bandNames)
    formattedTitle = [bandNames{bandIdx}];
    figure;
    bandPowers = squeeze(bandConditionPowers(bandIdx, :, :));
    
    % Calculate the maximum value for y-axis limit
    maxYLim = max(bandPowers, [], 'all');
    
    % Create subplots for each condition
    for c = 1:length(baseDirs)
        subplot(1, length(baseDirs), c);
        bar(bandPowers(c, :));
        title(allConditions{c});
        xticks(1:3);
        xticklabels(sleepStateNames);
        ylabel('Raw Power');
        ylim([0, maxYLim]); % Set the same y-axis limit for all subplots
    end
    sgtitle(formattedTitle);
    fullPath = fullfile(saveDir, [formattedTitle, '.png']);
    saveas(gcf, fullPath);
end