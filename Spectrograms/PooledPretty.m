%% IN PROGRESS
% TODO:
% 1) fix error where it displays many ZT hours even across only a few 5 minute
% bins
% 2) fix outliers
% log transform the power
% 3) change colormap to parula

%% Gathering and aggregating data
% Define directories for different experimental conditions
baseDirs = {'/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/300Lux', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/1000LuxWk1', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/1000LuxWk4'};

% Set the channels to include in the analysis
channelsToAnalyze = 80; % Set this to desired channel(s) for analysis

% Initialize a structure to hold the results for each condition
results = struct();

stdDevThreshold = 2; % for eliminating based on stdev
threshold = 10; % for eliminating 60 and 120 Hz bands

for b = 1:length(baseDirs)
    baseFolder = baseDirs{b};
    [~, condition] = fileparts(baseFolder); % Extract condition name (e.g., '300Lux')
    [~, animalName] = fileparts(fileparts(baseFolder)); % Extract animal name (e.g., 'Canute')

    % Modify condition name to be a valid field name
    validCondition = ['Cond', condition];

    % Initialize variables for the condition
    all5MinPower = [];
    all5MinTimestampsZT = [];
    all5MinTimestampsEST = [];
    channels = [];

    % Get a list of all subfolders in the base folder
    subFolders = dir(baseFolder);
    subFolders = subFolders([subFolders.isdir]);  % Keep only directories
    subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % Remove '.' and '..'

    for k = 1:length(subFolders)
        currentSubFolder = fullfile(baseFolder, subFolders(k).name);
        matFiles = dir(fullfile(currentSubFolder, '*.eegstates.mat'));

        for m = 1:length(matFiles)
            rootPath = fullfile(matFiles(m).folder, matFiles(m).name);
            data = load(rootPath);

            % Extract start time from the subfolder name
            [~, folderName, ~] = fileparts(currentSubFolder);
            startDatetimeEST = datetime(folderName(8:end), 'InputFormat', 'yyMMdd_HHmmss', 'TimeZone', 'America/New_York');
            
            % Determine the DST offset
            isDST = isdst(startDatetimeEST);
            lightsOnHour = isDST * 6 + ~isDST * 5;

            % Get the channel numbers
            channels = zeros(1, length(data.StateInfo.fspec));
            for i = 1:length(channels)
                channels(1, i) = data.StateInfo.fspec{1,i}.info.Ch;
            end

            % Logical index for the channels to analyze
            channelIndicesToAnalyze = ismember(channels, channelsToAnalyze);

            % Loop over the specified channels
            for i = find(channelIndicesToAnalyze)
                % Extract relevant fields
                spec = data.StateInfo.fspec{1, i}.spec;  % Spectrogram data
                freqs = data.StateInfo.fspec{1, i}.fo;   % Frequencies in Hz
                times = data.StateInfo.fspec{1, i}.to;   % Time points in seconds

                % Zero out or attenuate power at 60 Hz and 120 Hz
                clean_spec = spec;
                clean_spec(:, ismember(round(freqs), [60, 120])) = NaN;

                % Adjust times based on lights on hour for ZT time
                initialDatetimeZT = startDatetimeEST - hours(lightsOnHour);
                adjustedDatetimesZT = initialDatetimeZT + seconds(times);
                adjustedDatetimesEST = startDatetimeEST + seconds(times);

                % Initialize variables for 5-minute bins
                fiveMinBinSize = 300; % 5 minutes in seconds
                numOf5MinBins = ceil(max(times) / fiveMinBinSize);
                fiveMinPower = zeros(numOf5MinBins, length(freqs));
                timestampsZT = NaT(numOf5MinBins, 1, 'TimeZone', 'America/New_York');
                timestampsEST = NaT(numOf5MinBins, 1, 'TimeZone', 'America/New_York');

                % Bin power in each 5-minute interval
                for bin = 1:numOf5MinBins
                    indices = (times >= (bin-1) * fiveMinBinSize) & (times < bin * fiveMinBinSize);
                    if any(indices)
                        currentSpec = clean_spec(indices, :);
                        fiveMinPower(bin, :) = mean(currentSpec, 1, 'omitnan');
                        timestampsZT(bin) = initialDatetimeZT + seconds((bin-1) * fiveMinBinSize);
                        timestampsEST(bin) = startDatetimeEST + seconds((bin-1) * fiveMinBinSize);
                    end
                end

                % Concatenate results across files for the condition
                all5MinPower = [all5MinPower; fiveMinPower];
                all5MinTimestampsZT = [all5MinTimestampsZT; timestampsZT];
                all5MinTimestampsEST = [all5MinTimestampsEST; timestampsEST];
            end
        end
    end

    % Store results for the condition
    results.(validCondition).FiveMinPower = all5MinPower;
    results.(validCondition).TimestampsZT = all5MinTimestampsZT;
    results.(validCondition).TimestampsEST = all5MinTimestampsEST;
    results.(validCondition).Freqs = freqs;
    results.(validCondition).Channels = channels(channelIndicesToAnalyze);
end

%% saving data
save('/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/metrics.mat', 'results');


%% plotting data
conditions = fieldnames(results);

for c = 1:length(conditions)
    condition = conditions{c};
    condData = results.(condition);
    
    % Extract data for the current condition
    powerData = condData.FiveMinPower;
    timestampsZT = condData.TimestampsZT;
    freqs = condData.Freqs;
    
    % Find indices of frequencies to be set to NaN
    noiseBands = (freqs >= 50 & freqs <= 70) | (freqs >= 110 & freqs <= 130);
    powerData(:, noiseBands) = NaN;  % Set noise ranges to NaN

    % Determine the start day (Day 1) and split by day
    days = unique(dateshift(timestampsZT, 'start', 'day'));
    
    for d = 1:length(days)
        dayStart = days(d);
        dayEnd = dayStart + duration(24, 0, 0);  % End of the current day
        
        % Index for the data of the current day
        dayIndices = (timestampsZT >= dayStart) & (timestampsZT < dayEnd);
        
        if any(dayIndices)
            dayPowerData = powerData(dayIndices, :);
            dayTimestampsZT = timestampsZT(dayIndices);
            
            % Calculate the number of data points for this day
            numDataPoints = sum(dayIndices);
            
            % Determine ZT hours for x-axis labels as a difference in time
            ztHours = hour(dayTimestampsZT);
            
            % Create spectrogram plot for the current day
            figure;
            imagesc(ztHours, freqs, dayPowerData');
            axis xy;
            xlabel('ZT Hour', 'FontSize', 14);
            ylabel('Frequency (Hz)', 'FontSize', 14);
            title(['Spectrogram - ', condition, ' - Day ', num2str(d), ' (N = ', num2str(numDataPoints), ')'], 'FontSize', 16);
            colorbar;
            colormap('jet');
            ylim([min(freqs) max(freqs)]);
            set(gca, 'FontSize', 12);
            xticks(0:23);
            xticklabels(0:23);  % Label x-axis with ZT hours 0-23
            grid on;

            hTitle = get(gca, 'Title'); % Get the title handle
            titleText = get(hTitle, 'String'); % Extract the title string
            formattedTitle = regexprep(titleText, '[^\w]', '_'); % Replace invalid characters with underscores
            saveDir = '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute';

            fullPath = fullfile(saveDir, [formattedTitle, '.png']);
            saveas(gcf, fullPath); % Save the current figure


        end
    end
end