% Define the base directories for each condition
baseDirs = {'/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/300Lux', ...
            '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/1000LuxWk1', ...
            '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/1000LuxWk4'};

% Initialize the main structure
HaraldCombined = struct();
HaraldCombined.MetaData = struct();
HaraldCombined.MetaData.Ch = 80; % Channel of interest
HaraldCombined.MetaData.Rat = 'Harald';
HaraldCombined.MetaData.fo = []; % Initialize MetaData.fo as an empty array

% Define each condition's data structure
conditions = {'300Lux', '1000LuxWk1', '1000LuxWk4'};
validConditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
for i = 1:length(conditions)
    condition = conditions{i};
    validCondition = validConditions{i};
    
    % Initialize fields for each condition
    HaraldCombined.(validCondition) = struct();
    HaraldCombined.(validCondition).Datetime = datetime.empty();
    HaraldCombined.(validCondition).ZT_Datetime = [];
    HaraldCombined.(validCondition).SleepState = [];
    HaraldCombined.(validCondition).FrequencyPower = {}; % Cell array for power spectra

    % List all subdirectories for each condition
    dayDirs = dir(baseDirs{i});
    dayDirs = dayDirs([dayDirs.isdir]);
    dayDirs = dayDirs(~ismember({dayDirs.name}, {'.', '..'}));

    % Process each day's data
    for j = 1:length(dayDirs)
        dayDir = fullfile(baseDirs{i}, dayDirs(j).name);

        % Parse the start datetime from the folder name
        folderName = dayDirs(j).name;
        tokens = regexp(folderName, '_(\d{6})_(\d{6})$', 'tokens');
        if isempty(tokens)
            warning("Folder name does not match expected format: %s", folderName);
            continue;
        end
        dateStr = tokens{1}{1};
        timeStr = tokens{1}{2};
        startDatetime = datetime(['20', dateStr, timeStr], 'InputFormat', 'yyyyMMddHHmmss', 'TimeZone', 'America/New_York');

        % Load and process SleepState file
        sleepFile = fullfile(dayDir, strcat(folderName(1:end-7), '.SleepState.states.mat'));
        if exist(sleepFile, 'file')
            [appropriateTimestamps, ZTData] = getDataFromMatFile(sleepFile);
            
            % Append to pooled data
            HaraldCombined.(validCondition).Datetime = [HaraldCombined.(validCondition).Datetime; appropriateTimestamps];
            HaraldCombined.(validCondition).ZT_Datetime = [HaraldCombined.(validCondition).ZT_Datetime; ZTData.ZT_Datetime];
            HaraldCombined.(validCondition).SleepState = [HaraldCombined.(validCondition).SleepState; ZTData.Sleep_State];
        else
            warning("SleepState file not found: %s", sleepFile);
        end
        
        % Load and process EEG states file
        eegFile = fullfile(dayDir, strcat(folderName(1:end-7), '.eegstates.mat'));
        if exist(eegFile, 'file')
            eegData = load(eegFile);
            chIdx = find(eegData.StateInfo.Chs == 80, 1); % Find index of channel 80
            
            if ~isempty(chIdx)
                % Retrieve spectrogram data
                spec = eegData.StateInfo.fspec{1, chIdx}.spec; % TxF
                fo = eegData.StateInfo.fspec{1, chIdx}.fo; % Fx1

                % Ensure MetaData has frequencies
                if isempty(HaraldCombined.MetaData.fo)
                    HaraldCombined.MetaData.fo = fo; % Save frequencies once
                end

                % Append power data only if timestamps and states are valid
                numRecords = numel(appropriateTimestamps);
                for t = 1:min(size(spec, 1), numRecords) % Choose the minimum
                    HaraldCombined.(validCondition).FrequencyPower{end+1, 1} = spec(t, :)';
                end
            else
                warning("Channel 80 not found in EEG data for folder: %s", folderName);
            end
        else
            warning("EEG file not found: %s", eegFile);
        end
    end
    
    % Ensure consistency across all data fields
    numEntries = length(HaraldCombined.(validCondition).Datetime);
    if length(HaraldCombined.(validCondition).FrequencyPower) > numEntries
        HaraldCombined.(validCondition).FrequencyPower = HaraldCombined.(validCondition).FrequencyPower(1:numEntries);
    end
end

%% Normalization
% Exclude frequencies from analysis
exclude_ranges = [55, 65; 115, 125];

% Frequency values from metadata
frequencies = HaraldCombined.MetaData.fo;

% Create an exclusion mask
exclude_mask = false(size(frequencies));
for i = 1:size(exclude_ranges, 1)
    exclude_mask = exclude_mask | (frequencies >= exclude_ranges(i, 1) & frequencies <= exclude_ranges(i, 2));
end

% Indices of frequencies to include
include_indices = find(~exclude_mask);
numIncludedFrequencies = length(include_indices);

% Define sleep states
sleepStates = [1, 3, 5]; % 1 = WAKE, 3 = NREM, 5 = REM

% Use the last 4 days of data from the 300Lux condition as the normalization baseline
luxField = 'Cond_300Lux';
data300Lux = HaraldCombined.(luxField);

% Calculate timeframe
datetime_list = data300Lux.Datetime;
endTime = datetime_list(end);
startTime = endTime - days(4);

% Last 4 days' data
last4DaysMask = (datetime_list >= startTime) & (datetime_list <= endTime);
last4DaysSleepState = data300Lux.SleepState(last4DaysMask);
last4DaysFrequencyPower = data300Lux.FrequencyPower(last4DaysMask);

% Initialize mean and std holders for included frequencies and sleep states
mean_sleepState = zeros(numIncludedFrequencies, length(sleepStates));
std_sleepState = zeros(numIncludedFrequencies, length(sleepStates));

% Compute fractional power and the mean and std dev per frequency per sleep state
fractionalPower_last4Days = cell(size(last4DaysFrequencyPower));

for idx = 1:length(last4DaysFrequencyPower)
    powerArray = last4DaysFrequencyPower{idx};
    filteredPowerArray = powerArray(include_indices); % Apply frequency filter
    totalPower = sum(filteredPowerArray);
    fractionalPower_last4Days{idx} = filteredPowerArray / totalPower; % Fractional power
end

for sleepStateIdx = 1:length(sleepStates)
    currentState = sleepStates(sleepStateIdx);
    currentMask = last4DaysSleepState == currentState;
    currentFractionalPower = fractionalPower_last4Days(currentMask);

    for freqIdx = 1:numIncludedFrequencies
        freqData = cell2mat(cellfun(@(x) x(freqIdx), currentFractionalPower, 'UniformOutput', false));
        mean_sleepState(freqIdx, sleepStateIdx) = mean(freqData, 'omitnan');
        std_sleepState(freqIdx, sleepStateIdx) = std(freqData, 'omitnan');
    end
end

% Normalize each condition based on fractional power
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};

for condIdx = 1:length(conditions)
    cond = conditions{condIdx};
    cond_data = HaraldCombined.(cond);
    zscoreFreqPower = cell(size(cond_data.FrequencyPower));
    
    for idx = 1:length(cond_data.FrequencyPower)
        state = cond_data.SleepState(idx);
        powerArray = cond_data.FrequencyPower{idx};
        filteredPowerArray = powerArray(include_indices); % Apply frequency filter
        totalPower = sum(filteredPowerArray);
        fractionalPower = filteredPowerArray / totalPower;

        % Find the corresponding sleep state index
        sleepStateIdx = find(sleepStates == state, 1);

        % Safety check for unexpected sleep states
        if isempty(sleepStateIdx)
            warning('Unexpected sleep state encountered: %d', state);
            zscoreFreqPower{idx} = NaN(numIncludedFrequencies, 1);
            continue;
        end
        
        mean_usage = mean_sleepState(:, sleepStateIdx);
        std_usage = std_sleepState(:, sleepStateIdx);

        % Compute Z-score
        zscoreFreqPower{idx} = (fractionalPower - mean_usage) ./ std_usage;
        outlier_mask = abs(zscoreFreqPower{idx}) <= 2;
        zscoreFreqPower{idx}(~outlier_mask) = NaN;
    end
    
    HaraldCombined.(cond).ZscoredFrequencyPower = zscoreFreqPower;
end

%% Save Processed Data
save('HaraldCombinedNormalized.mat', 'HaraldCombined', '-v7.3');

%% functions
function isDst = isDST(timestamp)
    % DST starts on the second Sunday in March and ends on the first Sunday in November
    startDST = datetime(timestamp.Year, 3, 8) + days(7 - weekday(datetime(timestamp.Year, 3, 8), 'dayofweek'));
    endDST = datetime(timestamp.Year, 11, 1) + days(7 - weekday(datetime(timestamp.Year, 11, 1), 'dayofweek'));
    isDst = (timestamp >= startDST) && (timestamp < endDST);
end

function [appropriateTimestamps, ZTData] = getDataFromMatFile(fullFilePath)
    % Load the .mat file
    data = load(fullFilePath);
    sleepStatesStruct = data.SleepState;

    % Access the 'timestamps' and 'states' fields within 'idx'
    if isfield(sleepStatesStruct, 'idx') && isfield(sleepStatesStruct.idx, 'timestamps') && isfield(sleepStatesStruct.idx, 'states')
        timestamps = sleepStatesStruct.idx.timestamps;
        sleepStates = sleepStatesStruct.idx.states;
    else
        error('The .mat file does not have the expected structure with idx.timestamps and idx.states');
    end

    % Extract the folder name from the full file path
    [folderPath, ~, ~] = fileparts(fullFilePath);
    [~, folderName] = fileparts(folderPath);

    % Extract the date and time from the folder name
    tokens = regexp(folderName, '_(\d{6})_(\d{6})$', 'tokens');
    if isempty(tokens)
        error('The folder name does not match the expected format Animal_YYMMDD_HHMMSS');
    end
    dateStr = tokens{1}{1};
    timeStr = tokens{1}{2};

    % Combine date and time strings to create a datetime object
    startDateStr = ['20' dateStr]; % Assuming the date is in the format YYMMDD
    startTimeStr = timeStr;
    startDatetime = datetime([startDateStr, startTimeStr], 'InputFormat', 'yyyyMMddHHmmss');

    % Calculate appropriate timestamps by adding 'timestamps' to 'startDatetime'
    appropriateTimestamps = startDatetime + seconds(timestamps);

    % Calculate ZTDatetime by adjusting for DST
    ZTDatetimes = appropriateTimestamps;

    for i = 1:length(ZTDatetimes)
        if isDST(appropriateTimestamps(i))
            ZTDatetimes(i) = appropriateTimestamps(i) - hours(6);
        else
            ZTDatetimes(i) = appropriateTimestamps(i) - hours(5);
        end
    end

    % Combine the data into a table
    ZTData = table(appropriateTimestamps, ZTDatetimes, sleepStates, ...
                   'VariableNames', {'Timestamp', 'ZT_Datetime', 'Sleep_State'});
end
