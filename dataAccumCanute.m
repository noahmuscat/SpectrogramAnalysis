% Define the base directories for each condition
baseDirs = {'/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Canute/300Lux', ...
            '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Canute/1000LuxWk1', ...
            '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Canute/1000LuxWk4'};

% Initialize the main structure
CanuteCombined = struct();
CanuteCombined.MetaData = struct();
CanuteCombined.MetaData.Ch = 80; % Channel of interest
CanuteCombined.MetaData.Rat = 'Canute';
CanuteCombined.MetaData.fo = []; % Initialize MetaData.fo as an empty array

% Define each condition's data structure
conditions = {'300Lux', '1000LuxWk1', '1000LuxWk4'};
validConditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
for i = 1:length(conditions)
    condition = conditions{i};
    validCondition = validConditions{i};
    
    % Initialize fields for each condition
    CanuteCombined.(validCondition) = struct();
    CanuteCombined.(validCondition).Datetime = datetime.empty();
    CanuteCombined.(validCondition).ZT_time = [];
    CanuteCombined.(validCondition).SleepState = [];
    CanuteCombined.(validCondition).FrequencyPower = {}; % Cell array for power spectra

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
            CanuteCombined.(validCondition).Datetime = [CanuteCombined.(validCondition).Datetime; appropriateTimestamps];
            CanuteCombined.(validCondition).ZT_time = [CanuteCombined.(validCondition).ZT_time; ZTData.ZT_Hour];
            CanuteCombined.(validCondition).SleepState = [CanuteCombined.(validCondition).SleepState; ZTData.Sleep_State];
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
                if isempty(CanuteCombined.MetaData.fo)
                    CanuteCombined.MetaData.fo = fo; % Save frequencies once
                end

                % Append power data only if timestamps and states are valid
                numRecords = numel(appropriateTimestamps);
                for t = 1:min(size(spec, 1), numRecords) % Choose the minimum
                    CanuteCombined.(validCondition).FrequencyPower{end+1, 1} = spec(t, :)';
                end
            else
                warning("Channel 80 not found in EEG data for folder: %s", folderName);
            end
        else
            warning("EEG file not found: %s", eegFile);
        end
    end
    
    % Ensure consistency across all data fields
    numEntries = length(CanuteCombined.(validCondition).Datetime);
    if length(CanuteCombined.(validCondition).FrequencyPower) > numEntries
        CanuteCombined.(validCondition).FrequencyPower = CanuteCombined.(validCondition).FrequencyPower(1:numEntries);
    end
end

%% Normalization
numStdDevs = 4;  % Number of standard deviations for outlier detection
luxField = 'Cond_300Lux';
data300Lux = CanuteCombined.(luxField);

% Calculate timeframe
datetime_list = data300Lux.Datetime;
endTime = datetime_list(end);
startTime = endTime - days(4);

% Last 4 days' data
last4DaysMask = (datetime_list >= startTime) & (datetime_list <= endTime);
last4DaysZT = data300Lux.ZT_time(last4DaysMask);
last4DaysFrequencyPower = data300Lux.FrequencyPower(last4DaysMask);

% Mean and STD holders
mean_day = zeros(size(CanuteCombined.MetaData.fo));
std_day = zeros(size(CanuteCombined.MetaData.fo));
mean_night = zeros(size(CanuteCombined.MetaData.fo));
std_night = zeros(size(CanuteCombined.MetaData.fo));

% Compute mean and std dev per frequency
for freqIdx = 1:length(CanuteCombined.MetaData.fo)
    dayMask = (last4DaysZT >= 0) & (last4DaysZT <= 11);
    dayData = cell2mat(cellfun(@(x) x(freqIdx), last4DaysFrequencyPower(dayMask), 'UniformOutput', false));
    mean_day(freqIdx) = mean(dayData, 'omitnan');
    std_day(freqIdx) = std(dayData, 'omitnan');
    
    nightMask = (last4DaysZT >= 12) & (last4DaysZT <= 23);
    nightData = cell2mat(cellfun(@(x) x(freqIdx), last4DaysFrequencyPower(nightMask), 'UniformOutput', false));
    mean_night(freqIdx) = mean(nightData, 'omitnan');
    std_night(freqIdx) = std(nightData, 'omitnan');
end

% Define normalization function with outlier removal
function z_power = zscore_with_outlier_removal(power, mean_val, std_val, numStdDevs)
    z_power = (power - mean_val) ./ std_val;
    outlier_mask = abs(z_power) <= numStdDevs;
    z_power(~outlier_mask) = NaN;  % Remove outliers by setting them to NaN
end

% Normalize each condition
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
for cond = 1:length(conditions)
    cond_data = CanuteCombined.(conditions{cond});
    zscoreFreqPower = cell(size(cond_data.FrequencyPower));
    
    for idx = 1:length(cond_data.FrequencyPower)
        isDay = (cond_data.ZT_time(idx) >= 0 && cond_data.ZT_time(idx) <= 11);
        mean_usage = isDay * mean_day + ~isDay * mean_night;
        std_usage = isDay * std_day + ~isDay * std_night;
        
        powerArray = cond_data.FrequencyPower{idx};
        zscoreFreqPower{idx} = zscore_with_outlier_removal(powerArray, mean_usage, std_usage, numStdDevs);
    end
    
    CanuteCombined.(conditions{cond}).ZscoredFrequencyPower = zscoreFreqPower;
end

%% Save Processed Data
save('CanuteCombinedDataOutliersRemoved.mat', 'CanuteCombined', '-v7.3');

%% functions
function ZTHour = calculateZT(timestamp)
    baseHour = 5; % Default start for ZT 0 is 5 AM
    
    if isDST(timestamp)
        baseHour = 6; % Adjust for DST
    end
    
    % Calculate the hour of the day since base hour and adjust to ZT hour
    hourOfDay = hour(timestamp);
    minuteOfDay = minute(timestamp) / 60;
    
    elapsedHours = hourOfDay + minuteOfDay - baseHour;
    if elapsedHours < 0
        elapsedHours = elapsedHours + 24; % Wrap around midnight
    end
    
    ZTHour = floor(elapsedHours);
end

function isDst = isDST(timestamp)
    % DST starts on the second Sunday in March and ends on the first Sunday in November
    startDST = datetime(timestamp.Year, 3, 8) + days(7 - weekday(datetime(timestamp.Year, 3, 8), 'dayofweek'));
    endDST = datetime(timestamp.Year, 11, 1) + days(7 - weekday(datetime(timestamp.Year, 11, 1), 'dayofweek'));
    isDst = (timestamp >= startDST) && (timestamp < endDST);
end

function [appropriateTimestamps, ZTData] = getDataFromMatFile(fullFilePath)
    % Extract the folder name from the full file path
    [folderPath, ~, ~] = fileparts(fullFilePath);

    % Extract the base folder name
    [~, folderName] = fileparts(folderPath);

    % Extract the animal name and initial timestamp from the folder name
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

    % Generate appropriate timestamps by adding 'timestamps' to 'startDatetime'
    appropriateTimestamps = startDatetime + seconds(timestamps);

    % Create an empty table
    ZTHours = zeros(length(timestamps), 1);
    sleepStatesColumn = sleepStates;

    % Calculate ZT hours
    for i = 1:length(timestamps)
        ZTHours(i) = calculateZT(appropriateTimestamps(i));
    end

    % Combine the data into a table
    ZTData = table(appropriateTimestamps, ZTHours, sleepStatesColumn, ...
                   'VariableNames', {'Timestamp', 'ZT_Hour', 'Sleep_State'});
end