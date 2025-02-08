% Define the base directories for each condition
baseDirs = {'/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/300Lux', ...
            '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/1000LuxWk1', ...
            '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/1000LuxWk4'};

% Initialize the main structure
dataStruct = struct();
dataStruct.MetaData = struct();
dataStruct.MetaData.Ch = 80; % Channel of interest
dataStruct.MetaData.Rat = 'Canute';
dataStruct.MetaData.fo = []; % Initialize MetaData.fo as an empty array

% Define each condition's data structure
conditions = {'300Lux', '1000LuxWk1', '1000LuxWk4'};
validConditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
for i = 1:length(conditions)
    condition = conditions{i};
    validCondition = validConditions{i};
    
    % Initialize fields for each condition
    dataStruct.(validCondition) = struct();
    dataStruct.(validCondition).Datetime = datetime.empty();
    dataStruct.(validCondition).ZT_time = [];
    dataStruct.(validCondition).SleepState = [];
    dataStruct.(validCondition).FrequencyPower = {}; % Cell array for power spectra

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
            dataStruct.(validCondition).Datetime = [dataStruct.(validCondition).Datetime; appropriateTimestamps];
            dataStruct.(validCondition).ZT_time = [dataStruct.(validCondition).ZT_time; ZTData.ZT_Hour];
            dataStruct.(validCondition).SleepState = [dataStruct.(validCondition).SleepState; ZTData.Sleep_State];
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
                if isempty(dataStruct.MetaData.fo)
                    dataStruct.MetaData.fo = fo; % Save frequencies once
                end

                % Append power data only if timestamps and states are valid
                numRecords = numel(appropriateTimestamps);
                for t = 1:min(size(spec, 1), numRecords) % Choose the minimum
                    dataStruct.(validCondition).FrequencyPower{end+1, 1} = spec(t, :)';
                end
            else
                warning("Channel 80 not found in EEG data for folder: %s", folderName);
            end
        else
            warning("EEG file not found: %s", eegFile);
        end
    end
    
    % Ensure consistency across all data fields
    numEntries = length(dataStruct.(validCondition).Datetime);
    if length(dataStruct.(validCondition).FrequencyPower) > numEntries
        dataStruct.(validCondition).FrequencyPower = dataStruct.(validCondition).FrequencyPower(1:numEntries);
    end
end

%% Save the structured data to a .mat file (v7.3 cuz its too big)
save('Canute_combined_data.mat', 'dataStruct', '-v7.3');

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