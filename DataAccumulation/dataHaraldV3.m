%% Step 0: Preparing the data

% Define the base directories for each condition
baseDirs = {
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/300Lux', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/1000LuxWk1', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/1000LuxWk4'
};

% Initialize the main structure
HaraldV3Combined80 = struct();
HaraldV3Combined80.MetaData = struct();
HaraldV3Combined80.MetaData.Ch = 80; % Channel of interest
HaraldV3Combined80.MetaData.Rat = 'Harald';
HaraldV3Combined80.MetaData.fo = []; % Initialize MetaData.fo as an empty array

% Define each condition's data structure
conditions = {'300Lux', '1000LuxWk1', '1000LuxWk4'};
validConditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};

%% Step 1: Data wrangling and fractional power calculation
% This step collects all of the data and calculates the fractional power of
% each time point. It also removes 55-65 Hz and 115-125 Hz to eliminate
% electrical outlet noise. This is done before fractional power calculation
% so that noise does not contribute to the calculation.

% Iterate through each condition
for i = 1:length(conditions)
    condition = conditions{i};
    validCondition = validConditions{i};

    % Initialize fields for each condition
    HaraldV3Combined80.(validCondition) = struct();
    HaraldV3Combined80.(validCondition).ZT_Datetime = [];
    HaraldV3Combined80.(validCondition).SleepState = [];
    HaraldV3Combined80.(validCondition).FractionalPower = {}; % Cell array for fractional power spectra

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
            [~, ZTData] = getDataFromMatFile(sleepFile);

            % Append to pooled data
            HaraldV3Combined80.(validCondition).ZT_Datetime = [HaraldV3Combined80.(validCondition).ZT_Datetime; ZTData.ZT_Datetime];
            HaraldV3Combined80.(validCondition).SleepState = [HaraldV3Combined80.(validCondition).SleepState; ZTData.Sleep_State];
        else
            warning("SleepState file not found: %s", sleepFile);
        end

        % Load and process EEG states file
        eegFile = fullfile(dayDir, strcat(folderName(1:end-7), '.eegstates.mat'));
        if exist(eegFile, 'file')
            eegData = load(eegFile);
            chIdx = find(eegData.StateInfo.Chs == 80, 1); % Find index of channel 80
            
            if ~isempty(chIdx)
                % Retrieve spectrogram data and frequency info
                spec = eegData.StateInfo.fspec{1, chIdx}.spec; % TxF
                fo = eegData.StateInfo.fspec{1, chIdx}.fo; % Fx1

                if isempty(HaraldV3Combined80.MetaData.fo)
                    HaraldV3Combined80.MetaData.fo = fo; % Save frequencies once
                end

                % Identify frequencies to blank based on electrical noise
                noiseIndices = (fo >= 55 & fo <= 65) | (fo >= 115 & fo <= 125);

                % Calculate fractional power for each time point
                numRecords = numel(ZTData.ZT_Datetime);
                for t = 1:min(size(spec, 1), numRecords) % Choose the minimum

                    % Get the power spectrum for the current time point
                    currentPower = spec(t, :)';

                    % Set noise frequencies to NaN
                    currentPower(noiseIndices) = NaN;

                    % Calculate total power (ignoring NaNs with 'omitnan')
                    totalPower = sum(currentPower, 'omitnan');

                    % Calculate fractional power
                    if totalPower > 0
                        fractionalPower = currentPower / totalPower;
                    else
                        fractionalPower = NaN(size(currentPower));
                    end

                    % Append to data structure
                    HaraldV3Combined80.(validCondition).FractionalPower{end+1, 1} = fractionalPower;
                end
            else
                warning("Channel 80 not found in EEG data for folder: %s", folderName);
            end
        else
            warning("EEG file not found: %s", eegFile);
        end
    end

    % Ensure consistency across all data fields
    numEntries = length(HaraldV3Combined80.(validCondition).ZT_Datetime);
    if length(HaraldV3Combined80.(validCondition).FractionalPower) > numEntries
        HaraldV3Combined80.(validCondition).FractionalPower = HaraldV3Combined80.(validCondition).FractionalPower(1:numEntries);
    end
end

%% Step 2: Outlier removal per condition, sleep state, ZT hour, and frequency

% Define z-score threshold
zScoreThreshold = 2;

for i = 1:length(conditions)
    validCondition = validConditions{i};

    % Check if ZT_Datetime is a datetime array
    if isfield(HaraldV3Combined80.(validCondition), 'ZT_Datetime') && ...
       isa(HaraldV3Combined80.(validCondition).ZT_Datetime, 'datetime') && ...
       ~isempty(HaraldV3Combined80.(validCondition).ZT_Datetime)
        
        % Extract ZT hour integers
        ZT_Hour = hour(HaraldV3Combined80.(validCondition).ZT_Datetime);
        
    else
        error('ZT_Datetime is either missing, not a datetime array, or empty in %s', validCondition);
    end
    
    % Initialize cell array for cleaned data
    HaraldV3Combined80.(validCondition).CleanedFractionalPower = cell(size(HaraldV3Combined80.(validCondition).FractionalPower));

    % Obtain sleep state and fractional powers
    sleepStates = HaraldV3Combined80.(validCondition).SleepState;
    allFractionalPowers = cat(2, HaraldV3Combined80.(validCondition).FractionalPower{:});

    % Define unique sleep states and ZT hours
    uniqueStates = unique(sleepStates);
    uniqueHours = 0:23;

    % Iteration over sleep states and ZT hours
    for s = 1:length(uniqueStates)
        state = uniqueStates(s);

        for ZT_HourValue = uniqueHours
            % Find data indices matching current state and hour
            stateHourIndices = find(sleepStates == state & ZT_Hour == ZT_HourValue);
            if isempty(stateHourIndices)
                continue;
            end

            % Gather power data for this state and ZT hour
            stateHourFractionalPowers = allFractionalPowers(:, stateHourIndices);

            % Calculate mean and standard deviation for each frequency
            meanPower = mean(stateHourFractionalPowers, 2, 'omitnan');
            stdPower = std(stateHourFractionalPowers, 0, 2, 'omitnan');

            for idx = stateHourIndices'
                currentFractionalPower = HaraldV3Combined80.(validCondition).FractionalPower{idx};

                % Calculate z-scores for each frequency
                zScores = (currentFractionalPower - meanPower) ./ stdPower;

                % Identify outliers based on z-score
                isOutlier = abs(zScores) > zScoreThreshold;

                % Exclude outliers by setting to NaN
                cleanedPower = currentFractionalPower;
                cleanedPower(isOutlier) = NaN;

                % Store cleaned fractional power
                HaraldV3Combined80.(validCondition).CleanedFractionalPower{idx} = cleanedPower;
            end
        end
    end
end

%% Step 3: Normalize to 300Lux Baseline per sleep state, ZT hour, and frequency

% Define sleep states
sleepStates = [1, 3, 5]; % 1 = WAKE, 3 = NREM, 5 = REM

% Use the last 4 days of data from the 300Lux condition as the normalization baseline
luxField = 'Cond_300Lux';
data300Lux = HaraldV3Combined80.(luxField);

% Calculate the timeframe for the last 4 days
datetime_list = data300Lux.ZT_Datetime;
endTime = datetime_list(end);
startTime = endTime - days(4);

% Mask for the last 4 days' data
last4DaysMask = (datetime_list >= startTime) & (datetime_list <= endTime);

% Extract last 4 days' data
last4DaysSleepState = data300Lux.SleepState(last4DaysMask);
last4DaysFractionalPower = data300Lux.CleanedFractionalPower(last4DaysMask);
last4DaysZT_Hour = hour(datetime_list(last4DaysMask));

% Initialize holders for mean and std for each frequency, sleep state, and hour
numFrequencies = length(HaraldV3Combined80.MetaData.fo);
mean_baseline = NaN(numFrequencies, length(sleepStates), 24);
std_baseline = NaN(numFrequencies, length(sleepStates), 24);

% Calculate baseline mean and std
for s = 1:length(sleepStates)
    state = sleepStates(s);

    for hr = uniqueHours
        % Find corresponding indices
        stateHourIndices = find(last4DaysSleepState == state & last4DaysZT_Hour == hr);
        if isempty(stateHourIndices)
            continue;
        end

        % Gather power data for this sleep state and hour
        stateHourFractionalPowers = cat(2, last4DaysFractionalPower{stateHourIndices});

        % Calculate mean and std
        mean_baseline(:, s, hr + 1) = mean(stateHourFractionalPowers, 2, 'omitnan');
        std_baseline(:, s, hr + 1) = std(stateHourFractionalPowers, 0, 2, 'omitnan');
    end
end

% Normalize each condition's data
for i = 1:length(conditions)
    validCondition = validConditions{i};
    conditionData = HaraldV3Combined80.(validCondition);

    % Initialize z-score storage
    conditionData.ZScoreFractionalPower = cell(size(conditionData.CleanedFractionalPower));

    % Extract ZT hour for each data point
    ZT_Hour = hour(conditionData.ZT_Datetime);

    % Iterate over each sleep state and ZT hour
    for s = 1:length(sleepStates)
        state = sleepStates(s);

        for hrs = uniqueHours
            stateHourIndices = find(conditionData.SleepState == state & ZT_Hour == hrs);
            if isempty(stateHourIndices)
                continue;
            end

            % Retrieve baseline mean and std for the specific state and hour
            baselineMean = mean_baseline(:, s, hrs + 1);
            baselineStd = std_baseline(:, s, hrs + 1);

            % Iterate over each index
            for idx = stateHourIndices'
                currentFractionalPower = conditionData.CleanedFractionalPower{idx};

                % Calculate z-scores
                zScores = (currentFractionalPower - baselineMean) ./ baselineStd;

                % Store z-scores
                conditionData.ZScoreFractionalPower{idx} = zScores;
            end
        end
    end

    % Update the main structure with the z-score data
    HaraldV3Combined80.(validCondition).ZScoreFractionalPower = conditionData.ZScoreFractionalPower;
end

%% Save processed data
save('HaraldV3Ch80.mat', 'HaraldV3Combined80', '-v7.3');

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
