%% Step 0: Preparing the data

% Define the base directories for each condition
baseDirs = {
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/300Lux', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/1000LuxWk1', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/1000LuxWk4'
};

% Initialize the main structure
HaraldV2Combined80 = struct();
HaraldV2Combined80.MetaData = struct();
HaraldV2Combined80.MetaData.Ch = 80; % Channel of interest
HaraldV2Combined80.MetaData.Rat = 'Harald';
HaraldV2Combined80.MetaData.fo = []; % Initialize MetaData.fo as an empty array

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
    HaraldV2Combined80.(validCondition) = struct();
    HaraldV2Combined80.(validCondition).ZT_Datetime = [];
    HaraldV2Combined80.(validCondition).SleepState = [];
    HaraldV2Combined80.(validCondition).FractionalPower = {}; % Cell array for fractional power spectra

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
            HaraldV2Combined80.(validCondition).ZT_Datetime = [HaraldV2Combined80.(validCondition).ZT_Datetime; ZTData.ZT_Datetime];
            HaraldV2Combined80.(validCondition).SleepState = [HaraldV2Combined80.(validCondition).SleepState; ZTData.Sleep_State];
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

                if isempty(HaraldV2Combined80.MetaData.fo)
                    HaraldV2Combined80.MetaData.fo = fo; % Save frequencies once
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
                    HaraldV2Combined80.(validCondition).FractionalPower{end+1, 1} = fractionalPower;
                end
            else
                warning("Channel 80 not found in EEG data for folder: %s", folderName);
            end
        else
            warning("EEG file not found: %s", eegFile);
        end
    end

    % Ensure consistency across all data fields
    numEntries = length(HaraldV2Combined80.(validCondition).ZT_Datetime);
    if length(HaraldV2Combined80.(validCondition).FractionalPower) > numEntries
        HaraldV2Combined80.(validCondition).FractionalPower = HaraldV2Combined80.(validCondition).FractionalPower(1:numEntries);
    end
end

%% Step 2: Outlier removal 
% This part of the code uses the fractional powers and removes outliers
% based on condition, sleep state, and frequency

% Define z-score threshold
zScoreThreshold = 2;

for i = 1:length(conditions)
    validCondition = validConditions{i};

    % Initialize cell array for cleaned data
    HaraldV2Combined80.(validCondition).CleanedFractionalPower = cell(size(HaraldV2Combined80.(validCondition).FractionalPower));

    % Obtain sleep state and fractional powers
    sleepStates = HaraldV2Combined80.(validCondition).SleepState;
    allFractionalPowers = cat(2, HaraldV2Combined80.(validCondition).FractionalPower{:}); % FxT

    % Unique sleep states for separation
    uniqueStates = unique(sleepStates);

    for s = 1:length(uniqueStates)
        state = uniqueStates(s);

        % Select indices for the current sleep state
        stateIndices = find(sleepStates == state);
        if isempty(stateIndices)
            continue;
        end

        % Gather power data for only this sleep state
        stateFractionalPowers = allFractionalPowers(:, stateIndices);

        % Calculate mean and standard deviation for each frequency in this state
        meanPower = mean(stateFractionalPowers, 2, 'omitnan');
        stdPower = std(stateFractionalPowers, 0, 2, 'omitnan');

        % Iterate over each time point in this state
        for idx = stateIndices'
            currentFractionalPower = HaraldV2Combined80.(validCondition).FractionalPower{idx};

            % Calculate z-score for each frequency
            zScores = (currentFractionalPower - meanPower) ./ stdPower;

            % Identify outliers based on z-score
            isOutlier = abs(zScores) > zScoreThreshold;

            % Exclude outliers by setting to NaN, retaining structure
            cleanedPower = currentFractionalPower;
            cleanedPower(isOutlier) = NaN;

            % Store cleaned fractional power back in the original index
            HaraldV2Combined80.(validCondition).CleanedFractionalPower{idx} = cleanedPower;
        end
    end
end

%% Step 3: Normalize to 300Lux Baseline

% Define sleep states
sleepStates = [1, 3, 5]; % 1 = WAKE, 3 = NREM, 5 = REM

% Use the last 4 days of data from the 300Lux condition as the normalization baseline
luxField = 'Cond_300Lux';
data300Lux = HaraldV2Combined80.(luxField);

% Calculate the timeframe for the last 4 days
datetime_list = data300Lux.ZT_Datetime; % Using ZT_Datetime, assuming you have updated it
endTime = datetime_list(end);
startTime = endTime - days(4);

% Mask for the last 4 days' data
last4DaysMask = (datetime_list >= startTime) & (datetime_list <= endTime);

% Extract last 4 days' data
last4DaysSleepState = data300Lux.SleepState(last4DaysMask);
last4DaysFractionalPower = data300Lux.CleanedFractionalPower(last4DaysMask);

% Initialize holders for mean and std for each frequency and sleep state
numFrequencies = length(HaraldV2Combined80.MetaData.fo);
mean_sleepState = zeros(numFrequencies, length(sleepStates));
std_sleepState = zeros(numFrequencies, length(sleepStates));

% Calculate baseline mean and std for each sleep state
for s = 1:length(sleepStates)
    state = sleepStates(s);
    
    % Select data specific to the sleep state
    stateMask = (last4DaysSleepState == state);
    stateFractionalPower = cat(2, last4DaysFractionalPower{stateMask});
    
    % Calculate mean and std, omitting NaNs
    mean_sleepState(:, s) = mean(stateFractionalPower, 2, 'omitnan');
    std_sleepState(:, s) = std(stateFractionalPower, 0, 2, 'omitnan');
end

% Calculate z-scores across all conditions based on baseline statistics
for i = 1:length(conditions)
    validCondition = validConditions{i};
    conditionData = HaraldV2Combined80.(validCondition);

    % Initialize z-score storage
    conditionData.ZScoreFractionalPower = cell(size(conditionData.CleanedFractionalPower));

    % Iterate over each sleep state
    for s = 1:length(sleepStates)
        state = sleepStates(s);
        stateIndices = (conditionData.SleepState == state);

        % Retrieve baseline mean and std
        baselineMean = mean_sleepState(:, s);
        baselineStd = std_sleepState(:, s);

        % Iterate over each time point index for the current sleep state
        for idx = find(stateIndices)'
            currentFractionalPower = conditionData.CleanedFractionalPower{idx};

            % Calculate z-scores for current time point
            zScores = (currentFractionalPower - baselineMean) ./ baselineStd;

            % Store z-scores
            conditionData.ZScoreFractionalPower{idx} = zScores;
        end
    end
    
    % Update the main structure with the z-score data
    HaraldV2Combined80.(validCondition).ZScoreFractionalPower = conditionData.ZScoreFractionalPower;
end

%% Save processed data
save('HaraldV2Ch80.mat', 'HaraldV2Combined80', '-v7.3');

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
