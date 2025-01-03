%% Makes pretty, heat mapped spectrograms
% file comes from TheStateEditor

% Define the base folder for the 300 Lux lighting condition
baseFolder = '/data/Jeremy/Sleepscoring_Data_Noah/Canute/300Lux';
[~, condition] = fileparts(baseFolder);
[~, animalName] = fileparts(fileparts(baseFolder)); % Get animal name (Canute)

% Get a list of all subfolders in the base folder
subFolders = dir(baseFolder);
subFolders = subFolders([subFolders.isdir]);  % Keep only directories
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % Remove '.' and '..'

% Optional outlier removal setting
stdDevThreshold = 2; % Set the number of standard deviations to use for outlier removal; set to empty ([]) to disable

% Initialize variables for pooled data
pooledSpec = [];
pooledIndices = [];
numChannels = 0;
firstRun = true;

% Iterate through all subfolders and process the data
for k = 1:length(subFolders)
    currentSubFolder = fullfile(baseFolder, subFolders(k).name);
    matFile = dir(fullfile(currentSubFolder, '*.eegstates.mat'));
    
    if ~isempty(matFile)
        rootPath = fullfile(currentSubFolder, matFile.name);
        data = load(rootPath);

        % Extract start time from the file path
        [~, folderName, ~] = fileparts(fileparts(rootPath));
        startDatetime = datetime(folderName(8:end), 'InputFormat', 'yyMMdd_HHmmss', 'TimeZone', 'America/New_York');

        % Determine the DST offset
        isDST = isdst(startDatetime);
        if isDST
            lightsOnHour = 6;
        else
            lightsOnHour = 5;
        end

        % Get the number of channels
        numChannels = length(data.StateInfo.fspec);

        % Loop over the specified channels and process data
        for i = 1:numChannels
            % Extract relevant fields
            spec = data.StateInfo.fspec{1, i}.spec; % Spectrogram data (power or magnitude)
            times = data.StateInfo.fspec{1, i}.to; % Time points in seconds
            freqs = data.StateInfo.fspec{1, i}.fo; % Frequencies in Hz

            % Adjust times based on the initial start time and lights on hour
            initialDatetime = startDatetime - hours(lightsOnHour);
            adjustedDatetimes = initialDatetime + seconds(times);
            
            % Determine ZT hours for each data point
            adjustedHours = hour(adjustedDatetimes); % Extract the hour part of the adjusted times

            if firstRun
                % Initialize pooled data structures for each channel
                pooledSpec = zeros(24, length(freqs), numChannels); % Accumulate binned spectrogram
                pooledIndices = zeros(24, numChannels); % Count of data points in each bin
                firstRun = false;
            end

            % Average power for each frequency in each ZT hour bin
            for zt = 0:23
                indices = (adjustedHours == zt);
                if any(indices)
                    pooledSpec(zt + 1, :, i) = pooledSpec(zt + 1, :, i) + sum(spec(indices, :), 1);
                    pooledIndices(zt + 1, i) = pooledIndices(zt + 1, i) + sum(indices);
                end
            end
        end
    end
end

% Calculate the average spectrogram for pooled data
for i = 1:numChannels
    for zt = 0:23
        if pooledIndices(zt + 1, i) > 0
            pooledSpec(zt + 1, :, i) = pooledSpec(zt + 1, :, i) / pooledIndices(zt + 1, i);
        end
    end
end

% Remove outliers if threshold is set
if ~isempty(stdDevThreshold)
    for i = 1:numChannels
        meanSpec = mean(pooledSpec(:, :, i), 1);
        stdSpec = std(pooledSpec(:, :, i), 0, 1);

        for zt = 0:23
            if pooledIndices(zt + 1, i) > 0
                outliers = abs(pooledSpec(zt + 1, :, i) - meanSpec) > (stdDevThreshold * stdSpec);
                pooledSpec(zt + 1, outliers, i) = meanSpec(outliers);  % Replace outliers with mean
            end
        end
    end
end

% Plotting pooled spectrogram data
for i = 1:numChannels
    channel = data.StateInfo.fspec{1, i}.info.Ch; % Channel number
    binnedSpec = pooledSpec(:, :, i); % Extract binned spectrogram

    % Plotting spectrogram (binned per hour)
    zt_labels = cellstr(num2str((0:23)')); % Create ZT labels
    figure;
    imagesc(0:23, freqs, binnedSpec'); % Convert ZT hours to the x-axis
    axis xy;
    xlabel('Zeitgeber Time (ZT)');
    xticks(0:23); % ZT 0 to 23
    xticklabels(zt_labels);
    ylabel('Frequency (Hz)');
    title(sprintf('%s %s - Spectrogram - Channel %d', animalName, condition, channel));
    colorbar;
    colormap('jet');
    set(gca, 'FontSize', 14);
    grid on;
end

%% function definitions
function specs = SaveSpectrogramsAsStruct(StateInfo)
    specs = struct('spec', {}, 'freqs', {}, 'times', {}, 'info', {});

    for a = 1:length(StateInfo.fspec)
        specs(a).spec = StateInfo.fspec{a}.spec;
        specs(a).freqs = StateInfo.fspec{a}.fo;
        specs(a).times = StateInfo.fspec{a}.to;    
        specs(a).info = StateInfo.fspec{a}.info;
    end
end