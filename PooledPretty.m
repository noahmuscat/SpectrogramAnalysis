%% Makes pretty, heat mapped spectrograms for All Conditions
%baseDirs = {'/data/Jeremy/Sleepscoring_Data_Noah/Canute/300Lux', ...
%            '/data/Jeremy/Sleepscoring_Data_Noah/Canute/1000LuxWk1', ...
%            '/data/Jeremy/Sleepscoring_Data_Noah/Canute/1000LuxWk4'};

% Define directories for different experimental conditions
baseDirs = {'/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/300Lux', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/1000LuxWk1', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/1000LuxWk4'};

% Initialize a structure to hold the results for each condition
results = struct();

for b = 1:length(baseDirs)
    baseFolder = baseDirs{b};
    [~, condition] = fileparts(baseFolder); % Extract condition name (e.g., '300Lux')
    [~, animalName] = fileparts(fileparts(baseFolder)); % Extract animal name (e.g., 'Canute')

    % Modify condition name to be valid field name
    validCondition = ['Cond', condition];

    % Initialize variables for pooled data for this condition
    % Using cell arrays to handle different sizes of data per channel
    pooledBinnedSpec = {}; 
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

            % Get channel numbers
            channels = zeros(1, length(data.StateInfo.fspec));
            for i = 1:length(channels)
                channels(1, i) = data.StateInfo.fspec{1,i}.info.Ch;
            end

            % Initialize storage for this file's binned data
            binnedSpecs = cell(1, length(channels)); % Channels
            nBins = cell(1, length(channels)); % Channels

            % Loop over the specified channels
            for i = 1:length(channels)
                % Extract relevant fields
                spec = data.StateInfo.fspec{1, i}.spec; % Spectrogram data (power or magnitude)
                times = data.StateInfo.fspec{1, i}.to; % Time points in seconds
                freqs = data.StateInfo.fspec{1, i}.fo; % Frequencies in Hz

                % Adjust times based on the initial start time and lights on hour
                initialDatetime = startDatetime - hours(lightsOnHour);
                adjustedDatetimes = initialDatetime + seconds(times);

                % Determine ZT hours for each data point
                adjustedHours = hour(adjustedDatetimes); % Extract the hour part of the adjusted times

                % Initialize variables for binning
                numOfHours = 24; % ZT 0 to 23
                binnedSpec = zeros(numOfHours, length(freqs)); % Rows for each ZT hour, columns for frequencies
                nBin = zeros(1, numOfHours); % Count of data points in each bin

                % Average power for each frequency in each ZT hour bin
                for zt = 0:23
                    indices = (adjustedHours == zt);
                    if any(indices)
                        binnedSpec(zt + 1, :) = mean(spec(indices, :), 1);
                        nBin(zt + 1) = sum(indices);
                    end
                end

                % If pooledBinnedSpec for a channel is empty, initialize it with the current binnedSpec
                if isempty(pooledBinnedSpec)
                    pooledBinnedSpec{i} = binnedSpec;
                else
                    if length(pooledBinnedSpec) < i || isempty(pooledBinnedSpec{i})
                        pooledBinnedSpec{i} = binnedSpec;
                    else
                        % Aggregate current file's data into the pooled data for the channel
                        pooledBinnedSpec{i} = pooledBinnedSpec{i} + binnedSpec;
                    end
                end

                % If nBins for a channel is empty, initialize it with the current nBin
                if length(nBins) < i || isempty(nBins{i})
                    nBins{i} = nBin;
                else
                    nBins{i} = nBins{i} + nBin;
                end
            end
        end
    end

    % No normalization, just accumulate data directly

    % Store results for the condition
    results.(validCondition).BinnedSpec = pooledBinnedSpec;
    results.(validCondition).Freqs = freqs;
    results.(validCondition).Channels = channels;

    %% Plotting spectrogram for the condition
    numOfHours = 24;
    zt_labels = cellstr(num2str((0:numOfHours-1)')); % Create ZT labels
    for i = 1:length(channels)
        figure;
        imagesc(0:23, freqs, squeeze(pooledBinnedSpec{i})'); % Convert ZT hours to x-axis
        axis xy;
        xlabel('Zeitgeber Time (ZT)');
        xticks(0:23); % ZT 0 to 23
        xticklabels(zt_labels);
        ylabel('Frequency (Hz)');
        title(['Spectrogram - Channel ', num2str(channels(i)), ' - ', condition]);
        colorbar;
        colormap('jet'); 
        set(gca, 'FontSize', 14);
        grid on;
    end
end

%% Comparisons Across Conditions
conditions = fieldnames(results);
comparisonMetrics = {'BinnedSpec'};
numOfHours = 24;
zt_labels = cellstr(num2str((0:numOfHours-1)'));

% Plot comparisons
for metric = 1:length(comparisonMetrics)
    metricName = comparisonMetrics{metric};
    for i = 1:length(results.(conditions{1}).Channels)
        figure;
        hold on;
        for c = 1:length(conditions)
            condition = conditions{c};
            channels = results.(condition).Channels;
            plotData = squeeze(results.(condition).(metricName){i});
            plot(0:numOfHours-1, mean(plotData, 2), 'DisplayName', [condition ' - Channel ' num2str(channels(i))]); % Mean across frequencies for each ZT hour
        end
        hold off;
        ylabel('Average Power');
        xlabel('Zeitgeber Time (ZT)');
        title(['Comparison of ', metricName, ' Across Conditions - Channel ', num2str(channels(i))]);
        xticks(0:23);
        xticklabels(zt_labels);
        legend show;
        set(gca, 'FontSize', 14);
        grid on;
    end
end