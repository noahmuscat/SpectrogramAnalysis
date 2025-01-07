%% Makes pretty, heat mapped spectrograms for All Conditions
baseDirs = {'/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/SleepStuff/SleepScoringFilesScatha/Canute/300Lux', ...
            '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/SleepStuff/SleepScoringFilesScatha/Canute/1000LuxWk1', ...
            '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/SleepStuff/SleepScoringFilesScatha/Canute/1000LuxWk4'};
        
% Initialize a structure to hold the results for each condition
results = struct();

for b = 1:length(baseDirs)
    baseDir = baseDirs{b};
    [~, condition] = fileparts(baseDir); % Extract condition name (e.g., '300Lux')

    % Get a list of all .eegstates.mat files in the baseDir
    matFiles = dir(fullfile(baseDir, '*.eegstates.mat'));

    % Initialize variables for pooled data for this condition
    pooledBinnedSpec = [];
    pooledNBins = [];
    numChannels = [];

    for k = 1:length(matFiles)
        rootPath = fullfile(matFiles(k).folder, matFiles(k).name);
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

        % Loop over the specified channels
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

            % Initialize variables for binning
            numOfHours = 24; % ZT 0 to 23
            binnedSpec = zeros(numOfHours, length(freqs)); % Rows for each ZT hour, columns for frequencies
            nBins = zeros(1, numOfHours); % Count of data points in each bin

            % Average power for each frequency in each ZT hour bin
            for zt = 0:(numOfHours - 1)
                indices = (adjustedHours == zt);
                if any(indices)
                    binnedSpec(zt + 1, :) = mean(spec(indices, :), 1);
                    nBins(zt + 1) = sum(indices);
                end
            end

            % Aggregate data across files
            if isempty(pooledBinnedSpec)
                pooledBinnedSpec = binnedSpec;
                pooledNBins = nBins;
            else
                pooledBinnedSpec = pooledBinnedSpec + binnedSpec;
                pooledNBins = pooledNBins + nBins;
            end
        end
    end

    % Normalize the pooled data
    for zt = 1:numOfHours
        if pooledNBins(zt) > 0
            pooledBinnedSpec(zt, :) = pooledBinnedSpec(zt, :) / pooledNBins(zt);
        end
    end

    % Store results for the condition
    results.(condition).BinnedSpec = pooledBinnedSpec;
    results.(condition).Freqs = freqs;
    results.(condition).NumChannels = numChannels;

    %% Plotting spectrogram for the condition
    zt_labels = cellstr(num2str((0:23)')); % Create ZT labels
    for i = 1:numChannels
        figure;
        imagesc(0:numOfHours-1, freqs, pooledBinnedSpec'); % Convert ZT hours to x-axis
        axis xy;
        xlabel('Zeitgeber Time (ZT)');
        xticks(0:23); % ZT 0 to 23
        xticklabels(zt_labels);
        ylabel('Frequency (Hz)');
        title(['Spectrogram - Channel ', num2str(i), ' - ', condition]);
        colorbar;
        colormap('jet'); 
        set(gca, 'FontSize', 14);
        grid on;
    end
end

%% Comparisons Across Conditions
conditions = fieldnames(results);
comparisonMetrics = {'BinnedSpec'};
states = conditions;

% Plot comparisons
for metric = 1:length(comparisonMetrics)
    metricName = comparisonMetrics{metric};
    for i = 1:numChannels
        figure;
        hold on;
        for c = 1:length(conditions)
            condition = conditions{c};
            plotData = results.(condition).(metricName);
            plot(0:numOfHours-1, mean(plotData, 2), 'DisplayName', condition); % Mean across frequencies for each ZT hour
        end
        hold off;
        ylabel('Average Power');
        xlabel('Zeitgeber Time (ZT)');
        title(['Comparison of ', metricName, ' Across Conditions - Channel ', num2str(i)]);
        xticks(0:23);
        xticklabels(zt_labels);
        legend show;
        set(gca, 'FontSize', 14);
        grid on;
    end
end