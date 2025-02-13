%% Combined Script for Histograms and Heat Maps of Spectrograms

% Define channels to analyze (can use [80, 112] for both)
channelsToAnalyze = 80;
threshold = 10; % for getting rid of 60 and 120 Hz (10 indicates getting rid of 50-70 and 110-130)

% Define directories for different experimental conditions
baseDirs = {'/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/300Lux', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/1000LuxWk1', ...
    '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute/1000LuxWk4'};

% Initialize a structure to hold the results for each condition
resultsSpectrogram = struct();

for b = 1:length(baseDirs)
    baseFolder = baseDirs{b};
    [~, condition] = fileparts(baseFolder); % Extract condition name (e.g., '300Lux')
    [~, animalName] = fileparts(fileparts(baseFolder)); % Extract animal name (e.g., 'Canute')

    % Modify condition name to be valid field name
    validCondition = ['Cond', condition];

    % Initialize variables for pooled data for this condition
    pooledBands = [];
    pooledEpochs = [];
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
            data = load(rootPath);  % Load the .mat file

            % Extract start time
            [~, folderName, ~] = fileparts(fileparts(rootPath));
            startDatetime = datetime(folderName(8:end), 'InputFormat', 'yyMMdd_HHmmss', 'TimeZone', 'America/New_York');

            % Determine the DST offset
            isDST = isdst(startDatetime);
            if isDST
                lightsOnHour = 6;
            else
                lightsOnHour = 5;
            end

            % Get channel numbers and filter by desired channels
            channels = zeros(1, length(data.StateInfo.fspec));
            for i = 1:length(channels)
                channels(1, i) = data.StateInfo.fspec{1,i}.info.Ch;
            end

            % Select only the specified channels to analyze
            channelIndicesToAnalyze = ismember(channels, channelsToAnalyze);

            % Initialize storage for this file's binned data
            binnedSpecs = cell(1, sum(channelIndicesToAnalyze)); % Channels to analyze
            nBins = cell(1, sum(channelIndicesToAnalyze));

            % Convert StateInfo to struct array specs
            specs = SaveSpectrogramsAsStruct(data.StateInfo);

            % Process the spectrogram data
            [bands, epochs] = PowerFreqFromSpecFreqInator(specs, startDatetime);

            % Aggregate data for histogram-based analysis
            pooledBands = aggregateBands(pooledBands, bands);
            pooledEpochs = aggregateEpochs(pooledEpochs, epochs);

            % Process heat map-based data for selected channels
            count = 1;
            for i = find(channelIndicesToAnalyze)
                % Extract relevant fields from the spectrogram data
                spec = data.StateInfo.fspec{1, i}.spec;  % Spectrogram data (power or magnitude)
                times = data.StateInfo.fspec{1, i}.to;   % Time points in seconds
                freqs = data.StateInfo.fspec{1, i}.fo;   % Frequencies in Hz
                
                % Find indices of frequencies close to 60 Hz and 120 Hz
                freq_indices_60Hz = find(abs(freqs - 60) <= threshold);
                freq_indices_120Hz = find(abs(freqs - 120) <= threshold);
                
                % Create a copy of the spectral data to modify
                clean_spec = spec;
                
                % Zero out or attenuate power at 60 Hz
                clean_spec(:, freq_indices_60Hz) = NaN;
                
                % Zero out or attenuate power at 120 Hz
                clean_spec(:, freq_indices_120Hz) = NaN;

                % Adjust times based on the initial start time and lights on hour
                initialDatetime = startDatetime - hours(lightsOnHour);
                adjustedDatetimes = initialDatetime + seconds(times);

                % Determine ZT hours for each data point
                adjustedHours = hour(adjustedDatetimes);  % Extract the hour part of the adjusted times

                % Initialize variables for binning
                numOfHours = 24;  % ZT 0 to 23
                binnedSpec = zeros(numOfHours, length(freqs));  % Rows for each ZT hour, columns for frequencies
                nBin = zeros(1, numOfHours);  % Count of data points in each bin

                % Average power for each frequency in each ZT hour bin
                for zt = 0:23
                    indices = (adjustedHours == zt);
                    if any(indices)
                        binnedSpec(zt + 1, :) = mean(clean_spec(indices, :), 1);
                        nBin(zt + 1) = sum(indices);
                    end
                end

                % If pooledBinnedSpec for a channel is empty, initialize it with the current binnedSpec
                if isempty(pooledBinnedSpec)
                    pooledBinnedSpec{count} = binnedSpec;
                else
                    if length(pooledBinnedSpec) < count || isempty(pooledBinnedSpec{count})
                        pooledBinnedSpec{count} = binnedSpec;
                    else
                        % Aggregate current file's data into the pooled data for the channel
                        pooledBinnedSpec{count} = pooledBinnedSpec{count} + binnedSpec;
                    end
                end

                % If nBins for a channel is empty, initialize it with the current nBin
                if length(nBins) < count || isempty(nBins{count})
                    nBins{count} = nBin;
                else
                    nBins{count} = nBins{count} + nBin;
                end

                count = count + 1;
            end
        end
    end

    % Store results for the condition
    resultsSpectrogram.(validCondition).Bands = pooledBands;
    resultsSpectrogram.(validCondition).Epochs = pooledEpochs;
    resultsSpectrogram.(validCondition).BinnedSpec = pooledBinnedSpec;
    resultsSpectrogram.(validCondition).Freqs = freqs;
    resultsSpectrogram.(validCondition).Channels = channels(channelIndicesToAnalyze);

    %% Plotting histogram-based results for the condition
    plotPowerVectors(pooledBands, pooledEpochs.HourlyBinIndices, channels(channelIndicesToAnalyze), condition, channelsToAnalyze);
    plotPercentOscillatoryPower(pooledBands, pooledEpochs.HourlyBinIndices, channels(channelIndicesToAnalyze), condition, channelsToAnalyze);

    %% Plotting heat map-based spectrogram for the condition
    numOfHours = 24;
    zt_labels = cellstr(num2str((0:numOfHours-1)'));  % Create ZT labels
    count = 1;
    for i = find(channelIndicesToAnalyze)
        figure;
        imagesc(0:23, freqs, squeeze(pooledBinnedSpec{count})');  % Convert ZT hours to x-axis
        axis xy;
        xlabel('Zeitgeber Time (ZT)');
        xticks(0:23);  % ZT 0 to 23
        xticklabels(zt_labels);
        ylabel('Frequency (Hz)');
        title(['Spectrogram - Channel ', num2str(channels(i)), ' - ', condition]);
        colorbar;
        colormap('jet'); 
        set(gca, 'FontSize', 14);
        grid on;
        count = count + 1;
    end
end

%% saving .mat
matFileName = 'spectrogramMetrics.mat';
matFolderPath = '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/SampleFiles/Canute';
matFilePath = fullfile(matFolderPath, matFileName);
save(matFilePath, "resultsSpectrogram");

%% Comparisons Across Conditions
conditions = fieldnames(resultsSpectrogram);
bandnames = {resultsSpectrogram.(conditions{1}).Bands.name};
numOfHours = 24;
zt_labels = cellstr(num2str((0:numOfHours-1)'));

% Iterate over specified channels to analyze
for i = 1:length(resultsSpectrogram.(conditions{1}).Channels)
    channel = resultsSpectrogram.(conditions{1}).Channels(i);

    if ismember(channel, channelsToAnalyze)  % Check if the channel is in the specified channelsToAnalyze
        % Histogram-based power vectors comparison
        figure;
        sgtitle(['Power Across Conditions - Channel ', num2str(channel)]); 
        for b = 1:length(bandnames)
            subplot(4, 3, b);
            hold on;
            for c = 1:length(conditions)
                condition = conditions{c};
                channels = resultsSpectrogram.(condition).Channels;
                bands = resultsSpectrogram.(condition).Bands;

                % Find the index of the current channel within this condition
                channelIndex = find(channels == channel, 1);
                if ~isempty(channelIndex)
                    powerData = bands(b).powervectors.HourlyBin(channelIndex, :);
                    plotData = cellfun(@mean, powerData);
                    plot(0:numOfHours-1, plotData, 'DisplayName', condition);
                end
            end
            hold off;
            xlabel('Zeitgeber Time (ZT)');
            ylabel('Power');
            title(bandnames{b});
            legend show;
        end

        % Histogram-based percent oscillatory power comparison
        figure;
        sgtitle(['Percent Oscillatory Power Across Conditions - Channel ', num2str(channel)]);
        for b = 1:length(bandnames)
            subplot(4, 3, b);
            hold on;
            for c = 1:length(conditions)
                condition = conditions{c};
                channels = resultsSpectrogram.(condition).Channels;
                bands = resultsSpectrogram.(condition).Bands;

                totalPower = zeros(numOfHours, 1);
                percentagePower = zeros(numOfHours, 1);
                % Identify the correct indices to refer to selected channel
                channelIndex = find(channels == channel, 1);
                if ~isempty(channelIndex)
                    for zt = 0:(numOfHours-1)
                        if ~isempty(bands(b).powervectors.HourlyBin{channelIndex, zt+1})
                            totalPower(zt+1) = totalPower(zt+1) + mean(cellfun(@mean, bands(b).powervectors.HourlyBin(:, zt+1)));
                        end
                    end

                    for zt = 0:(numOfHours-1)
                        if totalPower(zt+1) > 0
                            if ~isempty(bands(b).powervectors.HourlyBin{channelIndex, zt+1})
                                percentagePower(zt+1) = 100 * mean(cellfun(@mean, bands(b).powervectors.HourlyBin(channelIndex, zt+1))) / totalPower(zt+1);
                            end
                        end
                    end
                    plot(0:numOfHours-1, percentagePower, 'DisplayName', num2str(channel));
                end
            end
            hold off;
            xlabel('Zeitgeber Time (ZT)');
            ylabel('Percent Power (%)');
            title(bandnames{b});
            legend show;
        end

        % Heat map-based spectrogram comparison
        figure;
        hold on;
        for c = 1:length(conditions)
            condition = conditions{c};
            channels = resultsSpectrogram.(condition).Channels;

            % Locate the correct index for selecting the current channel
            channelIndex = find(channels == channel, 1);
            if ~isempty(channelIndex)
                plotData = squeeze(resultsSpectrogram.(condition).BinnedSpec{channelIndex});
                plot(0:numOfHours-1, mean(plotData, 2), 'DisplayName', num2str(channel));  % Mean across frequencies for each ZT hour
            end
        end
        hold off;
        ylabel('Average Power');
        xlabel('Zeitgeber Time (ZT)');
        title(['Comparison of BinnedSpec Across Conditions - Channel ', num2str(channel)]);
        xticks(0:23);
        xticklabels(zt_labels);
        legend show;
        set(gca, 'FontSize', 14);
        grid on;
    end
end

%% Sub-functions

function specs = SaveSpectrogramsAsStruct(StateInfo)
    specs = struct('spec', {}, 'freqs', {}, 'times', {});

    for a = 1:length(StateInfo.fspec)
        specs(a).spec = StateInfo.fspec{a}.spec;
        specs(a).freqs = StateInfo.fspec{a}.fo;
        specs(a).times = StateInfo.fspec{a}.to;    
    end
end

function [bands, epochs] = PowerFreqFromSpecFreqInator(specs, startDatetime)
    if ~exist('startDatetime', 'var') || isempty(startDatetime)
        error('startDatetime is required as an input.');
    end
    
    isDST = isdst(startDatetime);
    if isDST
        lightsOnHour = 6;
    else
        lightsOnHour = 5;
    end

    if ~isstruct(specs) || isempty(specs)
        error('Input specs must be a non-empty struct array.');
    end

    requiredFields = {'spec', 'freqs', 'times'};
    for i = 1:length(specs)
        if ~all(isfield(specs(i), requiredFields))
            error('Each element of specs must contain fields: spec, freqs, and times.');
        end
    end

    bands = initializeBands();

    initialDatetime = startDatetime + hours(lightsOnHour);
    timestamps = specs(1).times;
    adjustedDatetimes = initialDatetime + seconds(timestamps);
    adjustedHours = hour(adjustedDatetimes);

    numOfHours = 24;
    epochs.HourlyBinIndices = cell(1, numOfHours);
    for zt = 0:(numOfHours - 1)
        epochs.HourlyBinIndices{zt + 1} = find(adjustedHours == zt);
    end

    for b = 1:length(bands)
        bands(b).powervectors.All = cell(length(specs), 1);
        bands(b).powervectors.HourlyBin = cell(length(specs), numOfHours);
    end

    for b = 1:length(bands)
        tbandname = bands(b).name;
        tband = bands(b);
        tband.freqidxs = find(specs(1).freqs >= tband.startstopfreq(1) & specs(1).freqs < tband.startstopfreq(2));

        if isempty(tband.freqidxs) && ~strcmp(tbandname, 'thetaratio')
            warning(['No frequencies found in the range for band: ', tbandname]);
            continue;
        end
        
        for a = 1:length(specs)
            if strcmp(tbandname, 'thetaratio')
                [delta_power, spindle_power, theta_power] = getBandPowers(bands, a);
                if isempty(delta_power) || isempty(spindle_power) || isempty(theta_power)
                    tband.powervectors.All{a} = [];
                else
                    tband.powervectors.All{a} = theta_power ./ (delta_power + spindle_power);
                end
            else
                specData = specs(a).spec(:, tband.freqidxs);
                tband.powervectors.All{a} = mean(specData, 2);
            end

            for zt = 0:(numOfHours - 1)
                indices = epochs.HourlyBinIndices{zt + 1};
                if ~isempty(indices)
                    tband.powervectors.HourlyBin{a, zt + 1} = mean(tband.powervectors.All{a}(indices));
                end
            end
        end
        
        fields = fieldnames(tband);
        for f = 1:length(fields)
            bands(b).(fields{f}) = tband.(fields{f});
        end
    end
end

function bands = initializeBands()
    bandNames = {'delta', 'theta', 'spindle', 'lowbeta', 'highbeta', 'lowgamma', 'midgamma', 'highgamma', 'ripple', 'thetaratio'};
    startStopFreqs = {[0.5 4], [5 10], [11 19], [20 30], [30 40], [40 60], [60 100], [100 140], [140 180], [0 0]};
    
    bands = struct('name', {}, 'startstopfreq', {}, 'freqidxs', {}, 'powervectors', {});
    for i = 1:length(bandNames)
        bands(i).name = bandNames{i};
        bands(i).startstopfreq = startStopFreqs{i};
        bands(i).freqidxs = [];
        bands(i).powervectors.All = {};
        bands(i).powervectors.HourlyBin = cell(1, 24);
    end
end

function [delta_power, spindle_power, theta_power] = getBandPowers(bands, channel)
    delta_power = bands(1).powervectors.All{channel};
    spindle_power = bands(3).powervectors.All{channel};
    theta_power = bands(2).powervectors.All{channel};
end

function plotPowerVectors(bands, hourlyBinIndices, channels, condition, channelsToAnalyze)
    bandnames = {bands.name};
    numOfHours = 24;
    zt_labels = cellstr(num2str((0:numOfHours-1)'));

    for a = 1:length(channels)
        channel = channels(a);
        
        if ismember(channel, channelsToAnalyze)  % Check if the channel is in the specified channelsToAnalyze
            figure;
            sgtitle(['Power Across Bands - Channel ', num2str(channel), ' - ', condition]);
            plotcounter = 0;
            for b = 1:length(bandnames)
                tbandname = bandnames{b};
                tband = bands(b);
                plotcounter = plotcounter + 1;
                subplot(4, 3, plotcounter);
                hourlyPower = zeros(1, numOfHours);
                for zt = 0:(numOfHours-1)
                    if ~isempty(hourlyBinIndices{zt+1})
                        hourlyPower(zt+1) = mean(tband.powervectors.HourlyBin{a, zt+1});
                    end
                end
                plot(0:numOfHours-1, hourlyPower);
                addShadedAreaToPlotZT24Hour();
                xlabel('Zeitgeber Time (ZT)');
                xticks(0:numOfHours-1);
                xticklabels(zt_labels);
                ylabel('Power');
                title(tbandname);
                axis tight;
            end
        end
    end
end

function plotPercentOscillatoryPower(bands, hourlyBinIndices, channels, condition, channelsToAnalyze)
    bandnames = {bands.name};
    numOfHours = 24;
    zt_labels = cellstr(num2str((0:numOfHours-1)'));

    for a = 1:length(channels)
        channel = channels(a);
        
        if ismember(channel, channelsToAnalyze)  % Check if the channel is in the specified channelsToAnalyze
            figure;
            sgtitle(['% Oscillatory Power - Channel ', num2str(channel), ' - ', condition]);
            plotcounter = 0;

            totalPower = zeros(numOfHours, 1);
            for b = 1:length(bands)
                for zt = 0:(numOfHours-1)
                    if ~isempty(hourlyBinIndices{zt+1})
                        totalPower(zt+1) = totalPower(zt+1) + mean(bands(b).powervectors.HourlyBin{a, zt+1});
                    end
                end
            end
            
            for b = 1:length(bandnames)
                tbandname = bandnames{b};
                tband = bands(b);
                plotcounter = plotcounter + 1;
                subplot(4, 3, plotcounter);
                percentPower = zeros(1, numOfHours);
                for zt = 0:(numOfHours-1)
                    if ~isempty(hourlyBinIndices{zt+1}) && totalPower(zt+1) > 0
                        percentPower(zt+1) = 100 * mean(tband.powervectors.HourlyBin{a, zt+1}) / totalPower(zt+1);
                    end
                end
                plot(0:numOfHours-1, percentPower);
                addShadedAreaToPlotZT24Hour();
                xlabel('Zeitgeber Time (ZT)');
                xticks(0:numOfHours-1);
                xticklabels(zt_labels);
                ylabel('Percent Power (%)');
                title(tbandname);
                axis tight;
            end
        end
    end
end

function addShadedAreaToPlotZT24Hour()
    hold on;
    x_shaded = [12, 24, 24, 12];
    y_lim = ylim;
    y_shaded = [y_lim(1), y_lim(1), y_lim(2), y_lim(2)];
    
    fill_color = [0.7, 0.7, 0.7]; % Light gray color for the shading
    h = fill(x_shaded, y_shaded, fill_color, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(h, 'bottom');
    
    xlabel('Hour of Day (ZT Time)');
    ylabel('Sum of PixelDifference');
    xlim([-0.5, 23.5]);
    xticks(0:23);
    xtickangle(0);
    hold off;
end

function aggregatedBands = aggregateBands(pooledBands, newBands)
    if isempty(pooledBands)  % If the first set of bands
        aggregatedBands = newBands;
        return;
    end

    for b = 1:length(pooledBands)
        for h = 1:length(newBands(b).powervectors.HourlyBin)
            if isempty(pooledBands(b).powervectors.HourlyBin{h})
                pooledBands(b).powervectors.HourlyBin{h} = newBands(b).powervectors.HourlyBin{h};
            else
                pooledBands(b).powervectors.HourlyBin{h} = ...
                    pooledBands(b).powervectors.HourlyBin{h} + newBands(b).powervectors.HourlyBin{h};
            end
        end
    end
    aggregatedBands = pooledBands;
end

function aggregatedEpochs = aggregateEpochs(pooledEpochs, newEpochs)
    if isempty(pooledEpochs)
        aggregatedEpochs = newEpochs;
        return;
    end

    fields = fieldnames(newEpochs);
    for i = 1:length(fields)
        pooledEpochs.(fields{i}) = [pooledEpochs.(fields{i}), newEpochs.(fields{i})];
    end
    aggregatedEpochs = pooledEpochs;
end

