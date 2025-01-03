%% Makes ugly, histogram based spectrograms
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
pooledSpecs = struct('spec', {}, 'freqs', {}, 'times', {}, 'info', {});
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

        % Convert StateInfo to struct array specs
        specs = SaveSpectrogramsAsStruct(data.StateInfo);

        % Pool the specs data
        pooledSpecs = aggregateSpecs(pooledSpecs, specs, startDatetime, firstRun);
        firstRun = false;
    end
end

% Get channel numbers
if ~isempty(pooledSpecs)
    channels = zeros(1, length(pooledSpecs));
    for i = 1:length(channels)
        channels(1, i) = pooledSpecs(i).info.Ch;
    end

    % Process the pooled spectrogram data
    PowerFreqFromSpecFreqInator(pooledSpecs, channels, stdDevThreshold);
end

%% functions
function specs = SaveSpectrogramsAsStruct(StateInfo)
    specs = struct('spec', {}, 'freqs', {}, 'times', {}, 'info', {});

    for a = 1:length(StateInfo.fspec)
        specs(a).spec = StateInfo.fspec{a}.spec;
        specs(a).freqs = StateInfo.fspec{a}.fo;
        specs(a).times = StateInfo.fspec{a}.to;    
        specs(a).info = StateInfo.fspec{a}.info;
    end
end

function pooledSpecs = aggregateSpecs(pooledSpecs, specs, startDatetime, firstRun)
    % Determine DST and lightsOnHour
    isDST = isdst(startDatetime);
    if isDST
        lightsOnHour = 6;
    else
        lightsOnHour = 5;
    end

    % Adjust times based on the initial start time and lights on hour
    initialDatetime = startDatetime - hours(lightsOnHour);
    for a = 1:length(specs)
        adjustedDatetimes = initialDatetime + seconds(specs(a).times);
        
        % Determine ZT hours for each data point
        adjustedHours = hour(adjustedDatetimes); % Extract the hour part of the adjusted times

        if firstRun
            pooledSpecs(a).spec = zeros(24, length(specs(a).freqs));
            pooledSpecs(a).times = (0:23)'; % Zeitgeber time hours
            pooledSpecs(a).freqs = specs(a).freqs;
            pooledSpecs(a).info = specs(a).info;
            pooledSpecs(a).count = zeros(24, 1); % Count of data points in each bin
        end

        % Average power for each frequency in each ZT hour bin
        for zt = 0:23
            indices = (adjustedHours == zt);
            if any(indices)
                pooledSpecs(a).spec(zt + 1, :) = pooledSpecs(a).spec(zt + 1, :) + sum(specs(a).spec(indices, :), 1);
                pooledSpecs(a).count(zt + 1) = pooledSpecs(a).count(zt + 1) + sum(indices);
            end
        end
    end
end

function [bands, epochs] = PowerFreqFromSpecFreqInator(pooledSpecs, channels, stdDevThreshold)
    % Validate input structure
    if ~isstruct(pooledSpecs) || isempty(pooledSpecs)
        error('Input pooledSpecs must be a non-empty struct array.');
    end
    
    % Ensure channel numbers are available
    if isempty(channels)
        if isfield(pooledSpecs, 'info') && isfield(pooledSpecs(1).info, 'Ch')
            channels = arrayfun(@(x) x.info.Ch, pooledSpecs);
        else
            error('Channels information is required.');
        end
    end

    numOfHours = 24;

    % Calculate average spectrogram for pooled data
    for i = 1:length(pooledSpecs)
        for zt = 0:23
            if pooledSpecs(i).count(zt + 1) > 0
                pooledSpecs(i).spec(zt + 1, :) = pooledSpecs(i).spec(zt + 1, :) / pooledSpecs(i).count(zt + 1);
            end
        end
    end

    % Remove outliers based on specified standard deviation threshold
    if ~isempty(stdDevThreshold)
        for i = 1:length(pooledSpecs)
            meanSpec = mean(pooledSpecs(i).spec, 1);
            stdSpec = std(pooledSpecs(i).spec, 0, 1);

            for zt = 0:23
                if pooledSpecs(i).count(zt + 1) > 0
                    outliers = abs(pooledSpecs(i).spec(zt + 1, :) - meanSpec) > (stdDevThreshold * stdSpec);
                    pooledSpecs(i).spec(zt + 1, outliers) = meanSpec(outliers);  % Replace outliers with mean
                end
            end
        end
    end

    % Proceed with existing function logic
    bands = initializeBands();
    epochs.HourlyBinIndices = cell(1, numOfHours);
    for zt = 0:numOfHours-1
        epochs.HourlyBinIndices{zt+1} = find(pooledSpecs(1).times == zt);
    end

    % Initialize power vectors for each band
    for b = 1:length(bands)
        bands(b).powervectors.All = cell(length(pooledSpecs), 1);
        bands(b).powervectors.HourlyBin = cell(length(pooledSpecs), numOfHours);
    end

    % Extracting powers for each frequency band
    for b = 1:length(bands)
        tbandname = bands(b).name;
        tband = bands(b);
        tband.freqidxs = find(pooledSpecs(1).freqs >= tband.startstopfreq(1) & pooledSpecs(1).freqs < tband.startstopfreq(2));

        if isempty(tband.freqidxs) && ~strcmp(tbandname, 'thetaratio')
            warning(['No frequencies found in the range for band: ', tbandname]);
            continue;
        end
        
        % Loop over pooled specs
        for a = 1:length(pooledSpecs)
            if strcmp(tbandname, 'thetaratio')
                [delta_power, spindle_power, theta_power] = getBandPowers(bands, a);
                if isempty(delta_power) || isempty(spindle_power) || isempty(theta_power)
                    tband.powervectors.All{a} = [];
                else
                    tband.powervectors.All{a} = theta_power ./ (delta_power + spindle_power);
                end
            else
                specData = pooledSpecs(a).spec(:, tband.freqidxs);
                tband.powervectors.All{a} = mean(specData, 2);
            end

            for zt = 0:numOfHours - 1
                indices = epochs.HourlyBinIndices{zt + 1};
                if ~isempty(indices)
                    tband.powervectors.HourlyBin{a, zt + 1} = mean(tband.powervectors.All{a}(indices));
                end
            end
        end
        
        % Ensuring struct consistency
        fields = fieldnames(tband);
        for field = fields'
            bands(b).(field{1}) = tband.(field{1});
        end
    end

    % Plotting with hourly bins for each channel
    plotPowerVectors(pooledSpecs, bands, epochs.HourlyBinIndices, channels);
    plotPercentOscillatoryPower(pooledSpecs, bands, epochs.HourlyBinIndices, channels);
end

function bands = initializeBands()
    % Initialize frequency bands
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

function plotPowerVectors(pooledSpecs, bands, hourlyBinIndices, channels)
    bandnames = {bands.name};
    numOfHours = 24;
    zt_labels = cellstr(num2str((0:numOfHours-1)')); % Create ZT hour labels
    
    % Create a figure for each pooled spec
    for a = 1:length(pooledSpecs)
        figure;
        sgtitle(['Power Across Bands - Channel ', num2str(channels(a))]); % Using sgtitle for a single title for the figure
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

function plotPercentOscillatoryPower(pooledSpecs, bands, hourlyBinIndices, channels)
    bandnames = {bands.name};
    numOfHours = 24;
    zt_labels = cellstr(num2str((0:numOfHours-1)')); % Create ZT hour labels
    
    % Initialize figure for each pooled spec
    for a = 1:length(pooledSpecs)
        figure;
        sgtitle(['% Oscillatory Power - Channel ', num2str(channels(a))]); % Using sgtitle for a single title for the figure
        plotcounter = 0;

        % Calculate total power across all bands for each hour
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

            % Calculate percentage power for each band for each hour
            percentPower = zeros(1, numOfHours);
            for zt = 0:(numOfHours-1)
                if ~isempty(hourlyBinIndices{zt+1}) && totalPower(zt+1) > 0
                    percentPower(zt+1) = 100 * mean(tband.powervectors.HourlyBin{a, zt+1}) / totalPower(zt+1);
                end
            end

            % Plot percentage power
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

function addShadedAreaToPlotZT24Hour()
    hold on;
    % Define x and y coordinates for the shaded area (from t=12 to t=24)
    x_shaded = [12, 24, 24, 12];
    y_lim = ylim;
    y_shaded = [y_lim(1), y_lim(1), y_lim(2), y_lim(2)];

    fill_color = [0.7, 0.7, 0.7]; % Light gray color for the shading

    % Add shaded areas to the plot with 'HandleVisibility', 'off' to exclude from the legend
    h = fill(x_shaded, y_shaded, fill_color, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % Use uistack to send the shaded area to the back
    uistack(h, 'bottom');

    % Additional Plot settings
    xlabel('Hour of Day (ZT Time)');
    ylabel('Sum of PixelDifference');
    xlim([-0.5, 23.5]);
    xticks(0:23);
    xtickangle(0);
    
    hold off;
end