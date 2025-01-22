%% Makes ugly, histogram-based spectrograms for All Conditions
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
    pooledBands = [];
    pooledEpochs = [];
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

            % Convert StateInfo to struct array specs
            specs = SaveSpectrogramsAsStruct(data.StateInfo);

            % Get channel numbers
            channels = zeros(1, length(data.StateInfo.fspec));
            for i = 1:length(channels)
                channels(1, i) = data.StateInfo.fspec{1,i}.info.Ch;
            end

            % Process the spectrogram data
            [bands, epochs] = PowerFreqFromSpecFreqInator(specs, startDatetime);

            % Aggregate data
            pooledBands = aggregateBands(pooledBands, bands);
            pooledEpochs = aggregateEpochs(pooledEpochs, epochs);
        end
    end

    % Store results for the condition
    results.(validCondition).Bands = pooledBands;
    results.(validCondition).Epochs = pooledEpochs;
    results.(validCondition).Channels = channels;

    %% Plotting for individuals
    plotPowerVectors(specs, pooledBands, pooledEpochs.HourlyBinIndices, channels, condition);
    plotPercentOscillatoryPower(specs, pooledBands, pooledEpochs.HourlyBinIndices, channels, condition);
end

%% Comparisons Across Conditions
conditions = fieldnames(results);
numOfHours = 24;
bandnames = {results.(conditions{1}).Bands.name};

% Plot comparisons for power vectors
for i = 1:length(results.(conditions{1}).Channels)
    figure;
    sgtitle(['Power Across Conditions - Channel ', num2str(results.(conditions{1}).Channels(i))]); 
    for b = 1:length(bandnames)
        subplot(4, 3, b);
        hold on;
        for c = 1:length(conditions)
            condition = conditions{c};
            channels = results.(condition).Channels;
            bands = results.(condition).Bands;

            % Assuming this channel for comparison
            if i <= length(channels)
                powerData = bands(b).powervectors.HourlyBin(i, :);
                plotData = cellfun(@mean, powerData);
                plot(0:numOfHours-1, plotData, 'DisplayName', [condition ' - Channel ' num2str(channels(i))]);
            end
        end
        hold off;
        xlabel('Zeitgeber Time (ZT)');
        ylabel('Power');
        title(bandnames{b});
        legend show;
    end
end

% Plot comparisons for percent oscillatory power
for i = 1:length(results.(conditions{1}).Channels)
    figure;
    sgtitle(['Percent Oscillatory Power Across Conditions - Channel ', num2str(results.(conditions{1}).Channels(i))]);
    for b = 1:length(bandnames)
        subplot(4, 3, b);
        hold on;
        for c = 1:length(conditions)
            condition = conditions{c};
            channels = results.(condition).Channels;
            bands = results.(condition).Bands;

            totalPower = zeros(numOfHours, 1);
            percentagePower = zeros(numOfHours, 1);
            for zt = 0:(numOfHours-1)
                if ~isempty(bands(b).powervectors.HourlyBin{i, zt+1})
                    totalPower(zt+1) = totalPower(zt+1) + mean(cellfun(@mean, bands(b).powervectors.HourlyBin(:, zt+1)));
                end
            end

            for zt = 0:(numOfHours-1)
                if totalPower(zt+1) > 0
                    if ~isempty(bands(b).powervectors.HourlyBin{i, zt+1})
                        percentagePower(zt+1) = 100 * mean(cellfun(@mean, bands(b).powervectors.HourlyBin(i, zt+1))) / totalPower(zt+1);
                    end
                end
            end
            plot(0:numOfHours-1, percentagePower, 'DisplayName', [condition ' - Channel ' num2str(channels(i))]);
        end
        hold off;
        xlabel('Zeitgeber Time (ZT)');
        ylabel('Percent Power (%)');
        title(bandnames{b});
        legend show;
    end
end

%% functions
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

function plotPowerVectors(specs, bands, hourlyBinIndices, channels, condition)
    bandnames = {bands.name};
    numOfHours = 24;
    zt_labels = cellstr(num2str((0:numOfHours-1)'));

    for a = 1:length(specs)
        figure;
        sgtitle(['Power Across Bands - Channel ', num2str(channels(a)), ' - ', condition]);
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

function plotPercentOscillatoryPower(specs, bands, hourlyBinIndices, channels, condition)
    bandnames = {bands.name};
    numOfHours = 24;
    zt_labels = cellstr(num2str((0:numOfHours-1)'));

    for a = 1:length(specs)
        figure;
        sgtitle(['% Oscillatory Power - Channel ', num2str(channels(a)), ' - ', condition]);
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
    if isempty(pooledBands)
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

