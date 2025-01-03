%% Makes ugly, histogram based spectrograms
% file comes from TheStateEditor

% Load the data
rootPath = '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/SleepStuff/SleepScoringFilesScatha/Canute_231210_041514/Canute_231210.eegstates.mat';
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
PowerFreqFromSpecFreqInator(specs, startDatetime, channels);

%% functions
function specs = SaveSpectrogramsAsStruct(StateInfo)
    specs = struct('spec', {}, 'freqs', {}, 'times', {});

    for a = 1:length(StateInfo.fspec)
        specs(a).spec = StateInfo.fspec{a}.spec;
        specs(a).freqs = StateInfo.fspec{a}.fo;
        specs(a).times = StateInfo.fspec{a}.to;    
    end
end

function [bands, epochs] = PowerFreqFromSpecFreqInator(specs, startDatetime, channels)
    % Default value for startDatetime
    if ~exist('startDatetime', 'var') || isempty(startDatetime)
        error('startDatetime is required as an input.');
    end
    
    % Determine DST and lightsOnHour
    isDST = isdst(startDatetime);
    if isDST
        lightsOnHour = 6;
    else
        lightsOnHour = 5;
    end

    % Validate input structure
    if ~isstruct(specs) || isempty(specs)
        error('Input specs must be a non-empty struct array.');
    end

    requiredFields = {'spec', 'freqs', 'times'};
    for i = 1:length(specs)
        if ~all(isfield(specs(i), requiredFields))
            error('Each element of specs must contain fields: spec, freqs, and times.');
        end
    end

    % Defining frequency bands
    bands = initializeBands();

    % Setting up time intervals
    initialDatetime = startDatetime + hours(lightsOnHour);
    timestamps = specs(1).times;
    adjustedDatetimes = initialDatetime + seconds(timestamps);
    adjustedHours = hour(adjustedDatetimes);

    % Prepare variables for binning
    numOfHours = 24;
    epochs.HourlyBinIndices = cell(1, numOfHours);
    for zt = 0:(numOfHours - 1)
        epochs.HourlyBinIndices{zt + 1} = find(adjustedHours == zt);
    end

    % Initialize power vectors for each band
    for b = 1:length(bands)
        bands(b).powervectors.All = cell(length(specs), 1);
        bands(b).powervectors.HourlyBin = cell(length(specs), numOfHours);
    end

    % Extracting powers for each frequency band
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
        
        % Avoiding subscripted assignment issue by ensuring bands struct consistency
        fields = fieldnames(tband);
        for field = fields'
            bands(b).(field{1}) = tband.(field{1});
        end
    end

    % Plotting with hourly bins
    plotPowerVectors(specs, bands, epochs.HourlyBinIndices, channels);
    plotPercentOscillatoryPower(specs, bands, epochs.HourlyBinIndices, channels);
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

function plotPowerVectors(specs, bands, hourlyBinIndices, channels)
    bandnames = {bands.name};
    numOfHours = 24;
    zt_labels = cellstr(num2str((0:numOfHours-1)')); % Create ZT hour labels
    
    for a = 1:length(specs)
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

function plotPercentOscillatoryPower(specs, bands, hourlyBinIndices, channels)
    bandnames = {bands.name};
    numOfHours = 24;
    zt_labels = cellstr(num2str((0:numOfHours-1)')); % Create ZT hour labels
    
    for a = 1:length(specs)
        % Initialize figure
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