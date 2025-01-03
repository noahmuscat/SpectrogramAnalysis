%% Makes pretty, heat mapped spectrograms
% file comes from TheStateEditor

% Load the data
rootPath = '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/SleepStuff/SleepScoringFilesScatha/Canute_231210_041514/Canute_231210.eegstates.mat';
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

    channel = data.StateInfo.fspec{1, i}.info.Ch;

    % Plotting spectrogram (binned per hour)
    zt_labels = cellstr(num2str((0:23)')); % Create ZT labels
    figure;
    imagesc(0:numOfHours-1, freqs, binnedSpec'); % Convert ZT hours to the x-axis
    axis xy;
    xlabel('Zeitgeber Time (ZT)');
    xticks(0:23); % ZT 0 to 23
    xticklabels(zt_labels);
    ylabel('Frequency (Hz)');
    title(['Spectrogram - Channel ', num2str(channel)]);
    colorbar;
    colormap('jet'); 
    set(gca, 'FontSize', 14);
    grid on;
end