load('/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/HaraldCombinedDataOutliersRemoved.mat');

% Metadata: fo contains the frequency values
fo = HaraldCombined.MetaData.fo;

% Define the frequency range to plot (0 to 40 Hz)
freqIndices = fo >= 0 & fo <= 40;
filteredFrequencies = fo(freqIndices);

% Conditions to loop through
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};

for iCond = 1:length(conditions)
    % Get the current condition
    condName = conditions{iCond};
    condData = HaraldCombined.(condName);

    % Extract ZscoredFrequencyPower and ZT_time information
    zscoredPower = condData.ZscoredFrequencyPower;
    ZT_time = condData.ZT_time;
    
    % Define hourly bins
    binEdges = floor(min(ZT_time)):1:ceil(max(ZT_time)); % 1-hour bins

    % Preallocate for binned power
    numBins = length(binEdges) - 1;
    binnedPower = zeros(sum(freqIndices), numBins);

    % Bin data based on ZT_time
    for bin = 1:numBins
        % Find indices within the current bin
        binInds = ZT_time >= binEdges(bin) & ZT_time < binEdges(bin + 1);
        
        % Average power for the current bin
        if any(binInds)
            % Filter frequencies to 0-40 Hz range and average
            binnedPower(:, bin) = mean(cell2mat(zscoredPower(binInds)'), 2);
            binnedPower(:, bin) = binnedPower(freqIndices, bin);
        end
    end

    % Plot the spectrogram
    figure;
    imagesc(binEdges(1:end-1) + 0.5, filteredFrequencies, binnedPower); % Offset bins for label centering
    axis xy; colormap(jet); colorbar;
    xlabel('ZT Time (hours)');
    ylabel('Frequency (Hz)');
    title(sprintf('Spectrogram for %s (0-40 Hz)', condName));
end