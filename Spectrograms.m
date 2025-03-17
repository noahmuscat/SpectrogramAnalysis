saveDir = '/Users/noahmuscat/Desktop';
dataStruct = HaraldV2Combined80;

% Get frequency values from MetaData within the specified range
fo = dataStruct.MetaData.fo;
validFreqIdx = (fo >= 0 & fo < 40);
frequencies = fo(validFreqIdx);

% Specify conditions and labels for plotting
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionLabels = {'300 Lux', '1000 Lux Week 1', '1000 Lux Week 4'};

% Define sleep states and their labels
sleepStates = [1, 3, 5];
sleepStateLabels = {'Wake', 'NREM', 'REM'};

% Define aggregation time interval (5 minutes)
aggInterval = minutes(5);

% Precompute global min and max
globalMin = Inf;
globalMax = -Inf;

% First pass to determine global min and max across all data
for stateIdx = 1:length(sleepStates)
    currentState = sleepStates(stateIdx);

    for conditionIdx = 1:length(conditions)
        condition = dataStruct.(conditions{conditionIdx});
        ztDatetime = condition.ZT_Datetime;
        sleepState = condition.SleepState;
        zscoredFrequencyPower = condition.ZScoreFractionalPower;

        % Logical indexing for current sleep state
        stateIndices = sleepState == currentState;

        % Extract relevant z-scored power and corresponding times
        currentTimes = ztDatetime(stateIndices);
        zPowerState = zscoredFrequencyPower(stateIndices);

        % Prepare array to hold accumulated power spectra
        numPoints = numel(zPowerState);
        powerAccum = zeros(numPoints, length(frequencies));

        for n = 1:numPoints
            powerAccum(n, :) = zPowerState{n}(validFreqIdx)';
        end

        % Sort data
        [sortedTimes, sortIdx] = sort(currentTimes);
        sortedPower = powerAccum(sortIdx, :);

        % Bin data
        binEdges = min(sortedTimes):aggInterval:max(sortedTimes);
        numBins = length(binEdges) - 1;
        binnedPower = zeros(numBins, length(frequencies));

        for b = 1:numBins
            inBin = sortedTimes >= binEdges(b) & sortedTimes < binEdges(b+1);
            if any(inBin)
                binnedPower(b, :) = mean(sortedPower(inBin, :), 1, 'omitnan');
            else
                binnedPower(b, :) = NaN;
            end
        end

        % Update global min and max
        globalMin = min(globalMin, min(binnedPower, [], 'all', 'omitnan'));
        globalMax = max(globalMax, max(binnedPower, [], 'all', 'omitnan'));
    end
end

% Create plots using the determined global min/max for color scaling
for stateIdx = 1:length(sleepStates)
    currentState = sleepStates(stateIdx);
    
    figure;
    for conditionIdx = 1:length(conditions)
        condition = dataStruct.(conditions{conditionIdx});
        ztDatetime = condition.ZT_Datetime;
        sleepState = condition.SleepState;
        zscoredFrequencyPower = condition.ZScoreFractionalPower;
        
        % Logical indexing for current sleep state
        stateIndices = sleepState == currentState;

        % Extract relevant z-scored power and corresponding times
        currentTimes = ztDatetime(stateIndices);
        zPowerState = zscoredFrequencyPower(stateIndices);

        % Prepare array to hold accumulated power spectra
        numPoints = numel(zPowerState);
        powerAccum = zeros(numPoints, length(frequencies));

        for n = 1:numPoints
            powerAccum(n, :) = zPowerState{n}(validFreqIdx)';
        end

        % Sort data
        [sortedTimes, sortIdx] = sort(currentTimes);
        sortedPower = powerAccum(sortIdx, :);

        % Bin data
        binEdges = min(sortedTimes):aggInterval:max(sortedTimes);
        numBins = length(binEdges) - 1;
        binnedPower = zeros(numBins, length(frequencies));

        for b = 1:numBins
            inBin = sortedTimes >= binEdges(b) & sortedTimes < binEdges(b+1);
            if any(inBin)
                binnedPower(b, :) = mean(sortedPower(inBin, :), 1, 'omitnan');
            else
                binnedPower(b, :) = NaN;
            end
        end

        % Plot for the current condition
        subplot(1, length(conditions), conditionIdx);
        imagesc(datetime(binEdges(1:end-1)), frequencies, binnedPower');
        colorbar;
        colormap(jet);
        clim([globalMin, globalMax]);
        axis xy;
        xlabel('Time');
        ylabel('Frequency (Hz)');
        title(conditionLabels{conditionIdx});
    end
    
    sgtitle(['Z-scored Frequency Power - ', sleepStateLabels{stateIdx}]);

    % Save the plot for each sleep state
    saveas(gcf, fullfile(saveDir, sprintf('CanuteV2Ch80Spectrogram_%s.png', sleepStateLabels{stateIdx})));
end