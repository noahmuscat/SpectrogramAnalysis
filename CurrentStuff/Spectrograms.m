saveDir = '/Users/noahmuscat/Desktop';
dataStruct = HaraldV2Combined80;

% Get frequency values from MetaData within the specified range
fo = dataStruct.MetaData.fo;
%validFreqIdx = (fo >= 0 & fo < 40); 
validFreqIdx = fo >= 40 & fo < 200 & ...
              ~(fo >= 55 & fo <= 65) & ~(fo >= 115 & fo <= 125);
frequencies = fo(validFreqIdx);

% Specify conditions and labels for plotting
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionLabels = {'300 Lux', '1000 Lux Week 1', '1000 Lux Week 4'};

% Define sleep states and their labels
sleepStates = [1, 3, 5];
sleepStateLabels = {'Wake', 'NREM', 'REM'};

% Define aggregation time interval (5 minutes)
aggInterval = minutes(5);

% Create spectrograms 
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

        % Sort and accumulate data into the ZT hour structure
        ztHours = hour(currentTimes) + minute(currentTimes)/60;
        binEdges = 0:(5/60):24;
        numBins = length(binEdges) - 1;
        binnedPower = zeros(numBins, length(frequencies));

        for b = 1:numBins
            inBin = ztHours >= binEdges(b) & ztHours < binEdges(b+1);
            if any(inBin)
                binnedPower(b, :) = mean(powerAccum(inBin, :), 1, 'omitnan');
            else
                binnedPower(b, :) = NaN;
            end
        end

        % Plot for the current condition
        subplot(1, length(conditions), conditionIdx);
        imagesc(binEdges(1:end-1), frequencies, binnedPower');
        colorbar;
        colormap(parula);
        clim([-1.5 1.5]);
        axis xy;
        xlabel('ZT Time (hours)');
        ylabel('Z-Scored Frequency (Hz)');
        title(conditionLabels{conditionIdx});
        xlim([0, 24]); % Restrict x-axis to ZT 0-23
    end
    
    sgtitle(['Z-scored Frequency Power - ', sleepStateLabels{stateIdx}]);

    % Save the plot for each sleep state
    saveas(gcf, fullfile(saveDir, sprintf('HaraldV2Ch80SpectrogramHighFreqs_%s.png', sleepStateLabels{stateIdx})));
end