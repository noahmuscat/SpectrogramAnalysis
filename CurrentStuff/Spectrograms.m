% Define directory and data structure
saveDir = '/Users/noahmuscat/Desktop';
dataStruct = HaraldV3Combined80;

% Frequency values
fo = dataStruct.MetaData.fo;
validFreqIdx = (fo >= 0 & fo < 40); 
%validFreqIdx = fo >= 40 & fo < 200 & ~(fo >= 55 & fo <= 65) & ~(fo >= 115 & fo <= 125);
frequencies = fo(validFreqIdx);

% Conditions and labels
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionLabels = {'300 Lux', '1000 Lux Week 1', '1000 Lux Week 4'};

% Sleep states and labels
sleepStates = [1, 3, 5];
sleepStateLabels = {'Wake', 'NREM', 'REM'};

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

        % Bin data according to the ZT hour structure
        ztHours = hour(currentTimes);
        binEdges = 0:24; % Bin edges for ZT hours 0 to 23
        numBins = length(binEdges) - 1;
        binnedPower = zeros(numBins, length(frequencies));

        for b = 1:numBins
            inBin = ztHours == (b - 1);
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
        xlim([0, 23]); % Restrict x-axis to complete ZT hours
    end
    
    sgtitle(['Z-scored Frequency Power - ', sleepStateLabels{stateIdx}]);

    % Save the plot for each sleep state
    saveas(gcf, fullfile(saveDir, sprintf('HaraldV3Ch80SpectrogramLowFreqs_%s.png', sleepStateLabels{stateIdx})));
end