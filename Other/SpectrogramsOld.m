saveDir = '/Users/noahmuscat/Desktop';

dataStruct = HaraldV2Combined80;

% Get frequency values from MetaData within the specified range
fo = dataStruct.MetaData.fo;
validFreqIdx = (fo >= 0 & fo < 40);
%validFreqIdx = fo >= 40 & fo < 200 & ...
%               ~(fo >= 55 & fo <= 65) & ~(fo >= 115 & fo <= 125);
frequencies = fo(validFreqIdx);

% Specify conditions and labels for plotting
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionLabels = {'300 Lux', '1000 Lux Week 1', '1000 Lux Week 4'};

% Define sleep states and their labels
sleepStates = [1, 3, 5];
sleepStateLabels = {'Wake', 'NREM', 'REM'};

% Iterate over each sleep state to compute global min and max for plotting
globalMin = zeros(1, length(sleepStates));
globalMax = zeros(1, length(sleepStates));

% Initialize global min and max with extremes
globalMin(:) = Inf;
globalMax(:) = -Inf;

for stateIdx = 1:length(sleepStates)
    avgZscoredPowerOverall = []; % Store all hourly averages for current sleep state
    
    currentState = sleepStates(stateIdx);
    
    for conditionIdx = 1:length(conditions)
        condition = dataStruct.(conditions{conditionIdx});
        ztDatetime = condition.ZT_Datetime;
        sleepState = condition.SleepState;
        zscoredFrequencyPower = condition.ZScoreFractionalPower;

        % Convert the ZT datetime to ZT hours
        ztHours = hour(ztDatetime);
        uniqueHours = 0:23;

        % Initialize storage for average z-scored power per ZT hour for the current sleep state
        avgZscoredPowerPerHour = zeros(length(uniqueHours), length(frequencies));

        % Logical indexing for the current sleep state
        stateIndices = (sleepState == currentState);

        % Iterate over each ZT hour
        for h = uniqueHours
            % Logical indexing for the current ZT hour within the current sleep state
            hourIndices = stateIndices & (ztHours == h);

            % Extract z-scored frequency power data for the current hour and sleep state
            hourStateZscoredPower = zscoredFrequencyPower(hourIndices);

            % Accumulate power data for averaging
            if ~isempty(hourStateZscoredPower)
                hourPowerAccum = zeros(length(hourStateZscoredPower), length(frequencies));
                for n = 1:length(hourStateZscoredPower)
                    freqPower = hourStateZscoredPower{n};
                    % Apply the valid frequency index
                    hourPowerAccum(n, :) = freqPower(validFreqIdx)';
                end
                avgZscoredPowerPerHour(h + 1, :) = mean(hourPowerAccum, 1, 'omitnan');
            end
        end

        % Store the computed averages for global min/max calculation
        avgZscoredPowerOverall = [avgZscoredPowerOverall; avgZscoredPowerPerHour];
    end
    
    % Calculate global min and max from the hourly binned averages
    globalMin(stateIdx) = min(avgZscoredPowerOverall(:), [], 'omitnan');
    globalMax(stateIdx) = max(avgZscoredPowerOverall(:), [], 'omitnan');
end

% Iterate over each sleep state to plot
for stateIdx = 1:length(sleepStates)
    currentState = sleepStates(stateIdx);

    % Create a new figure for the current sleep state
    figure;
    
    % Iterate over each condition to fill subplots
    for conditionIdx = 1:length(conditions)
        condition = dataStruct.(conditions{conditionIdx});
        ztDatetime = condition.ZT_Datetime;
        sleepState = condition.SleepState;
        zscoredFrequencyPower = condition.ZScoreFractionalPower;

        % Convert the ZT datetime to ZT hours
        ztHours = hour(ztDatetime);
        uniqueHours = 0:23;

        % Initialize storage for average z-scored power per ZT hour for the current sleep state
        avgZscoredPowerPerHour = zeros(length(uniqueHours), length(frequencies));

        % Logical indexing for the current sleep state
        stateIndices = (sleepState == currentState);

        % Iterate over each ZT hour
        for h = uniqueHours
            % Logical indexing for the current ZT hour within the current sleep state
            hourIndices = stateIndices & (ztHours == h);

            % Extract z-scored frequency power data for the current hour and sleep state
            hourStateZscoredPower = zscoredFrequencyPower(hourIndices);

            % Accumulate power data for averaging
            if ~isempty(hourStateZscoredPower)
                hourPowerAccum = zeros(length(hourStateZscoredPower), length(frequencies));
                for n = 1:length(hourStateZscoredPower)
                    freqPower = hourStateZscoredPower{n};
                    % Apply the valid frequency index
                    hourPowerAccum(n, :) = freqPower(validFreqIdx)';
                end
                avgZscoredPowerPerHour(h + 1, :) = mean(hourPowerAccum, 1, 'omitnan');
            end
        end

        % Create a subplot for each condition
        subplot(1, length(conditions), conditionIdx);
        imagesc(uniqueHours, frequencies, avgZscoredPowerPerHour');
        colorbar;
        clim([globalMin(stateIdx), globalMax(stateIdx)]); % Set the color axis for global min/max
        axis xy;
        xlabel('ZT Hour');
        ylabel('Frequency (Hz)');
        title(conditionLabels{conditionIdx});
        set(gca, 'XTick', 0:1:23);
    end
    
    sgtitle(['Z-scored Frequency Power - ', sleepStateLabels{stateIdx}]);

    %saveas(gcf, fullfile(saveDir, sprintf('CanuteV2Ch80SpectrogramLowFreqs_%s.png', sleepStateLabels{stateIdx})));
end