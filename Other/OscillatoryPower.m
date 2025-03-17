saveDir = '/Users/noahmuscat/Desktop';

% Define frequency bands and their corresponding names
bandNames = {'delta', 'theta', 'spindle', 'lowbeta', 'highbeta', 'lowgamma', 'midgamma', 'highgamma', 'ripple'};
startStopFreqs = {[0.5 4], [5 10], [11 19], [20 30], [30 40], [40 60], [60 100], [100 140], [140 180]};

% Get the frequency values from MetaData and apply range exclusions
fo = HaraldCombined.MetaData.fo;
selectedFreqIdx = (fo >= 0.5 & fo <= 180) & ~(fo >= 55 & fo <= 65) & ~(fo >= 115 & fo <= 125);

% Define the conditions and their labels
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionLabels = {'300 Lux', '1000 Lux Wk1', '1000 Lux Wk4'};

% Define sleep states and labels
sleepStates = [1, 3, 5];  % Wake, Non-REM, REM
stateLabels = {'Wake', 'NREM', 'REM'};

% Loop over each sleep state to create separate figures
for sleepStateIdx = 1:length(sleepStates)
    state = sleepStates(sleepStateIdx);
    
    % Create a new figure for the current sleep state
    figure;
    sgtitle(sprintf('Percentage Oscillatory Power During %s', stateLabels{sleepStateIdx}));
    
    % Loop over each condition to create subplots
    for condIdx = 1:length(conditions)
        condition = HaraldCombined.(conditions{condIdx});

        % Logical indexing for the current sleep state
        stateIndices = (condition.SleepState == state);

        % Filter FrequencyPower by the current sleep state
        stateFrequencyPower = condition.FrequencyPower(stateIndices);

        % Initialize storage for average power within each band
        avgBandPower = zeros(length(bandNames), 1);

        % Calculate mean power in each band for the selected sleep state
        for b = 1:length(bandNames)
            bandRange = startStopFreqs{b};
            bandIdx = find(fo >= bandRange(1) & fo <= bandRange(2));
            
            if isempty(bandIdx)
                avgBandPower(b) = 0;
            else
                % Exclude frequencies remove if selectedFreqIdx is false
                validBandIdx = bandIdx(selectedFreqIdx(bandIdx));
                bandPower = zeros(length(stateFrequencyPower), 1);
                
                for n = 1:length(stateFrequencyPower)
                    freqPower = stateFrequencyPower{n}; % 491x1 double
                    bandPower(n) = mean(freqPower(validBandIdx));  % Mean power for this band
                end

                % Average power across all time points for the current sleep state
                avgBandPower(b) = mean(bandPower);
            end
        end
        
        % Calculate percentage power for each band relative to total band power
        bandPercentages = (avgBandPower / sum(avgBandPower)) * 100;

        % Create a subplot
        subplot(1, 3, condIdx);
        bar(bandPercentages);
        set(gca, 'XTickLabel', bandNames);
        xlabel('Frequency Band');
        ylabel('Percentage Power (%)');
        title(conditionLabels{condIdx});
        ylim([0 50]);
        grid on;
    end

    %saveas(gcf, fullfile(saveDir, sprintf('HaraldOscillatoryPower_%s.png', stateLabels{sleepStateIdx})));
end