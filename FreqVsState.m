% Load the .mat file
% load('/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/DataFiles/Harald/HaraldCombinedDataOutliersRemoved.mat');
saveDir = '/Users/noahmuscat/Desktop';

% Define sleep state values and labels
sleepStates = {'WAKE', 'NREM', 'REM'};
stateValues = [1, 3, 5];

% Define conditions and their colors
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionTitles = {'300 Lux', '1000 Lux Week 1', '1000 Lux Week 4'};
colors = {'b', 'g', 'r'}; % Blue, Green, Red for the different conditions

% Get frequency indices excluding 55-65 Hz and 115-125 Hz
fo = HaraldCombined.MetaData.fo;
%freqIndices = ~(fo >= 55 & fo <= 65) & ~(fo >= 115 & fo <= 125);
freqIndices = fo <= 50; 

% Find global min and max power across all conditions, sleep states, and time periods
globalMin = inf;
globalMax = -inf;

for s = 1:length(sleepStates)
    for timePeriod = 1:2
        for c = 1:length(conditions)
            condition = conditions{c};
            
            hourRange = [0, 12];
            if timePeriod == 2
                hourRange = [12, 24];
            end
            
            % Identify respective indices
            idx = HaraldCombined.(condition).ZT_Datetime.Hour >= hourRange(1) & ...
                  HaraldCombined.(condition).ZT_Datetime.Hour < hourRange(2) & ...
                  HaraldCombined.(condition).SleepState == stateValues(s);

            if any(idx)
                % Gather power data and compute mean
                powerData = cat(2, HaraldCombined.(condition).FrequencyPower{idx});
                meanPower = mean(powerData, 2, 'omitnan');
                
                % Update global min and max
                currentMin = min(meanPower(freqIndices));
                currentMax = max(meanPower(freqIndices));
                globalMin = min(globalMin, currentMin);
                globalMax = max(globalMax, currentMax);
            end
        end
    end
end

% Plot each sleep state
for s = 1:length(sleepStates)
    figure('Name', sprintf('%s - Frequency Power', sleepStates{s}), 'NumberTitle', 'off');
    
    for timePeriod = 1:2
        timeTitle = 'Day';
        hourRange = [0, 12];
        if timePeriod == 2
            timeTitle = 'Night';
            hourRange = [12, 24];
        end
        
        % Plot each condition
        for c = 1:length(conditions)
            condition = conditions{c};
            color = colors{c};
            
            idx = HaraldCombined.(condition).ZT_Datetime.Hour >= hourRange(1) & ...
                  HaraldCombined.(condition).ZT_Datetime.Hour < hourRange(2) & ...
                  HaraldCombined.(condition).SleepState == stateValues(s);

            if any(idx)
                powerData = cat(2, HaraldCombined.(condition).FrequencyPower{idx});
                meanPower = mean(powerData, 2, 'omitnan');
                
                subplot(1, 2, timePeriod);
                hold on;
                plot(fo(freqIndices), meanPower(freqIndices), 'Color', color, 'DisplayName', conditionTitles{c});
                xlabel('Frequency (Hz)');
                ylabel('Power');
                title(sprintf('%s - %s', sleepStates{s}, timeTitle));
                ylim([globalMin, globalMax]);  % Set consistent y-axis limits
            end
        end
        hold off;
        legend('show', 'Location', 'northeast');
    end
    
    sgtitle(sprintf('%s - Frequency Power', sleepStates{s}));
    fullPath = fullfile(saveDir, sprintf('%s_HaraldFrequencyPower.png', sleepStates{s}));
    %saveas(gcf, fullPath);
end