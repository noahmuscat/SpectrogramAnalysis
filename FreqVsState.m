% Load the .mat file
load('HaraldCombinedDataOutliersRemoved.mat');

%% 
saveDir = '/Users/noahmuscat/Desktop';

% Define sleep state values and labels
sleepStates = {'WAKE', 'NREM', 'REM'};
stateValues = [1, 3, 5];

% Define conditions
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionTitles = {'300 Lux', '1000 Lux Week 1', '1000 Lux Week 4'};

% Get frequency indices from 0 to 40 Hz
freqIndices = HaraldCombined.MetaData.fo <= 40;
upperLimit = 50; % Remove any powers above this limit

% Loop through each condition
for c = 1:length(conditions)
    condition = conditions{c};
    figure('Name', conditionTitles{c}, 'NumberTitle', 'off'); % Create a new figure for each condition
    
    % Loop through each sleep state
    for s = 1:length(sleepStates)
        state = stateValues(s);
        
        % Identify day and night indices separately
        dayIdx = HaraldCombined.(condition).ZT_time >= 0 & HaraldCombined.(condition).ZT_time < 12 & HaraldCombined.(condition).SleepState == state;
        nightIdx = HaraldCombined.(condition).ZT_time >= 12 & HaraldCombined.(condition).ZT_time < 24 & HaraldCombined.(condition).SleepState == state;

        % Ensure there are data points
        if any(dayIdx)
            % Gather power data for day, eliminating powers above the set limit
            dayPowerData = cat(2, HaraldCombined.(condition).FrequencyPower{dayIdx});
            dayPowerData(dayPowerData > upperLimit) = NaN;

            % Compute mean power for the day
            dayPower = mean(dayPowerData, 2, 'omitnan');

            % Plot for day
            subplot(3, 2, 2*s-1);
            plot(HaraldCombined.MetaData.fo(freqIndices), dayPower(freqIndices));
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title(sprintf('%s - Day', sleepStates{s}));
        end

        if any(nightIdx)
            % Gather power data for night, eliminating powers above the set limit
            nightPowerData = cat(2, HaraldCombined.(condition).FrequencyPower{nightIdx});
            nightPowerData(nightPowerData > upperLimit) = NaN;

            % Compute mean power for the night
            nightPower = mean(nightPowerData, 2, 'omitnan');

            % Plot for night
            subplot(3, 2, 2*s);
            plot(HaraldCombined.MetaData.fo(freqIndices), nightPower(freqIndices));
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title(sprintf('%s - Night', sleepStates{s}));
        end
    end
    
    % Save the figure
    sgtitle(conditionTitles{c});
    fullPath = fullfile(saveDir, [conditionTitles{c}, '.png']);
    saveas(gcf, fullPath);
end