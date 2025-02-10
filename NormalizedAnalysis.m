%% get data
load('/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/EphysAnalysis/CanuteCombinedData.mat');
%% Spectrograms
% Create spectrograms for each condition
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionNames = {'300Lux', '1000LuxWk1', '1000LuxWk4'};
freqs = CanuteCombined.MetaData.fo;

% Find the indices of frequencies below 50 Hz
below50HzIndex = freqs < 50;
freqsBelow50Hz = freqs(below50HzIndex);

% Set up the figure
figure;

% Iterate over each condition to plot the spectrogram
for condIdx = 1:length(conditions)
    cond_data = CanuteCombined.(conditions{condIdx});
    zscoreFreqPower = cond_data.ZscoredFrequencyPower;
    zt_time = cond_data.ZT_time;
    
    % Initialize a matrix to accumulate frequency power for each ZT hour
    powerMatrix = zeros(sum(below50HzIndex), 24); % Using only frequencies below 50 Hz
    
    % Accumulator for the number of entries per ZT hour for averaging
    countMatrix = zeros(sum(below50HzIndex), 24);
    
    % Accumulate power data
    for idx = 1:length(zscoreFreqPower)
        hour = zt_time(idx) + 1; % Convert ZT time to 1-based indexing
        powerMatrix(:, hour) = powerMatrix(:, hour) + zscoreFreqPower{idx}(below50HzIndex);
        countMatrix(:, hour) = countMatrix(:, hour) + 1;
    end
    
    % Average the accumulated power by dividing by the count
    averagedPowerMatrix = powerMatrix ./ countMatrix;
    
    % Plot
    subplot(3, 1, condIdx); % subplot for each condition
    imagesc(0:23, freqsBelow50Hz, averagedPowerMatrix);
    set(gca, 'YDir', 'normal');
    colormap('parula');
    colorbar;
    xlabel('Zeitgeber Time (ZT)');
    ylabel('Frequency (Hz)');
    title(['Spectrogram for ' conditionNames{condIdx} ' (<50 Hz)']);
end

% Adjust figure
sgtitle('Z-Scored Frequency Power Spectrograms Below 50 Hz by Condition');

%% Line plots per band

% Define frequency bands and their corresponding names
bandNames = {'delta', 'theta', 'spindle', 'lowbeta', 'highbeta'};
startStopFreqs = {[0.5, 4], [5, 10], [11, 19], [20, 30], [30, 40]};

% Store color and line styles for different conditions
lineStyles = {'-', '--', ':'};
conditionNames = {'300Lux', '1000LuxWk1', '1000LuxWk4'};
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
freqs = CanuteCombined.MetaData.fo;

% Initiate figure
figure;

% Iterate over each frequency band
for bandIdx = 1:length(bandNames)
    
    % Find indices within this frequency band
    band = startStopFreqs{bandIdx};
    bandMask = (freqs >= band(1)) & (freqs <= band(2));
    
    % Plot setup
    subplot(length(bandNames), 1, bandIdx); % One subplot per band
    hold on;
    
    % Iterate over each condition
    for condIdx = 1:length(conditions)
        cond_data = CanuteCombined.(conditions{condIdx});
        zscoreFreqPower = cond_data.ZscoredFrequencyPower;
        
        % Initialize matrix to store mean power in each ZT hour
        meanPowerByZT = zeros(1, 24);
        countByZT = zeros(1, 24);
        
        % Accumulate power data for each ZT hour
        for idx = 1:length(zscoreFreqPower)
            hour = cond_data.ZT_time(idx) + 1; % Convert ZT time to 1-based indexing
            bandPower = mean(zscoreFreqPower{idx}(bandMask), 'omitnan');
            if ~isnan(bandPower)
                meanPowerByZT(hour) = meanPowerByZT(hour) + bandPower;
                countByZT(hour) = countByZT(hour) + 1;
            end
        end
        
        % Average the accumulated power
        averagedPower = meanPowerByZT ./ countByZT;
        
        % Plot the line
        plot(0:23, averagedPower, lineStyles{condIdx}, 'LineWidth', 1.5, 'DisplayName', conditionNames{condIdx});
    end
    
    % Customize plot
    xlabel('Zeitgeber Time (ZT)');
    ylabel('Normalized Power');
    title([bandNames{bandIdx} ' Band Power']);
    legend('show');
    grid on;
    hold off;
end

% Adjust figure
sgtitle('Normalized Power by Frequency Band and Condition');