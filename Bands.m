% Define save directory
saveDir = '/Users/noahmuscat/Desktop';
dataStruct = HaraldV2Combined80;

% Defining frequency bands
bands = struct();
bands.delta.startstopfreq = [0.5 4];
bands.theta.startstopfreq = [5 10];
bands.spindle.startstopfreq = [11 19];
bands.lowbeta.startstopfreq = [20 30];
bands.highbeta.startstopfreq = [30 40];
bands.lowgamma.startstopfreq = [40 60];
bands.midgamma.startstopfreq = [60 100];
bands.highgamma.startstopfreq = [100 140];
bands.ripple.startstopfreq = [140 180];

bandNames = fieldnames(bands);

% Specify conditions and labels for plotting
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionLabels = {'300 Lux', '1000 Lux Week 1', '1000 Lux Week 4'};

% Define sleep states and their labels
sleepStates = [1, 3, 5];
sleepStateLabels = {'Wake', 'NREM', 'REM'};

% Get frequency values and set valid indices
fo = dataStruct.MetaData.fo;

% Create indices for each band
bandIndices = cell(size(bandNames));
for b = 1:length(bandNames)
    bandFreq = bands.(bandNames{b}).startstopfreq;
    bandIndices{b} = fo >= bandFreq(1) & fo <= bandFreq(2);
end

% Calculate mean Z-score power in each band for each condition and state
for stateIdx = 1:length(sleepStates)
    currentState = sleepStates(stateIdx);
    
    % Prepare data for bar plotting
    meanPowers = zeros(length(bandNames), length(conditions));
    
    for conditionIdx = 1:length(conditions)
        condition = dataStruct.(conditions{conditionIdx});
        sleepState = condition.SleepState == currentState;
        zscoredFrequencyPower = condition.ZScoreFractionalPower;
        
        % Filter zscoredFrequencyPower for the current sleep state
        zPowerState = zscoredFrequencyPower(sleepState);
        
        % Precompute band power sums
        bandPowerSums = zeros(length(bandNames), numel(zPowerState));

        for n = 1:numel(zPowerState)
            currentPower = zPowerState{n};
            for b = 1:length(bandNames)
                bandPowerSums(b, n) = mean(currentPower(bandIndices{b}), 'omitnan');
            end
        end
        meanPowers(:, conditionIdx) = mean(bandPowerSums, 2, 'omitnan');
    end
    
    % Plotting the bar graph for the current sleep state
    figure;
    bar(meanPowers);
    set(gca, 'XTickLabel', bandNames);
    ylabel('Mean Z-Score Fractional Power');
    title(['Mean Z-Scored Power - ', sleepStateLabels{stateIdx}]);
    legend(conditionLabels, 'Location', 'best');
    
    % Save each figure
    saveas(gcf, fullfile(saveDir, sprintf('HaraldBars_%s.png', sleepStateLabels{stateIdx})));
end

%% Combining Animals
% Define save directory
saveDir = '/Users/noahmuscat/Desktop';

% Load data for both animals
haraldData = HaraldV2Combined80;
canuteData = CanuteV2Combined80;

% Defining frequency bands
bands = struct();
bands.delta.startstopfreq = [0.5 4];
bands.theta.startstopfreq = [5 10];
bands.spindle.startstopfreq = [11 19];
bands.lowbeta.startstopfreq = [20 30];
bands.highbeta.startstopfreq = [30 40];
bands.lowgamma.startstopfreq = [40 60];
bands.midgamma.startstopfreq = [60 100];
bands.highgamma.startstopfreq = [100 140];
bands.ripple.startstopfreq = [140 180];

bandNames = fieldnames(bands);

% Specify conditions and labels for plotting
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionLabels = {'300 Lux', '1000 Lux Week 1', '1000 Lux Week 4'};

% Define sleep states and their labels
sleepStates = [1, 3, 5];
sleepStateLabels = {'Wake', 'NREM', 'REM'};

% Get frequency values and set valid indices for the bands
fo = haraldData.MetaData.fo; % Assuming both animals have identical frequency metadata

% Create indices for each band
bandIndices = cell(size(bandNames));
for b = 1:length(bandNames)
    bandFreq = bands.(bandNames{b}).startstopfreq;
    bandIndices{b} = fo >= bandFreq(1) & fo <= bandFreq(2);
end

% Define the function for calculating mean powers.
% Note: Pass sleepStates array to avoid scope issue.
function meanPowers = calculateMeanPowers(dataStruct, sleepStateIdx, conditions, bandNames, bandIndices, sleepStates)
    sleepState = sleepStates(sleepStateIdx);
    meanPowers = zeros(length(bandNames), length(conditions));
    
    for conditionIdx = 1:length(conditions)
        condition = dataStruct.(conditions{conditionIdx});
        sleepStateMask = condition.SleepState == sleepState;
        zscoredFrequencyPower = condition.ZScoreFractionalPower;
        
        % Extract z-scored powers for current sleep state
        zPowerState = zscoredFrequencyPower(sleepStateMask);
        
        % Precompute band power sums
        bandPowerSums = zeros(length(bandNames), numel(zPowerState));

        for n = 1:numel(zPowerState)
            currentPower = zPowerState{n};
            for b = 1:length(bandNames)
                bandPowerSums(b, n) = mean(currentPower(bandIndices{b}), 'omitnan');
            end
        end
        meanPowers(:, conditionIdx) = mean(bandPowerSums, 2, 'omitnan');
    end
end

% Calculate mean Z-score power in each band for each condition and state
for stateIdx = 1:length(sleepStates)
    currentStateLabel = sleepStateLabels{stateIdx};
    
    % Calculate power for both animals
    haraldMeanPowers = calculateMeanPowers(haraldData, stateIdx, conditions, bandNames, bandIndices, sleepStates);
    canuteMeanPowers = calculateMeanPowers(canuteData, stateIdx, conditions, bandNames, bandIndices, sleepStates);

    % Plotting
    figure;
    for animalIdx = 1:2
        subplot(1, 2, animalIdx);
        if animalIdx == 1
            bar(haraldMeanPowers);
            title(['Harald - ', currentStateLabel]);
        else
            bar(canuteMeanPowers);
            title(['Canute - ', currentStateLabel]);
        end
        set(gca, 'XTickLabel', bandNames);
        ylabel('Mean Z-Score Fractional Power');
        legend(conditionLabels, 'Location', 'best');
    end
    
    % Annotate figure title and save
    sgtitle(['Mean Z-Scored Power - ', currentStateLabel]);
    saveas(gcf, fullfile(saveDir, sprintf('MeanZPower_%s.png', currentStateLabel)));
end