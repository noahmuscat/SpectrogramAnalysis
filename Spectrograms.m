saveDir = '/Users/noahmuscat/Desktop';

% Get the frequency values from MetaData within the specified range
fo = HaraldCombined.MetaData.fo;
%validFreqIdx = (fo >= 0 & fo < 40);
validFreqIdx = fo >= 40 & fo < 200 & ...
               ~(fo >= 55 & fo <= 65) & ~(fo >= 115 & fo <= 125);
frequencies = fo(validFreqIdx);

% Specify the conditions and their labels for plotting
conditions = {'Cond_300Lux', 'Cond_1000LuxWk1', 'Cond_1000LuxWk4'};
conditionLabels = {'300 Lux', '1000 Lux Week 1', '1000 Lux Week 4'};

% Gather spectrograms data for all conditions
spectrogramData = cell(length(conditions), 1);

% Iterate over each condition
for conditionIdx = 1:length(conditions)
    condition = HaraldCombined.(conditions{conditionIdx});
    ztDatetime = condition.ZT_Datetime;
    zscoredFrequencyPower = condition.ZscoredFrequencyPower;

    % Convert the ztDatetime to ZT hours
    ztHours = hour(ztDatetime);
    uniqueHours = 0:23;

    % Initialize storage for average z-scored power per ZT hour
    avgZscoredPowerPerHour = zeros(length(uniqueHours), length(frequencies));

    % Iterate over each ZT hour
    for h = uniqueHours
        % Logical indexing for the current ZT hour
        hourIndices = (ztHours == h);

        % Extract z-scored frequency power data for the current hour
        hourZscoredPower = zscoredFrequencyPower(hourIndices);

        % Calculate average z-scored power across all records in the hour
        if ~isempty(hourZscoredPower)
            hourPowerAccum = zeros(length(hourZscoredPower), length(frequencies));
            for n = 1:length(hourZscoredPower)
                freqPower = hourZscoredPower{n}; % 491x1 double
                hourPowerAccum(n, :) = freqPower(validFreqIdx)';
            end
            avgZscoredPowerPerHour(h + 1, :) = mean(hourPowerAccum, 1, 'omitnan'); % Handle NaNs if present
        end
    end

    % Store spectrogram data
    spectrogramData{conditionIdx} = avgZscoredPowerPerHour;
end

% Plot each condition using a consistent color scale
for conditionIdx = 1:length(conditions)
    figure;
    imagesc(uniqueHours, frequencies, spectrogramData{conditionIdx}');
    colorbar
    axis xy;
    xlabel('ZT Hour');
    ylabel('Frequency (Hz)');
    title(['Z-scored Frequency Power per ZT Hour - ', conditionLabels{conditionIdx}]);
    set(gca, 'XTick', 0:1:23, 'XTickLabel', 0:1:23);
    saveas(gcf, fullfile(saveDir, sprintf('HaraldSpectrogramHighFreqs_%s.png', conditionLabels{conditionIdx})));

end