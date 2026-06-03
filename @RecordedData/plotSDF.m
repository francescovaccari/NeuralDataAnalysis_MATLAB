function plotSDF(obj, neus2consider, conds2consider, order, displayErrorFactor)
    % PLOTSDF - Plot spike density function with condition comparison
    %
    % INPUTS:
    %   obj                     - RecordedData object
    %   neus2consider           - Neuron indices (numeric vector)
    %   conds2consider          - Condition indices (numeric vector)
    %   order                   - "preference" or "condition" (string)
    %   displayErrorFactor      - Error bar scaling factor (numeric scalar)
    %
    % OUTPUT:
    %   Figure with spike density function for each condition
    %
    % DESCRIPTION:
    %   Visualizes spike density function (mean firing rate) across multiple neurons
    %   and conditions. Can sort conditions by neural preference or keep original order.

    preprocessedData = [];

    preprocessedData.avg = obj.TensorWithFixBinCondAvg(neus2consider,conds2consider,:);
    if numel(neus2consider) == 1
        preprocessedData.std = std(obj.TensorWithFixBin(neus2consider,conds2consider,:,:), 0, 3, 'omitnan');
    else
        preprocessedData.std = squeeze(std(obj.TensorWithFixBin(neus2consider,conds2consider,:,:), 0, 3, 'omitnan'));
    end

    maxx = max(preprocessedData.avg, [], [3], 'omitnan');

    if strcmp(order, "preference")
        reorderedAvg = nan(size(preprocessedData.avg));
        reorderedStd = nan(size(preprocessedData.std));
        [~, sortIdx] = sort(maxx, 2, 'descend');

        for neu = 1:numel(neus2consider)
            reorderedAvg(neu, :, :) = preprocessedData.avg(neu, sortIdx(neu, :), :);
            reorderedStd(neu, :, :) = preprocessedData.std(neu, sortIdx(neu, :), :);
        end

        legendLabels = compose('Preference %d', 1:numel(conds2consider));
    elseif strcmp(order, "condition")
        reorderedAvg = preprocessedData.avg;
        reorderedStd = preprocessedData.std;
        legendLabels = compose('Condition %d', conds2consider);
    end

    figure;
    hold on;
    i = 0;
    cmap = colormap(jet(numel(conds2consider)));
    for cond = 1:numel(conds2consider)
        i = i+1;
        if numel(neus2consider) == 1
            avgTrace = squeeze(reorderedAvg(:, cond, :));
            errTrace = squeeze(reorderedStd(:, cond, :));
        else
            avgMatrix = squeeze(reorderedAvg(:, cond, :));
            stdMatrix = squeeze(reorderedStd(:, cond, :));
            avgTrace = mean(avgMatrix, 1, 'omitnan')';
            errTrace = mean(stdMatrix, 1, 'omitnan')';
        end

        shaded_areas( ...
            [1:numel(avgTrace)]', ...
            avgTrace(:), ...
            errTrace(:)*displayErrorFactor, ...
            'FaceColor', cmap(i, :), ...
            'DisplayName', char(legendLabels(i)));
    end
    xline(squeeze(mean(obj.MarkerTensorForFixBin(neus2consider,conds2consider,:,:),[1 2 3],'omitmissing')))
    hold off;
    legend show;
    xlabel('Time (s)');
    xticks(linspace(1,numel(avgTrace),11))
    window = obj.TensorWithFixBinInfo.window;
    xticklabels(linspace(window(1),window(3),11))
    ylabel('A.U.');
    if strcmp(order, "preference")
        title('Spike Density Function ordered by preference');
    else
        title(['Spike Density Function for conds ' num2str(conds2consider)]);
    end
end
