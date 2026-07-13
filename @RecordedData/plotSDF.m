function plotSDF(obj, neus2consider, conds2consider, window, binWidth, order, displayErrorFactor, varargin)
    % PLOTSDF - Plot spike density function with condition comparison
    %
    % SYNTAX
    %   obj.plotSDF(neus2consider, conds2consider, window, binWidth, order, displayErrorFactor)
    %   obj.plotSDF(..., 'markerNames', names)
    %   obj.plotSDF(..., 'smoothWindow', win)
    %
    % INPUTS:
    %   obj                 - RecordedData object
    %   neus2consider       - Neuron indices to include (numeric/logical vector)
    %   conds2consider      - Condition indices to include (numeric/logical vector)
    %   window              - [start(s), alignmentMarker, end(s)] (length-3 vector)
    %   binWidth            - [binSize, stepSize] in seconds (length-2 vector;
    %                         set stepSize = 0 for non-overlapping bins)
    %   order               - "preference" | "condition"
    %   displayErrorFactor  - Std scaling factor for shaded error band (numeric scalar)
    %
    % OPTIONAL NAME-VALUE PAIRS:
    %   'markerNames'  - Cell array of strings labelling each marker line.
    %                    Must have one entry per marker.
    %   'smoothWindow' - Integer window length in bins for Gaussian smoothing.
    %                    If omitted, no smoothing is applied.
    %
    % OUTPUT:
    %   Figure with spike density function for each condition.
    %
    % DESCRIPTION:
    %   Internally calls prepareTensorWithFixBin (and optionally refineNeuralData
    %   for smoothing) on a temporary copy of the object so the caller's data
    %   are never modified.  Conditions can be sorted by neural preference
    %   (peak firing rate) or kept in their original order.

    % ---- Parse optional arguments ----------------------------------------
    p = inputParser;
    validMarkerNames = @(x) isempty(x) || ...
        (iscell(x) && all(cellfun(@(s) ischar(s) || (isstring(s) && isscalar(s)), x(:))));
    addParameter(p, 'markerNames',  {},  validMarkerNames);
    addParameter(p, 'smoothWindow', [],  @(x) isempty(x) || ...
        (isnumeric(x) && isscalar(x) && x == round(x) && x > 0));
    parse(p, varargin{:});

    markerNames  = p.Results.markerNames;
    smoothWindow = p.Results.smoothWindow;

    if ~isempty(markerNames)
        markerNames = cellfun(@char, markerNames(:)', 'UniformOutput', false);
    end

    % ---- Build tensor (temporary local copy — caller object unchanged) ---
    tmpObj = obj.prepareTensorWithFixBin(window, binWidth, neus2consider, conds2consider);

    % ---- Optional Gaussian smoothing via refineNeuralData ---------------
    if ~isempty(smoothWindow)
        smoothOpt.smoothType   = 'gaussian';
        smoothOpt.smoothWindow = smoothWindow;
        tmpObj = tmpObj.refineNeuralData('smooth', smoothOpt, true);
    end

    % ---- Extract data ----------------------------------------------------
    nNeus  = size(tmpObj.TensorWithFixBin, 1);
    nConds = size(tmpObj.TensorWithFixBin, 2);
    nTime  = size(tmpObj.TensorWithFixBin, 4);

    avgTensor = reshape(mean(tmpObj.TensorWithFixBin, 3, 'omitnan'), nNeus, nConds, nTime);
    stdTensor = reshape(std( tmpObj.TensorWithFixBin, 0, 3, 'omitnan'), nNeus, nConds, nTime);

    if isfield(tmpObj.TensorWithFixBinInfo, 'conds2consider') && ...
            numel(tmpObj.TensorWithFixBinInfo.conds2consider) == nConds
        condLabels = tmpObj.TensorWithFixBinInfo.conds2consider;
    else
        condLabels = 1:nConds;
    end

    % ---- Condition ordering ----------------------------------------------
    maxx = max(avgTensor, [], 3, 'omitnan');

    if strcmp(order, "preference")
        reorderedAvg = nan(size(avgTensor));
        reorderedStd = nan(size(stdTensor));
        [~, sortIdx] = sort(maxx, 2, 'descend');
        for neu = 1:nNeus
            reorderedAvg(neu, :, :) = avgTensor(neu, sortIdx(neu, :), :);
            reorderedStd(neu, :, :) = stdTensor(neu, sortIdx(neu, :), :);
        end
        legendLabels = compose('Preference %d', 1:nConds);
    elseif strcmp(order, "condition")
        reorderedAvg = avgTensor;
        reorderedStd = stdTensor;
        legendLabels = compose('Condition %d', condLabels);
    else
        error('plotSDF: unknown order ''%s''. Choose ''preference'' or ''condition''.', order);
    end

    % ---- Build time axis (seconds) ---------------------------------------
    if numel(binWidth) > 1 && binWidth(2) > 0
        timeStep = binWidth(2);
    else
        timeStep = binWidth(1);
    end
    timeAxis = window(1) + (0:nTime-1) * timeStep + binWidth(1)/2;

    % ---- Plot ------------------------------------------------------------
    figure;
    hold on;
    cmap = colormap(jet(nConds));

    for cond = 1:nConds
        if nNeus == 1
            avgTrace = squeeze(reorderedAvg(1, cond, :));
            errTrace = squeeze(reorderedStd(1, cond, :));
        else
            avgTrace = mean(squeeze(reorderedAvg(:, cond, :)), 1, 'omitnan')';
            errTrace = mean(squeeze(reorderedStd(:, cond, :)), 1, 'omitnan')';
        end

        shaded_areas( ...
            timeAxis(:), ...
            avgTrace(:), ...
            errTrace(:) * displayErrorFactor, ...
            'FaceColor',   cmap(cond, :), ...
            'DisplayName', char(legendLabels(cond)));
    end

    % ---- Marker lines (bin positions → time) ----------------------------
    markerBins  = squeeze(mean(tmpObj.MarkerTensorForFixBin, [1 2 3], 'omitmissing'));
    markerBins  = markerBins(:)';
    markerTimes = interp1(1:nTime, timeAxis, markerBins, 'linear', NaN);
    markerIsPlotted = isfinite(markerTimes);

    if isempty(markerNames)
        if any(markerIsPlotted)
            xline(markerTimes(markerIsPlotted));
        end
    else
        if numel(markerNames) ~= numel(markerTimes)
            error('RecordedData:InvalidMarkerNames', ...
                'markerNames must contain one name per marker (%d expected).', ...
                numel(markerTimes));
        end
        if any(markerIsPlotted)
            xline(markerTimes(markerIsPlotted), '-', markerNames(markerIsPlotted), ...
                'LabelOrientation',       'aligned', ...
                'LabelHorizontalAlignment', 'left');
        end
    end

    hold off;
    legend show;
    xlabel('Time (s)');
    ylabel('Firing Rate (sp/s)');
    neuStr  = sprintf('neus %d-%d (n=%d)', min(neus2consider), max(neus2consider), numel(neus2consider));
    condStr = sprintf('conds [%s]', num2str(condLabels));
    winStr  = sprintf('win [%.3g %.3g] s', window(1), window(3));
    binStr  = sprintf('bin %.0f/%.0f ms', binWidth(1)*1000, binWidth(2)*1000);
    if ~isempty(smoothWindow)
        infoStr = [neuStr ' | ' condStr ' | ' winStr ' | ' binStr ' | smooth ' num2str(smoothWindow) ' bins'];
    else
        infoStr = [neuStr ' | ' condStr ' | ' winStr ' | ' binStr];
    end
    if strcmp(order, "preference")
        title({'SDF — ordered by preference', infoStr});
    else
        title({'SDF — condition order', infoStr});
    end
end
