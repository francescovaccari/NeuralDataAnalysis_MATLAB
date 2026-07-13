function anovaResults = computeAnovaOnCond(obj, neus2consider, conds2consider, refWindow, anaWindow, binWidth, varargin)
    % COMPUTEANOVAONCOND - Compute ANOVA comparing baseline and analysis windows
    %
    % INPUTS:
    %   obj                     - RecordedData object
    %   neus2consider           - Neuron indices (numeric vector)
    %   conds2consider          - Condition indices (numeric vector)
    %   refWindow               - Reference [start(s), marker, end(s)] (numeric vector, length 3)
    %   anaWindow               - Analysis [start(s), marker, end(s)] (numeric vector, length 3)
    %   binWidth                - [bin_size, step] in seconds (numeric vector, length 2)
    %   'markerNames'           - Cell array of marker names
    %   'display'               - Logical flag to show the summary figures (default: true)
    %   'factorNames'           - Names for epoch and condition factors. If 'epoch'
    %                             is omitted, it is added as the first factor.
    %
    % OUTPUT:
    %   anovaResults            - Struct containing:
    %                             - p2way: (neurons x bins x 3) p-values [epoch, condition, epoch x condition]
    %                             - p1way: (neurons x bins) p-values [condition]
    %                             - sign2way: (bins x 3) % neurons with p<0.05 (2-way)
    %                             - sign1way: (bins) % neurons with p<0.05 (1-way)
    %
    % DESCRIPTION:
    %   Two-way ANOVA (epoch x condition) comparing baseline vs. analysis window.
    %   One-way ANOVA (condition only) within analysis window.

    p = inputParser;
    validMarkerNames = @(names) isempty(names) || ...
        (iscell(names) && all(cellfun(@(name) ischar(name) || ...
        (isstring(name) && isscalar(name)), names(:))));
    validDisplay = @(value) (islogical(value) && isscalar(value)) || ...
        (isnumeric(value) && isscalar(value) && ismember(value, [0 1]));
    validFactorNames = @(names) isempty(names) || ischar(names) || isstring(names) || ...
        (iscell(names) && all(cellfun(@(name) ischar(name) || ...
        (isstring(name) && isscalar(name)), names(:))));
    addParameter(p, 'markerNames', {}, validMarkerNames);
    addParameter(p, 'display', true, validDisplay);
    addParameter(p, 'factorNames', {}, validFactorNames);
    parse(p, varargin{:});

    markerNames = p.Results.markerNames;
    if isempty(markerNames)
        markerNames = {};
    else
        markerNames = cellfun(@char, markerNames(:)', 'UniformOutput', false);
    end
    displayFigures = logical(p.Results.display);
    varNames = normalizeFactorNames(p.Results.factorNames);
    effectNames2way = makeEffectNames([1 0; 0 1; 1 1], varNames);
    effectNames1way = varNames(2);

    if (refWindow(3)-refWindow(1)) ~= binWidth(1)
        disp("NOTE: The binWidth chosen for the analysis window doesn't match the size of the reference window")
    end

    % Call prepareTensorWithFixBin to prepare reference data
    refObj = obj.prepareTensorWithFixBin(refWindow, [refWindow(3)-refWindow(1) 0], 0, neus2consider, conds2consider);
    refData = refObj.TensorWithFixBin;

    % Call prepareTensorWithFixBin to prepare data to be analyzed
    anaObj = obj.prepareTensorWithFixBin(anaWindow, binWidth, 0, neus2consider, conds2consider);
    anaData = anaObj.TensorWithFixBin;

    nNeus = size(anaData, 1);
    nConds = size(anaData, 2);
    nTime = size(anaData, 4);

    p2way = nan(nNeus,nTime,3);
    p1way = nan(nNeus,nTime,1);

    for neu = 1:nNeus
        for t = 1:nTime
            refY = refData(neu,:,:);
            anaY = anaData(neu,:,:,t);
            conds = nan(size(anaY));
            for cond = 1:nConds
                conds(1,cond,:) = cond;
            end

            % 2-way ANOVA (epoch, condition, epoch*condition)
            y = [refY(:); anaY(:)];
            g1 = [zeros(numel(refY),1); ones(numel(anaY),1)];
            g2 = [conds(:); conds(:)];
            p2way(neu,t,:) = anovan(y,{g1, g2}, ...
                'model','interaction', ...
                'varnames', varNames, ...
                'display', 'off')';

            % 1-way ANOVA (condition)
            y = anaY(:);
            g = conds(:);
            p1way(neu,t) = anovan(y, g, ...
                'varnames', varNames(2), ...
                'display', 'off')';
        end
    end

    anovaResults.p2way = p2way;
    anovaResults.p1way = p1way;
    anovaResults.sign2way = nan(nTime,3);
    anovaResults.sign1way = nan(nTime,1);
    anovaResults.effectNames2way = effectNames2way;
    anovaResults.effectNames1way = effectNames1way;
    anovaResults.factorNames = varNames;

    for t = 1:nTime
        for vvar = 1:3
            anovaResults.sign2way(t,vvar) = nnz(p2way(:,t,vvar)<0.05)/numel(p2way(:,t,vvar))*100;
        end
        anovaResults.sign1way(t) = nnz(p1way(:,t)<0.05)/numel(p1way(:,t))*100;
    end

    timeAxis = linspace(anaWindow(1), anaWindow(3), nTime);
    markerPositions = squeeze(mean(anaObj.MarkerTensorForFixBin,[1 2 3],'omitmissing'));
    markerPositions = markerPositions(:)';
    markerTimes = anaWindow(1) + markerPositions .* ((anaWindow(3)-anaWindow(1))/nTime);
    markerIsPlotted = isfinite(markerTimes) & markerTimes>=anaWindow(1) & markerTimes<=anaWindow(3);

    neuCondStr = sprintf('neus %d-%d (n=%d) | conds [%s]', min(neus2consider), max(neus2consider), numel(neus2consider), num2str(conds2consider));
    if displayFigures
        figure
        plot(timeAxis, anovaResults.sign2way, 'LineWidth', 3)
        hold on
        plotMarkerLines();
        ylabel("% of units modulated")
        xlabel('Time (s)')
        title({["2-way ANOVA | ref: " num2str(refWindow) " | ana: " num2str(anaWindow) ...
            " | bin: " num2str(binWidth(1)) "s step: " num2str(binWidth(2)) "s"], neuCondStr})
        legend(effectNames2way, 'Interpreter', 'none')

        figure
        plot(timeAxis, anovaResults.sign1way, 'LineWidth', 3)
        hold on
        plotMarkerLines();
        ylabel("% of units modulated")
        xlabel('Time (s)')
        title({["1-way ANOVA | ana: " num2str(anaWindow) ...
            " | bin: " num2str(binWidth(1)) "s step: " num2str(binWidth(2)) "s"], neuCondStr})
        legend(effectNames1way, 'Interpreter', 'none')
    end

    function names = makeEffectNames(effectTerms, factorNames)
        names = cell(1, size(effectTerms, 1));
        for termIdx = 1:size(effectTerms, 1)
            factorIdx = find(effectTerms(termIdx,:));
            names{termIdx} = strjoin(factorNames(factorIdx), '*');
        end
    end

    function names = normalizeFactorNames(inputNames)
        if isempty(inputNames)
            names = {'epoch', 'condition'};
            return
        end

        if ischar(inputNames)
            names = {inputNames};
        elseif isstring(inputNames)
            names = cellstr(inputNames(:))';
        elseif iscell(inputNames)
            names = cellfun(@char, inputNames(:)', 'UniformOutput', false);
        else
            error('factorNames must be a cell array, string array, or character vector.');
        end

        if any(cellfun(@isempty, names))
            error('factorNames cannot contain empty names.');
        end

        epochIdx = find(strcmpi(names, 'epoch'));
        if isempty(epochIdx)
            names = [{'epoch'}, names];
        else
            names(epochIdx) = [];
            names = [{'epoch'}, names];
        end

        if numel(names) ~= 2
            error(['factorNames must contain one name for each ANOVA factor. ' ...
                'Expected 2 names including epoch, or 1 condition-factor name without epoch.']);
        end
    end

    function plotMarkerLines()
        if isempty(markerTimes) || ~any(markerIsPlotted)
            return
        end

        plottedMarkerIdx = find(markerIsPlotted);
        if isempty(markerNames)
            xline(markerTimes(plottedMarkerIdx));
            return
        end

        if numel(markerNames) == numel(markerTimes)
            plottedMarkerNames = markerNames(plottedMarkerIdx);
        elseif numel(markerNames) == numel(plottedMarkerIdx)
            plottedMarkerNames = markerNames;
        else
            error(['markerNames must contain either one name for each marker ' ...
                'or one name for each plotted marker.']);
        end

        xline(markerTimes(plottedMarkerIdx), '-', plottedMarkerNames, ...
            'LabelOrientation', 'aligned', ...
            'LabelHorizontalAlignment', 'left');
    end
end
