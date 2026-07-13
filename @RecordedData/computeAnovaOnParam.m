function anovaResults = computeAnovaOnParam(obj, neus2consider, conds2consider, params, refWindow, anaWindow, binWidth, varargin)
    % COMPUTEANOVAONPARAM - Compute ANOVA with epoch and task-parameter factors
    %
    % INPUTS:
    %   obj                     - RecordedData object
    %   neus2consider           - Neuron indices (numeric vector)
    %   conds2consider          - Condition indices (numeric vector)
    %   params                  - Parameter matrix (conditions x parameters)
    %   refWindow               - Reference [start(s), marker, end(s)] (numeric vector, length 3)
    %   anaWindow               - Analysis [start(s), marker, end(s)] (numeric vector, length 3)
    %   binWidth                - [bin_size, step] in seconds (numeric vector, length 2)
    %   'markerNames'           - Cell array of marker names
    %   'display'               - Logical flag to show the summary figures (default: true)
    %   'factorNames'           - Names for epoch and parameter factors. If 'epoch'
    %                             is omitted, it is added as the first factor.
    %
    % OUTPUT:
    %   anovaResults            - Struct containing two ANOVA results:
    %                             - withEpoch: epoch + task-parameter factors
    %                             - withoutEpoch: task-parameter factors only

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

    if abs((refWindow(3)-refWindow(1)) - binWidth(1)) > 0.001
        disp("NOTE: The binWidth chosen for the analysis window doesn't match the size of the reference window")
    end
    
    nConds = numel(conds2consider);

    if isvector(params)
        params = params(:);
    end

    if size(params, 1) ~= nConds
        error(['params must have one row for each selected condition in the prepared tensor.' ...
             ' Expected %d rows from conds2consider but got %d.'], nConds, size(params, 1));
    end

    nParams = size(params, 2);
    varNames = normalizeFactorNames(p.Results.factorNames, nParams);
    paramVarNames = varNames(2:end);
    
    refObj = obj.prepareTensorWithFixBin(refWindow, [refWindow(3)-refWindow(1) 0], neus2consider, conds2consider);
    refObj = refObj.splitTensorWithParams(params, true);
    refData = refObj.TensorWithParams;

    anaObj = obj.prepareTensorWithFixBin(anaWindow, binWidth, neus2consider, conds2consider);
    anaObj = anaObj.splitTensorWithParams(params, true);
    anaData = anaObj.TensorWithParams;

    nNeus = size(anaData, 1);
    nTrials = size(anaData, nParams+2);
    nTime = size(anaData, nParams+3);
    
    expectedTrials = size(anaObj.TensorWithFixBin, 3);
    expectedTime = size(anaObj.TensorWithFixBin, 4);
    if nTrials ~= expectedTrials || nTime ~= expectedTime
        error(['TensorWithParams must have dimensions neurons x param1 x ... x trials x time. ' ...
            'Expected trials in dimension %d and time in dimension %d.'], nParams+2, nParams+3);
    end

    nFactors = nParams + 1; % epoch + task parameters

    pWithEpoch = [];
    pWithoutEpoch = [];
    termsWithEpoch = [];
    termsWithoutEpoch = [];

    paramLevels = anaObj.TensorWithParamsInfo.paramLevels;
    paramDims = cellfun(@numel, paramLevels);
    nParamCombinations = prod(paramDims);
    paramGrids = cell(1, nParams);
    [paramGrids{:}] = ndgrid(paramLevels{:});

    % Match MATLAB's column-major order for TensorWithParams(:):
    % param1 varies fastest, then param2, ..., then trials.
    paramGroups = cell(1, nParams);
    for paramIdx = 1:nParams
        paramGroups{paramIdx} = repmat(paramGrids{paramIdx}(:), nTrials, 1);
    end

    nObsPerEpoch = nParamCombinations*nTrials;
    epochGroup = [zeros(nObsPerEpoch,1); ones(nObsPerEpoch,1)];
    groupsWithEpoch = cell(1, nFactors);
    groupsWithEpoch{1} = epochGroup;
    groupsWithoutEpoch = cell(1, nParams);
    for paramIdx = 1:nParams
        groupsWithEpoch{paramIdx+1} = [paramGroups{paramIdx}; paramGroups{paramIdx}];
        groupsWithoutEpoch{paramIdx} = paramGroups{paramIdx};
    end

    paramSubs = repmat({':'}, 1, nParams);
    for neu = 1:nNeus
        for t = 1:nTime
            refY = reshape(refData(neu, paramSubs{:}, :, :), [], 1);
            anaY = reshape(anaData(neu, paramSubs{:}, :, t), [], 1);
            y = [refY; anaY];

            valid = isfinite(y);
            validGroups = cellfun(@(g) g(valid), groupsWithEpoch, 'UniformOutput', false);

            [currentP,~,~,currentTerms] = anovan(y(valid), validGroups, ...
                'model','interaction', ...
                'varnames', varNames, ...
                'display', 'off');

            if isempty(pWithEpoch)
                pWithEpoch = nan(nNeus, nTime, numel(currentP));
                termsWithEpoch = currentTerms;
            end
            pWithEpoch(neu,t,:) = reshape(currentP, 1, 1, []);

            valid = isfinite(anaY);
            validGroups = cellfun(@(g) g(valid), groupsWithoutEpoch, 'UniformOutput', false);

            [currentP,~,~,currentTerms] = anovan(anaY(valid), validGroups, ...
                'model','interaction', ...
                'varnames', paramVarNames, ...
                'display', 'off');

            if isempty(pWithoutEpoch)
                pWithoutEpoch = nan(nNeus, nTime, numel(currentP));
                termsWithoutEpoch = currentTerms;
            end
            pWithoutEpoch(neu,t,:) = reshape(currentP, 1, 1, []);
        end
    end

    effectNamesWithEpoch = makeEffectNames(termsWithEpoch, varNames);
    effectNamesWithoutEpoch = makeEffectNames(termsWithoutEpoch, paramVarNames);
    signWithEpoch = computeSignPercent(pWithEpoch);
    signWithoutEpoch = computeSignPercent(pWithoutEpoch);

    anovaResults.withEpoch.p = pWithEpoch;
    anovaResults.withEpoch.sign = signWithEpoch;
    anovaResults.withEpoch.effectNames = effectNamesWithEpoch;
    anovaResults.withEpoch.terms = termsWithEpoch;
    anovaResults.withEpoch.factorNames = varNames;
    anovaResults.withEpoch.params = params;

    anovaResults.withoutEpoch.p = pWithoutEpoch;
    anovaResults.withoutEpoch.sign = signWithoutEpoch;
    anovaResults.withoutEpoch.effectNames = effectNamesWithoutEpoch;
    anovaResults.withoutEpoch.terms = termsWithoutEpoch;
    anovaResults.withoutEpoch.factorNames = paramVarNames;
    anovaResults.withoutEpoch.params = params;

    anovaResults.pWithEpoch = pWithEpoch;
    anovaResults.signWithEpoch = signWithEpoch;
    anovaResults.effectNamesWithEpoch = effectNamesWithEpoch;
    anovaResults.termsWithEpoch = termsWithEpoch;
    anovaResults.factorNamesWithEpoch = varNames;
    anovaResults.pWithoutEpoch = pWithoutEpoch;
    anovaResults.signWithoutEpoch = signWithoutEpoch;
    anovaResults.effectNamesWithoutEpoch = effectNamesWithoutEpoch;
    anovaResults.termsWithoutEpoch = termsWithoutEpoch;
    anovaResults.factorNamesWithoutEpoch = paramVarNames;

    % Backward-compatible aliases for the epoch-inclusive result.
    anovaResults.p = pWithEpoch;
    anovaResults.sign = signWithEpoch;
    anovaResults.effectNames = effectNamesWithEpoch;
    anovaResults.terms = termsWithEpoch;
    anovaResults.params = params;
    anovaResults.factorNames = varNames;

    if binWidth(2) > 0
        timeStep = binWidth(2);
    else
        timeStep = binWidth(1);
    end
    timeAxis = anaWindow(1) + (0:nTime-1)*timeStep;
    markerPositions = squeeze(mean(anaObj.MarkerTensorForFixBin,[1 2 3],'omitmissing'));
    markerPositions = markerPositions(:)';
    markerTimes = anaWindow(1) + markerPositions*timeStep;
    markerIsPlotted = isfinite(markerTimes) & markerTimes>=anaWindow(1) & markerTimes<=anaWindow(3);

    anovaResults.time = timeAxis;
    anovaResults.markerTimes = markerTimes(markerIsPlotted);
    anovaResults.withEpoch.time = timeAxis;
    anovaResults.withEpoch.markerTimes = markerTimes(markerIsPlotted);
    anovaResults.withoutEpoch.time = timeAxis;
    anovaResults.withoutEpoch.markerTimes = markerTimes(markerIsPlotted);

    neuCondStr = sprintf('neus %d-%d (n=%d) | conds [%s]', min(neus2consider), max(neus2consider), numel(neus2consider), num2str(conds2consider));
    if displayFigures
        figure
        plot(timeAxis, anovaResults.withEpoch.sign, 'LineWidth', 3)
        hold on
        plotMarkerLines();
        ylabel("% of units modulated")
        xlabel('Time (s)')
        title({["ANOVA with epoch | ref: " num2str(refWindow) " | ana: " num2str(anaWindow) ...
            " | bin: " num2str(binWidth(1)) "s step: " num2str(binWidth(2)) "s"], neuCondStr})
        legend(effectNamesWithEpoch, 'Interpreter', 'none')

        figure
        plot(timeAxis, anovaResults.withoutEpoch.sign, 'LineWidth', 3)
        hold on
        plotMarkerLines();
        ylabel("% of units modulated")
        xlabel('Time (s)')
        title({["ANOVA without epoch | ana: " num2str(anaWindow) ...
            " | bin: " num2str(binWidth(1)) "s step: " num2str(binWidth(2)) "s"], neuCondStr})
        legend(effectNamesWithoutEpoch, 'Interpreter', 'none')
    end

    function signPercent = computeSignPercent(pValues)
        signPercent = nan(nTime, size(pValues, 3));
        for timeIdx = 1:nTime
            for effectIdx = 1:size(pValues, 3)
                signPercent(timeIdx,effectIdx) = nnz(pValues(:,timeIdx,effectIdx)<0.05)/ ...
                    numel(pValues(:,timeIdx,effectIdx))*100;
            end
        end
    end

    function names = makeEffectNames(effectTerms, factorNames)
        names = cell(1, size(effectTerms, 1));
        for termIdx = 1:size(effectTerms, 1)
            factorIdx = find(effectTerms(termIdx,:));
            names{termIdx} = strjoin(factorNames(factorIdx), '*');
        end
    end

    function names = normalizeFactorNames(inputNames, nTaskParams)
        defaultParamNames = cellstr(compose('param%d', 1:nTaskParams));
        if isempty(inputNames)
            names = [{'epoch'}, defaultParamNames(:)'];
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

        expectedNames = nTaskParams + 1;
        if numel(names) ~= expectedNames
            error(['factorNames must contain one name for each ANOVA factor. ' ...
                'Expected %d names including epoch, or %d task-parameter names without epoch.'], ...
                expectedNames, nTaskParams);
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
