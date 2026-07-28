function obj = splitTensorWithParams(obj, params, fixedBinFlag)
    % SPLITTENSORWITHPARAMS - Prepare multi-dimensional tensor split by task parameters
    %
    % INPUTS:
    %   obj                 - RecordedData object
    %   params              - Parameter matrix (conditions x parameters)
    %   fixedBinFlag        - true for fixed bins, false for variable bins (logical)
    %
    % OUTPUT:
    %   obj                 - Updated RecordedData object with newly filled:
    %            .TensorWithParams        - Trial-by-trial tensor (neurons x param1 x param2 x ... x trials x time)
    %            .TensorWithParamsCondAvg - Condition-averaged tensor (neurons x param1 x param2 x ... x time)
    %            .TensorWithParamsInfo    - Struct with metadata (params, fixedBinFlag, etc.)
    %
    % DESCRIPTION:
    %   Restructures trial-by-trial tensor into multi-dimensional format where each
    %   parameter gets its own dimension. Also stores the condition-averaged version
    %   used by dPCA.

    if isempty(params) || ~ismatrix(params)
        error('params must be a non-empty column vector or matrix');
    end

    if fixedBinFlag
        if isempty(obj.TensorWithFixBin)
            error('Call prepareTensorWithFixBin before splitTensorWithParams.');
        end
        sourceTensor = obj.TensorWithFixBin;
    else
        if isempty(obj.TensorWithVariableBin)
            error('Call prepareTensorWithVariableBin before splitTensorWithParams.');
        end
        sourceTensor = obj.TensorWithVariableBin;
    end

    if fixedBinFlag
        sourceTS = obj.BinTimestampsForFixBin;
    else
        sourceTS = obj.BinTimestampsForVariableBin;
    end

    nPreparedConds = size(sourceTensor, 2);
    if isvector(params) && nPreparedConds > 1
        params = params(:);
    end

    if size(params, 1) ~= nPreparedConds
        error('params must have one row for each condition in the prepared tensor.');
    end

    nParams = size(params, 2);
    if size(unique(params, 'rows'), 1) ~= size(params, 1)
        warning(['params contains repeated values. Duplicate rows will map multiple conditions to the same TensorWithParams cell.']);
    end

    % Initialize TensorWithParamsInfo struct
    obj.TensorWithParamsInfo = struct();
    obj.TensorWithParamsInfo.params = params;
    obj.TensorWithParamsInfo.fixedBinFlag = fixedBinFlag;

    if fixedBinFlag && ~isempty(obj.MarkerTensorForFixBin)
        obj.TensorWithParamsInfo.MarkersFixedBin = squeeze(mean(obj.MarkerTensorForFixBin,[1 2 3],'omitnan'));
    end

    if fixedBinFlag
        obj.TensorWithParamsInfo.binWidth = obj.TensorWithFixBinInfo.binWidth;
    else
        obj.TensorWithParamsInfo.expectedBinWidth = obj.TensorWithVariableBinInfo.expectedBinWidth;
        if isfield(obj.TensorWithVariableBinInfo, 'nBins')
            obj.TensorWithParamsInfo.nBins = obj.TensorWithVariableBinInfo.nBins;
        end
    end

    paramLevels = cell(1, nParams);
    paramLevelIdx = nan(size(params));
    paramDims = nan(1, nParams);
    for p = 1:nParams
        paramLevels{p} = unique(params(:, p), 'stable');
        [~, paramLevelIdx(:, p)] = ismember(params(:, p), paramLevels{p});
        paramDims(p) = numel(paramLevels{p});
    end

    obj.TensorWithParamsInfo.paramLevels = paramLevels;

    nNeurons = size(sourceTensor, 1);
    nTrials  = size(sourceTensor, 3);
    nTime    = size(sourceTensor, 4);
    nPreparedConds2 = size(sourceTensor, 2);
    nParamCombinations = prod(paramDims);

    % Compute linear index of each condition's parameter combination.
    % sub2ind requires at least 2 size elements, so handle single-param case.
    comboLinearIdx = zeros(nPreparedConds2, 1);
    for cond = 1:nPreparedConds2
        subs = num2cell(paramLevelIdx(cond, :));
        if nParams == 1
            comboLinearIdx(cond) = subs{1};
        else
            comboLinearIdx(cond) = sub2ind(paramDims, subs{:});
        end
    end

    % Count how many conditions map to each combination to size the trial dim.
    comboCounts  = histcounts(comboLinearIdx, 0.5:(nParamCombinations + 0.5));
    maxMergedTrials = max(comboCounts) * nTrials;

    TensorWithParams    = nan([nNeurons, paramDims, maxMergedTrials, nTime]);
    TensorWithParamsCondAvg = nan([nNeurons, paramDims, nTime]);
    BinTimestampsForParams = nan([nNeurons, paramDims, maxMergedTrials, nTime]);

    % Track how many trials have been filled for each parameter combination.
    trialFill = zeros(1, nParamCombinations);

    % Populate tensor, stacking trials from conditions with the same params.
    for cond = 1:nPreparedConds2
        linIdx     = comboLinearIdx(cond);
        startTrial = trialFill(linIdx) + 1;
        endTrial   = startTrial + nTrials - 1;
        trialFill(linIdx) = endTrial;

        indices        = cell(1, nParams + 3);
        indices{1}     = ':';  % neurons
        for p = 1:nParams
            indices{p+1} = paramLevelIdx(cond, p);
        end
        indices{nParams+2} = startTrial:endTrial;  % merged trial range
        indices{end}   = ':';  % time

        TensorWithParams(indices{:}) = reshape(sourceTensor(:, cond, :, :), ...
            [nNeurons, ones(1, nParams), nTrials, nTime]);
        if ~isempty(sourceTS)
            BinTimestampsForParams(indices{:}) = reshape(sourceTS(:, cond, :, :), ...
                [nNeurons, ones(1, nParams), nTrials, nTime]);
        end
    end

    % Condition average: mean over the merged trials dimension.
    TensorWithParamsCondAvg = mean(TensorWithParams, nParams + 2, 'omitnan');
    TensorWithParamsCondAvg = reshape(TensorWithParamsCondAvg, [nNeurons, paramDims, nTime]);

    % Compute trial counts: N x param1 x param2 x ...
    % For each (neuron, param combination), count non-NaN trials from first time bin.
    idx = repmat({':'}, 1, ndims(TensorWithParams));
    idx{end} = 1;  % first time point
    firstTimeSlice = TensorWithParams(idx{:});  % [nNeurons, paramDims, maxMergedTrials]
    TrialNumWithParams = sum(~isnan(firstTimeSlice), nParams + 2);  % [nNeurons, paramDims, 1]
    TrialNumWithParams = reshape(TrialNumWithParams, [nNeurons, paramDims]);  % [nNeurons, paramDims]

    obj.TensorWithParams = TensorWithParams;
    obj.TensorWithParamsCondAvg = TensorWithParamsCondAvg;
    obj.TrialNumWithParams = TrialNumWithParams;
    obj.BinTimestampsForParams = BinTimestampsForParams;
end
