function [tensorOut, infoOut] = addVariable2Tensor(obj, tensorIn, var, requiredTS, opt)
%ADDVARIABLE2TENSOR  Append a preprocessed variable as a new slice in the last dimension.
%
% [tensorOut, infoOut] = obj.addVariable2Tensor(tensorIn, var, requiredTS)
% [tensorOut, infoOut] = obj.addVariable2Tensor(tensorIn, var, requiredTS, opt)
%
% Designed to be called iteratively to build a regressor tensor.  Each call
% appends one variable (after optional derivative / lag / filter processing)
% as a new slice along the last dimension of tensorIn.
%
% INPUTS
%   tensorIn   – Current regressor tensor  [d1 × … × dN × K].
%                Pass [] on the first call; tensorOut will be [d1 × … × dN × 1].
%
%   var        (required) Variable source.  Either:
%                • {'FR','FixBin'}      – uses obj.TensorWithFixBin
%                • {'FR','VariableBin'} – uses obj.TensorWithVariableBin
%                • {'TC','xCursor'}     – x coordinate from obj.TC (cell mode)
%                • {'TC','yCursor'}     – y coordinate from obj.TC (cell mode)
%                • {'KIN',…}            – not yet implemented
%                • {'EY',…}             – not yet implemented
%                • numeric tensor       – opt.timeDim and opt.varTS are required
%                • cell array           – {nNeu×nCond}{1×nTrials}[1×T];
%                                         opt.varTS must match the outer shape.
%
%   requiredTS (required) Numeric tensor whose non-timeDim dimensions match
%                the variable.  Its size along timeDim defines the output sample
%                count and its values give the target timestamps for resampling.
%                Resampling is performed per slice as:
%                  interp1(varTS_slice, sig_slice, requiredTS_slice,'linear','extrap')
%
%   opt        (optional) Options struct:
%
%     .varTS      Timestamp source for the variable.
%                 Required for numeric tensor and cell array inputs.
%                 Auto-inferred for named sources:
%                   {'FR','FixBin'}      → obj.BinTimestampsForFixBin
%                   {'FR','VariableBin'} → obj.BinTimestampsForVariableBin
%                   {'TC',…}            → obj.TCTS
%                 Must be a numeric tensor (same shape as var) for tensor-mode
%                 sources, or a cell array {nNeu×nCond}{1×nTrials}[1×T] for
%                 cell-mode sources (TC, cell array).
%
%     .timeDim   Positive integer: index of the time dimension inside the variable
%                 (not tensorIn).  Auto-set to 4 for all named sources and cell
%                 arrays (assumed format: neu × cond × trial × time).
%                 Required for numeric tensor inputs.
%
%     .reuseAcrossUnits  Optional logical flag.
%                 For cell-mode sources (TC, EY, KIN, cell arrays) — default true:
%                   true  → process only first neuron's data, copy across all
%                            neurons.  A warning is issued when
%                            obj.simultaneousRecording is false.
%                   false → process every neuron independently.  A warning is
%                            issued when obj.simultaneousRecording is true.
%                 For FR: always false; setting true raises an error.
%                 For numeric tensors: if true, opt.unitDim is required and data
%                   is processed once and copied across the unit dimension.
%
%     .unitDim   Positive integer: neuron dimension index in a numeric variable
%                 tensor.  Required when reuseAcrossUnits = true for numeric
%                 tensor inputs.  Not needed for cell-mode inputs (neuron is
%                 always the first outer dimension).
%
%     .der        Derivative order (non-negative integer, default 0).
%                 Each order applies central_diff once (output length preserved).
%
%     .lag        Lag in samples (integer, default 0).
%                 Positive = delay (shift right); negative = anticipate (shift left).
%                 Boundary samples are padded with the mean of the valid portion.
%
%     .filter     Convolution kernel for conv(…,'same') (default []).
%                 A warning is issued if sum(filter) ≠ 1.
%
%     .varName    Label string for bookkeeping (default '').
%
%     .infoIn     Info struct from a previous call (default struct()).
%
% OUTPUTS
%   tensorOut – [d1 × … × dN × K+1].
%   infoOut   – Copy of opt.infoIn with a new field 'n<K+1>' storing all settings
%               used for this addition.  First variable → 'n1', second → 'n2', etc.
%
% EXAMPLE
%   % Build a regressor tensor from FR (FixBin) and TC x-coordinate.
%   opt1.varName = 'FiringRate';
%   [T, info] = rd.addVariable2Tensor([], {'FR','FixBin'}, requiredTS, opt1);
%
%   opt2.der     = 1;          % first derivative
%   opt2.lag     = 2;          % delay by 2 bins
%   opt2.varName = 'xVelocity';
%   opt2.infoIn  = info;
%   [T, info] = rd.addVariable2Tensor(T, {'TC','xCursor'}, requiredTS, opt2);

    %% 1 ── Validate inputs and apply opt defaults ──────────────────────────

    if nargin < 5 || isempty(opt)
        opt = struct();
    end
    if ~isstruct(opt)
        error('RecordedData:addVariable2Tensor:invalidOpt', ...
            'opt must be a struct.');
    end
    if isempty(var)
        error('RecordedData:addVariable2Tensor:missingVar', ...
            'var is required and cannot be empty.');
    end
    if isempty(requiredTS)
        error('RecordedData:addVariable2Tensor:missingRequiredTS', ...
            'requiredTS is required and cannot be empty.');
    end

    % timeDim may be auto-set later depending on var type (Section 3)
    if ~isfield(opt, 'timeDim') || isempty(opt.timeDim)
        timeDim = [];
    else
        if ~isscalar(opt.timeDim) || opt.timeDim < 1 || opt.timeDim ~= floor(opt.timeDim)
            error('RecordedData:addVariable2Tensor:invalidTimeDim', ...
                'opt.timeDim must be a positive integer scalar.');
        end
        timeDim = opt.timeDim;
    end

    if ~isfield(opt, 'der')    || isempty(opt.der),    opt.der    = 0;        end
    if ~isfield(opt, 'lag')    || isempty(opt.lag),    opt.lag    = 0;        end
    if ~isfield(opt, 'filter'),                         opt.filter = [];       end
    if ~isfield(opt, 'varName'),                        opt.varName = '';      end
    if ~isfield(opt, 'infoIn') || isempty(opt.infoIn), opt.infoIn = struct(); end

    if ~isfield(opt, 'unitDim') || isempty(opt.unitDim)
        opt.unitDim = [];
    else
        if ~isnumeric(opt.unitDim) || ~isscalar(opt.unitDim) || ...
                ~isfinite(opt.unitDim) || opt.unitDim < 1 || ...
                opt.unitDim ~= floor(opt.unitDim)
            error('RecordedData:addVariable2Tensor:invalidUnitDim', ...
                'opt.unitDim must be a positive integer scalar.');
        end
    end

    reuseOptionSet = isfield(opt, 'reuseAcrossUnits') && ~isempty(opt.reuseAcrossUnits);
    if ~reuseOptionSet
        opt.reuseAcrossUnits = [];
    elseif ~(islogical(opt.reuseAcrossUnits) || isnumeric(opt.reuseAcrossUnits)) || ...
            ~isscalar(opt.reuseAcrossUnits) || ...
            ~ismember(double(opt.reuseAcrossUnits), [0, 1])
        error('RecordedData:addVariable2Tensor:invalidReuseAcrossUnits', ...
            'opt.reuseAcrossUnits must be a logical scalar.');
    else
        opt.reuseAcrossUnits = logical(opt.reuseAcrossUnits);
    end

    if ~isempty(opt.filter) && abs(sum(opt.filter) - 1) > 1e-6
        warning('RecordedData:addVariable2Tensor:filterPower', ...
            'opt.filter sums to %.6f (not 1). Proceeding anyway.', sum(opt.filter));
    end

    %% 2 ── First-call flag (tensorIn shape resolved in Section 4) ──────────
    %   nDataDims is determined from the variable in Section 4 so that we
    %   can detect and correct the K=1 trailing-singleton issue:
    %   when K=1, MATLAB's size() silently drops the last dimension, making
    %   numel(size(tensorIn)) equal nDataDims instead of nDataDims+1.

    isFirstCall = isempty(tensorIn);

    %% 3 ── Extract variable from var ──────────────────────────────────────
    %   Two operating modes:
    %     Tensor mode – var and varTS are numeric tensors.
    %     Cell mode   – var and varTS are {nNeu×nCond}{1×nTrials}[1×T] cell arrays.

    isCellMode = false;
    varSource  = '';

    if iscell(var) && ~isempty(var) && (ischar(var{1}) || isstring(var{1}))
        % Named source specification: {'FR','...'}, {'TC','...'}, etc.
        varSource = upper(char(var{1}));
        switch varSource

            case 'FR'
                if numel(var) < 2
                    error('RecordedData:addVariable2Tensor:missingSpec', ...
                        'Specify {''FR'',''FixBin''} or {''FR'',''VariableBin''}.');
                end
                switch lower(char(var{2}))
                    case 'fixbin'
                        if isempty(obj.TensorWithFixBin)
                            error('RecordedData:addVariable2Tensor:emptyData', ...
                                'obj.TensorWithFixBin is empty.');
                        end
                        varTensor = obj.TensorWithFixBin;
                        if ~isfield(opt, 'varTS') || isempty(opt.varTS)
                            opt.varTS = obj.BinTimestampsForFixBin;
                        end
                    case 'variablebin'
                        if isempty(obj.TensorWithVariableBin)
                            error('RecordedData:addVariable2Tensor:emptyData', ...
                                'obj.TensorWithVariableBin is empty.');
                        end
                        varTensor = obj.TensorWithVariableBin;
                        if ~isfield(opt, 'varTS') || isempty(opt.varTS)
                            opt.varTS = obj.BinTimestampsForVariableBin;
                        end
                    otherwise
                        error('RecordedData:addVariable2Tensor:unknownBinType', ...
                            'Unknown FR bin type "%s". Use ''FixBin'' or ''VariableBin''.', ...
                            char(var{2}));
                end
                if isempty(timeDim), timeDim = 4; end

            case 'TC'
                if numel(var) < 2
                    error('RecordedData:addVariable2Tensor:missingSpec', ...
                        'Specify {''TC'',''xCursor''} or {''TC'',''yCursor''}.');
                end
                if isempty(obj.TC)
                    error('RecordedData:addVariable2Tensor:emptyData', ...
                        'obj.TC is empty.');
                end
                switch lower(char(var{2}))
                    case 'xcursor', tcField = 'xCursor';
                    case 'ycursor', tcField = 'yCursor';
                    otherwise
                        error('RecordedData:addVariable2Tensor:unknownCoord', ...
                            'Unknown TC coordinate "%s". Use ''xCursor'' or ''yCursor''.', ...
                            char(var{2}));
                end
                isCellMode = true;
                varCell = cellfun( ...
                    @(condCell) cellfun(@(s) s.(tcField)(:)', condCell, ...
                        'UniformOutput', false), ...
                    obj.TC, 'UniformOutput', false);
                if ~isfield(opt, 'varTS') || isempty(opt.varTS)
                    if isempty(obj.TCTS)
                        error('RecordedData:addVariable2Tensor:emptyData', ...
                            'obj.TCTS is empty; provide opt.varTS manually.');
                    end
                    opt.varTS = obj.TCTS;
                end
                if isempty(timeDim), timeDim = 4; end

            case {'KIN', 'EY'}
                error('RecordedData:addVariable2Tensor:notImplemented', ...
                    '"%s" support is not yet implemented.', varSource);

            otherwise
                error('RecordedData:addVariable2Tensor:unknownSource', ...
                    'Unknown source "%s". Use ''FR'', ''TC'', ''KIN'' or ''EY''.', varSource);
        end

    elseif iscell(var)
        % Direct cell array input: {nNeu×nCond}{1×nTrials}[1×T]
        isCellMode = true;
        varCell    = var;
        if isempty(timeDim), timeDim = 4; end
        if ~isfield(opt, 'varTS') || isempty(opt.varTS)
            error('RecordedData:addVariable2Tensor:missingVarTS', ...
                'opt.varTS is required for cell array inputs.');
        end

    elseif isnumeric(var)
        varTensor = var;
        if isempty(timeDim)
            error('RecordedData:addVariable2Tensor:missingTimeDim', ...
                'opt.timeDim is required for numeric tensor inputs.');
        end
        if ~isfield(opt, 'varTS') || isempty(opt.varTS)
            error('RecordedData:addVariable2Tensor:missingVarTS', ...
                'opt.varTS is required for numeric tensor inputs.');
        end

    else
        error('RecordedData:addVariable2Tensor:invalidVar', ...
            'var must be a named source spec, cell array, or numeric tensor.');
    end

    %% 4 ── Resolve tensor shape; correct for K=1 trailing-singleton ─────────
    %   nDataDims is derived from the variable here (after extraction) so it
    %   is available for the tensorIn dimension check that follows.

    if isCellMode
        nDataDims = numel(size(requiredTS));
    else
        nDataDims = numel(size(varTensor));
    end

    if isFirstCall
        nSlices = 0;
        nDimsIn = nDataDims + 1;
        if isCellMode
            sizeIn = [size(requiredTS), 1];   % placeholder
        else
            sizeIn = [size(varTensor), 1];    % placeholder: K = 1 after first addition
        end
    else
        sizeIn_raw  = size(tensorIn);
        nDimsIn_raw = numel(sizeIn_raw);

        if nDimsIn_raw == nDataDims
            % K=1: MATLAB's size() dropped the trailing singleton — restore it
            nDimsIn = nDataDims + 1;
            nSlices  = 1;
            sizeIn   = [sizeIn_raw, 1];
        elseif nDimsIn_raw == nDataDims + 1
            nDimsIn = nDimsIn_raw;
            nSlices  = sizeIn_raw(end);
            sizeIn   = sizeIn_raw;
        else
            error('RecordedData:addVariable2Tensor:tensorInDimMismatch', ...
                ['tensorIn has %d dimension(s) but the variable has %d; ' ...
                 'expected %d (data dims, K=1 dropped) or %d (data dims + K).'], ...
                nDimsIn_raw, nDataDims, nDataDims, nDataDims + 1);
        end
    end

    %% 5 ── Validate variable shape against tensorIn ───────────────────────

    if isCellMode
        % Cell mode: reference shape is requiredTS; validate varTS cell layout
        sizeVar     = size(requiredTS);
        cellOuterSz = size(varCell);
        if ~iscell(opt.varTS) || ~isequal(size(opt.varTS), cellOuterSz)
            error('RecordedData:addVariable2Tensor:varTSshapeMismatch', ...
                ['In cell mode opt.varTS must be a cell array with outer shape ' ...
                 '[%s] matching var.'], num2str(cellOuterSz));
        end
    else
        sizeVar  = size(varTensor);
        nDimsVar = numel(sizeVar);

        if nDimsVar ~= nDataDims
            error('RecordedData:addVariable2Tensor:dimMismatch', ...
                ['var has %d dimension(s); expected %d ' ...
                 '(one fewer than tensorIn).'], nDimsVar, nDataDims);
        end

        for d = 1:nDataDims
            if d ~= timeDim && sizeVar(d) ~= sizeIn(d)
                error('RecordedData:addVariable2Tensor:sizeMismatch', ...
                    ['var size at dimension %d is %d, but tensorIn ' ...
                     'size at dimension %d is %d.'], d, sizeVar(d), d, sizeIn(d));
            end
        end

        if ~isequal(size(opt.varTS), sizeVar)
            error('RecordedData:addVariable2Tensor:varTSshapeMismatch', ...
                'opt.varTS must have the same shape as var.');
        end
        if timeDim > nDataDims
            error('RecordedData:addVariable2Tensor:invalidTimeDim', ...
                'opt.timeDim (%d) exceeds the number of data dimensions (%d).', ...
                timeDim, nDataDims);
        end
        if ~isempty(opt.unitDim) && opt.unitDim > nDataDims
            error('RecordedData:addVariable2Tensor:invalidUnitDim', ...
                'opt.unitDim (%d) exceeds the number of data dimensions (%d).', ...
                opt.unitDim, nDataDims);
        end
    end

    %% 6 ── Determine reuseAcrossUnits and reuseDim ────────────────────────

    if isCellMode
        % Case A: TC, EY, KIN, or direct cell array.
        % Default reuseAcrossUnits = true; warn when it contradicts simultaneousRecording.
        if reuseOptionSet && opt.reuseAcrossUnits
            reuseAcrossUnits = true;
            if ~obj.simultaneousRecording
                warning('RecordedData:addVariable2Tensor:reuseNotSimultaneous', ...
                    ['Re-using data across units even though ' ...
                     'obj.simultaneousRecording is false because ' ...
                     'opt.reuseAcrossUnits was set to true. Every unit ' ...
                     'will receive the same data as the first one.']);
            end
        elseif reuseOptionSet && ~opt.reuseAcrossUnits
            reuseAcrossUnits = false;
            if obj.simultaneousRecording
                warning('RecordedData:addVariable2Tensor:noReuseSimultaneous', ...
                    ['Not re-using data across units even though ' ...
                     'obj.simultaneousRecording is true because ' ...
                     'opt.reuseAcrossUnits was set to false.']);
            end
        else
            % Not set: default true for cell-mode sources
            reuseAcrossUnits = true;
            if ~obj.simultaneousRecording
                warning('RecordedData:addVariable2Tensor:reuseNotSimultaneous', ...
                    ['Re-using data across units even though ' ...
                     'obj.simultaneousRecording is false (default for cell-mode ' ...
                     'sources). Every unit will receive the same data as the first one.']);
            end
        end
        reuseDim = 1;   % neuron is always the first outer dimension for cell arrays

    elseif strcmpi(varSource, 'FR')
        % Case B: FR — always process each neuron independently.
        if reuseOptionSet && opt.reuseAcrossUnits
            error('RecordedData:addVariable2Tensor:invalidReuseForFR', ...
                ['FR variables are always processed independently per neuron; ' ...
                 'opt.reuseAcrossUnits cannot be set to true.']);
        end
        reuseAcrossUnits = false;
        reuseDim         = [];

    else
        % Case C: numeric tensor.
        if reuseOptionSet && opt.reuseAcrossUnits
            if isempty(opt.unitDim)
                error('RecordedData:addVariable2Tensor:missingUnitDim', ...
                    ['opt.unitDim is required when opt.reuseAcrossUnits = true ' ...
                     'for numeric tensor inputs.']);
            end
            if opt.unitDim == timeDim
                error('RecordedData:addVariable2Tensor:unitDimIsTimeDim', ...
                    'opt.unitDim must identify a non-time dimension.');
            end
            reuseAcrossUnits = true;
            reuseDim         = opt.unitDim;
        else
            reuseAcrossUnits = false;
            reuseDim         = [];
        end
    end

    % Safety: if requiredTS varies across reuseDim, independent resampling
    % is needed per unit — silently disable reuse.
    if reuseAcrossUnits && ~isRepeatedAlongDim(requiredTS, reuseDim)
        reuseAcrossUnits = false;
    end

    %% 7 ── Validate requiredTS ────────────────────────────────────────────

    sizeReqTS = size(requiredTS);
    if numel(sizeReqTS) ~= nDataDims
        error('RecordedData:addVariable2Tensor:requiredTSdimMismatch', ...
            'requiredTS must have %d dimension(s) to match the variable.', nDataDims);
    end
    for d = 1:nDataDims
        if d ~= timeDim && sizeReqTS(d) ~= sizeVar(d)
            error('RecordedData:addVariable2Tensor:requiredTSsizeMismatch', ...
                ['requiredTS size at dimension %d is %d, but the variable ' ...
                 'size at that dimension is %d.'], d, sizeReqTS(d), sizeVar(d));
        end
    end
    if ~isFirstCall && sizeReqTS(timeDim) ~= sizeIn(timeDim)
        error('RecordedData:addVariable2Tensor:requiredTSlengthMismatch', ...
            ['requiredTS has %d sample(s) along timeDim, but tensorIn ' ...
             'has %d. All added variables must share the same target sample count.'], ...
            sizeReqTS(timeDim), sizeIn(timeDim));
    end

    %% 8 ── Determine output length and shape from requiredTS ──────────────

    targetLen = sizeReqTS(timeDim);
    sizeOut   = sizeVar;
    sizeOut(timeDim) = targetLen;

    %% 9 ── Iterate over non-signal indices; resample and process each 1-D signal ──

    newSlice  = nan(sizeOut);
    otherDims = [1:timeDim-1, timeDim+1:nDataDims];

    if isCellMode
        % Cell mode: varCell{outer...}{t}[1×T] – trials stored as nested cells.
        if timeDim == 1
            error('RecordedData:addVariable2Tensor:cellModeTimeDim', ...
                'Cell mode requires timeDim >= 2 because trials occupy timeDim - 1.');
        end
        trialDim      = timeDim - 1;
        cellOuterDims = [1:trialDim-1, timeDim+1:nDataDims];
        cellOuterSz   = size(varCell);

        % MATLAB omits trailing singleton dimensions from size(). Restore
        % them so outSubs has one subscript per outer tensor dimension.
        if numel(cellOuterSz) < numel(cellOuterDims)
            cellOuterSz(end+1:numel(cellOuterDims)) = 1;
        elseif numel(cellOuterSz) > numel(cellOuterDims)
            if any(cellOuterSz(numel(cellOuterDims)+1:end) ~= 1)
                error('RecordedData:addVariable2Tensor:cellOuterShapeMismatch', ...
                    'Cell-array outer dimensions do not match the variable tensor.');
            end
            cellOuterSz = cellOuterSz(1:numel(cellOuterDims));
        end
        if ~isequal(cellOuterSz, sizeVar(cellOuterDims))
            error('RecordedData:addVariable2Tensor:cellOuterShapeMismatch', ...
                'Cell-array outer dimensions do not match the variable tensor.');
        end

        processOuterSz = cellOuterSz;
        if reuseAcrossUnits
            % reuseDim == 1 always for cell mode; locate it in cellOuterDims
            unitOuterPos = find(cellOuterDims == reuseDim, 1, 'first');
            if isempty(unitOuterPos)
                error('RecordedData:addVariable2Tensor:invalidUnitDim', ...
                    'The cell input must have neuron as its first outer dimension.');
            end
            processOuterSz(unitOuterPos) = 1;
        end
        nOuterCombos = prod(processOuterSz);

        for outIdx = 1:nOuterCombos
            outSubs     = ind2subVec(processOuterSz, outIdx);
            outSubsCell = num2cell(outSubs);

            trialCell = varCell{outSubsCell{:}};
            tsCell    = opt.varTS{outSubsCell{:}};

            for t = 1:numel(trialCell)
                sig    = trialCell{t}(:);
                tsOrig = tsCell{t}(:);

                % Build index into newSlice / requiredTS
                idxOut    = cell(1, nDataDims);
                outerRank = 0;
                for d = 1:nDataDims
                    if d == timeDim
                        idxOut{d} = ':';
                    elseif d == trialDim   % trial dimension
                        idxOut{d} = t;
                    else
                        outerRank = outerRank + 1;
                        idxOut{d} = outSubs(outerRank);
                    end
                end

                tsTarget  = requiredTS(idxOut{:});
                validMask = ~isnan(tsOrig) & ~isnan(sig);
                if sum(validMask) >= 2
                    sigResampled = interp1(tsOrig(validMask), sig(validMask), ...
                        tsTarget(:), 'linear', 'extrap');
                else
                    sigResampled = nan(numel(tsTarget), 1);
                end
                processed = processSignal(sigResampled, opt);
                if reuseAcrossUnits
                    for unitIdx = 1:sizeOut(reuseDim)
                        idxOut{reuseDim} = unitIdx;
                        newSlice(idxOut{:}) = processed;
                    end
                else
                    newSlice(idxOut{:}) = processed;
                end
            end
        end

    elseif isempty(otherDims)   %% TENSOR MODE — single time dimension

        tsOrig1d  = opt.varTS(:);
        tsTarget  = requiredTS(:);
        sig1d     = varTensor(:);
        valid1d   = ~isnan(tsOrig1d) & ~isnan(sig1d);
        if sum(valid1d) >= 2
            sigResampled = interp1(tsOrig1d(valid1d), sig1d(valid1d), tsTarget, 'linear', 'extrap');
        else
            sigResampled = nan(numel(tsTarget), 1);
        end
        newSlice(:) = processSignal(sigResampled, opt);

    else                    %% TENSOR MODE — multi-dimensional

        % Omit the unit dimension from the loop when one shared pass suffices.
        processDims = otherDims;
        if reuseAcrossUnits
            processDims(processDims == reuseDim) = [];
        end
        processSizes = sizeVar(processDims);
        nProcess     = prod(processSizes);

        for linIdx = 1:nProcess
            % Map linear index to the dimensions that still need processing.
            subProcess = ind2subVec(processSizes, linIdx);

            % Build full indexing cell (colon for timeDim, scalar for the rest)
            idxCell = cell(1, nDataDims);
            pCount  = 0;
            for d = 1:nDataDims
                if d == timeDim
                    idxCell{d} = ':';
                elseif reuseAcrossUnits && d == reuseDim
                    idxCell{d} = 1;
                else
                    pCount     = pCount + 1;
                    idxCell{d} = subProcess(pCount);
                end
            end

            tsOrig    = opt.varTS(idxCell{:});
            tsTarget  = requiredTS(idxCell{:});
            sig       = varTensor(idxCell{:});
            tsOrig    = tsOrig(:);
            tsTarget  = tsTarget(:);
            sig       = sig(:);
            validMask = ~isnan(tsOrig) & ~isnan(sig);
            if sum(validMask) >= 2
                sigResampled = interp1(tsOrig(validMask), sig(validMask), tsTarget, 'linear', 'extrap');
            else
                sigResampled = nan(numel(tsTarget), 1);
            end
            processed = processSignal(sigResampled, opt);
            if reuseAcrossUnits
                for unitIdx = 1:sizeOut(reuseDim)
                    idxCell{reuseDim} = unitIdx;
                    newSlice(idxCell{:}) = processed;
                end
            else
                newSlice(idxCell{:}) = processed;
            end
        end
    end

    %% 10 ── Concatenate new slice to form tensorOut ───────────────────────

    newSliceExpanded = reshape(newSlice, [sizeOut, 1]);

    if isFirstCall
        tensorOut = newSliceExpanded;
    else
        tensorOut = cat(nDimsIn, tensorIn, newSliceExpanded);
    end

    %% 11 ── Build infoOut ─────────────────────────────────────────────────

    infoOut                = opt.infoIn;
    fieldName              = sprintf('n%d', nSlices + 1);
    entry.var              = var;
    entry.requiredTS       = requiredTS;
    entry.varTS            = opt.varTS;
    entry.timeDim          = timeDim;
    entry.unitDim          = opt.unitDim;
    entry.reuseAcrossUnits = reuseAcrossUnits;
    entry.der              = opt.der;
    entry.lag              = opt.lag;
    entry.filter           = opt.filter;
    entry.varName          = opt.varName;
    infoOut.(fieldName)    = entry;

end   % addVariable2Tensor

%% =========================================================================
%  LOCAL HELPER FUNCTIONS
%% =========================================================================

function sig = processSignal(sig, opt)
%PROCESSSIGNAL  Apply derivative → lag → filter to a column vector.

    % ---- Derivative (iterated central differences) ----
    for k = 1:opt.der                      %#ok<FXUP>
        sig = central_diff(sig);
    end

    % ---- Lag ----
    if opt.lag ~= 0
        n   = numel(sig);
        lag = opt.lag;

        if lag > 0
            % Delay: shift right; pad left with mean of the retained samples
            nValid = max(n - lag, 0);
            if nValid == 0
                sig = mean(sig) * ones(n, 1);
            else
                padVal = mean(sig(1:nValid));
                sig    = [padVal * ones(lag, 1); sig(1:nValid)];
            end
        else
            % Anticipate: shift left; pad right with mean of the retained samples
            lag    = abs(lag);
            nValid = max(n - lag, 0);
            if nValid == 0
                sig = mean(sig) * ones(n, 1);
            else
                padVal = mean(sig(lag + 1:end));
                sig    = [sig(lag + 1:end); padVal * ones(lag, 1)];
            end
        end
    end

    % ---- Filter ----
    if ~isempty(opt.filter)
        sig = conv(sig(:)', opt.filter(:)', 'same')';
    end

end

% -------------------------------------------------------------------------

function Tout = interpAlongDim(Tin, dim, targetLen)
%INTERPALONGDIM  Linearly resample a tensor along one dimension.
%   Uses a normalised [0,1] time axis.  For TC data, replace with a
%   timestamp-aware version using obj.TCTS.

    sz      = size(Tin);
    origLen = sz(dim);
    nD      = numel(sz);

    % Permute target dimension to front for vectorised interp1
    permOrd  = [dim, 1:dim-1, dim+1:nD];
    Tp       = permute(Tin, permOrd);
    szP      = sz(permOrd);
    nRest    = prod(szP(2:end));   % product of empty = 1 when nD == 1

    tOrig    = linspace(0, 1, origLen);
    tNew     = linspace(0, 1, targetLen);

    Tm       = reshape(Tp, origLen, nRest);
    TmNew    = interp1(tOrig, Tm, tNew, 'linear');   % [targetLen × nRest]

    szPNew    = szP;
    szPNew(1) = targetLen;
    Tout      = ipermute(reshape(TmNew, szPNew), permOrd);

end

% -------------------------------------------------------------------------

function tf = isRepeatedAlongDim(data, dim)
%ISREPEATEDALONGDIM  True when every slice along dim equals slice one.

    nSlices = size(data, dim);
    if nSlices <= 1
        tf = true;
        return
    end

    nDims = max(ndims(data), dim);
    idx   = repmat({':'}, 1, nDims);
    idx{dim} = 1;
    reference = data(idx{:});

    tf = true;
    for sliceIdx = 2:nSlices
        idx{dim} = sliceIdx;
        if ~isequaln(data(idx{:}), reference)
            tf = false;
            return
        end
    end

end

% -------------------------------------------------------------------------

function sub = ind2subVec(sz, linIdx)
%IND2SUBVEC  Convert a column-major linear index to a subscript vector.

    nD  = numel(sz);
    sub = zeros(1, nD);
    r   = linIdx - 1;
    for d = 1:nD
        sub(d) = mod(r, sz(d)) + 1;
        r      = floor(r / sz(d));
    end

end
