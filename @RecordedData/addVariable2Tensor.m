function [tensorOut, infoOut] = addVariable2Tensor(obj, tensorIn, opt)
%ADDVARIABLE2TENSOR  Append a preprocessed variable as a new slice in the last dimension.
%
% [tensorOut, infoOut] = obj.addVariable2Tensor(tensorIn, opt)
%
% Designed to be called iteratively to build a regressor tensor.  Each call
% appends one variable (after optional derivative / lag / filter processing)
% as a new slice along the last dimension of tensorIn.
%
% INPUTS
%   tensorIn  – Current regressor tensor  [d1 × … × dN × K].
%               Pass [] on the first call; tensorOut will be [d1 × … × dN × 1].
%   opt       – Options struct with the following fields:
%
%     .timeDim  (required) Positive integer: index of the time/signal
%                 dimension inside the variable tensor (not tensorIn).
%                 Example: if the variable is neu×cond×trial×time, set
%                 timeDim = 4 to process each time trace independently.
%
%     .var        Variable source.  Either:
%                   • {'FR','FixBin'}      – obj.TensorWithFixBin 
%                   • {'FR','VariableBin'} – obj.TensorWithVariableBin 
%                   • {'TC','xCursor'}     – x coordinate from obj.TC
%                   • {'TC','yCursor'}     – y coordinate from obj.TC 
%                   • {'KIN',…}            – not yet implemented
%                   • {'EY',…}             – not yet implemented
%                   • numeric tensor       – if tensorIn is not empty, .var must have exactly ndims(tensorIn)-1 dimensions
%                   • cell array           – {nOuter}{1×nTrials}[1×T]. NOTE: opt.varTS must be a cell array with the same outer shape, see below.
%                    
%     .varTS      (required) Numeric tensor with the same shape OR cell array as the variable
%                 extracted from opt.var.  Each element holds the
%                 timestamp (in any consistent unit).
%
%     .requiredTS (required) Numeric tensor whose non-timeDim dimensions
%                 match the variable tensor.  Its size along timeDim defines
%                 the output sample count and its values give the target
%                 timestamps for resampling.  Different non-signal slices may
%                 carry different timestamp sequences (e.g. variable trial
%                 durations or mixed sample rates).  Resampling is performed
%                 per slice as:
%                   interp1(varTS_slice, sig_slice, requiredTS_slice,
%                           ''linear'', ''extrap'')
%
%     .unitDim  Optional positive integer: index of the unit/neuron
%                 dimension inside a numeric variable tensor (not tensorIn).
%                 Required when reusing a numeric tensor across units.
%
%     .reuseAcrossUnits  Optional logical flag. For simultaneous cell inputs,
%                 true processes only the first neuron and copies the result
%                 across the fixed neuron axis. For numeric tensors, explicit
%                 true enables the same optimization using .unitDim. The
%                 default enables it for simultaneous cell inputs and
%                 disables it for FR and numeric tensors. Set it for a numeric
%                 tensor only when its values and .varTS are identical across
%                 units.
%
%     .der        Derivative order (non-negative integer, default 0).
%                 Each order applies central_diff once (output length is
%                 preserved).
%
%     .lag        Lag in samples (integer, default 0).
%                 Positive = delay (signal shifted right).
%                 Negative = anticipate (signal shifted left).
%                 Boundary samples are padded with the mean of the valid
%                 (non-padded) portion.
%
%     .filter     Convolution kernel vector for conv(…,'same') (default []).
%                 If empty, no filtering is applied.  A warning is issued if
%                 sum(filter) ≠ 1.
%
%     .varName    Optional label string for bookkeeping (default '').
%
%     .infoIn     Info struct returned by a previous call (default struct()).
%
% OUTPUTS
%   tensorOut – [d1 × … × dN × K+1].
%   infoOut   – Copy of opt.infoIn with a new field 'n<K+1>' that stores
%               all opt settings used for this addition.
%               First variable → field 'n1', second → 'n2', etc.
%
% EXAMPLE
%   % Build a regressor tensor from firing rate (FixBin) and TC x-coordinate
%   opt1.timeDim = 4;
%   opt1.var       = {'FR','FixBin'};
%   opt1.varName   = 'FiringRate';
%   [T, info] = recordedData.addVariable2Tensor([], opt1);
%
%   opt2.timeDim = 4;
%   opt2.var       = {'TC','xCursor'};
%   opt2.der       = 1;          % first derivative
%   opt2.lag       = 2;          % delay by 2 bins
%   opt2.varName   = 'xVelocity';
%   opt2.infoIn    = info;
%   [T, info] = recordedData.addVariable2Tensor(T, opt2);

    %% 1 ── Validate and apply defaults ────────────────────────────────────

    if ~isstruct(opt) || ~isfield(opt, 'timeDim') || isempty(opt.timeDim)
        error('RecordedData:addVariable2Tensor:missingtimeDim', ...
            'opt.timeDim is required and cannot be empty.');
    end
    timeDim = opt.timeDim;
    if ~isscalar(timeDim) || timeDim < 1 || timeDim ~= floor(timeDim)
        error('RecordedData:addVariable2Tensor:invalidtimeDim', ...
            'opt.timeDim must be a positive integer scalar.');
    end

    if ~isfield(opt, 'der') || isempty(opt.der)
        opt.der = 0;
    end
    if ~isfield(opt, 'lag') || isempty(opt.lag)
        opt.lag = 0;
    end
    if ~isfield(opt, 'filter')
        opt.filter = [];
    end
    if ~isfield(opt, 'varName')
        opt.varName = '';
    end
    if ~isfield(opt, 'infoIn') || isempty(opt.infoIn)
        opt.infoIn = struct();
    end
    if ~isfield(opt, 'unitDim') || isempty(opt.unitDim)
        opt.unitDim = [];
    elseif ~isnumeric(opt.unitDim) || ~isscalar(opt.unitDim) || ...
            ~isfinite(opt.unitDim) || opt.unitDim < 1 || ...
            opt.unitDim ~= floor(opt.unitDim)
        error('RecordedData:addVariable2Tensor:invalidUnitDim', ...
            'opt.unitDim must be a positive integer scalar.');
    end
    reuseOptionSet = isfield(opt, 'reuseAcrossUnits') && ...
        ~isempty(opt.reuseAcrossUnits);
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
    if ~isfield(opt, 'var')
        error('RecordedData:addVariable2Tensor:missingVar', ...
            'opt.var is required.');
    end

    if ~isfield(opt, 'varTS') || isempty(opt.varTS)
        error('RecordedData:addVariable2Tensor:missingVarTS', ...
            'opt.varTS is required and cannot be empty.');
    end
    if ~isfield(opt, 'requiredTS') || isempty(opt.requiredTS)
        error('RecordedData:addVariable2Tensor:missingRequiredTS', ...
            'opt.requiredTS is required and cannot be empty.');
    end

    if ~isempty(opt.filter) && abs(sum(opt.filter) - 1) > 1e-6
        warning('RecordedData:addVariable2Tensor:filterPower', ...
            'opt.filter sums to %.6f (not 1). Proceeding anyway.', ...
            sum(opt.filter));
    end

    %% 2 ── Determine current tensor shape and slice count ─────────────────
    %   numel(size(...)) is used throughout because ndims() drops trailing
    %   singleton dimensions (which occur when K = 1).

    isFirstCall = isempty(tensorIn);

    if ~isFirstCall
        sizeIn    = size(tensorIn);
        nDimsIn   = numel(sizeIn);   % total dims, including the K last dim
        nDataDims = nDimsIn - 1;
        nSlices   = sizeIn(end);     % current K
    end

    %% 3 ── Extract variable tensor from opt.var ───────────────────────────
    %   Two operating modes:
    %     Tensor mode – opt.var and opt.varTS are numeric tensors.
    %     Cell mode   – opt.var and opt.varTS are {nUnit x nCond}{1×nTrials}[1×T]
    %                   cell arrays; no large tensor is pre-allocated.

    isCellMode = false;
    varSource  = '';

    if iscell(opt.var) && ~isempty(opt.var) && (ischar(opt.var{1}) || isstring(opt.var{1}))
        varSource = upper(opt.var{1});
        switch varSource

            case 'FR'
                if numel(opt.var) < 2
                    error('RecordedData:addVariable2Tensor:missingSpec', ...
                        'Specify {''FR'',''FixBin''} or {''FR'',''VariableBin''}.');
                end
                switch lower(opt.var{2})
                    case 'fixbin'
                        if isempty(obj.TensorWithFixBin)
                            error('RecordedData:addVariable2Tensor:emptyData', ...
                                'obj.TensorWithFixBin is empty.');
                        end
                        varTensor = obj.TensorWithFixBin;
                    case 'variablebin'
                        if isempty(obj.TensorWithVariableBin)
                            error('RecordedData:addVariable2Tensor:emptyData', ...
                                'obj.TensorWithVariableBin is empty.');
                        end
                        varTensor = obj.TensorWithVariableBin;
                    otherwise
                        error('RecordedData:addVariable2Tensor:unknownBinType', ...
                            'Unknown FR bin type "%s". Use ''FixBin'' or ''VariableBin''.', ...
                            opt.var{2});
                end

            case 'TC'
                if numel(opt.var) < 2
                    error('RecordedData:addVariable2Tensor:missingSpec', ...
                        'Specify {''TC'',''xCursor''} or {''TC'',''yCursor''}.');
                end
                if isempty(obj.TC)
                    error('RecordedData:addVariable2Tensor:emptyData', ...
                        'obj.TC is empty.');
                end
                switch lower(opt.var{2})
                    case 'xcursor', coordIdx = 1; tcField = 'xCursor';
                    case 'ycursor', coordIdx = 2; tcField = 'yCursor';
                    otherwise
                        error('RecordedData:addVariable2Tensor:unknownCoord', ...
                            'Unknown TC coordinate "%s". Use ''xCursor'' or ''yCursor''.', ...
                            opt.var{2});
                end

                if iscell(obj.TC)
                    % TC is {nNeu × nCond}{1 × nTrials} struct cell array.
                    % Use cell mode – no large tensor pre-allocation needed.
                    isCellMode = true;
                    varCell = cellfun( ...
                        @(condCell) cellfun(@(s) s.(tcField)(:)', condCell, ...
                            'UniformOutput', false), ...
                        obj.TC, 'UniformOutput', false);
                    % opt.varTS is expected to be TCTS with a matching cell layout
                else
                    % TC is a numeric tensor with shape [..., time_tc, 2]
                    tcSz    = size(obj.TC);
                    tcNd    = numel(tcSz);
                    idxCell = repmat({':'}, 1, tcNd);
                    idxCell{end} = coordIdx;
                    varTensor = obj.TC(idxCell{:});
                end

            case {'KIN', 'EY'}
                error('RecordedData:addVariable2Tensor:notImplemented', ...
                    '"%s" support is not yet implemented.', varSource);

            otherwise
                error('RecordedData:addVariable2Tensor:unknownSource', ...
                    'Unknown source "%s". Use ''FR'', ''TC'', ''KIN'' or ''EY''.', ...
                    varSource);
        end

    elseif iscell(opt.var)
        % Direct cell array data: {nOuter...}{1×nTrials}[1×T]
        isCellMode = true;
        varCell    = opt.var;

    elseif isnumeric(opt.var)
        varTensor = opt.var;

    else
        error('RecordedData:addVariable2Tensor:invalidVar', ...
            'opt.var must be a text spec, cell array, or numeric tensor.');
    end

    %% 4 ── Resolve shape on first call (empty tensorIn) ───────────────────

    if isFirstCall
        if isCellMode
            % In cell mode the output shape is governed by requiredTS
            nDataDims = numel(size(opt.requiredTS));
            nSlices   = 0;
            nDimsIn   = nDataDims + 1;
            sizeIn    = [size(opt.requiredTS), 1];   % placeholder
        else
            sizeVar   = size(varTensor);
            nDataDims = numel(sizeVar);
            nSlices   = 0;
            nDimsIn   = nDataDims + 1;
            sizeIn    = [sizeVar, 1];   % placeholder: K = 1 after first addition
        end
    end

    %% 5 ── Validate variable against tensorIn ────────────────────────────

    if isCellMode
        % Cell mode: reference shape is requiredTS; validate varTS cell layout
        sizeVar     = size(opt.requiredTS);
        cellOuterSz = size(varCell);
        if ~iscell(opt.varTS) || ~isequal(size(opt.varTS), cellOuterSz)
            error('RecordedData:addVariable2Tensor:varTSshapeMismatch', ...
                ['In cell mode opt.varTS must be a cell array with outer shape ' ...
                 '[%s] matching opt.var.'], num2str(cellOuterSz));
        end
    else
        sizeVar  = size(varTensor);
        nDimsVar = numel(sizeVar);

        if nDimsVar ~= nDataDims
            error('RecordedData:addVariable2Tensor:dimMismatch', ...
                ['opt.var has %d dimension(s); expected %d ' ...
                 '(one fewer than tensorIn).'], nDimsVar, nDataDims);
        end

        for d = 1:nDataDims
            if d ~= timeDim && sizeVar(d) ~= sizeIn(d)
                error('RecordedData:addVariable2Tensor:sizeMismatch', ...
                    ['opt.var size at dimension %d is %d, but tensorIn ' ...
                     'size at dimension %d is %d.'], d, sizeVar(d), d, sizeIn(d));
            end
        end

        if ~isequal(size(opt.varTS), sizeVar)
            error('RecordedData:addVariable2Tensor:varTSshapeMismatch', ...
                'opt.varTS must have the same shape as the variable tensor.');
        end
    end

    if timeDim > nDataDims
        error('RecordedData:addVariable2Tensor:invalidtimeDim', ...
            'opt.timeDim (%d) exceeds the number of data dimensions (%d).', ...
            timeDim, nDataDims);
    end
    if ~isempty(opt.unitDim) && opt.unitDim > nDataDims
        error('RecordedData:addVariable2Tensor:invalidUnitDim', ...
            'opt.unitDim (%d) exceeds the number of data dimensions (%d).', ...
            opt.unitDim, nDataDims);
    end

    % Cell/behavioral inputs use the fixed first outer dimension as the
    % neuron axis. Numeric tensors require an explicit unitDim opt-in.
    if isCellMode
        if reuseOptionSet && strcmpi(varSource, 'FR') && opt.reuseAcrossUnits
            error('RecordedData:addVariable2Tensor:invalidReuseForFR', ...
                'FR variables are always processed independently per neuron.');
        end
        if reuseOptionSet
            reuseAcrossUnits = obj.simultaneousRecording && opt.reuseAcrossUnits;
        else
            reuseAcrossUnits = obj.simultaneousRecording;
        end
        reuseDim = 1;
    elseif strcmpi(varSource, 'FR')
        if reuseOptionSet && opt.reuseAcrossUnits
            error('RecordedData:addVariable2Tensor:invalidReuseForFR', ...
                'FR variables are always processed independently per neuron.');
        end
        reuseAcrossUnits = false;
        reuseDim = [];
    else
        if reuseOptionSet
            reuseAcrossUnits = opt.reuseAcrossUnits;
        else
            reuseAcrossUnits = false;
        end
        reuseDim = opt.unitDim;
    end
    if reuseAcrossUnits
        if isempty(reuseDim)
            error('RecordedData:addVariable2Tensor:missingUnitDim', ...
                'opt.unitDim is required when reusing a numeric tensor across units.');
        end
        if reuseDim == timeDim
            error('RecordedData:addVariable2Tensor:unitDimIsTimeDim', ...
                'opt.unitDim must identify a non-time dimension.');
        end
    end

    sizeReqTS = size(opt.requiredTS);
    if numel(sizeReqTS) ~= nDataDims
        error('RecordedData:addVariable2Tensor:requiredTSdimMismatch', ...
            'opt.requiredTS must have %d dimension(s) to match the variable tensor.', ...
            nDataDims);
    end
    for d = 1:nDataDims
        if d ~= timeDim && sizeReqTS(d) ~= sizeVar(d)
            error('RecordedData:addVariable2Tensor:requiredTSsizeMismatch', ...
                ['opt.requiredTS size at dimension %d is %d, but the variable ' ...
                 'tensor size at that dimension is %d.'], d, sizeReqTS(d), sizeVar(d));
        end
    end
    if ~isFirstCall && sizeReqTS(timeDim) ~= sizeIn(timeDim)
        error('RecordedData:addVariable2Tensor:requiredTSlengthMismatch', ...
            ['opt.requiredTS has %d sample(s) along timeDim, but tensorIn ' ...
             'has %d. All added variables must share the same target sample count.'], ...
            sizeReqTS(timeDim), sizeIn(timeDim));
    end
    if reuseAcrossUnits && ~isRepeatedAlongDim(opt.requiredTS, reuseDim)
        % A shared source still needs a shared target grid. Otherwise each
        % unit must be resampled independently to preserve the old result.
        reuseAcrossUnits = false;
    end

    %% 6 ── Determine output length and shape from requiredTS ──────────────

    targetLen = sizeReqTS(timeDim);
    sizeOut   = sizeVar;
    sizeOut(timeDim) = targetLen;

    %% 7 ── Iterate over non-signal indices; resample and process each 1-D signal ──

    newSlice  = nan(sizeOut);
    otherDims = [1:timeDim-1, timeDim+1:nDataDims];

    if isCellMode
        % Cell mode: varCell{outer...}{t}(1:T) — no large tensor pre-allocated.
        % Tensor dim (timeDim-1) maps to the trial index t within each cell.
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
                    elseif d == timeDim - 1   % trial dimension
                        idxOut{d} = t;
                    else
                        outerRank = outerRank + 1;
                        idxOut{d} = outSubs(outerRank);
                    end
                end

                tsTarget  = opt.requiredTS(idxOut{:});
                validMask = ~isnan(tsOrig) & ~isnan(sig);
                if sum(validMask) >= 2
                    sigResampled = interp1(tsOrig(validMask), sig(validMask), tsTarget(:), 'linear', 'extrap');
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

    elseif isempty(otherDims)   %% TENSOR MODE! but only one dimension (timeDim) exists

        tsOrig1d  = opt.varTS(:);
        tsTarget  = opt.requiredTS(:);
        sig1d     = varTensor(:);
        valid1d   = ~isnan(tsOrig1d) & ~isnan(sig1d);
        if sum(valid1d) >= 2
            sigResampled = interp1(tsOrig1d(valid1d), sig1d(valid1d), tsTarget, 'linear', 'extrap');
        else
            sigResampled = nan(numel(tsTarget), 1);
        end
        newSlice(:)  = processSignal(sigResampled, opt);

    else                    %% TENSOR MODE! but multi-D
        % Keep the full output shape, but omit the unit dimension from the
        % expensive processing loop when one shared representative suffices.
        processDims = otherDims;
        if reuseAcrossUnits
            processDims(processDims == reuseDim) = [];
        end
        processSizes = sizeVar(processDims);
        nProcess     = prod(processSizes);


        for linIdx = 1:nProcess
            % Map linear index to the dimensions that still need processing.
            subProcess = ind2subVec(processSizes, linIdx);

            % Build full indexing cell (colon for timeDim, scalar for rest)
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
            tsTarget  = opt.requiredTS(idxCell{:});
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

    %% 8 ── Concatenate new slice to form tensorOut ────────────────────────

    newSliceExpanded = reshape(newSlice, [sizeOut, 1]);

    if isFirstCall
        tensorOut = newSliceExpanded;
    else
        tensorOut = cat(nDimsIn, tensorIn, newSliceExpanded);
    end

    %% 9 ── Build infoOut ──────────────────────────────────────────────────

    infoOut             = opt.infoIn;
    fieldName           = sprintf('n%d', nSlices + 1);
    entry.timeDim     = opt.timeDim;
    entry.var           = opt.var;
    entry.varTS         = opt.varTS;
    entry.requiredTS    = opt.requiredTS;
    entry.unitDim       = opt.unitDim;
    entry.reuseAcrossUnits = reuseAcrossUnits;
    entry.der           = opt.der;
    entry.lag           = opt.lag;
    entry.filter        = opt.filter;
    entry.varName       = opt.varName;
    infoOut.(fieldName) = entry;

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
