function obj = refineNeuralData(obj, refType, refOpt, fixedFlag)
% REFINENEURALDATA  Apply normalization, subtraction, or smoothing to neural tensors.
%
% -------------------------------------------------------------------------
% SYNTAX
%   obj = obj.refineNeuralData(refType, refOpt, fixedFlag)
%
% -------------------------------------------------------------------------
% INPUTS
%   refType    - Operation to apply (string):
%                  'normalize' | 'subtract' | 'smooth'
%
%   refOpt     - Struct with operation-specific options (see below).
%
%   fixedFlag  - true  → operate on TensorWithFixBin / TensorWithFixBinCondAvg
%                false → operate on TensorWithVariableBin / TensorWithVariableBinCondAvg
%
% -------------------------------------------------------------------------
% OPERATION OPTIONS
%
%   refType = 'normalize'
%     refOpt.normType  - 'z-score' | 'range' | 'soft-norm'
%
%     normType = 'z-score'
%       Z-scores each neuron individually using its own mean and std.
%       The trial tensor and the condAvg tensor are z-scored independently
%       (each using their own mean/std).
%
%     normType = 'range'
%       refOpt.normRange  - [mi ma]  target range
%       Scales each neuron's activity to [mi, ma] (min/max computed
%       independently for the trial tensor and for the condAvg tensor).
%
%     normType = 'soft-norm'
%       refOpt.normConstant           - number or 'mean'
%                                       If 'mean', the mean firing rate of
%                                       the neuron is used as constant.
%       refOpt.normConstantMultiplier - (optional) scalar multiplier applied
%                                       to the constant (default: 1).
%       Formula: B = A / (range(A) + constant * multiplier)
%       Range, constant, and multiplier are computed independently for the
%       trial tensor and for the condAvg tensor.
%
%   refType = 'subtract'
%     refOpt.value  - 'mean' | 'median' | scalar | (nNeu x 1) numeric vector
%
%       'mean'   : subtracts each neuron's mean, computed independently
%                  from the trial tensor and from the condAvg tensor.
%       'median' : same as 'mean' but using the median.
%       scalar   : subtracts the same constant from every neuron in both
%                  tensors.
%       vector   : (nNeu x 1) pre-computed per-neuron values. Entry n is
%                  subtracted from all entries of neuron n in both the
%                  trial tensor and the condAvg tensor. Use this form to
%                  apply externally computed statistics (e.g. a mean
%                  derived from a separate training set) consistently
%                  across multiple objects.
%
%   refType = 'smooth'
%     refOpt.smoothType    - 'movmean' | 'movmedian' | 'gaussian'
%     refOpt.smoothWindow  - integer window length (samples)
%       Calls MATLAB's smoothdata along the time dimension.
%       Applied per neuron per trial independently via smoothdata's
%       element-wise behaviour along the time axis.
%       NaN-padded entries are ignored ('omitnan').
%
% -------------------------------------------------------------------------
% OUTPUT
%   obj  - Updated RecordedData object (tensors replaced in place).
%
% -------------------------------------------------------------------------
% NOTES
%   - All per-neuron statistics are computed ignoring NaN values (used as
%     padding for variable-length trial data).
%   - The trial tensor and the condition-averaged tensor are treated
%     independently: statistics are derived separately from each tensor and
%     applied to it. This preserves the internal consistency of each
%     representation.
%
% -------------------------------------------------------------------------
% EXAMPLES
%   % Z-score fixed-bin tensors
%   opt.normType = 'z-score';
%   obj = obj.refineNeuralData('normalize', opt, true);
%
%   % Soft-normalize variable-bin tensors with mean constant × 0.5
%   opt.normType               = 'soft-norm';
%   opt.normConstant           = 'mean';
%   opt.normConstantMultiplier = 0.5;
%   obj = obj.refineNeuralData('normalize', opt, false);
%
%   % Subtract baseline mean (fixed-bin)
%   opt.value = 'mean';
%   obj = obj.refineNeuralData('subtract', opt, true);
%
%   % Subtract a pre-computed per-neuron mean vector (e.g. from training data)
%   %   meanVec is nNeu x 1 and was computed externally from a training set
%   opt.value = meanVec;
%   trainObj = trainObj.refineNeuralData('subtract', opt, true);
%   testObj  = testObj.refineNeuralData('subtract', opt, true);
%
%   % Gaussian smoothing (variable-bin)
%   opt.smoothType   = 'gaussian';
%   opt.smoothWindow = 5;
%   obj = obj.refineNeuralData('smooth', opt, false);
%
% See also: smoothdata, prepareTensorWithFixBin, prepareTensorWithVariableBin

% ---- Select tensors based on fixedFlag -----------------------------------
if fixedFlag
    T    = obj.TensorWithFixBin;         % nNeu x nCond x nTrials x nTime
    Tavg = obj.TensorWithFixBinCondAvg;  % nNeu x nCond x nTime
    Info = obj.TensorWithFixBinInfo;
else
    T    = obj.TensorWithVariableBin;         % nNeu x nCond x nTrials x nTime
    Tavg = obj.TensorWithVariableBinCondAvg;  % nNeu x nCond x nTime
    Info = obj.TensorWithVariableBinInfo;
end

if isempty(T) && isempty(Tavg)
    warning('refineNeuralData: selected tensors are empty. Run prepareTensor first.');
    return
end

nNeu       = size(T,    1);
timeDimT   = ndims(T);     % 4 for trial tensor
timeDimAvg = ndims(Tavg);  % 3 for condAvg tensor

% ---- Check for repeated preprocessing step ------------------------------
typeKey = lower(refType);
if isfield(Info, typeKey)
    switch typeKey
        case 'normalize'
            warning('refineNeuralData: Normalization already performed. Another normalization step could lead to unwanted results.');
        case 'subtract'
            warning('refineNeuralData: Subtraction already performed. Another subtraction step could lead to unwanted results.');
        case 'smooth'
            warning('refineNeuralData: Smoothing already performed. Another smoothing step could lead to unwanted results.');
    end
end

% ---- Dispatch ------------------------------------------------------------
switch typeKey

    % ======================================================================
    case 'normalize'
        normType = refOpt.normType;

        for n = 1 : nNeu

            % Slices for this neuron (kept as arrays, not squeezed,
            % so assignment back preserves dimensions)
            neuT   = T(n, :, :, :);
            neuAvg = Tavg(n, :, :);

            switch lower(normType)

                % ----------------------------------------------------------
                case 'z-score'
                    % --- trial tensor ---
                    vT  = neuT(:);   vT  = vT(~isnan(vT));
                    muT = mean(vT);  sdT = std(vT);
                    if sdT == 0,  sdT = 1;  end
                    T(n,:,:,:) = (neuT - muT) ./ sdT;

                    % --- condAvg tensor ---
                    vA   = neuAvg(:); vA  = vA(~isnan(vA));
                    muA  = mean(vA);  sdA = std(vA);
                    if sdA == 0, sdA = 1; end
                    Tavg(n,:,:) = (neuAvg - muA) ./ sdA;

                % ----------------------------------------------------------
                case 'range'
                    mi = refOpt.normRange(1);
                    ma = refOpt.normRange(2);

                    % --- trial tensor ---
                    vT   = neuT(:);  vT = vT(~isnan(vT));
                    minT = min(vT);  maxT = max(vT);
                    rngT = maxT - minT;
                    if rngT == 0, rngT = 1; end
                    T(n,:,:,:) = mi + (neuT - minT) ./ rngT .* (ma - mi);

                    % --- condAvg tensor ---
                    vA   = neuAvg(:); vA = vA(~isnan(vA));
                    minA = min(vA);   maxA = max(vA);
                    rngA = maxA - minA;
                    if rngA == 0, rngA = 1; end
                    Tavg(n,:,:) = mi + (neuAvg - minA) ./ rngA .* (ma - mi);

                % ----------------------------------------------------------
                case 'soft-norm'
                    % Default multiplier
                    if isfield(refOpt, 'normConstantMultiplier')
                        multiplier = refOpt.normConstantMultiplier;
                    else
                        multiplier = 1;
                    end

                    % --- trial tensor ---
                    vT   = neuT(:); vT = vT(~isnan(vT));
                    rngT = max(vT) - min(vT);
                    if ischar(refOpt.normConstant) || isstring(refOpt.normConstant)
                        cT = mean(vT);
                    else
                        cT = refOpt.normConstant;
                    end
                    T(n,:,:,:) = neuT ./ (rngT + cT * multiplier);

                    % --- condAvg tensor ---
                    vA   = neuAvg(:); vA = vA(~isnan(vA));
                    rngA = max(vA) - min(vA);
                    if ischar(refOpt.normConstant) || isstring(refOpt.normConstant)
                        cA = mean(vA);
                    else
                        cA = refOpt.normConstant;
                    end
                    Tavg(n,:,:) = neuAvg ./ (rngA + cA * multiplier);

                otherwise
                    error('refineNeuralData: unknown normType ''%s''. Choose ''z-score'', ''range'', or ''soft-norm''.', normType);
            end
        end

    % ======================================================================
    case 'subtract'
        val = refOpt.value;

        for n = 1 : nNeu
            neuT   = T(n, :, :, :);
            neuAvg = Tavg(n, :, :);

            if ischar(val) || isstring(val)
                switch lower(val)
                    case 'mean'
                        vT  = neuT(:);   vT  = vT(~isnan(vT));   subT  = mean(vT);
                        vA  = neuAvg(:); vA  = vA(~isnan(vA));   subA  = mean(vA);
                    case 'median'
                        vT  = neuT(:);   vT  = vT(~isnan(vT));   subT  = median(vT);
                        vA  = neuAvg(:); vA  = vA(~isnan(vA));   subA  = median(vA);
                    otherwise
                        error('refineNeuralData: unknown subtract value ''%s''. Choose ''mean'', ''median'', or a scalar.', val);
                end
            else
                % Scalar constant or per-neuron vector
                if isscalar(val)
                    subT = val;
                    subA = val;
                else
                    subT = val(n);
                    subA = val(n);
                end
            end

            T(n,:,:,:) = neuT   - subT;
            Tavg(n,:,:) = neuAvg - subA;
        end

    % ======================================================================
    case 'smooth'
        switch lower(refOpt.smoothType)
            case 'movmean'
                method = 'movmean';
            case 'movmedian'
                method = 'movmedian';
            case 'gaussian'
                method = 'gaussian';
            otherwise
                error('refineNeuralData: unknown smoothType ''%s''. Choose ''movmean'', ''movmedian'', or ''gaussian''.', refOpt.smoothType);
        end

        smoothWindow = refOpt.smoothWindow;

        % smoothdata operates along the time dimension on the full tensor.
        % Each (neuron, condition, trial) fibre is smoothed independently,
        % matching the "per neuron per trial" requirement.
        % 'omitnan' prevents NaN padding from contaminating valid bins.
        if ~isempty(T)
            T    = smoothdata(T,    timeDimT,   method, smoothWindow, 'omitnan');
        end
        if ~isempty(Tavg)
            Tavg = smoothdata(Tavg, timeDimAvg, method, smoothWindow, 'omitnan');
        end

    otherwise
        error('refineNeuralData: unknown refType ''%s''. Choose ''normalize'', ''subtract'', or ''smooth''.', refType);
end

% ---- Record preprocessing step in Info ----------------------------------
Info.(typeKey) = refOpt;

% ---- Write back to object ------------------------------------------------
if fixedFlag
    obj.TensorWithFixBin        = T;
    obj.TensorWithFixBinCondAvg = Tavg;
    obj.TensorWithFixBinInfo    = Info;
else
    obj.TensorWithVariableBin        = T;
    obj.TensorWithVariableBinCondAvg = Tavg;
    obj.TensorWithVariableBinInfo    = Info;
end

end
