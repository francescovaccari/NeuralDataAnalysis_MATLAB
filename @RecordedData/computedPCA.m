function [dPCAmodel, dPCAprojTrain, dPCAprojTest] = computedPCA(obj, dPCAOpt, processingOpt, trainOpt, testOpt, visOpt)
% COMPUTEDPCA  Fit dPCA on training data and project train/test sets.
%
% -------------------------------------------------------------------------
% SYNTAX
%   [dPCAmodel, dPCAprojTrain, dPCAprojTest] = obj.computedPCA( ...
%       dPCAOpt, processingOpt, trainOpt, testOpt)
%
%   [dPCAmodel, dPCAprojTrain, dPCAprojTest] = obj.computedPCA( ...
%       dPCAOpt, processingOpt, trainOpt, testOpt, visOpt)
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS
%   dPCAOpt        - Struct controlling the dPCA decomposition:
%     .combinedParams          - cell array of cell arrays specifying which
%                                marginalizations to group together (see dpca.m).
%                                E.g. {{1,[1 3]},{2,[2 3]},{3},{[1 2],[1 2 3]}}
%     .numComps                - number of dPCA components to extract (scalar)
%     .simultaneousRecording   - logical; true if all neurons were recorded
%                                simultaneously (same trial counts), false for
%                                sequentially recorded neurons
%     .procedure               - 'simple'               → standard dPCA (no regularization)
%                                'regularized'           → dPCA with optimal lambda and
%                                                          noise-covariance (dpca_optimizeLambda +
%                                                          dpca_getNoiseCovariance)
%                                'regularized&noiseCovar'→ same as 'regularized' but also
%                                                          estimates signal variance in
%                                                          dpca_explainedVariance (passing
%                                                          Cnoise and numOfTrials)
%
%   processingOpt  - Struct controlling neuron selection and binning:
%     .neus2consider       - (1 x N) neuron indices to include
%     .fixedFlag           - true  → fixed-bin pipeline
%                            false → variable-bin pipeline
%     if fixedFlag = true:
%       .binWidth           - (1 x 2) [binSize, stepSize] in seconds.
%                             Set stepSize = 0 for non-overlapping bins.
%     if fixedFlag = false:
%       .expectedBinWidth        - expected bin width (s) for prepareTensorWithVariableBin
%     .refType  (optional)  - 'normalize' or 'smooth'. Applied to both train and
%                             test objects before mean-centering. 'subtract' is
%                             forbidden here (use the built-in mean-centering step).
%     .refOpt   (optional)  - Struct passed verbatim to refineNeuralData as refOpt.
%                             Required when .refType is provided.
%
%   trainOpt       - Struct defining the training set:
%     .conds2consider4train - condition indices for the training tensor
%     .trainParams          - (nConds x nParams) parameter matrix
%     .trainWindow          - analysis window (same format as
%                             prepareTensorWithFixBin / prepareTensorWithVariableBin)
%
%   testOpt        - Struct defining the test set:
%     .conds2consider4test  - condition indices for the test tensor
%     .testParams           - (nConds x nParams) parameter matrix
%     .testWindow           - analysis window (same format as trainWindow)
%
% -------------------------------------------------------------------------
% OPTIONAL INPUTS
%   visOpt         - Struct controlling output visualisation (all optional):

%     .marginalizationNames   - cell array of marginalization names
%     .marginalizationColours - (nMarg x 3) colour matrix
%     .legendSubplot          - subplot index for the legend panel
%                               (default: 16)
%
% -------------------------------------------------------------------------
% OUTPUTS
%   dPCAmodel      - Struct with the dPCA model fitted on the training set:
%     .W            - decoder matrix  (neurons x numComps)
%     .V            - encoder matrix  (neurons x numComps)
%     .whichMarg    - (1 x numComps) marginalization assignment per component
%
%   dPCAprojTrain  - Struct with the training-set projection:
%     .explVar      - explained-variance struct from dpca_explainedVariance
%
%   dPCAprojTest   - Struct with the test-set projection:
%     .explVar      - explained-variance struct from dpca_explainedVariance
%
% -------------------------------------------------------------------------
% PROCESSING PIPELINE
%   The following steps are applied internally to both train and test sets,
%   in this order:
%     1. prepareTensor*         Bin spikes into the requested time window.
%     2. refineNeuralData       (optional) Normalize or smooth, controlled
%                               by processingOpt.refType / .refOpt.
%                               'subtract' is not allowed at this stage.
%     3. Mean-centering         Subtract the per-neuron global mean, derived
%                               from the concatenated train+test tensor.
%                               'simple'      → mean of the condAvg tensor.
%                               'regularized' → mean of the full trial tensor.
%     4. splitTensorWithParams  Reorganise the tensor by task parameters.
%     5. dpca / dpca_plot       Decompose and visualise.
%
% -------------------------------------------------------------------------
% EXAMPLE  (fixed-bin, two parameters)
%   dpOpt.combinedParams = {{1,[1 3]},{2,[2 3]},{3},{[1 2],[1 2 3]}};
%   dpOpt.numComps       = 20;
%
%   procOpt.fixedFlag      = true;
%   procOpt.neus2consider  = 1:96;
%   procOpt.binWidth       = [0.05 0.01];
%   % optional: soft-normalize before mean-centering
%   procOpt.refType        = 'normalize';
%   procOpt.refOpt.normType = 'soft-norm';
%   procOpt.refOpt.normConstant = 'mean';
%
%   trOpt.conds2consider4train = 1:8;
%   trOpt.trainParams          = conditions.Target6(1:8);
%   trOpt.trainWindow          = [-0.5 2 1.5];
%
%   teOpt.conds2consider4test  = 1:8;
%   teOpt.testParams           = conditions.Obstacle6(1:8);
%   teOpt.testWindow           = [-0.5 2 1.5];
%
%   vOpt.marginalizationNames = {'target','time','target x time'};
%   vOpt.legendSubplot        = 16;
%
%   [mdl, projTr, projTe] = obj.computedPCA(dpOpt, procOpt, trOpt, teOpt, vOpt);
%
% See also: dpca, dpca_explainedVariance, dpca_plot,
%           prepareTensorWithFixBin, prepareTensorWithVariableBin,
%           splitTensorWithParams

% ---- Unpack dPCAOpt ------------------------------------------------------
combinedParams         = dPCAOpt.combinedParams;
numComps               = dPCAOpt.numComps;
simultaneousRecording = dPCAOpt.simultaneousRecording;
procedure              = dPCAOpt.procedure;

% ---- Unpack processingOpt ------------------------------------------------
fixedFlag     = processingOpt.fixedFlag;
neus2consider = processingOpt.neus2consider;

if fixedFlag
    binWidth = processingOpt.binWidth;
else
    expectedBinWidth       = processingOpt.expectedBinWidth;
end

% ---- Unpack trainOpt / testOpt -------------------------------------------
conds2consider4train = trainOpt.conds2consider4train;
trainParams          = trainOpt.trainParams;
trainWindow          = trainOpt.trainWindow;

conds2consider4test  = testOpt.conds2consider4test;
testParams           = testOpt.testParams;
testWindow           = testOpt.testWindow;

% ---- Unpack visOpt (optional) -------------------------------------------
if nargin < 6 || isempty(visOpt)
    visOpt = struct();
end

nMarg = numel(combinedParams);

if isfield(visOpt, 'marginalizationColours') && ~isempty(visOpt.marginalizationColours)
    margColours = visOpt.marginalizationColours;
else
    margColours = jet(nMarg);
end

if isfield(visOpt, 'marginalizationNames') && ~isempty(visOpt.marginalizationNames)
    margNames = visOpt.marginalizationNames;
else
    margNames = arrayfun(@(i) ['MargNum' num2str(i)], 1:nMarg, 'UniformOutput', false);
end

legendSubplot = getOptField(visOpt, 'legendSubplot', 16);

% ---- Build tensors -------------------------------------------------------
if fixedFlag
    trainObj = obj.prepareTensorWithFixBin(trainWindow, binWidth, neus2consider, conds2consider4train);

    testObj = obj.prepareTensorWithFixBin(testWindow, binWidth, neus2consider, conds2consider4test);
else
    trainObj = obj.prepareTensorWithVariableBin(trainWindow, expectedBinWidth, neus2consider, conds2consider4train);
    
    testObj = obj.prepareTensorWithVariableBin(testWindow, expectedBinWidth, neus2consider, conds2consider4test);
end

% ---- Normalize neural data (optional)  ---------------------------------------
if isfield(processingOpt, 'refType') && ~isempty(processingOpt.refType)
    if strcmpi(processingOpt.refType, 'subtract')
        error(['computedPCA: ''subtract'' is not permitted in processingOpt.refType. ' ...
               'The method will automatically mean-centering the data. You can use ''normalize'' or ''smooth'' instead.']);
    end
    trainObj = trainObj.refineNeuralData(processingOpt.refType, processingOpt.refOpt, fixedFlag);
    testObj  = testObj.refineNeuralData(processingOpt.refType, processingOpt.refOpt, fixedFlag);
end

% ---- Mean centering -------------------------------------------------------
% Per-neuron mean is computed from the combined train+test tensor (ignoring NaN)
% and subtracted from both objects to keep centering consistent across sets.
% 'simple' derives the mean from the condAvg tensor; 'regularized' from the
% full trial tensor.
refOpt_mc = struct();
if fixedFlag
    switch procedure
        case 'simple'
            refOpt_mc.value = mean(trainObj.TensorWithFixBinCondAvg, [2, 3], 'omitnan');  % nNeu x 1
            trainObj = trainObj.refineNeuralData('subtract', refOpt_mc, fixedFlag);
            testObj  = testObj.refineNeuralData('subtract', refOpt_mc, fixedFlag);

        case {'regularized', 'regularized&noiseCovar'}
            refOpt_mc.value = mean(trainObj.TensorWithFixBin, [2, 3, 4], 'omitnan');  % nNeu x 1
            trainObj = trainObj.refineNeuralData('subtract', refOpt_mc, fixedFlag);
            testObj  = testObj.refineNeuralData('subtract', refOpt_mc, fixedFlag);
    end
else
    switch procedure
        case 'simple'
            refOpt_mc.value = mean(trainObj.TensorWithVariableBinCondAvg, [2, 3], 'omitnan');  % nNeu x 1
            trainObj = trainObj.refineNeuralData('subtract', refOpt_mc, fixedFlag);
            testObj  = testObj.refineNeuralData('subtract', refOpt_mc, fixedFlag);

        case {'regularized', 'regularized&noiseCovar'}
            refOpt_mc.value = mean(trainObj.TensorWithVariableBin, [2, 3, 4], 'omitnan');  % nNeu x 1
            trainObj = trainObj.refineNeuralData('subtract', refOpt_mc, fixedFlag);
            testObj  = testObj.refineNeuralData('subtract', refOpt_mc, fixedFlag);
    end
end

% ---- Split tensors by parameters ----------------------------------------
trainObj = trainObj.splitTensorWithParams(trainParams, fixedFlag);
testObj  = testObj.splitTensorWithParams(testParams, fixedFlag);

XtrainCondAvg = trainObj.TensorWithParamsCondAvg;
XtestCondAvg  = testObj.TensorWithParamsCondAvg;
trainInfo     = trainObj.TensorWithParamsInfo;
testInfo      = testObj.TensorWithParamsInfo;

nTimeTrain = size(XtrainCondAvg, ndims(XtrainCondAvg));
nTimeTest  = size(XtestCondAvg,  ndims(XtestCondAvg));

% ---- Build time axes and event markers -----------------------------------
if fixedFlag
    [timeTrain, timeEventsTrain] = buildTimeOptions(trainInfo, nTimeTrain, true,  trainWindow, binWidth);
    [timeTest,  timeEventsTest]  = buildTimeOptions(testInfo,  nTimeTest,  true,  testWindow,  binWidth);
else
    [timeTrain, timeEventsTrain] = buildTimeOptions(trainInfo, nTimeTrain, false, [], []);
    [timeTest,  timeEventsTest]  = buildTimeOptions(testInfo,  nTimeTest,  false, [], []);
end

% ---- Run dPCA on training data -------------------------------------------
switch procedure
    case 'simple'
        [W, V, whichMarg] = dpca(XtrainCondAvg, numComps, ...
            'combinedParams', combinedParams);
        Cnoise        = [];
        trainTrialNum = [];
        optimalLambda = [];

    case {'regularized', 'regularized&noiseCovar'}
        trainTrialNum    = trainObj.TrialNumWithParams;
        % TensorWithParams is N x [params] x maxTrials x T.
        % dPCA functions expect N x [params] x T x maxTrials → swap last two dims.
        trainFiringRates = trainObj.TensorWithParams;
        nTFDims          = ndims(trainFiringRates);
        trainFiringRates = permute(trainFiringRates, [1:nTFDims-2, nTFDims, nTFDims-1]);

        optimalLambda = dpca_optimizeLambda(XtrainCondAvg, trainFiringRates, trainTrialNum, ...
            'combinedParams',  combinedParams, ...
            'simultaneous',    simultaneousRecording, ...
            'numRep',          10);

        Cnoise = dpca_getNoiseCovariance(XtrainCondAvg, trainFiringRates, trainTrialNum, ...
            'simultaneous', simultaneousRecording);

        [W, V, whichMarg] = dpca(XtrainCondAvg, numComps, ...
            'combinedParams', combinedParams, ...
            'lambda',         optimalLambda, ...
            'Cnoise',         Cnoise);

    otherwise
        error('computedPCA: unknown procedure ''%s''. Choose ''simple'', ''regularized'', or ''regularized&noiseCovar''.', procedure);
end

% % Allow caller to override the marginalization assignments
% if isfield(visOpt, 'whichMarg') && ~isempty(visOpt.whichMarg)
%     whichMarg = visOpt.whichMarg;
% end

% ---- Explained variance + plot — training set ----------------------------
if strcmp(procedure, 'regularized&noiseCovar')
    explVarTrain = dpca_explainedVariance(XtrainCondAvg, W, V, ...
        'combinedParams', combinedParams, ...
        'Cnoise',         Cnoise, ...
        'numOfTrials',    trainTrialNum);
else
    explVarTrain = dpca_explainedVariance(XtrainCondAvg, W, V, ...
        'combinedParams', combinedParams);
end

dpca_plot(XtrainCondAvg, W, V, @dpca_plot_default, ...
    'whichMarg',              whichMarg, ...
    'explainedVar',           explVarTrain, ...
    'marginalizationNames',   margNames, ...
    'marginalizationColours', margColours, ...
    'time',                   timeTrain, ...
    'timeEvents',             timeEventsTrain, ...
    'timeMarginalization',    ndims(XtrainCondAvg), ...
    'legendSubplot',          legendSubplot);

try
    infoStr = sprintf('neus %d-%d (n=%d) | train win [%.3g %.3g]s | %s | %d comps', ...
        min(neus2consider), max(neus2consider), numel(neus2consider), ...
        trainWindow(1), trainWindow(3), procedure, numComps);
    sgtitle({'dPCA — Training set', infoStr});
catch; end

% ---- Explained variance + plot — test set --------------------------------
explVarTest = dpca_explainedVariance(XtestCondAvg, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(XtestCondAvg, W, V, @dpca_plot_default, ...
    'whichMarg',              whichMarg, ...
    'explainedVar',           explVarTest, ...
    'marginalizationNames',   margNames, ...
    'marginalizationColours', margColours, ...
    'time',                   timeTest, ...
    'timeEvents',             timeEventsTest, ...
    'timeMarginalization',    ndims(XtestCondAvg), ...
    'legendSubplot',          legendSubplot);

try
    infoStr = sprintf('neus %d-%d (n=%d) | test win [%.3g %.3g]s | %s | %d comps', ...
        min(neus2consider), max(neus2consider), numel(neus2consider), ...
        testWindow(1), testWindow(3), procedure, numComps);
    sgtitle({'dPCA — Test set', infoStr});
catch; end

% ---- Pack outputs --------------------------------------------------------
dPCAmodel.W         = W;
dPCAmodel.V         = V;
dPCAmodel.whichMarg = whichMarg;
dPCAmodel.Cnoise     = Cnoise;
dPCAmodel.optimalLambda = optimalLambda;
dPCAmodel.simultaneousRecording = simultaneousRecording;
dPCAmodel.procedure = procedure;

dPCAprojTrain.explVar = explVarTrain;
dPCAprojTest.explVar  = explVarTest;

end

%% ---- LOCAL FUNCTIONS -------------------------------------------------------

function [timeAxis, timeEvents] = buildTimeOptions(info, nTime, fixedFlag, window, binWidth)
% Build a time axis and event-marker times for dpca_plot.
%
%  Fixed-bin  : timeAxis in seconds, aligned to window start, bin-centred.
%               timeEvents interpolated from MarkersFixedBin bin positions.
%  Variable-bin: timeAxis in seconds using expectedBinWidth as step.
%               timeEvents at epoch boundaries (cumsum of nBins × binWidth).

    if fixedFlag
        if numel(binWidth) > 1 && binWidth(2) > 0
            timeStep = binWidth(2);
        else
            timeStep = binWidth(1);
        end
        timeAxis = window(1) + (0:nTime-1) * timeStep + binWidth(1)/2;

        if isfield(info, 'MarkersFixedBin') && ~isempty(info.MarkersFixedBin)
            markerBins = info.MarkersFixedBin(:)';
            % Interpolate bin positions onto the time axis (out-of-range → NaN)
            timeEvents = interp1(1:nTime, timeAxis, markerBins, 'linear', NaN)';
            timeEvents = timeEvents(isfinite(timeEvents));
        else
            timeEvents = [];
        end
    else
        bw       = info.expectedBinWidth;
        timeAxis = (0:nTime-1) .* bw;

        if isfield(info, 'nBins') && ~isempty(info.nBins)
            timeEvents = cumsum(info.nBins(1:end-1))' .* bw;
        else
            timeEvents = [];
        end
    end
end

function val = getOptField(s, field, default)
% Return s.field if present, otherwise default.
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end