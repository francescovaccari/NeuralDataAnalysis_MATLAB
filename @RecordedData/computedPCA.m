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
%     .paramNames             - cell array of names for the task parameters
%                               (the columns of trainParams/testParams),
%                               used to label the legend panel in the
%                               summary figure instead of the generic
%                               'Param 1', 'Param 2', ... .
%     .markerNames            - cell array labelling the temporal markers
%                               drawn as vertical lines in the component
%                               plots (same convention as in
%                               neuralDecodingClassification). May either
%                               span all experiment markers (indexed the
%                               same way as AlignedMarkerTensorForFixBin / the
%                               variable-bin window markers) or have one
%                               entry per marker actually plotted within
%                               trainWindow/testWindow. Default: {} (lines
%                               are drawn unlabeled).
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
%   vOpt.paramNames           = {'target','obstacle'};
%   vOpt.markerNames          = {'go-cue', 'movement-onset'};
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
paramNames    = getOptField(visOpt, 'paramNames', {});
markerNames   = getOptField(visOpt, 'markerNames', {});

if ischar(markerNames) || isstring(markerNames)
    markerNames = cellstr(markerNames(:))';
elseif iscell(markerNames)
    markerNames = cellfun(@char, markerNames(:)', 'UniformOutput', false);
end

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
    timeTrain = trainObj.TensorWithFixBinInfo.relativeBinCenters;
    timeTest  = testObj.TensorWithFixBinInfo.relativeBinCenters;
    [timeEventsTrain, timeEventsIdxTrain] = trainObj.computeAlignedMarkerTimes(trainWindow);
    [timeEventsTest,  timeEventsIdxTest]  = testObj.computeAlignedMarkerTimes(testWindow);
else
    [timeTrain, timeEventsTrain, timeEventsIdxTrain] = buildTimeOptions(trainInfo, nTimeTrain, false, trainWindow, []);
    [timeTest,  timeEventsTest,  timeEventsIdxTest]  = buildTimeOptions(testInfo,  nTimeTest,  false, testWindow,  []);
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

%% Original plotting function as provided by Kobak2016
% dpca_plot(XtrainCondAvg, W, V, @dpca_plot_default, ...
%     'whichMarg',              whichMarg, ...
%     'explainedVar',           explVarTrain, ...
%     'marginalizationNames',   margNames, ...
%     'marginalizationColours', margColours, ...
%     'time',                   timeTrain, ...
%     'timeEvents',             timeEventsTrain, ...
%     'timeMarginalization',    ndims(XtrainCondAvg), ...
%     'legendSubplot',          legendSubplot);

%% Generalized plotting function for more parameters
dpca_plot_general(XtrainCondAvg, W, V, @dpca_plot_general_default, ...
    'whichMarg',              whichMarg, ...
    'explainedVar',           explVarTrain, ...
    'marginalizationNames',   margNames, ...
    'time',                   timeTrain, ...
    'timeEvents',             timeEventsTrain, ...
    'timeEventsIdx',          timeEventsIdxTrain, ...
    'timeMarginalization',    ndims(XtrainCondAvg), ...
    'legendSubplot',          legendSubplot, ...
    'paramNames',             paramNames, ...
    'markerNames',            markerNames, ...
    'uniqueYLim',             true);

try
    infoStr = sprintf('neus %d-%d (n=%d) | train win [%.3g %.3g]s | %s | %d comps', ...
        min(neus2consider), max(neus2consider), numel(neus2consider), ...
        trainWindow(1), trainWindow(3), procedure, numComps);
    sgtitle({'dPCA — Training set', infoStr});
catch; end

% ---- Explained variance + plot — test set --------------------------------
explVarTest = dpca_explainedVariance(XtestCondAvg, W, V, ...
    'combinedParams', combinedParams);

%% Original plotting function as provided by Kobak2016
%dpca_plot(XtestCondAvg, W, V, @dpca_plot_default, ...
   % 'whichMarg',              whichMarg, ...
    %'explainedVar',           explVarTest, ...
    %'marginalizationNames',   margNames, ...
    %'marginalizationColours', margColours, ...
    %'time',                   timeTest, ...
    %'timeEvents',             timeEventsTest, ...
    %'timeMarginalization',    ndims(XtestCondAvg), ...
    %'legendSubplot',          legendSubplot);

%% Generalized plotting function for more parameters
dpca_plot_general(XtestCondAvg, W, V, @dpca_plot_general_default, ...
    'whichMarg',              whichMarg, ...
    'explainedVar',           explVarTest, ...
    'marginalizationNames',   margNames, ...
    'time',                   timeTest, ...
    'timeEvents',             timeEventsTest, ...
    'timeEventsIdx',          timeEventsIdxTest, ...
    'timeMarginalization',    ndims(XtestCondAvg), ...
    'legendSubplot',          legendSubplot, ...
    'paramNames',             paramNames, ...
    'markerNames',            markerNames, ...
    'uniqueYLim',             true);

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
dPCAmodel.simultaneousRecording = obj.simultaneousRecording;
dPCAmodel.procedure = procedure;

dPCAprojTrain.explVar = explVarTrain;
dPCAprojTest.explVar  = explVarTest;

end

%% ---- LOCAL FUNCTIONS -------------------------------------------------------

function [timeAxis, timeEvents, timeEventsIdx] = buildTimeOptions(info, nTime, fixedFlag, window, binWidth)
% Build a time axis and event-marker times for dpca_plot (variable-bin only;
% the fixed-bin case is handled directly in the main function body via
% TensorWithFixBinInfo.relativeBinCenters and computeAlignedMarkerTimes).
%
%  Variable-bin: timeAxis in seconds using expectedBinWidth as step.
%               timeEvents at epoch boundaries (cumsum of nBins × binWidth).
%
%  timeEventsIdx holds, for each entry in timeEvents, the absolute
%  experiment-marker index it corresponds to: window(1) + 1:numel(nBins)-1,
%  since nBins only covers the markers within [window(1) window(2)]. Lets
%  visOpt.markerNames (which may span all experiment markers) be matched to
%  the markers actually plotted.

    if fixedFlag
        error('buildTimeOptions: the fixed-bin case is handled directly in computedPCA, not via this helper.');
    end

    bw       = info.expectedBinWidth;
    timeAxis = (0:nTime-1) .* bw;

    if isfield(info, 'nBins') && ~isempty(info.nBins)
        timeEvents = cumsum(info.nBins(1:end-1))' .* bw;
        if ~isempty(window)
            timeEventsIdx = window(1) + (1:numel(info.nBins)-1)';
        else
            timeEventsIdx = (1:numel(info.nBins)-1)';
        end
    else
        timeEvents = [];
        timeEventsIdx = [];
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