function [neuralDecodingScores, neuralDecodingData, neuralDecodingModels] = neuralDecodingClassification(obj, processingOpt, trainOpt, testOpt, decoderOpt, visOpt)
% NEURALDECODINGCLASSIFICATION  Decode task parameters from neural population activity.
%
%   For each cross-validation fold the method trains one classifier on all
%   time bins of the training set, then evaluates per-time-bin performance on
%   the training data, on the held-out (validation) fold, and on the
%   independent test set.  A random (chance) decoder is run in parallel so
%   that chance-level performance can be visualised with shaded error bands.
%
% -------------------------------------------------------------------------
% SYNTAX
%   [scores, data, models] = obj.neuralDecodingClassification( ...
%       processingOpt, trainOpt, testOpt)
%
%   [scores, data, models] = obj.neuralDecodingClassification( ...
%       processingOpt, trainOpt, testOpt, decoderOpt, visOpt)
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS
%   processingOpt  - Struct controlling neuron selection and binning:
%     .neus2consider       - (1 x N) neuron indices to include
%     .fixedFlag           - true  → fixed-bin pipeline
%                            false → variable-bin pipeline 
%     if fixedFlag = true:
%       .binWidth           - (1 x 2) [binSize, stepSize] in seconds.
%                             Set stepSize = 0 for non-overlapping bins.
%     if fixedFlag = false:
%       .expectedBinWidth   - expected bin width for prepareTensorWithVariableBin
%
%   trainOpt       - Struct defining the training set:
%     .conds2consider4train - condition indices for the training tensor
%     .trainParams          - (nConds x nParams) parameter matrix
%     .trainWindow          - analysis window; same format as
%                             prepareTensorWithFixBin / prepareTensorWithVariableBin
%
%   testOpt        - Struct defining the test set:
%     .conds2consider4test  - condition indices for the test tensor
%     .testParams           - (nConds x nParams) parameter matrix
%     .testWindow           - analysis window (same format as trainWindow)
%
% -------------------------------------------------------------------------
% OPTIONAL INPUTS
%   decoderOpt     - Struct controlling the classifier (all fields optional):
%     .CVKfold      - Number of CV folds (integer >= 2).  Default: 10.
%     .decoder      - 'naivebayes' (default) | 'lda' | 'svm'
%     .decoderOpts  - Struct forwarded as name-value pairs to the underlying
%                     fit* function (fitcnb / fitcdiscr / fitcsvm).
%                     Default: struct().
%
%   visOpt         - Struct controlling output visualisation (all optional):
%     .factorNames  - Cell array naming each task parameter
%                     (e.g. {'target', 'obstacle'}).  'epoch' is appended
%                     automatically.  Default: {'param1', 'param2', ...}
%     .display      - Metric to plot vs. time, or 'off' to suppress figures.
%                     Must match a metricNames entry (e.g. 'accuracy').
%                     Default: 'off'.
%     .markerNames  - Cell array labelling temporal markers drawn as vertical
%                     lines in the decoding plots.  Default: {}.
%
% -------------------------------------------------------------------------
% OUTPUTS
%   neuralDecodingScores  - Struct with decoding performance over time:
%     .factorNames        - {1 x nFactors} factor labels (params + 'epoch')
%     .metricNames        - {nMetrics x 1} metric labels
%     .train              - {CVKfold x nFactors x nTimeTrain} statsOfMeasure
%     .validation         - {CVKfold x nFactors x nTimeTrain} statsOfMeasure
%     .test               - {CVKfold x nFactors x nTimeTest}  statsOfMeasure
%     .chanceTrain        - same as .train  but for the random decoder
%     .chanceValidation   - same as .validation but for the random decoder
%     .chanceTest         - same as .test    but for the random decoder
%     .timeAxisTrain      - (1 x nTimeTrain) time axis in seconds
%     .timeAxisTest       - (1 x nTimeTest)  time axis in seconds
%     .markerTimesTrain   - marker times (s) within the training window
%     .markerTimesTest    - marker times (s) within the test window
%
%   neuralDecodingData    - Struct with the flattened feature matrices:
%     .Xtrain             - (nSamplesTrain x nNeurons)
%     .labelsTrain        - (nSamplesTrain x nFactors)  [param1 ... epoch]
%     .sampleInfoTrain    - table: time, trial, paramCombination, param1 ...
%     .Xtest              - (nSamplesTest  x nNeurons)
%     .labelsTest         - (nSamplesTest  x nFactors)
%     .sampleInfoTest     - table
%
%   neuralDecodingModels  - Struct with the trained classifiers:
%     .decoder            - name of the algorithm used
%     .decoderOpts        - options struct forwarded to fit*
%     .models             - {CVKfold x nFactors} cell array of fitted models
%     .cvpartition        - cvpartition object used for CV
%     .uniqueTrials       - trial IDs used for partitioning
%
% -------------------------------------------------------------------------
% EXAMPLE  (fixed-bin LDA, single target parameter)
%   procOpt.fixedFlag        = true;
%   procOpt.neus2consider    = 1:96;
%   procOpt.binWidth         = [0.05 0.01];
%
%   trOpt.conds2consider4train = 1:8;
%   trOpt.trainParams          = conditions.Target6(1:8);
%   trOpt.trainWindow          = [-0.5 1.0];
%
%   teOpt.conds2consider4test  = 1:8;
%   teOpt.testParams           = conditions.Obstacle6(1:8);
%   teOpt.testWindow           = [-0.5 1.0];
%
%   decOpt.decoder    = 'lda';
%   decOpt.CVKfold    = 5;
%
%   visOpt.factorNames = {'target'};
%   visOpt.display     = 'accuracy';
%   visOpt.markerNames = {'go-cue', 'movement-onset'};
%
%   [scores, data, mdls] = obj.neuralDecodingClassification( ...
%       procOpt, trOpt, teOpt, decOpt, visOpt);
%
% See also: prepareTensorWithFixBin, splitTensorWithParams, statsOfMeasure

% ---- Unpack processingOpt -----------------------------------------------
fixedFlag     = processingOpt.fixedFlag;
neus2consider = processingOpt.neus2consider;

if fixedFlag
    binWidth = processingOpt.binWidth;
else
    expectedBinWidth        = processingOpt.expectedBinWidth;
end

% ---- Unpack trainOpt / testOpt ------------------------------------------
conds2consider4train = trainOpt.conds2consider4train;
trainParams          = trainOpt.trainParams;
trainWindow          = trainOpt.trainWindow;

conds2consider4test  = testOpt.conds2consider4test;
testParams           = testOpt.testParams;
testWindow           = testOpt.testWindow;

% ---- Unpack decoderOpt (optional) ---------------------------------------
if nargin < 5 || isempty(decoderOpt)
    decoderOpt = struct();
end
CVKfold     = getOptField(decoderOpt, 'CVKfold',     10);
decoderName = getOptField(decoderOpt, 'decoder',     'naivebayes');
decoderOpts = getOptField(decoderOpt, 'decoderOpts', struct());

if ~isnumeric(CVKfold) || ~isscalar(CVKfold) || CVKfold < 2 || round(CVKfold) ~= CVKfold
    error('decoderOpt.CVKfold must be an integer >= 2.');
end
decoderName = validatestring(char(string(decoderName)), {'naivebayes','lda','svm'}, mfilename, 'decoderOpt.decoder');
if ~isstruct(decoderOpts), error('decoderOpt.decoderOpts must be a struct.'); end

% ---- Unpack visOpt (optional) -------------------------------------------
if nargin < 6 || isempty(visOpt)
    visOpt = struct();
end
factorNamesIn = getOptField(visOpt, 'factorNames', {});
displayMetric = char(string(getOptField(visOpt, 'display',     'off')));
markerNames   = getOptField(visOpt, 'markerNames', {});

if ischar(markerNames) || isstring(markerNames)
    markerNames = cellstr(markerNames(:))';
elseif iscell(markerNames)
    markerNames = cellfun(@char, markerNames(:)', 'UniformOutput', false);
end

% -------------------------------------------------------------------------
metricNames = ["true_positive"; "false_positive"; "false_negative"; ...
    "true_negative"; "precision"; "sensitivity"; "specificity"; ...
    "accuracy"; "F-measure"];

% ---- Build tensors -------------------------------------------------------
if fixedFlag
    trainObj = obj.prepareTensorWithFixBin(trainWindow, binWidth, neus2consider, conds2consider4train);
    trainObj = trainObj.splitTensorWithParams(trainParams, true);
    [Xtrain, labelsTrain, sampleInfoTrain, metaTrain] = tensorToRows(trainObj, trainWindow, binWidth, true);

    testObj = obj.prepareTensorWithFixBin(testWindow, binWidth, neus2consider, conds2consider4test);
    testObj = testObj.splitTensorWithParams(testParams, true);
    [Xtest, labelsTest, sampleInfoTest, metaTest] = tensorToRows(testObj, testWindow, binWidth, true);

    markerDataTrain = computeMarkerTimes(trainObj, trainWindow, binWidth, true);
    markerDataTest  = computeMarkerTimes(testObj,  testWindow,  binWidth, true);
    timeAxisTrain   = makeTimeAxis(trainWindow, binWidth, metaTrain.nTime);
    timeAxisTest    = makeTimeAxis(testWindow,  binWidth, metaTest.nTime);
    xLabel          = 'Time (s)';
else
    trainObj = obj.prepareTensorWithVariableBin(trainWindow, expectedBinWidth, neus2consider, conds2consider4train);
    trainObj = trainObj.splitTensorWithParams(trainParams, false);
    [Xtrain, labelsTrain, sampleInfoTrain, metaTrain] = tensorToRows(trainObj, [], [], false);

    testObj = obj.prepareTensorWithVariableBin(testWindow, expectedBinWidth, neus2consider, conds2consider4test);
    testObj = testObj.splitTensorWithParams(testParams, false);
    [Xtest, labelsTest, sampleInfoTest, metaTest] = tensorToRows(testObj, [], [], false);

    markerDataTrain = computeMarkerTimes(trainObj, [], [], false);
    markerDataTest  = computeMarkerTimes(testObj,  [], [], false);
    timeAxisTrain   = 1:metaTrain.nTime;
    timeAxisTest    = 1:metaTest.nTime;
    xLabel          = 'Bin index';
end

% nFactors = nParams + 1  (epoch occupies the last column)
nFactors   = size(labelsTrain, 2);
nParams    = nFactors - 1;
nTimeTrain = metaTrain.nTime;
nTimeTest  = metaTest.nTime;



if isempty(factorNamesIn)
    factorNamesIn = cellstr(compose('param%d', 1:nParams));
end
factorNamesIn = factorNamesIn(1:min(numel(factorNamesIn), nParams));
while numel(factorNamesIn) < nParams
    factorNamesIn{end+1} = sprintf('param%d', numel(factorNamesIn)+1);
end
allFactorNames = [factorNamesIn(:)', {'epoch'}];

% Trial-grouped CV partition on training trials
uniqueTrials = unique(sampleInfoTrain.trial, 'stable');
CVKfold = min(CVKfold, numel(uniqueTrials));
if CVKfold < 2, error('At least 2 unique training trials are required.'); end
cvp = cvpartition(numel(uniqueTrials), 'KFold', CVKfold);

rawResultsTrain  = cell(CVKfold, nFactors, nTimeTrain);
rawResultsVal    = cell(CVKfold, nFactors, nTimeTrain);
rawResultsTest   = cell(CVKfold, nFactors, nTimeTest);
chanceResultsTrain = cell(CVKfold, nFactors, nTimeTrain);
chanceResultsVal   = cell(CVKfold, nFactors, nTimeTrain);
chanceResultsTest  = cell(CVKfold, nFactors, nTimeTest);
modelsByFoldFactor = cell(CVKfold, nFactors);

for cv = 1:CVKfold
    disp(['Training, validating and testing on CV ' num2str(cv)])
    trTrials  = uniqueTrials(training(cvp, cv));
    valTrials = uniqueTrials(test(cvp, cv));
    isTrain   = ismember(sampleInfoTrain.trial, trTrials);
    isVal     = ismember(sampleInfoTrain.trial, valTrials);

    XtrainCV = Xtrain(isTrain, :);
    XvalCV   = Xtrain(isVal, :);
    timeTrain = sampleInfoTrain.time(isTrain);
    timeVal   = sampleInfoTrain.time(isVal);

    for factorIdx = 1:nFactors
        yTrain = labelsTrain(isTrain, factorIdx);
        yVal   = labelsTrain(isVal,   factorIdx);
        yTest  = labelsTest(:, factorIdx);

        if numel(unique(yTrain)) < 2, continue; end

        mdl = fitDecoder(XtrainCV, yTrain, decoderName, decoderOpts);
        modelsByFoldFactor{cv, factorIdx} = mdl;
        if isempty(mdl), continue; end

        predTrain = predictDecoder(mdl, XtrainCV);
        predVal   = predictDecoder(mdl, XvalCV);
        predTest  = predictDecoder(mdl, Xtest);
        classLabels = unique([yTrain; yVal; yTest], 'stable');

        % Random (chance) predictions sampled uniformly from classLabels
        nCls = numel(classLabels);
        randTrain = classLabels(randi(nCls, numel(yTrain), 1));
        randVal   = classLabels(randi(nCls, numel(yVal),   1));
        randTest  = classLabels(randi(nCls, numel(yTest),  1));

        for t = 1:nTimeTrain
            maskTr = timeTrain == t;
            if any(maskTr)
                C = confusionmat(yTrain(maskTr), predTrain(maskTr), 'Order', classLabels);
                rawResultsTrain{cv, factorIdx, t} = statsOfMeasure(C, false);
                Cr = confusionmat(yTrain(maskTr), randTrain(maskTr), 'Order', classLabels);
                chanceResultsTrain{cv, factorIdx, t} = statsOfMeasure(Cr, false);
            end
            maskV = timeVal == t;
            if any(maskV)
                C = confusionmat(yVal(maskV), predVal(maskV), 'Order', classLabels);
                rawResultsVal{cv, factorIdx, t} = statsOfMeasure(C, false);
                Cr = confusionmat(yVal(maskV), randVal(maskV), 'Order', classLabels);
                chanceResultsVal{cv, factorIdx, t} = statsOfMeasure(Cr, false);
            end
        end

        for t = 1:nTimeTest
            maskTe = sampleInfoTest.time == t;
            if any(maskTe)
                C = confusionmat(yTest(maskTe), predTest(maskTe), 'Order', classLabels);
                rawResultsTest{cv, factorIdx, t} = statsOfMeasure(C, false);
                Cr = confusionmat(yTest(maskTe), randTest(maskTe), 'Order', classLabels);
                chanceResultsTest{cv, factorIdx, t} = statsOfMeasure(Cr, false);
            end
        end
    end
end

neuralDecodingData = struct();
neuralDecodingData.Xtrain         = Xtrain;
neuralDecodingData.labelsTrain    = labelsTrain;
neuralDecodingData.sampleInfoTrain = sampleInfoTrain;
neuralDecodingData.Xtest          = Xtest;
neuralDecodingData.labelsTest     = labelsTest;
neuralDecodingData.sampleInfoTest = sampleInfoTest;

neuralDecodingModels = struct();
neuralDecodingModels.decoder     = decoderName;
neuralDecodingModels.decoderOpts = decoderOpts;
neuralDecodingModels.models      = modelsByFoldFactor;
neuralDecodingModels.cvpartition = cvp;
neuralDecodingModels.uniqueTrials = uniqueTrials;

neuralDecodingScores = struct();
neuralDecodingScores.factorNames     = allFactorNames;
neuralDecodingScores.metricNames     = metricNames;
neuralDecodingScores.train           = rawResultsTrain;
neuralDecodingScores.validation      = rawResultsVal;
neuralDecodingScores.test            = rawResultsTest;
neuralDecodingScores.chanceTrain     = chanceResultsTrain;
neuralDecodingScores.chanceValidation = chanceResultsVal;
neuralDecodingScores.chanceTest      = chanceResultsTest;
neuralDecodingScores.timeAxisTrain    = timeAxisTrain;
neuralDecodingScores.timeAxisTest     = timeAxisTest;
neuralDecodingScores.markerTimesTrain = markerDataTrain.plotted;
neuralDecodingScores.markerTimesTest  = markerDataTest.plotted;

if ~strcmpi(displayMetric, 'off')
    metricIdx = find(strcmpi(metricNames, displayMetric), 1);
    if ~isempty(metricIdx)
        infoStr = sprintf('neus %d-%d (n=%d) | train win [%.3g %.3g]s | test win [%.3g %.3g]s | %s %d-fold CV', ...
            min(neus2consider), max(neus2consider), numel(neus2consider), ...
            trainWindow(1), trainWindow(3), testWindow(1), testWindow(3), ...
            decoderName, CVKfold);
        plotDecodingCurves(rawResultsVal,  chanceResultsVal,  'Validation (Training Hold out)', timeAxisTrain, allFactorNames, displayMetric, metricIdx, CVKfold, markerDataTrain, markerNames, xLabel, infoStr);
        plotDecodingCurves(rawResultsTest, chanceResultsTest, 'Test',                           timeAxisTest,  allFactorNames, displayMetric, metricIdx, CVKfold, markerDataTest,  markerNames, xLabel, infoStr);
    end
end

end

%% ---- LOCAL FUNCTIONS -------------------------------------------------------

function [X, labels, sampleInfo, meta] = tensorToRows(dataObj, window, binWidth, fixedFlag)
% Convert TensorWithParams to a flat sample matrix.
%
% labels format: [param1, param2, ..., epoch]   epoch is the LAST column.
%
% For fixed-bin tensors (fixedFlag=true):
%   - time axis is in seconds; epoch boundaries come from MarkerTensorForFixBin.
% For variable-bin tensors (fixedFlag=false):
%   - time axis is bin index (1:totalBins); epoch boundaries are derived from
%     TensorWithVariableBinInfo.nBins stored by prepareTensorWithVariableBin.

    if nargin < 4, fixedFlag = true; end

    tensor = dataObj.TensorWithParams;
    paramLevels = dataObj.TensorWithParamsInfo.paramLevels;
    nParamsLocal = numel(paramLevels);
    paramDims = cellfun(@numel, paramLevels);
    nNeurons = size(tensor, 1);
    nTrials  = size(tensor, nParamsLocal + 2);
    nTimeLocal = size(tensor, nParamsLocal + 3);
    nParamCombinations = prod(paramDims);

    if fixedFlag
        % Time axis in seconds; epoch boundaries from MarkerTensorForFixBin
        timeAxis = makeTimeAxis(window, binWidth, nTimeLocal);
        if ~isempty(dataObj.MarkerTensorForFixBin)
            markerBins = squeeze(mean(dataObj.MarkerTensorForFixBin, [1 2 3], 'omitmissing'));
            markerBins = markerBins(:)';
            if binWidth(2) > 0, timeStep = binWidth(2); else, timeStep = binWidth(1); end
            epochBoundaries = window(1) + markerBins * timeStep;
            epochBoundaries = sort(epochBoundaries(isfinite(epochBoundaries)));
        else
            epochBoundaries = [];
        end
    else
        % Time axis as bin index; epoch boundaries from nBins
        timeAxis = 1:nTimeLocal;
        if isfield(dataObj.TensorWithVariableBinInfo, 'nBins')
            nBins = dataObj.TensorWithVariableBinInfo.nBins;
            % Boundaries at 0.5 above cumulative bin counts: bin k is epoch e
            % iff cumsum(nBins(1:e-1)) < k <= cumsum(nBins(1:e))
            epochBoundaries = cumsum(nBins(1:end-1)) + 0.5;
        else
            epochBoundaries = [];
        end
    end

    % Parameter combination grid
    paramGrids = cell(1, nParamsLocal);
    [paramGrids{:}] = ndgrid(paramLevels{:});
    paramValues = nan(nParamCombinations, nParamsLocal);
    for pIdx = 1:nParamsLocal
        paramValues(:, pIdx) = paramGrids{pIdx}(:);
    end

    maxRows = nParamCombinations * nTrials * nTimeLocal;
    X              = nan(maxRows, nNeurons);
    labels         = nan(maxRows, nParamsLocal + 1);  % +1 for epoch (last)
    sampleTimeIdx  = nan(maxRows, 1);
    trialIdx       = nan(maxRows, 1);
    paramComboIdx  = nan(maxRows, 1);
    sampleParamValues = nan(maxRows, nParamsLocal);

    % Precompute subscript indices for each parameter combination.
    % ind2sub requires at least 2 size elements, so handle single-param case.
    paramSubsAll = zeros(nParamCombinations, nParamsLocal);
    if nParamsLocal == 1
        paramSubsAll(:, 1) = (1:nParamCombinations)';
    else
        tmpCell = cell(1, nParamsLocal);
        for ci = 1:nParamCombinations
            [tmpCell{:}] = ind2sub(paramDims, ci);
            paramSubsAll(ci, :) = [tmpCell{:}];
        end
    end

    row = 0;
    for tIdx = 1:nTimeLocal
        epochLabel = sum(timeAxis(tIdx) > epochBoundaries) + 1;
        for trIdx = 1:nTrials
            for comboIdx = 1:nParamCombinations
                tensorIdx = [{':'}, num2cell(paramSubsAll(comboIdx, :)), {trIdx, tIdx}];
                sample = reshape(tensor(tensorIdx{:}), 1, nNeurons);

                if all(isfinite(sample))
                    row = row + 1;
                    X(row, :)      = sample;
                    labels(row, :) = [paramValues(comboIdx, :), epochLabel];
                    sampleTimeIdx(row) = tIdx;
                    trialIdx(row)      = trIdx;
                    paramComboIdx(row) = comboIdx;
                    sampleParamValues(row, :) = paramValues(comboIdx, :);
                end
            end
        end
    end

    X      = X(1:row, :);
    labels = labels(1:row, :);
    sampleTimeIdx     = sampleTimeIdx(1:row);
    trialIdx          = trialIdx(1:row);
    paramComboIdx     = paramComboIdx(1:row);
    sampleParamValues = sampleParamValues(1:row, :);

    sampleInfo = table(sampleTimeIdx, trialIdx, paramComboIdx, ...
        'VariableNames', {'time', 'trial', 'paramCombination'});
    for pIdx = 1:nParamsLocal
        sampleInfo.(sprintf('param%d', pIdx)) = sampleParamValues(:, pIdx);
    end

    meta = struct();
    meta.nNeurons   = nNeurons;
    meta.nParams    = nParamsLocal;
    meta.nMaxTrials = nTrials;
    meta.nTime      = nTimeLocal;
end

function model = fitDecoder(X, y, decoderName, decoderOpts)
% Fit a classifier. X is (samples x neurons), y is (samples x 1).
    y = y(:);
    if isempty(X) || isempty(y) || numel(unique(y)) < 2
        model = [];
        return
    end

    X = X + rand(size(X)) * 1e-8;  % tiny jitter to avoid zero-variance features
    opts = struct2pairs(decoderOpts);

    switch lower(decoderName)
        case 'naivebayes'
            model = fitcnb(X, y, opts{:});
        case 'lda'
            model = fitcdiscr(X, y, opts{:});
        case 'svm'
            if numel(unique(y)) <= 2
                model = fitcsvm(X, y, opts{:});
            else
                learner = templateSVM(opts{:});
                model = fitcecoc(X, y, 'Learners', learner);
            end
        otherwise
            error('Unsupported decoder: %s', decoderName);
    end
end

function yPred = predictDecoder(model, X)
% Apply a classifier. Returns NaNs if model is empty.
    if isempty(model)
        yPred = nan(size(X, 1), 1);
        return
    end
    X = X + rand(size(X)) * 1e-8;
    yPred = predict(model, X);
    yPred = yPred(:);
end

function timeAxis = makeTimeAxis(window, binWidth, nTimeLocal)
    if binWidth(2) > 0
        timeStep = binWidth(2);
    else
        timeStep = binWidth(1);
    end
    timeAxis = window(1) + (0:nTimeLocal-1) * timeStep + binWidth(1)/2;
end

function pairs = struct2pairs(s)
    names = fieldnames(s);
    if isempty(names)
        pairs = {};
        return
    end
    pairs = cell(1, 2*numel(names));
    for i = 1:numel(names)
        pairs{2*i-1} = names{i};
        pairs{2*i}   = s.(names{i});
    end
end

function plotDecodingCurves(rawResults, chanceResults, setLabel, timeAxis, factorNames, metricName, metricIdx, CVKfold, markerData, markerNames, xLabel, infoStr)
% Plot mean +/- std across CV folds for each factor.
% Chance level (gray shaded area) is overlaid when chanceResults is non-empty.
% Vertical lines mark temporal markers; labels come from markerNames.
% xLabel controls the x-axis label ('Time (s)' for fixed bin, 'Bin index' for variable bin).
    if nargin < 11 || isempty(xLabel), xLabel = 'Time (s)'; end
    nFactors = size(rawResults, 2);
    nTime    = size(rawResults, 3);
    meanVals   = nan(nTime, nFactors);
    stdVals    = nan(nTime, nFactors);
    meanChance = nan(nTime, nFactors);
    stdChance  = nan(nTime, nFactors);

    for f = 1:nFactors
        for t = 1:nTime
            vals  = nan(CVKfold, 1);
            cvals = nan(CVKfold, 1);
            for cv = 1:CVKfold
                s = rawResults{cv, f, t};
                if istable(s) && height(s) >= metricIdx
                    vals(cv) = s.macroAVG(metricIdx);
                end
                if ~isempty(chanceResults)
                    sc = chanceResults{cv, f, t};
                    if istable(sc) && height(sc) >= metricIdx
                        cvals(cv) = sc.macroAVG(metricIdx);
                    end
                end
            end
            vf = vals(isfinite(vals));
            if ~isempty(vf)
                meanVals(t, f) = mean(vf);
                stdVals(t, f)  = std(vf, 0);
            end
            cf = cvals(isfinite(cvals));
            if ~isempty(cf)
                meanChance(t, f) = mean(cf);
                stdChance(t, f)  = std(cf, 0);
            end
        end
    end

    figure;
    colors = lines(nFactors);
    for f = 1:nFactors
        % Chance level first (behind)
        if ~isempty(chanceResults)
            shaded_areas(timeAxis(:), meanChance(:,f), stdChance(:,f), ...
                'FaceColor', colors(f,:), 'FaceAlpha', 0.25, ...
                'LineWidth', 1, 'LineStyle', '--',...
                'DisplayName', [factorNames{f} ' (chance)']);
        end
        % Decoder performance
        shaded_areas(timeAxis(:), meanVals(:,f), stdVals(:,f), ...
            'FaceColor', colors(f,:), 'FaceAlpha', 0.5, ...
            'LineWidth', 2, 'DisplayName', factorNames{f});
    end
    hold on;
    plotMarkerLines(markerData, markerNames);
    legend('Interpreter', 'none');
    xlabel(xLabel);  ylabel(metricName);
    title({sprintf('Decoding %s — %s', setLabel, metricName), infoStr});
end

function mt = computeMarkerTimes(dataObj, window, binWidth, fixedFlag)
% Extract marker positions that fall within the analysis window.
%
% Fixed-bin: returns times in seconds aligned with makeTimeAxis output.
% Variable-bin: returns fractional bin positions (0.5-shifted) aligned with
%   the 1:totalBins time axis produced by tensorToRows.
    if nargin < 4, fixedFlag = true; end

    if fixedFlag
        if isempty(dataObj.MarkerTensorForFixBin)
            mt.all = []; mt.plotted = []; mt.plottedIdx = [];
            return
        end
        if binWidth(2) > 0, timeStep = binWidth(2); else, timeStep = binWidth(1); end
        markerBins = squeeze(mean(dataObj.MarkerTensorForFixBin, [1 2 3], 'omitmissing'));
        markerBins = markerBins(:)';
        allTimes = window(1) + markerBins * timeStep;
        mask = isfinite(allTimes) & allTimes >= window(1) & allTimes <= window(end);
        % For fixed bin, allTimes spans all experiment markers, so find(mask)
        % already equals the absolute experiment marker indices.
        absPlottedIdx = find(mask);
    else
        % Variable bin: marker positions are cumulative bin counts stored in nBins.
        % Shift by +0.5 so that position k falls between bin k and bin k+1 on the
        % 1-indexed time axis produced by tensorToRows.
        if ~isfield(dataObj.TensorWithVariableBinInfo, 'nBins')
            mt.all = []; mt.plotted = []; mt.plottedIdx = [];
            return
        end
        nBins     = dataObj.TensorWithVariableBinInfo.nBins;
        totalBins = sum(nBins);
        allTimes  = [0, cumsum(nBins)] + 0.5;  % length = nMarker (= nEpoch+1)
        % Keep only interior boundaries (exclude start=0.5 and end=totalBins+0.5)
        mask = allTimes > 0.5 & allTimes < totalBins + 0.5;
        % Translate local mask positions to absolute experiment marker indices
        % (marker k in the window corresponds to experiment marker windowStart+k-1)
        % so that markerNames(plottedIdx) selects the correct labels when
        % markerNames spans all experiment markers.
        windowStart       = dataObj.TensorWithVariableBinInfo.window(1);
        localPlottedIdx   = find(mask);
        absPlottedIdx     = windowStart + localPlottedIdx - 1;
    end

    mt.all        = allTimes;
    mt.plotted    = allTimes(mask);
    mt.plottedIdx = absPlottedIdx;
end

function plotMarkerLines(markerData, markerNames)
% Draw vertical dashed lines at marker times, with optional labels.
    if isempty(markerData) || isempty(markerData.plotted)
        return
    end
    plottedTimes = markerData.plotted;
    plottedIdx   = markerData.plottedIdx;
    allTimes     = markerData.all;

    if isempty(markerNames)
        xline(plottedTimes, 'k--');
        return
    end

    if numel(markerNames) == numel(allTimes)
        % markerNames has one entry per entry in allTimes (fixed-bin typical case)
        names = markerNames(plottedIdx);
    elseif ~isempty(plottedIdx) && numel(markerNames) >= max(plottedIdx)
        % plottedIdx holds absolute experiment marker indices
        % (variable-bin case, or fixed-bin when markerNames spans all experiment markers)
        names = markerNames(plottedIdx);
    elseif numel(markerNames) == numel(plottedTimes)
        % markerNames has one entry per plotted line
        names = markerNames;
    else
        xline(plottedTimes, 'k--');
        return
    end

    xline(plottedTimes, 'k--', names, ...
        'LabelOrientation', 'aligned', ...
        'LabelHorizontalAlignment', 'left');
end

function val = getOptField(s, field, default)
% Return s.field if it exists, otherwise return default.
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end

