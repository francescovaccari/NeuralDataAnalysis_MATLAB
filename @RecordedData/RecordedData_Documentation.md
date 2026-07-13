# RecordedData Class Documentation

Load, process, and analyse neural recording data (spike times, state markers, kinematics, eye tracking).

---

## Class Summary

```matlab
% Create from raw files
obj = RecordedData.makeCSMS(filenameSpikes, filenameEvents, eventType, stateSeqs, matchOpt, otherMarkers)

% Load from MAT file
obj = RecordedData.loadData(filename)

% (Optional) merge conditions
obj = obj.mergeConditionsCSMS(cond2merge)

% Prepare binned tensors
obj = obj.prepareTensorWithFixBin(window, binWidth, neus2consider, conds2consider)
obj = obj.prepareTensorWithVariableBin(window, expectedBinWidth, neus2consider, conds2consider)

% (Optional) refine data
obj = obj.refineNeuralData(refType, refOpt, fixedFlag)

% (Optional) split tensor by task parameters
obj = obj.splitTensorWithParams(params, fixedBinFlag)

% Visualise
obj.plotPSTH(neu2consider, conds2consider, window, bin, vars2plot)
obj.plotSDF(neus2consider, conds2consider, window, binWidth, order, displayErrorFactor)

% Analyse
anovaResults  = obj.computeAnovaOnCond(neus2consider, conds2consider, refWindow, anaWindow, binWidth)
anovaResults  = obj.computeAnovaOnParam(neus2consider, conds2consider, params, refWindow, anaWindow, binWidth)
[scores, data, models] = obj.neuralDecodingClassification(processingOpt, trainOpt, testOpt, decoderOpt, visOpt)
[model, projTrain, projTest] = obj.computedPCA(dPCAOpt, processingOpt, trainOpt, testOpt, visOpt)
```

---

## Properties

| Property | Dimensions | Description |
|---|---|---|
| `CS` | nNeu x nCond | Spike times (cell: 1 x nTrials per entry) |
| `MS` | nNeu x nCond | State marker times (cell) |
| `KIN`, `EY`, `TC`, `TCTS` | — | Kinematics, eye tracking, touchscreen (optional) |
| `TRIAL_OUTCOME` | — | Trial outcome metadata |
| `TensorWithFixBin` | nNeu x nCond x nTrials x nTime | Fixed-bin trial activity |
| `TensorWithFixBinCondAvg` | nNeu x nCond x nTime | Fixed-bin condition average |
| `TensorWithFixBinInfo` | struct | Window, binWidth, preprocessing log |
| `MarkerTensorForFixBin` | nNeu x nCond x nTrials x nMarker | Marker bin positions |
| `TrialNumWithFixBin` | nNeu x nCond | Valid trial counts |
| `TensorWithVariableBin` | nNeu x nCond x nTrials x nTime | Variable-bin trial activity |
| `TensorWithVariableBinCondAvg` | nNeu x nCond x nTime | Variable-bin condition average |
| `TensorWithVariableBinInfo` | struct | Window, binWidth, preprocessing log |
| `MarkerTensorForVariableBin` | nNeu x nCond x nTrials x nMarker | Marker bin positions |
| `TrialNumWithVariableBin` | nNeu x nCond | Valid trial counts |
| `TensorWithParams` | nNeu x p1 x p2 x ... x nTrials x nTime | Parameter-split trial tensor |
| `TensorWithParamsCondAvg` | nNeu x p1 x p2 x ... x nTime | Parameter-split condition average |
| `TensorWithParamsInfo` | struct | Metadata for parameter-split tensor |

---

## Static Methods — Object Construction

### `RecordedData.makeCSMS`

Create object from raw spike and behavioral event files.

```matlab
obj = RecordedData.makeCSMS(filenameSpikes, filenameEvents, eventType, stateSeqs, matchOpt, otherMarkers)
```

| Argument | Type | Description |
|---|---|---|
| `filenameSpikes` | string | Path to .mat spike file |
| `filenameEvents` | string | Path to behavioral events CSV (`state`, `condID`, `currentTime`, optionally `xCursor`, `yCursor`, `errorType`) |
| `eventType` | string | `'per_frame_events'` or `'simple_events'` |
| `stateSeqs` | cell / matrix | State sequences to extract (cell for `'strict'`, matrix for `'first&last'`) |
| `matchOpt` | string | `'strict'` — exact sequence; `'first&last'` — first-to-last transition |
| `otherMarkers` | numeric / `[]` | Additional state indices to track as markers |

---

### `RecordedData.loadData`

Load object from a MAT file containing `CS`, `MS`, and optionally `KIN`, `EY`, `TC`, `TCTS`.

```matlab
obj = RecordedData.loadData(filename)
```

---

## Methods — Condition Management

### `mergeConditionsCSMS`

Merge original conditions into a smaller set. Clears all tensors (must re-run `prepareTensor` afterwards).

```matlab
obj = obj.mergeConditionsCSMS(cond2merge)
```

`cond2merge` — cell array; each cell is a vector of original condition indices to collapse.

```matlab
% 8 conditions -> 4
obj = obj.mergeConditionsCSMS({[1 2], [3 4], [5 6], [7 8]});
```

---

## Methods — Tensor Preparation & Refinement

### `prepareTensorWithFixBin`

Bin spike times into a fixed-width tensor. Firing rates in sp/s.

```matlab
obj = obj.prepareTensorWithFixBin(window, binWidth, neus2consider, conds2consider)
```

| Argument | Description |
|---|---|
| `window` | `[start_s, alignmentMarkerIdx, end_s]` |
| `binWidth` | `[binSize_s, stepSize_s]` — set `stepSize = 0` for non-overlapping bins |
| `neus2consider` | Neuron indices (numeric or logical vector) |
| `conds2consider` | Condition indices (numeric or logical vector) |

Updates: `TensorWithFixBin`, `TensorWithFixBinCondAvg`, `MarkerTensorForFixBin`, `TensorWithFixBinInfo`, `TrialNumWithFixBin`.

```matlab
% 50 ms non-overlapping bins
obj = obj.prepareTensorWithFixBin([-0.5 6 2], [0.05 0], 1:96, 1:8);

% 100 ms bins with 25 ms step
obj = obj.prepareTensorWithFixBin([-0.5 6 2], [0.1 0.025], 1:96, 1:8);
```

---

### `prepareTensorWithVariableBin`

Bin spike times with epoch-adaptive bin widths. Epochs shorter than 50% of `expectedBinWidth` are merged with the next epoch.

```matlab
obj = obj.prepareTensorWithVariableBin(window, expectedBinWidth, neus2consider, conds2consider)
```

| Argument | Description |
|---|---|
| `window` | `[startMarkerIdx, endMarkerIdx]` |
| `expectedBinWidth` | Target bin duration in seconds |
| `neus2consider` | Neuron indices |
| `conds2consider` | Condition indices |

Updates: `TensorWithVariableBin`, `TensorWithVariableBinCondAvg`, `MarkerTensorForVariableBin`, `TensorWithVariableBinInfo`, `TrialNumWithVariableBin`.

```matlab
obj = obj.prepareTensorWithVariableBin([3 8], 0.05, 1:96, 1:8);
```

---

### `refineNeuralData`

Apply normalization, baseline subtraction, or smoothing to the selected tensor pair (trial + condAvg). Each tensor is processed independently per neuron. Appends a record of the operation to the corresponding `Info` struct and warns if the same operation type is applied twice.

```matlab
obj = obj.refineNeuralData(refType, refOpt, fixedFlag)
```

`fixedFlag` — `true`: fix-bin tensors; `false`: variable-bin tensors.

**`refType = 'normalize'`**

| `refOpt.normType` | Additional fields | Description |
|---|---|---|
| `'z-score'` | — | Per-neuron z-score |
| `'range'` | `normRange = [mi ma]` | Scale to `[mi, ma]` |
| `'soft-normalize'` | `normConstant` (number or `'mean'`), `normConstantMultiplier` (default 1) | `B = A / (range(A) + c*m)` |

**`refType = 'subtract'`**

| `refOpt.value` | Description |
|---|---|
| `'mean'` | Per-neuron mean, computed from the object's own tensor (NaN-ignored) |
| `'median'` | Per-neuron median, computed from the object's own tensor (NaN-ignored) |
| scalar | Same constant subtracted from every neuron in both tensors |
| `(nNeu × 1)` vector | Pre-computed per-neuron values; entry `n` is subtracted from all entries of neuron `n`. Use this to apply externally computed statistics (e.g. a mean derived from a training set) consistently across multiple objects. |

**`refType = 'smooth'`**

`refOpt.smoothType` — `'movmean'`, `'movmedian'`, or `'gaussian'`.
`refOpt.smoothWindow` — integer window length in bins.
Uses `smoothdata` along the time dimension with `'omitnan'`.

```matlab
% Z-score fixed-bin tensors
opt.normType = 'z-score';
obj = obj.refineNeuralData('normalize', opt, true);

% Subtract per-neuron mean (variable-bin)
opt.value = 'mean';
obj = obj.refineNeuralData('subtract', opt, false);

% Subtract a pre-computed per-neuron mean vector (e.g. from a training set)
% meanVec is nNeu x 1, computed externally from trainObj
opt.value = meanVec;
trainObj = trainObj.refineNeuralData('subtract', opt, true);
testObj  = testObj.refineNeuralData('subtract', opt, true);

% Gaussian smoothing, 5-bin window (fixed-bin)
opt.smoothType = 'gaussian';  opt.smoothWindow = 5;
obj = obj.refineNeuralData('smooth', opt, true);
```

---

### `splitTensorWithParams`

Restructure the condition-averaged tensor into a multi-dimensional format for dPCA.

```matlab
obj = obj.splitTensorWithParams(params, fixedBinFlag)
```

`params` — `(nConds x nParams)` matrix; each column is one task parameter; rows correspond to the conditions selected during `prepareTensor`.
`fixedBinFlag` — `true`: source is `TensorWithFixBin`; `false`: `TensorWithVariableBin`.

Updates: `TensorWithParams`, `TensorWithParamsCondAvg`, `TrialNumWithParams`, `TensorWithParamsInfo`.

```matlab
params = [conditions.Target(:), conditions.Obstacle(:)];  % 8 conds x 2 params
obj = obj.splitTensorWithParams(params, true);
```

---

## Methods — Visualisation

### `plotPSTH`

Peri-stimulus time histogram with raster and optional behavioral overlays. Calls `prepareTensorWithFixBin` internally.

```matlab
obj.plotPSTH(neu2consider, conds2consider, window, bin, vars2plot)
obj.plotPSTH(..., 'markerNames', names)
```

| Argument | Description |
|---|---|
| `neu2consider` | Single neuron index |
| `conds2consider` | Condition indices |
| `window` | `[start_s, alignmentMarkerIdx, end_s]` |
| `bin` | Bin size in seconds |
| `vars2plot` | `{}`, `{'EY'}`, `{'TC'}`, or `{'EY','TC'}` |
| `'markerNames'` | Cell array of strings, one per marker line (optional) |

```matlab
obj.plotPSTH(1, 1:4, [-0.5 6 2], 0.05, {'EY'}, 'markerNames', {'start','cue','go','reward'});
```

---

### `plotSDF`

Spike density function across multiple neurons and conditions. Calls `prepareTensorWithFixBin` (and optionally `refineNeuralData` for Gaussian smoothing) internally on a temporary copy — the caller object is not modified.

```matlab
obj.plotSDF(neus2consider, conds2consider, window, binWidth, order, displayErrorFactor)
obj.plotSDF(..., 'markerNames', names)
obj.plotSDF(..., 'smoothWindow', win)
```

| Argument | Description |
|---|---|
| `neus2consider` | Neuron indices |
| `conds2consider` | Condition indices |
| `window` | `[start_s, alignmentMarkerIdx, end_s]` |
| `binWidth` | `[binSize_s, stepSize_s]` |
| `order` | `"preference"` — sort by peak FR; `"condition"` — original order |
| `displayErrorFactor` | Std multiplier for the error band |
| `'markerNames'` | Cell array of marker label strings (optional) |
| `'smoothWindow'` | Integer — Gaussian smoothing window in bins (optional) |

```matlab
obj.plotSDF(1:96, 1:8, [-0.5 6 2], [0.05 0.025], "preference", 0.1, ...
    'smoothWindow', 5, 'markerNames', {'start','cue','go','reward'});
```

---

## Methods — Analysis

### `computeAnovaOnCond`

Two-way ANOVA (epoch x condition) and one-way ANOVA (condition) across time bins.

```matlab
anovaResults = obj.computeAnovaOnCond(neus2consider, conds2consider, refWindow, anaWindow, binWidth)
anovaResults = obj.computeAnovaOnCond(..., 'markerNames', names)
```

`refWindow` / `anaWindow` — `[start_s, markerIdx, end_s]`.
`binWidth` — `[binSize_s, stepSize_s]`.

Output `anovaResults`:
- `.p2way` — `(nNeu x nBins x 3)` p-values: epoch, condition, epoch x condition
- `.p1way` — `(nNeu x nBins)` p-values: condition
- `.sign2way`, `.sign1way` — percentage of neurons with p < 0.05 over time

```matlab
anovaResults = obj.computeAnovaOnCond(1:96, 1:8, [-0.25 2 0.25], [-0.5 6 2], [0.5 0.05]);
```

---

### `computeAnovaOnParam`

Two-way ANOVA with epoch and task-parameter factors.

```matlab
anovaResults = obj.computeAnovaOnParam(neus2consider, conds2consider, params, refWindow, anaWindow, binWidth)
anovaResults = obj.computeAnovaOnParam(..., 'markerNames', names, 'factorNames', fnames, 'display', tf)
```

`params` — `(nConds x nParams)` parameter matrix.
Output: `.withEpoch` and `.withoutEpoch` structs (same p-value fields as `computeAnovaOnCond`).

---

### `neuralDecodingClassification`

Time-resolved decoding with cross-validation and optional train/test split across task contexts.

```matlab
[scores, data, models] = obj.neuralDecodingClassification(processingOpt, trainOpt, testOpt, decoderOpt, visOpt)
```

**`processingOpt`**

| Field | Description |
|---|---|
| `neus2consider` | Neuron indices |
| `fixedFlag` | `true` — fixed-bin; `false` — variable-bin |
| `binWidth` | `[binSize_s, stepSize_s]` (fixed-bin only) |
| `expectedBinWidth` | Target bin duration in s (variable-bin only) |

**`trainOpt` / `testOpt`**

| Field | Description |
|---|---|
| `conds2consider4train` / `conds2consider4test` | Condition indices |
| `trainParams` / `testParams` | `(nConds x nParams)` parameter matrix |
| `trainWindow` / `testWindow` | `[start_s, markerIdx, end_s]` |

**`decoderOpt`** (optional)

| Field | Default | Description |
|---|---|---|
| `decoder` | `'naivebayes'` | `'naivebayes'`, `'lda'`, or `'svm'` |
| `CVKfold` | `10` | Number of CV folds |
| `decoderOpts` | `struct()` | Options forwarded to the underlying `fit*` function |

**`visOpt`** (optional): `factorNames`, `markerNames`, and `display` — one of:
`'accuracy'`, `'precision'`, `'sensitivity'`, `'specificity'`, `'F-measure'`,
`'true_positive'`, `'false_positive'`, `'true_negative'`, `'false_negative'`, or `'off'` (default).

Output `scores`: `.train`, `.validation`, `.test`, `.timeAxisTrain`, `.timeAxisTest`, `.factorNames`, `.metricNames`, and chance-level equivalents.

```matlab
procOpt.fixedFlag = true;  procOpt.neus2consider = 1:96;
trOpt.conds2consider4train = 1:8;  trOpt.trainParams = cond.Target(1:8);  trOpt.trainWindow = [-0.5 6 2];
teOpt.conds2consider4test  = 1:8;  teOpt.testParams  = cond.Obstacle(1:8); teOpt.testWindow  = [-0.5 6 2];
decOpt.decoder = 'lda';  decOpt.CVKfold = 5;
visOpt.display = 'accuracy';  visOpt.factorNames = {'target'};
[scores, data, mdls] = obj.neuralDecodingClassification(procOpt, trOpt, teOpt, decOpt, visOpt);
```

---

### `computedPCA`

Demixed PCA: fit on training set, project train and test sets, plot explained variance.

```matlab
[dPCAmodel, dPCAprojTrain, dPCAprojTest] = obj.computedPCA(dPCAOpt, processingOpt, trainOpt, testOpt)
[dPCAmodel, dPCAprojTrain, dPCAprojTest] = obj.computedPCA(dPCAOpt, processingOpt, trainOpt, testOpt, visOpt)
```

**`dPCAOpt`**

| Field | Description |
|---|---|
| `combinedParams` | Cell of cell arrays — marginalization groupings (see `dpca.m`) |
| `numComps` | Number of dPCA components |
| `simultaneousRecording` | `true` if all neurons recorded simultaneously |
| `procedure` | `'simple'`, `'regularized'`, or `'regularized&noiseCovar'` |

**`processingOpt`**

| Field | Required | Description |
|---|---|---|
| `neus2consider` | yes | Neuron indices |
| `fixedFlag` | yes | `true` — fixed-bin; `false` — variable-bin |
| `binWidth` | fixed-bin | `[binSize_s, stepSize_s]` |
| `expectedBinWidth` | variable-bin | Target bin duration in s |
| `refType` | no | `'normalize'` or `'smooth'` — applied to both sets **before** mean-centering. Passing `'subtract'` raises an error. |
| `refOpt` | if `refType` set | Options struct forwarded verbatim to `refineNeuralData`. |

> **Internal pipeline order:** `prepareTensor*` → `refineNeuralData` (optional) → mean-centering (per-neuron, from concatenated train+test) → `splitTensorWithParams` → `dpca`.

**`trainOpt` / `testOpt`**

| Field | Description |
|---|---|
| `conds2consider4train` / `conds2consider4test` | Condition indices |
| `trainParams` / `testParams` | `(nConds x nParams)` parameter matrix |
| `trainWindow` / `testWindow` | `[start_s, markerIdx, end_s]` |

**`visOpt`** (optional): `marginalizationNames`, `marginalizationColours`, `legendSubplot`.

Output `dPCAmodel`: `.W` (decoder), `.V` (encoder), `.whichMarg`, `.Cnoise`, `.optimalLambda`, `.procedure`.
Output `dPCAprojTrain` / `dPCAprojTest`: `.explVar` struct from `dpca_explainedVariance`.

```matlab
dpOpt.combinedParams = {{1,[1 3]},{2,[2 3]},{3},{[1 2],[1 2 3]}};
dpOpt.numComps = 20;  dpOpt.simultaneousRecording = false;  dpOpt.procedure = 'regularized';

procOpt.fixedFlag = true;  procOpt.neus2consider = 1:96;  procOpt.binWidth = [0.05 0.025];

trOpt.conds2consider4train = 1:8;  trOpt.trainParams = cond.Target(1:8);  trOpt.trainWindow = [-0.5 6 2];
teOpt.conds2consider4test  = 1:8;  teOpt.testParams  = cond.Obstacle(1:8); teOpt.testWindow  = [-0.5 6 2];

vOpt.marginalizationNames = {'target','time','target x time'};
[mdl, projTr, projTe] = obj.computedPCA(dpOpt, procOpt, trOpt, teOpt, vOpt);
```

---

## Typical Workflow

```matlab
% 1. Build object
obj = RecordedData.makeCSMS(spikeFile, eventFile, 'per_frame_events', stateSeqs, 'strict', []);

% 2. (Optional) collapse conditions
obj = obj.mergeConditionsCSMS({[1 2 3 4], [5 6 7 8]});

% 3. Bin and optionally refine
obj = obj.prepareTensorWithFixBin([-0.5 6 2], [0.05 0.025], 1:96, 1:8);
opt.normType = 'soft-normalize';  opt.normConstant = 'mean';
obj = obj.refineNeuralData('normalize', opt, true);

% 4. Visualise
obj.plotSDF(1:96, 1:8, [-0.5 6 2], [0.05 0.025], "preference", 0.1, 'smoothWindow', 5);
obj.plotPSTH(1, 1:8, [-0.5 6 2], 0.05, {});

% 5. Statistics
anova = obj.computeAnovaOnCond(1:96, 1:8, [-0.25 2 0.25], [-0.5 6 2], [0.5 0.05]);

% 6. Decoding
procOpt.fixedFlag = true;  procOpt.neus2consider = 1:96;
trOpt.conds2consider4train = 1:8;  trOpt.trainParams = cond.Target(1:8);  trOpt.trainWindow = [-0.5 6 2];
teOpt.conds2consider4test  = 1:8;  teOpt.testParams  = cond.Obstacle(1:8); teOpt.testWindow  = [-0.5 6 2];
[scores, ~, ~] = obj.neuralDecodingClassification(procOpt, trOpt, teOpt);

% 7. dPCA
obj = obj.splitTensorWithParams([cond.Target(1:8), cond.Obstacle(1:8)], true);
dpOpt.combinedParams = {{1,[1 3]},{2,[2 3]},{3},{[1 2],[1 2 3]}};
dpOpt.numComps = 20;  dpOpt.simultaneousRecording = false;  dpOpt.procedure = 'simple';
[mdl, projTr, projTe] = obj.computedPCA(dpOpt, procOpt, trOpt, teOpt);
```
