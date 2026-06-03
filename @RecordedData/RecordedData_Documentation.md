# RecordedData Class Documentation

## Description

Load and analyze neural recording data (spike times, state markers, kinematic data, eye tracking). Provides comprehensive tools for neural data binning, visualization, decoding, statistical analysis, and dimensionality reduction via dPCA.

---

## Quick Start

```matlab
% Load neural recording data from MAT file
recordedData = RecordedData.loadData('/path/to/recording.mat');

% OR create from OE and CSV files
recordedData = RecordedData.makeCSMS(filenameSpikes, filenameEvents, ...
    'per_frame_events', stateSeqs, 'strict', otherMarkers);

% Prepare tensor with fixed-width binning
recordedData = recordedData.prepareTensorWithFixBin([-2 6 +2], [0.05 0], 0);

% Visualize spike activity
recordedData.plotPSTH(1, [1 2 3], [-2 6 +2], 0.05, {'EY'});

% Decode conditions
scores = recordedData.decodeCond([1:10], [1:4], [-0.5 6 0.5]);
```

---

## Syntax

```matlab
% Create object (empty)
recordedData = RecordedData();

% Load from MAT file
recordedData = RecordedData.loadData(filename);

% Create from raw OE/CSV data
recordedData = RecordedData.makeCSMS(filenameSpikes, filenameEvents, ...
    eventType, stateSeqs, matchOpt, otherMarkers);
```

---

## Input Arguments (Static Methods)

### loadData

**`filename`** — Path to MAT file  
string | char

Path to a `.mat` file containing neural recording data. Expected variables: `CS`, `MS`, and optionally `KIN`, `EY`, `TC`, `TCTS`. Variables not found are set to empty.

### makeCSMS

**`filenameSpikes`** — Path to spike data file  
string | char

File containing spike times, typically a `.mat` file with spike data from spike sorting results.

**`filenameEvents`** — Path to behavioral events CSV file  
string | char

CSV file containing state transitions and behavioral markers with columns: `state`, `condID`, `currentTime`, optionally `xCursor`, `yCursor`, `errorType`.

**`eventType`** — Type of behavioral event encoding  
'per_frame_events' | 'simple_events'

Specifies how states are encoded: 'per_frame_events' for frame-by-frame state sequences.

**`stateSeqs`** — State sequence patterns to extract  
cell array | numeric array

Cell array of state sequences to match (for 'strict' matching) or matrix of start/end states (for 'first&last').

**`matchOpt`** — Matching strategy  
'strict' | 'first&last'

'strict': match exact state sequences; 'first&last': match trials starting with first state and ending with last state.

**`otherMarkers`** — Additional state markers to track  
numeric vector | []

Indices of states to extract as event markers beyond standard trial markers.

---

## Output Arguments

**`recordedData`** — RecordedData object

Object containing raw and pre-processed neural recording data with properties and methods for analysis.

---

## Properties

| Property | Type | Description |
|----------|------|-------------|
| **Data** | | |
| `Filename` | string | Path to loaded MAT file |
| `CS` | cell array | Spike times (neurons × conditions × trials) |
| `MS` | cell array | State marker times (neurons × conditions × trials) |
| `KIN` | cell array | Kinematic data (optional) |
| `EY` | cell array | Eye tracking data (optional) |
| `TC` | cell array | Touchscreen x,y coordinates (optional) |
| `TCTS` | cell array | Touchscreen coordinate timestamps (optional) |
| `TRIAL_OUTCOME` | cell array | Trial outcome metadata (optional) |
| **Preprocessed Data** | | |
| `TensorWithVariableBin` | numeric array | (neurons × conditions × trials × time) with adaptive binning |
| `TensorWithVariableBinCondAvg` | numeric array | (neurons × conditions × time) condition-averaged |
| `TensorWithFixBin` | numeric array | (neurons × conditions × trials × time) with fixed bins |
| `TensorWithFixBinCondAvg` | numeric array | (neurons × conditions × time) condition-averaged |
| `Tensor4dPCA` | numeric array | Multi-dimensional tensor for dPCA |
| **Metadata** | | |
| `TensorWithVariableBinInfo` | struct | Parameters and info for variable-bin tensor |
| `TensorWithFixBinInfo` | struct | Parameters and info for fixed-bin tensor |
| `Tensor4dPCAInfo` | struct | Parameters and info for dPCA tensor |
| `MarkerTensorForVariableBin` | numeric array | State marker positions in variable bins |
| `MarkerTensorForFixBin` | numeric array | State marker positions in fixed bins |

---

## Examples

### Example 1: Load and Plot PSTH

```matlab
% Load neural recording data
recordedData = RecordedData.loadData('/path/to/recording.mat');

% Create peri-stimulus time histogram for neuron 1, conditions 1-3
recordedData.plotPSTH(1, [1 2 3], [-2 6 +2], 0.05, {});
```

### Example 2: Prepare Tensors and Plot SDF

```matlab
% Prepare fixed-width binning: 50ms bins, no overlap
recordedData = recordedData.prepareTensorWithFixBin([-2 6 +2], [0.05 0], 0);

% Plot spike density function for neurons 1-5 across conditions 1-3
recordedData.plotSDF([1:5], [1:3], "condition", 0.5);
```

### Example 3: Condition Decoding

```matlab
% Decode which condition was presented using neural activity
scores = recordedData.decodeCond([1:10], [1:4], [-0.5 6 0.5]);
% scores contains confusion matrix and performance metrics
```

### Example 4: Merge Conditions and Analyze

```matlab
% Merge original conditions into new condition groups
recordedData = recordedData.mergeConditionsCSMS({[1 2], [3 4 5]});

% Now CS and MS have 2 conditions instead of 5
```

### Example 5: dPCA Analysis

```matlab
% Prepare variable-bin tensor
recordedData = recordedData.prepareTensorWithVariableBin([0 5], 0.05, 0);

% Prepare multi-dimensional tensor for dPCA
recordedData = recordedData.prepareTensor4dPCA({'Target', 'Obstacle'}, ...
    experiment.conditions, false, [1:20], [1:8]);

% Compute dPCA
dPCA_results = recordedData.computedPCA({[1 2]}, 10);
```

---

## Methods

### prepareTensorWithVariableBin

Prepare neural activity tensor with epoch-adaptive binning.

**Syntax:**
```matlab
recordedData = recordedData.prepareTensorWithVariableBin(window, expectedBinWidth, smoothGaussianStdInBin)
```

**Input Arguments:**
- `window` — Time window [initial_marker (index), final_marker (index)] (numeric vector)
  - First element: start marker index
  - Second element: end marker index
- `expectedBinWidth` — Target bin duration in seconds (numeric scalar)
- `smoothGaussianStdInBin` — Gaussian smoothing standard deviation in bins (numeric scalar ≥ 0)

**Output:**
- `recordedData` — Updated with:
  - `TensorWithVariableBin`: (neurons × conditions × trials × time) with trial-by-trial activity
  - `TensorWithVariableBinCondAvg`: (neurons × conditions × time) condition-averaged activity

**Description:**
Creates binned neural activity where bin widths adapt to epoch durations. Epochs are divided into bins based on expected bin width. If a bin becomes too small (< 50% of expected width), it is merged with the next epoch. Gaussian smoothing is optional for temporal smoothing.

**Remarks:**
- Automatically computes bin widths for each epoch based on marker timing
- Stores real bin widths in `TensorWithVariableBinInfo.realBinWidth`
- Uses progress bar for monitoring
- Firing rates are in spikes/second

**Example:**
```matlab
% Prepare tensor with 50ms target bins, no smoothing
recordedData = recordedData.prepareTensorWithVariableBin([0 5], 0.05, 0);

% Check output
sizeCheck = size(recordedData.TensorWithVariableBinCondAvg);  % (neurons, conditions, time)

% With Gaussian smoothing (25 bin std)
recordedData = recordedData.prepareTensorWithVariableBin([0 5], 0.05, 25);
```

---

### prepareTensorWithFixBin

Prepare neural activity tensor with fixed-width binning.

**Syntax:**
```matlab
recordedData = recordedData.prepareTensorWithFixBin(window, binWidth, smoothGaussianStdInBin)
```

**Input Arguments:**
- `window` — Time interval [start (in s), alignment_marker (index), end (in s)] (numeric vector with 3 elements)
- `binWidth` — [bin_width, overlap] in seconds (numeric vector)
  - `binWidth(1)`: size of each bin
  - `binWidth(2) > 0`: overlap width (for moving bins)
  - `binWidth(2) = 0`: non-overlapping bins
- `smoothGaussianStdInBin` — Gaussian smoothing standard deviation in bins (numeric scalar ≥ 0)

**Output:**
- `recordedData` — Updated with:
  - `TensorWithFixBin`: (neurons × conditions × trials × time) with trial-by-trial activity
  - `TensorWithFixBinCondAvg`: (neurons × conditions × time) condition-averaged activity

**Description:**
Creates binned neural activity with uniform bin widths across entire time window. Supports overlapping bins for intrinsic temporal smoothing. State markers are tracked in `MarkerTensorForFixBin`.

**Remarks:**
- Firing rates are in spikes/second
- Overlapping bins are computed using moving sum
- Uses progress bar for monitoring
- State markers are detected and stored for reference

**Example:**
```matlab
% 100ms bins, no overlap, 2s before to 2s after the 5th marker
recordedData = recordedData.prepareTensorWithFixBin([-2 5 +2], [0.1 0], 0);

% 100ms bins with 50ms overlap
recordedData = recordedData.prepareTensorWithFixBin([-2 6 +2], [0.1 0.05], 0);

% With Gaussian smoothing (5 bin std)
recordedData = recordedData.prepareTensorWithFixBin([-2 5 +2], [0.05 0], 5);
```

---

### plotPSTH

Plot peri-stimulus time histogram with raster plot and optional behavioral variables.

**Syntax:**
```matlab
recordedData.plotPSTH(neu2consider, conds2consider, window, bin, vars2plot)
```

**Input Arguments:**
- `neu2consider` — Single neuron index to plot (numeric scalar, positive integer)
- `conds2consider` — Vector of condition indices to plot (numeric vector)
- `window` — [onset, marker, offset] where:
  - `onset`: start time relative to alignment marker (seconds)
  - `marker`: marker index for alignment
  - `offset`: end time relative to alignment marker (seconds)
- `bin` — Bin size in seconds for histogram (numeric scalar, positive)
- `vars2plot` — Cell array of variable names to plot below raster (optional)
  - Examples: `{'EY'}`, `{'TC'}`, `{'EY', 'TC'}`
  - Leave empty `{}` for no additional variables

**Output:**
- Figure with:
  - Histogram (top): Mean firing rate per condition with error shading
  - Raster plot (middle): Individual spike times and state markers by trial
  - Variable plots (bottom, optional): Behavioral variables (eye position, touchscreen)

**Description:**
Comprehensive visualization of neural activity aligned to a specific marker. Displays mean firing rate histogram with raster plot showing all trials. Optional behavioral variables (eye tracking, touchscreen position) can be overlaid below.

**Remarks:**
- Grid layout optimized based on number of conditions
- Consistent y-axis scaling across all conditions (global max)
- State markers color-coded by marker index (HSV colormap)
- Vertical alignment line at marker time (t=0, dashed)
- Automatic tick spacing based on number of conditions

**Example:**
```matlab
% Plot neuron 1 for conditions 1, 2, 3
recordedData.plotPSTH(1, [1 2 3], [-2 6 +2], 0.05, {});

% With eye tracking variable
recordedData.plotPSTH(1, [1 2 3], [-2 6 +2], 0.05, {'EY'});

% With touchscreen coordinates
recordedData.plotPSTH(5, [1:4], [-1.5 3 +1.5], 0.02, {'TC'});

% With multiple variables
recordedData.plotPSTH(10, [1 2], [-2 5 +2], 0.05, {'EY', 'TC'});
```

---

### plotSDF

Plot spike density function with condition comparison and state markers.

**Syntax:**
```matlab
recordedData.plotSDF(neus2consider, conds2consider, order, displayErrorFactor)
```

**Input Arguments:**
- `neus2consider` — Neuron indices (numeric vector)
- `conds2consider` — Condition indices (numeric vector)
- `order` — Sorting order for display (string)
  - `"preference"`: sort by neural response magnitude (peak firing)
  - `"condition"`: keep original condition order
- `displayErrorFactor` — Error bar scaling factor (numeric scalar)
  - Error bars displayed as: `displayErrorFactor × standard_deviation`

**Output:**
- Figure with mean spike density function for each condition with error shading and state markers

**Description:**
Visualizes spike density function (mean firing rate) across multiple neurons and conditions. Can sort conditions by neural preference (highest response first) or keep original condition order. Displays mean ± error with state markers overlaid.

**Requirements:**
- Must call `prepareTensorWithFixBin` before using this method

**Remarks:**
- Uses fixed-bin tensor (`TensorWithFixBinCondAvg`)
- Color-coded by condition (jet colormap)
- State marker positions extracted from `MarkerTensorForFixBin`
- If multiple neurons selected, mean across neurons is plotted

**Example:**
```matlab
% Prepare fixed-bin tensor first
recordedData = recordedData.prepareTensorWithFixBin([-2 6 +2], [0.001 0], 25);

% Plot neurons 1-5 for conditions 1-3, ordered by preference
recordedData.plotSDF([1:5], [1:3], "preference", 0.5);

% Same but keep original condition order
recordedData.plotSDF([1:5], [1:3], "condition", 1);

% Single neuron, all conditions
recordedData.plotSDF(7, [1:4], "preference", 0.3);
```

---

### decodeCond

Decode experimental conditions from neural activity using Naïve Bayes classifier.

**Syntax:**
```matlab
decodingScores = recordedData.decodeCond(neus2consider, conds2consider, window)
```

**Input Arguments:**
- `neus2consider` — Neuron indices (numeric vector)
- `conds2consider` — Condition indices (numeric vector)
- `window` — [start (in s), marker_index, end (in s)] relative to marker (numeric vector)

**Output:**
- `decodingScores` — Struct with decoding performance metrics:
  - Confusion matrix displayed as figure
  - Struct contains accuracy, sensitivity, specificity, etc. (from `statsOfMeasure`)

**Description:**
Performs 10-fold cross-validation classification using Naïve Bayes decoder. Spike counts within the specified window are used as features. Displays confusion matrix and returns performance metrics.

**Remarks:**
- Uses 10-fold cross-validation
- Adds small random noise (10^-8) to avoid zero variance
- Spike counts converted to binary classification per trial
- Results depend on condition separation in neural space

**Example:**
```matlab
% Decode conditions 1-4 using neurons 1-10
% Window: 0.5s before to 0.5s after 6th marker
scores = recordedData.decodeCond([1:10], [1:4], [-0.5 6 0.5]);

% Can also decode subset
scores = recordedData.decodeCond([5 10 15], [2 3], [-1 3 +1]);
```

---

### computeContAnova

Compute ANOVA comparing baseline and analysis windows across conditions.

**Syntax:**
```matlab
anovaResults = recordedData.computeContAnova(neus2consider, conds2consider, refWindow, anaWindow, binWidth)
```

**Input Arguments:**
- `neus2consider` — Neuron indices (numeric vector)
- `conds2consider` — Condition indices (numeric vector)
- `refWindow` — Reference (baseline) interval [start, marker_index, end] (numeric vector)
- `anaWindow` — Analysis interval [start, marker_index, end] (numeric vector)
- `binWidth` — [bin_width, step] in seconds (numeric vector)
  - `binWidth(1)`: bin size
  - `binWidth(2)`: step between bins (overlap if < binWidth(1))

**Output:**
- `anovaResults` — Struct containing:
  - `p2way`: (neurons × bins × 3) p-values [epoch, condition, epoch×condition]
  - `p1way`: (neurons × bins) p-values [condition]
  - `sign2way`: (bins × 3) percentage of neurons with p < 0.05 (2-way)
  - `sign1way`: (bins) percentage of neurons with p < 0.05 (1-way)

**Description:**
Two-way ANOVA (epoch × condition) comparing baseline vs. analysis window. One-way ANOVA (condition only) within analysis window. Plots percentage of modulated units across time bins.

**Remarks:**
- Note: binWidth for analysis should typically match reference window size
- Generates two figures: one for 2-way ANOVA, one for 1-way ANOVA
- Percentage curves show temporal evolution of condition modulation

**Example:**
```matlab
% Compare 500ms baseline (centered on 2nd marker) 
% with 2-second analysis (centered on 6th marker)
% in 500ms bins with 50ms steps
anovaResults = recordedData.computeContAnova([1:20], [1:4], ...
    [-0.25 2 +0.25], [-2 6 +2], [0.5 0.05]);
```

---

### prepareTensor4dPCA

Prepare multi-dimensional tensor for dPCA analysis with parameter decomposition.

**Syntax:**
```matlab
recordedData = recordedData.prepareTensor4dPCA(paramNames, conditions, fixedBinFlag, neus2consider, conds2consider)
```

**Input Arguments:**
- `paramNames` — Parameter names to extract (cell array of strings)
  - Example: `{'Target6', 'Obstacle6'}`
  - Each unique value in the condition table for each parameter creates a dimension
- `conditions` — Condition table/struct with fields matching `paramNames`
  - Each row corresponds to a condition
- `fixedBinFlag` — Use fixed or variable bins (logical)
  - `true`: uses `TensorWithFixBinCondAvg` (requires prior `prepareTensorWithFixBin`)
  - `false`: uses `TensorWithVariableBinCondAvg` (requires prior `prepareTensorWithVariableBin`)
- `neus2consider` — Neuron indices (numeric vector)
- `conds2consider` — Condition indices (numeric vector)

**Output:**
- `recordedData` — Updated with:
  - `Tensor4dPCA`: multi-dimensional array (neurons × param1 × param2 × ... × time)
  - `Tensor4dPCAInfo`: struct with metadata (paramNames, params, bin info)

**Description:**
Restructures condition-averaged tensor into multi-dimensional format where each parameter gets its own dimension. This allows dPCA to decompose neural activity by parameter.

**Remarks:**
- Requires prior tensor preparation (`prepareTensorWithFixBin` or `prepareTensorWithVariableBin`)
- Parameters must be numeric or categorical
- Supports arbitrary number of parameters
- Output tensor dimensions: (neurons, unique_param1_values, unique_param2_values, ..., time_bins)

**Example:**
```matlab
% First prepare variable-bin tensor
recordedData = recordedData.prepareTensorWithVariableBin([0 5], 0.05, 0);

% Create experiment conditions
experiment.conditions.Target6 = [1 1 1 1 2 2 2 2];
experiment.conditions.Obstacle6 = [1 2 1 2 1 2 1 2];

% Prepare 4D tensor
recordedData = recordedData.prepareTensor4dPCA(...
    {'Target6', 'Obstacle6'}, ...
    experiment.conditions, ...
    false,  % use variable bins
    [1:20], % neurons 1-20
    [1:8]); % all 8 conditions
```

---

### computedPCA

Compute demixed principal component analysis (dPCA) for parametric decomposition.

**Syntax:**
```matlab
dPCAResults = recordedData.computedPCA(combinedParams, numComps, options)
```

**Input Arguments:**
- `combinedParams` — Parameter combinations for dPCA (cell array)
  - Example: `{[1 2]}` combines dimensions 1 and 2
  - Example: `{[1], [2], [1 2]}` three components: marginal for 1, for 2, combined
- `numComps` — Number of components to compute (numeric scalar)
- `options` — Optional analysis options (struct, optional)
  - `margColours`: colormap for parameter marginalizations
  - `margNames`: labels for parameter dimensions
  - `timeMarginalization`: which dimension to marginalize for time (default: last)
  - `time`: time axis values (numeric vector)
  - `timeEvents`: marker positions on time axis (numeric vector)

**Output:**
- `dPCAResults` — Struct containing:
  - `W`: decoder matrix weights (time × components)
  - `V`: encoder matrix weights (features × components)
  - `explVar`: explained variance per component
  - `whichMarg`: component marginalization assignments

**Description:**
Performs demixed PCA to decompose neural population activity into components that selectively represent each parameter. Generates visual plots of component trajectories in PC space.

**Requirements:**
- Must call `prepareTensor4dPCA` before using this method
- dPCA toolbox must be in path: `C:\Program Files\MATLAB\SuppFunctions\dPCA-master`

**Remarks:**
- Automatically detects time axis and marker positions from `Tensor4dPCAInfo`
- Plots show component trajectories for each parameter marginalization
- Explained variance indicates component importance

**Example:**
```matlab
% Prepare tensor first
recordedData = recordedData.prepareTensorWithVariableBin([0 5], 0.05, 0);
recordedData = recordedData.prepareTensor4dPCA({'Target', 'Obstacle'}, ...
    cond, false, [1:20], [1:8]);

% Compute dPCA with combined parameters
dPCA_results = recordedData.computedPCA({[1 2]}, 10);

% Or with more detailed options
options = struct();
options.margNames = {'Target', 'Obstacle', 'Mixed'};
options.margColours = [1 0 0; 0 1 0; 0 0 1];
dPCA_results = recordedData.computedPCA({[1], [2], [1 2]}, 15, options);
```

---

### mergeConditionsCSMS

Merge multiple original conditions into fewer collapsed conditions.

**Syntax:**
```matlab
recordedData = recordedData.mergeConditionsCSMS(cond2merge)
```

**Input Arguments:**
- `cond2merge` — Cell array specifying condition merging (cell array)
  - Each cell contains vector of original condition indices to merge
  - Example: `{[1 2], [3 4 5], [6], [7 8]}`
  - Result: 4 new conditions (merge of 1-2, 3-5, 6 alone, 7-8)

**Output:**
- `recordedData` — Modified object with:
  - `CS`, `MS`, `TC`, `TCTS`, `TRIAL_OUTCOME` reorganized into new condition structure
  - Trials sorted by onset time within each new condition

**Description:**
Restructures neural data by merging multiple conditions into fewer conditions. All trial data from merged conditions are concatenated and sorted by trial onset time. Automatically clears computed tensors since dimensions have changed.

**Remarks:**
- Input validation ensures:
  - All condition indices are within valid range
  - No duplicate indices across merged groups
  - Warning if some original conditions not included
- Trials sorted by first marker (onset) within each merged condition
- All tensor data cleared (must recompute: `prepareTensorWithFixBin`, etc.)

**Example:**
```matlab
% Original data: 8 conditions
% Merge into 4 new conditions
recordedData = recordedData.mergeConditionsCSMS({[1 2], [3 4], [5 6], [7 8]});

% Or merge into 2 conditions
recordedData = recordedData.mergeConditionsCSMS({[1:4], [5:8]});

% Keep one condition separate
recordedData = recordedData.mergeConditionsCSMS({[1 2 3], [4], [5 6 7 8]});
```

---

### loadData (Static)

Load neural recording data from a MAT file.

**Syntax:**
```matlab
recordedData = RecordedData.loadData(filename)
```

**Input Arguments:**
- `filename` — Path to MAT file (string | char)

**Output:**
- `recordedData` — RecordedData object with loaded data

**Description:**
Static method to create and populate a RecordedData object from a MAT file. Expected variables: `CS` (spike times), `MS` (state markers), and optionally `KIN`, `EY`, `TC`, `TCTS`. Missing variables are set to empty.

**Example:**
```matlab
% Load from file
recordedData = RecordedData.loadData('/path/to/mydata.mat');

% Check what was loaded
disp(size(recordedData.CS));   % neurons × conditions
```

---

### makeCSMS (Static)

Create RecordedData from spike times and behavioral event files.

**Syntax:**
```matlab
recordedData = RecordedData.makeCSMS(filenameSpikes, filenameEvents, eventType, stateSeqs, matchOpt, otherMarkers)
```

**Input Arguments:**
- `filenameSpikes` — Path to spike data file (string | char)
- `filenameEvents` — Path to behavioral events CSV file (string | char)
- `eventType` — Event encoding type: `'per_frame_events'` | `'simple_events'`
- `stateSeqs` — State sequences to extract (cell array for 'strict', matrix for 'first&last')
- `matchOpt` — Matching strategy: `'strict'` | `'first&last'`
- `otherMarkers` — Additional state indices to track (numeric vector | [])

**Output:**
- `recordedData` — RecordedData object with CS, MS, TC, TCTS, TRIAL_OUTCOME populated

**Description:**
Static method to create RecordedData from raw spike data and behavioral event files. Extracts trials matching specified state sequences and segments spike times and behavioral variables accordingly.

**Remarks:**
- CSV file expected columns: `state`, `condID`, `currentTime`, optionally `xCursor`, `yCursor`, `errorType`
- Supports two matching strategies:
  - `'strict'`: match exact state sequences
  - `'first&last'`: match first→last state transitions
- Automatically segments touchscreen data (TC) if available
- Tracks trial outcomes (error types, final states)

**Example:**
```matlab
% Create from raw data
stateSeqs = {[1 2 3 4], [1 2 4]};  % Two trial types
recordedData = RecordedData.makeCSMS(...
    '/path/to/spikes.mat', ...
    '/path/to/events.csv', ...
    'per_frame_events', ...
    stateSeqs, ...
    'strict', ...
    [2 3]);  % track states 2 and 3 as markers
```

---

## Tips and Best Practices

1. **Always prepare tensors before analysis**: Use `prepareTensorWithFixBin` or `prepareTensorWithVariableBin` before calling visualization or analysis methods.

2. **Choose binning strategy**:
   - **Variable binning**: Use for epoch-aligned analysis where bins match behavioral periods
   - **Fixed binning**: Use for continuous analysis or when uniform temporal resolution is needed

3. **Smoothing**: Add Gaussian smoothing (5-25 bins typical) to reduce noise in plots without losing temporal resolution

4. **Window specifications**:
   - Time values in seconds (negative for before marker)
   - Marker index starts at 1 (first marker in each trial)
   - Example: `[-1 3 +1]` means 1s before to 1s after marker 3

5. **dPCA requirements**:
   - Design experiments with multiple parameters (e.g., target location, obstacle presence)
   - Ensure balanced factorial combinations of parameters
   - May need 20+ neurons for meaningful results

6. **Merging conditions**: Use to simplify analysis (e.g., group similar trial types), but remember this loses detailed condition structure


---

## Complete Workflow

```matlab
clearvars; close all; clc;

% Load recording
recordedData = RecordedData('/path/to/neural_data.mat');

% Prepare fixed-bin tensor with specific parameters
window = [0 2 5];                       % 0 to 5 seconds relative to marker 2
binWidth = [0.1 0];                     % 100ms bins, no overlap
gaussianWidth = 1;                      % Gaussian smoothing
recordedData = recordedData.prepareTensorWithFixBin(window, binWidth, gaussianWidth);

% Plot spike density
recordedData.plotSDF([1:10], [1:4], "condition", 1.96);

% Decode conditions (requires spike counts, no tensor needed)
scores = recordedData.decodeCond([1:20], [1:8], [0 2 2]);

% ANOVA analysis (prepares its own tensors)
recordedData.computeContAnova([1:20], [1:8], [-2 2 0], [0 2 3], 0.1);

% dPCA analysis
recordedData = recordedData.prepareTensor4dPCA({'Target6', 'Obstacle6'}, ...
    experiment.conditions, true, [1:20], [1:8]);
dPCA_res = recordedData.computedPCA({[1 2]}, 10);
```

---

## Data Format Requirements

**MAT File Structure**

Must contain four variables as cells organized by (neuron, condition):

| Variable | Content | Format |
|----------|---------|--------|
| `CS` | Spike times | (n_neurons × n_conditions) cell, each cell: {trial} = spike times (seconds) |
| `MS` | State markers | (n_neurons × n_conditions) cell, each cell: {trial} = marker times (seconds) |
| `KIN` | Kinematics | (n_neurons × n_conditions) cell or [] if unavailable |
| `EY` | Eye tracking | (n_neurons × n_conditions) cell or [] if unavailable |

**Example CS structure:**
```matlab
CS{1,1}{1} = [0.12 0.34 0.56 0.89 ...]  % Neuron 1, Condition 1, Trial 1: spike times
CS{2,1}{1} = [0.15 0.37 0.59 ...]       % Neuron 2, Condition 1, Trial 1
```
