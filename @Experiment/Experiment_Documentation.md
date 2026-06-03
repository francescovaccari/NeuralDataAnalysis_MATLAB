# Experiment Class Documentation

## Description

Load, process, and organize experimental data from touch screen studies. Handles multiple file formats (CSV, JSON) and supports data segmentation, metadata building, and metric computation.

---

## Syntax

```matlab
experiment = Experiment(folder, ID)
```

---

## Input Arguments

**`folder`** — Path to data folder  
char | string

Path to the folder containing data files.

**`ID`** — Experiment identifier  
char | string | numeric scalar

Experiment ID. If numeric, automatically prepends 'ID' prefix (e.g., `2002` → `'ID2002'`).

---

## Output Arguments

**`experiment`** — Experiment object  
Experiment

Object containing loaded and processed experimental data. All available data files are automatically loaded at construction.

---

## Properties

| Property | Type | Description |
|----------|------|-------------|
| `Folder` | string | Path to data folder |
| `ID` | string | Experiment ID (e.g., 'ID2002') |
| `targets` | struct array | Target definitions |
| `conditions` | table | Experimental conditions |
| `obstacles` | struct array | Obstacle definitions |
| `taskSettings` | struct | Task configuration |
| `rawData` | table | Raw sensor data |
| `fixPts` | struct array | Fixation point definitions |
| `segmentedData` | cell | Segmented trials by condition |
| `metaData` | cell | Metadata for each trial |
| `epochDuration` | cell | State duration per trial |

---

## Examples

### Basic Initialization

```matlab
% Load experiment from string ID
experiment = Experiment('/path/to/data', 'ID2002');

% Load experiment from numeric ID
experiment = Experiment('/path/to/data', 2002);  % Auto-converted to 'ID2002'

% Verify loaded data
disp(height(experiment.conditions));  % Number of conditions
disp(numel(experiment.targets));      % Number of targets
```

---

## Methods

### LoadFixationPoints

Load or create fixation point data.

**Syntax:**
```matlab
experiment = experiment.LoadFixationPoints()
experiment = experiment.LoadFixationPoints(defaultX, defaultY, defaultSize)
```

**Input Arguments:**
- `defaultX`, `defaultY`, `defaultSize` — Default fixation point values (numeric)

**Output:**
- `experiment` — Updated Experiment object

**Example:**
```matlab
% Load from file
experiment = experiment.LoadFixationPoints();

% Create with default values
experiment = experiment.LoadFixationPoints(0, -100, 100);
```

---

### ComputeFixPtTargetDistance

Compute distance between fixation points and targets.

**Syntax:**
```matlab
experiment = experiment.ComputeFixPtTargetDistance()
experiment = experiment.ComputeFixPtTargetDistance(defaultX, defaultY, defaultSize)
```

**Description:**
Calculates Euclidean distance from each fixation point to its associated target. Results stored in `conditions.FixPtTargetDistance`.

**Output:**
- `experiment` — Updated with computed distances

**Example:**
```matlab
experiment = experiment.ComputeFixPtTargetDistance();
disp(experiment.conditions.FixPtTargetDistance);
```

---

### ComputeDifficulty

Compute task difficulty using Fitts' Law.

**Syntax:**
```matlab
experiment = experiment.ComputeDifficulty(method)
```

**Input Arguments:**
- `method` — "ClassicFittsLaw" or "AdjustedFittsLaw" (char | string)
  - Classic: ID = log₂(2D/W)
  - Adjusted: ID = log₂(D/W + 1)

**Requirements:**
- Must call `ComputeFixPtTargetDistance` first
- All target sizes must be positive

**Example:**
```matlab
experiment = experiment.ComputeFixPtTargetDistance();
experiment = experiment.ComputeDifficulty('ClassicFittsLaw');
disp(experiment.conditions.Difficulty);
```

---

### SegmentData

Segment raw data into trials based on state sequences.

**Syntax:**
```matlab
experiment = experiment.SegmentData(stateSequence, strictOpt)
```

**Input Arguments:**
- `stateSequence` — State sequence marking trial boundaries (numeric vector)
- `strictOpt` — Matching mode: "strict" | "first&last" (char | string)
  - `"strict"`: Exact complete state sequence
  - `"first&last"`: First and last states, intermediate variation allowed

**Output:**
- `experiment` — Updated with `segmentedData{condition}{trial}` structure

**Example:**
```matlab
% Segment by complete state sequence
experiment = experiment.SegmentData([1 2 5], 'strict');

% Access first trial of first condition
firstTrial = experiment.segmentedData{1}{1};
```

---

### BuildMetaData

Build metadata structure from segmented data.

**Syntax:**
```matlab
experiment = experiment.BuildMetaData()
```

**Description:**
Creates metadata structures containing trial information including target details and condition parameters. Results in `metaData{condition}{trial}` cell array.

**Example:**
```matlab
experiment = experiment.BuildMetaData();
metadata = experiment.metaData{1}{1};  % First trial, first condition
disp(metadata.target6);
```

---

### ComputeEpochDuration

Compute (combined) duration of specific states per trial.

**Syntax:**
```matlab
experiment = experiment.ComputeEpochDuration(states2comp)
```

**Input Arguments:**
- `states2comp` — States to compute duration for (numeric vector)

**Output:**
- `experiment` — Updated with `epochDuration{condition}` cell array containing total duration (sum of all specified states)

**Description:**
Calculates the total duration when any of the specified states were active. All states in `states2comp` are summed together into a single duration value per trial.

**Example:**
```matlab
% Compute combined duration for states 2 and 3 together
experiment = experiment.ComputeEpochDuration([2 3]);
dur = experiment.epochDuration{1};  % Total duration of states 2 AND 3 for each trial in first condition
```

---

### PlotFixPtsTargets

Visualize fixation points and targets.

**Syntax:**
```matlab
fig = experiment.PlotFixPtsTargets(conds2plot)
```

**Input Arguments:**
- `conds2plot` — Condition IDs to plot (numeric vector) or "all" (char | string)

**Output:**
- `fig` — Figure handle

**Plot Details:**
- Blue circles: Targets (labeled T{ID})
- Red diamonds: Fixation points (labeled F{ID})

**Example:**
```matlab
fig = experiment.PlotFixPtsTargets([1 2 3]);
fig = experiment.PlotFixPtsTargets("all");
```

---

### Complete Workflow

```matlab
clearvars; close all; clc;

% Initialize
experiment = Experiment(pwd, 2002);
experiment = experiment.LoadFixationPoints(0, -100, 100);

% Compute metrics
experiment = experiment.ComputeFixPtTargetDistance();
experiment = experiment.ComputeDifficulty("ClassicFittsLaw");

% Process data
stateSequence = [experiment.taskSettings.stateSequence(:); 100];
experiment = experiment.SegmentData(stateSequence, "first&last");
experiment = experiment.BuildMetaData();
experiment = experiment.ComputeEpochDuration([6 7]);

% Visualize
nTrialsByCond = cellfun(@numel, experiment.segmentedData);
fprintf('ID: %s | Conditions: %d | Targets: %d | Trials: %d\n', ...
    experiment.ID, height(experiment.conditions), ...
    numel(experiment.targets), sum(nTrialsByCond));

fig = experiment.PlotFixPtsTargets("all");
```

---

## File Format Requirements

**CSV Files (Semicolon-Delimited)**
- Header row with column names
- Semicolon-separated values
- Decimal separator: `.`

**JSON Files**
- Standard JSON format
- Top-level keys: `"targets"`, `"conditions"`, `"obstacles"`, `"taskSettings"`, `"fixationPoints"`
