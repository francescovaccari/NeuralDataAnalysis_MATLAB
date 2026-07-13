clear all;   % Clear all variables from workspace
close all;   % Close all open figure windows

% Add the analysis scripts folder (and subfolders) to the MATLAB path
addpath(genpath("C:\Users\frada\Desktop\WORK in PROGRESS\LAB\TouchScreen_project\Matlab Script\Analisi dati"));

% Root folder containing the experimental datasets
dataFolder = "C:\Users\frada\Desktop\WORK in PROGRESS\LAB\TouchScreen_project\Data\Datasets";

% Load pre-processed neural data (spike counts / continuous signals) for subject ID2055
% The .mat file contains sorted/unsorted neural recordings with trial-aligned data
recordedData = RecordedData.loadData(fullfile(dataFolder, "ID2055_CSMS_chs128_192_ave_bpass300_6.0k_Int16_r_4_unsorted.mat"));

% Create an Experiment object for subject 2055, which holds metadata about conditions/trials
experiment = Experiment(dataFolder, 2055);

% Load trial/condition information from disk into the experiment object
experiment = experiment.LoadExperiment();

nNeu      = size(recordedData.CS,1);   % Total number of recorded neurons
nOrigCond = size(recordedData.CS,2);   % Total number of original (pre-merge) conditions

% Identify neurons that have NO spike data in ANY condition and flag them for exclusion
neus2NOTconsider = [];
for neu = 1:nNeu
    for cond = 1:nOrigCond
        % If all trials for this neuron/condition are empty, mark the neuron
        if all(cellfun(@isempty, recordedData.CS{neu,cond}), 'all');
            neus2NOTconsider = [neus2NOTconsider; neu];
        end
    end
end
% Keep only neurons that have data in at least one condition
neus2consider = setdiff(1:nNeu, unique(neus2NOTconsider));


% Names of the trial event markers, in chronological order within a trial
markerNames = {'TrialStart','FIXon','TgObsOn','Early2LateDelay','GO','Early2LateRT','MOV','Early2LateMov','Touch','Rew','TrialEnd'};

%% Merge conditions
% Some conditions differ only in the delay duration (short vs. long delay).
% Here we collapse those pairs so that delay length is no longer a factor.

% Each cell contains the indices of conditions to be pooled together
cond2merge = {[1 5], [2 6], [3 7], [4 8], ... %merge conds that have only different delays
    [9 13], [10 14], [11 15], [12 16], ...
    [17 21], [18 22], [19 23], [20 24], ...
    [25 29], [26 30], [27 31], [28 32], ...
    [33 37], [34 38], [35 39], [36 40]};

% Remove the '_longDLY' suffix to get a clean list of unique condition type names
tmp = erase(experiment.conditions.Type, '_longDLY');
uniqueCondTypes = unique(tmp, 'stable');

% Merge spike and other data across the paired conditions in recordedData
recordedData = recordedData.mergeConditionsCSMS(cond2merge);

% Build a condition metadata table for the merged conditions (take metadata from the first condition of each pair)
for c = 1:numel(cond2merge)
    mergedConditions(c,:) = experiment.conditions(cond2merge{c}(1),:);
end

conditions = mergedConditions;
% Re-code obstacle labels: original codes 3→1 and 4→2 for a cleaner 2-level factor
conditions.Obstacle6(conditions.Obstacle6 == 3) = 1;
conditions.Obstacle6(conditions.Obstacle6 == 4) = 2;


% %% Plot PSTH  (Peri-Stimulus Time Histogram)
% Plots the trial-averaged firing rate of a single neuron around a task event.
%
neu2consider   = 61;       % Index of the neuron to plot
conds2consider = [15];     % Which merged condition(s) to include in the plot
window =  [-3 7 2];        % Analysis window: [time before marker, marker index, time after marker] (seconds)
bin = 0.04;                % Histogram bin width in seconds (40 ms)
vars2plot = {'TC'}; % Which signal to plot: 'TC' = touchscreen coordinates

recordedData.plotPSTH(neu2consider, conds2consider, window, bin, vars2plot, 'markerNames', markerNames)

%% dPCA  (demixed Principal Component Analysis)
% dPCA separates the population activity into components that are "demixed" by
% experimental parameter (target, obstacle, their interaction, condition-independent).

% Define which parameter combinations to marginalize over
% Each cell corresponds to one marginalization (e.g., {1,[1 3]} = target + target×time)
dPCAOpt.combinedParams = {{1,[1 3]},{2,[2 3]},{3},{[1 2],[1 2 3]}};
dPCAOpt.numComps       = 15;           % Number of dPCA components to extract per marginalization
dPCAOpt.simultaneousRecording = true;  % Data was recorded simultaneously (not across sessions)
dPCAOpt.procedure = 'simple';          % Use the standard dPCA fitting procedure

% --- Data pre-processing options ---
processingOpt.fixedFlag      = true;           % Use fixed-length time bins
processingOpt.neus2consider  = neus2consider;  % Only include neurons with valid spike data
processingOpt.binWidth       = [0.1 0.05];     % [bin size, step size] in seconds (fixed-bin mode)
% processingOpt.expectedBinWidth = 0.1;        % Alternative: expected bin width for variable-bin mode
processingOpt.smoothGaussianStdInBin = 0;      % Gaussian smoothing SD in bins (0 = no smoothing)
processingOpt.refType = 'normalize';           % Normalize firing rates before dPCA
processingOpt.refOpt.normType = 'soft-norm'; processingOpt.refOpt.normConstant = 'mean'; processingOpt.refOpt.normConstantMultiplier = 0.5; % extra normalization options for soft-normalization

% --- Training set options ---
trainOpt.conds2consider4train = [2 1 3 4] ;    % Condition used to fit the dPCA model
% Build a [nCond x nParam] matrix of task parameters (in this case, 2: target, obstacle) for each training condition
trainOpt.trainParams          = [conditions.Target6(trainOpt.conds2consider4train), conditions.Obstacle6(trainOpt.conds2consider4train)];
trainOpt.trainWindow          = [-0.5 7 0.75]; % Training time window: [pre, marker index, post] in seconds
% trainOpt.trainWindow          = [3 9];       % Alternative: first and last markers for the variable-bin mode

% --- Test set options --- NOTE: the test set can be different from the training set to assess generalization
testOpt.conds2consider4test  = [2 1 3 4] ;     % Condition used to project/test the model ()
testOpt.testParams           = [conditions.Target6(testOpt.conds2consider4test), conditions.Obstacle6(testOpt.conds2consider4test)];
testOpt.testWindow           = [-1 7 1];        % Test time window: [pre, marker index, post] in seconds
% testOpt.testWindow           = [1 11];        % Alternative: first and last markers for the variable-bin mode

% --- Visualization options ---
vOpt.marginalizationNames = {'TG', 'OBS', 'Condition-independent', 'TG/OBS Interaction'}; % Labels for each dPCA axis
vOpt.legendSubplot        = 16;  % Which subplot to place the legend in

% Run dPCA: returns the model, projections on training data, and projections on test data
[mdl, projTr, projTe] = recordedData.computedPCA(dPCAOpt, processingOpt, trainOpt, testOpt, vOpt);

%% Decoding PLANNING PHASE
% Train a classifier on one set of conditions (obstacle-present) and test it
% on another set (no-obstacle) to assess how well the neural population encodes
% target and obstacle information during the planning period.

processingOpt.fixedFlag        = true;           % Use fixed-length time bins
processingOpt.neus2consider    = neus2consider;  % Restrict to neurons with valid data
processingOpt.binWidth         = [0.1 0.05];     % [bin size, step size] in seconds
%processingOpt.expectedBinWidth = 0.1;           % Alternative for variable-bin mode
processingOpt.smoothGaussianStdInBin = 0;        % No Gaussian smoothing

% --- Training set options  ---
trainOpt.conds2consider4train = [6 5 7 8];       % Condition indices for training
trainOpt.trainParams          = [conditions.Target6(trainOpt.conds2consider4train), conditions.Obstacle6(trainOpt.conds2consider4train)];
trainOpt.trainWindow          = [-1 3 1];         % Window around marker 3 (TgObsOn): 1 s before to 1 s after
%trainOpt.trainWindow = [6 7];                    % Alternative: first and last markers for the variable-bin mode

% --- Test set options  --- NOTE: the test set can be different from the training set to assess generalization
testOpt.conds2consider4test  = [2 1 3 4];
testOpt.testParams           = [conditions.Target6(testOpt.conds2consider4test), conditions.Obstacle6(testOpt.conds2consider4test)];
testOpt.testWindow           = [-1 3 1];          % Same window as training
%testOpt.testWindow  = [1 11];                    % Alternative: first and last markers for the variable-bin mode

% --- Classifier options ---
decoderOpt.decoder    = 'naivebayes'; % Naïve Bayes classifier
decoderOpt.CVKfold    = 5;            % 5-fold cross-validation

% --- Display options ---
visOpt.factorNames = {'target','obs'};  % Labels for the two decoded factors
visOpt.display     = 'accuracy';        % Show classification accuracy (not confusion matrix)
visOpt.markerNames = markerNames;       % Event markers for the time-axis labels

% Run decoding: returns time-resolved accuracy scores, binned data, and trained models
[scores, data, mdls] = recordedData.neuralDecodingClassification( ...
    processingOpt, trainOpt, testOpt, decoderOpt, visOpt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SDF  (Spike Density Function)
% Plots smoothed firing-rate traces for each neuron across selected conditions.
% Neurons are sorted either by their preferred condition or by condition index.
window = [-3 7 2];          % Analysis window: [pre, marker index, post] in seconds
binWidth = [0.05 0];        % [bin size, step (0 = non-overlapping bins)] in seconds
gaussianWidth = 5;          % Width of the Gaussian smoothing kernel (in bins)
order = "preference";       % Sort neurons by preferred condition ("preference") or condition number ("condition")
displayErrorFactor = 0.05;  % Scaling factor for the shaded error region (e.g. SEM)
neus2consider = neus2consider;   % Neurons to include in the plot
conds2consider = [6 5 7 8];      % Conditions to display (obstacle-present group)

recordedData.plotSDF(neus2consider, conds2consider, window, binWidth, order, displayErrorFactor, ...
    'smoothWindow', 5, 'markerNames', markerNames);

%% ANOVA su params  (ANOVA on task parameters)
% Runs a sliding-window ANOVA for each neuron to identify time bins where
% firing rate is significantly modulated by obstacle, target, or their interaction.
refWindow = [-1 7 1];       % Reference window for baseline normalization: [pre, marker, post]
anaWindow = [-1 7 1];       % Analysis window (same as reference here)
binWidth = [0.10 0]         % Bin width in seconds (100 ms non-overlapping bins)
neus2consider = neus2consider;    % Neurons to analyze
conds2consider = [6 5 7 8];      % Conditions to include (obstacle-present group)
factorNames = {'epoch','obstacle','target'};  % Labels for the ANOVA factors
% Build parameter matrix: each row is a condition, columns are obstacle and target codes
params = [conditions.Obstacle6(conds2consider), conditions.Target6(conds2consider)];

% Run time-resolved ANOVA; results contain p-values and effect sizes per neuron per bin
anovaResults = recordedData.computeAnovaOnParam(neus2consider, conds2consider, params, refWindow, anaWindow, binWidth, ...
    'markerNames', markerNames, 'factorNames', factorNames);

%% Decoding  (all neurons, variable-bin mode)
% Same classification analysis as above but using ALL neurons (no pre-filtering)
% and variable-length bins aligned to trial events.

processingOpt.fixedFlag        = false;    % Variable-bin mode: bins are aligned to trial markers
processingOpt.neus2consider    = 1:nNeu;   % Use every recorded neuron (no exclusion)
% processingOpt.binWidth         = [0.1 0]; % Fixed-bin alternative (commented out)
processingOpt.expectedBinWidth = 0.1;      % Target bin duration in seconds for variable-bin mode
processingOpt.smoothGaussianStdInBin = 0;  % No smoothing

% --- Training: obstacle-present conditions ---
trainOpt.conds2consider4train = [6 5 7 8];
trainOpt.trainParams          = [conditions.Obstacle6(trainOpt.conds2consider4train), conditions.Target6(trainOpt.conds2consider4train)];
% trainOpt.trainWindow          = [-0.05 3 1.05];  % Fixed-bin alternative
trainOpt.trainWindow = [6 7];   % Variable-bin range: bins 6–7 (relative to marker alignment)

% --- Testing: same conditions as training (cross-validated within the same set) ---
testOpt.conds2consider4test  = [6 5 7 8];
testOpt.testParams           = [conditions.Obstacle6(testOpt.conds2consider4test), conditions.Target6(testOpt.conds2consider4test)];
% testOpt.testWindow           = [-1 3 1];   % Fixed-bin alternative
testOpt.testWindow  = [1 11];  % Variable-bin range: bins 1–11

% --- Classifier options ---
decoderOpt.decoder    = 'lda';  % Linear Discriminant Analysis (better for high-dimensional data)
decoderOpt.CVKfold    = 5;      % 5-fold cross-validation

% --- Display options ---
visOpt.factorNames = {'obs', 'target'}; % Decoded factors: obstacle and target
visOpt.display     = 'accuracy';        % Plot classification accuracy over time
visOpt.markerNames = markerNames;

[scores, data, mdls] = recordedData.neuralDecodingClassification( ...
    processingOpt, trainOpt, testOpt, decoderOpt, visOpt);


%% DPCA  (second dPCA run — all neurons, variable-bin mode)
% Repeats the dPCA analysis using all neurons and variable-length bins
% to check how the demixed components look without pre-filtering neurons.

dPCAOpt.combinedParams = {{1,[1 3]},{2,[2 3]},{3},{[1 2],[1 2 3]}}; % Marginalization structure (same as above)
dPCAOpt.numComps       = 15;           % Number of components per marginalization
dPCAOpt.simultaneousRecording = true;  % Simultaneous recording assumption
dPCAOpt.procedure = 'simple';          % Standard fitting procedure

processingOpt.fixedFlag      = false;    % Variable-bin mode
processingOpt.neus2consider  = 1:nNeu;   % Include all neurons
% processingOpt.binWidth       = [0.1 0.05]; % Fixed-bin alternative
processingOpt.expectedBinWidth = 0.1;    % Target bin width in seconds
processingOpt.smoothGaussianStdInBin = 0; % No smoothing

% --- Training ---
trainOpt.conds2consider4train = [6 5 7 8] ;
trainOpt.trainParams          = [conditions.Target6(trainOpt.conds2consider4train), conditions.Obstacle6(trainOpt.conds2consider4train)];
% trainOpt.trainWindow          = [-1.5 7 1.5]; % Fixed-bin alternative
trainOpt.trainWindow          = [3 9];   % Variable-bin range for training

% --- Testing ---
testOpt.conds2consider4test  = [6 5 7 8] ;
testOpt.testParams           = [conditions.Target6(testOpt.conds2consider4test), conditions.Obstacle6(testOpt.conds2consider4test)];
% testOpt.testWindow           = [-1.5 7 1.5]; % Fixed-bin alternative
testOpt.testWindow           = [1 11];   % Variable-bin range for testing

vOpt.marginalizationNames = {'TG', 'OBS', 'Condition-independent', 'TG/OBS Interaction'};
vOpt.legendSubplot        = 16;

% Run dPCA: returns model, training projections, and test projections
[mdl, projTr, projTe] = recordedData.computedPCA(dPCAOpt, processingOpt, trainOpt, testOpt, vOpt);


smooth = 0;
window = [-1 3 1];
bin = [0.10 0.025];
fixedBinFlag = true;
neus2consider = 1:nNeu;
conds2consider = [6 5 7 8] ; %[2 1 3 4];  %[18 17 19 20]; %[2 1 3 4];

recordedData = recordedData.prepareTensorWithFixBin(window,bin, smooth, neus2consider, conds2consider);
%
% % window = [3 10];
% % expectedBinWidth = 0.05;
% % smoothGaussianStdinBin = 4; % Define standard deviation for Gaussian smoothing
% % fixedBinFlag = false;
% % recordedData = recordedData.prepareTensorWithVariableBin(window, expectedBinWidth, smoothGaussianStdinBin, neus2consider, conds2consider);
%
% conditions = mergedConditions;
% conditions.Obstacle6(conditions.Obstacle6 == 3) = 1;
% conditions.Obstacle6(conditions.Obstacle6 == 4) = 2;
%
% recordedData = recordedData.splitTensorWithParams([conditions.Target6(conds2consider), conditions.Obstacle6(conds2consider)],fixedBinFlag);
%
% % recordedData = recordedData.prepareTensor4dPCA(paramNames, conditions, fixedBinFlag, neus2consider, conds2consider);
%
% % TENSOR 'NEU x TG x OBS x TIME'
% combinedParams = { {1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]} };
% options.margNames = {'TG', 'OBS', 'Condition-independent', 'TG/OBS Interaction'};
% options.margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
%
% numComps = 10;
%
% dPCAResults = recordedData.computedPCA(combinedParams, numComps, options);
