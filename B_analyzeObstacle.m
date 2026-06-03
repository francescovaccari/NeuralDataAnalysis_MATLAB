clear all;
close all;

addpath(genpath("C:\Users\frada\Desktop\WORK in PROGRESS\LAB\TouchScreen_project\Matlab Script\Analisi dati"));

recordedData = RecordedData.loadData("C:\Users\frada\Desktop\WORK in PROGRESS\LAB\TouchScreen_project\Data\Datasets\CS_MS_TC_mef29_ID2020_st0_100_ch128_192.mat");

experiment = Experiment("C:\Users\frada\Desktop\WORK in PROGRESS\LAB\TouchScreen_project\Data\Datasets", 2020);

experiment = experiment.LoadExperiment();

nNeu = size(recordedData.CS,1);
nOrigCond = size(recordedData.CS,2);

offsetNeu = 192;

%% Merge conditions

cond2merge = {[1 5], [2 6], [3 7], [4 8], ... %merge conds that have only different delays
    [9 13], [10 14], [11 15], [12 16], ...
    [17 21], [18 22], [19 23], [20 24], ...
    [25 29], [26 30], [27 31], [28 32], ...
    [33 37], [34 38], [35 39], [36 40]};

tmp = erase(experiment.conditions.Type, '_longDLY');
uniqueCondTypes = unique(tmp, 'stable');

recordedData = recordedData.mergeConditionsCSMS(cond2merge);

for c = 1:numel(cond2merge)
mergedConditions(c,:) = experiment.conditions(cond2merge{c}(1),:); 
end

%% Plot PSTH

% neu2consider = 255-offsetNeu; % Neuron to consider for PSTH
% conds2consider = [6 5 7 8];   % Conditions to consider
% window = [-2 7 0.5];      % Time window for PSTH [onset index_of_the_marker offset] (in seconds) 
% bin = 0.04;                 % Bin size for PSTH (in seconds)
% vars2plot = {}; % Variables to plot
% 
% recordedData.plotPSTH(neu2consider, conds2consider, window, bin, vars2plot)


% %% SDF
% window = [-2 7 1];
% binWidth = [0.05 0.025];
% gaussianWidth = 5;
% order = "preference";
% displayErrorFactor = 0.05;
% neus2consider = 1:nNeu;
% conds2consider = [6 5 7 8]; 
% 
% recordedData = recordedData.prepareTensorWithFixBin(window,binWidth, gaussianWidth);
% recordedData.plotSDF(neus2consider, conds2consider, order, displayErrorFactor)

% %% Condition Decoding
% window = [-0.5 6 0.5];
% neus2consider = 1:nNeu;
% conds2consider = 9:12;
% decodingScores = recordedData.decodeCond(neus2consider, conds2consider, window);

% %% ANOVA (da sistemare)
% refWindow = [0 1 0.5];
% anaWindow = [-2 6 1];
% binWidth = [0.50 0.05]
% neus2consider = 1:nNeu;
% conds2consider = 9:12;
% 
% anovaResults = recordedData.computeContAnova(neus2consider, conds2consider, refWindow, anaWindow, binWidth);

%% DPCA
smooth = 0;
window = [-0.5 4 0.5];
bin = [0.050 0.025];
fixedBinFlag = true;
recordedData = recordedData.prepareTensorWithFixBin(window,bin, smooth);

% window = [3 10];
% expectedBinWidth = 0.05;
% smoothGaussianStdinBin = 4; % Define standard deviation for Gaussian smoothing
% fixedBinFlag = false;
% recordedData = recordedData.prepareTensorWithVariableBin(window, expectedBinWidth, smoothGaussianStdinBin);


paramNames = {'Target6', 'Obstacle6'}; 

neus2consider = [221 234 235	236	237	238	239	249	255 210	219	228	230	233	250]-offsetNeu; %1:nNeu;
conds2consider = [6 5 7 8]; 

conditions = mergedConditions;
conditions.Obstacle6(conditions.Obstacle6 == 3) = 1;
conditions.Obstacle6(conditions.Obstacle6 == 4) = 2;
recordedData = recordedData.prepareTensor4dPCA(paramNames, conditions, fixedBinFlag, neus2consider, conds2consider);

% TENSOR 'NEU x TG x OBS x TIME'
combinedParams = { {1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]} };
options.margNames = {'TG', 'OBS', 'Condition-independent', 'TG/OBS Interaction'};
% options.margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

numComps = 10;

dPCAResults = recordedData.computedPCA(combinedParams, numComps, options);