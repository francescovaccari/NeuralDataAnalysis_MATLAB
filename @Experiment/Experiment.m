%% Experiment - Class for loading and processing experimental data
%
% This class handles loading, processing, and organizing experimental data
% from CSV and JSON files, including conditions, raw data, targets, obstacles,
% task settings, and fixation points. It supports data segmentation, metadata
% building, and computation of derived metrics like difficulty and epoch duration.
%
% The class maintains internal experiment state and supports method chaining
% for streamlined data processing workflows.
%
% See also: Experiment.Experiment

classdef Experiment
    properties (SetAccess = private)
        Folder (1,1) string      % Path to folder containing data files
        ID (1,1) string          % Experiment ID (e.g., 'ID2002')
        
        targets                  % Loaded target definitions
        conditions               % Loaded experimental conditions
        obstacles                % Loaded obstacle definitions
        taskSettings             % Loaded task settings
        rawData                  % Loaded raw recorded data
        fixPts                   % Loaded or generated fixation points
        segmentedData            % Segmented trials by condition
        metaData                 % Trial metadata built from conditions and targets
        epochDuration            % Computed epoch durations by condition and trial
    end

    methods
        %% Constructor
        % Syntax: experiment = Experiment(folder, ID)
        % 
        % Initializes an Experiment object and automatically loads all available data
        % using LoadExperiment() implicitly.
        %
        % Inputs:
        %   folder - (char/string) Path to the data folder
        %   ID     - (char/string or numeric) Experiment ID. If numeric, 
        %            automatically prepends 'ID' prefix.
        %
        % Output:
        %   experiment - (Experiment) Initialized experiment object
        %
        % Error cases:
        %   - If folder does not exist
        %   - If ID is empty or invalid
        %   - If fewer than 2 arguments provided
        %
        % Example:
        %   experiment = Experiment('/path/to/data', 'ID2002')
        %   experiment = Experiment('/path/to/data', 2002)  % Automatically becomes 'ID2002'
        %
        % See also: normalizeID, LoadExperiment
        
        function obj = Experiment(folder, ID)
            if nargin ~= 2
                error('Experiment:InvalidConstructor', ...
                    'Use: experiment = Experiment(folder, ID).');
            end

            folder = string(folder);
            if ~isscalar(folder) || strlength(strtrim(folder)) == 0
                error('Experiment:InvalidFolder', ...
                    'folder must be a non-empty text value.');
            end
            folder = strtrim(folder);

            if ~isfolder(folder)
                error('Experiment:FolderNotFound', ...
                    'Folder not found: "%s".', folder);
            end

            obj.Folder = folder;
            obj.ID = obj.normalizeID(ID);
            
            % Automatically load all experiment data
            obj = obj.LoadExperiment();
        end

        %% LoadExperiment - Load all available experiment data
        % Convenience method that loads all experiment components (targets, conditions,
        % obstacles, task settings, raw data, fixation points). Missing files are
        % silently skipped.
        %
        % Syntax:
        %   experiment = experiment.LoadExperiment()
        %
        % Output:
        %   experiment - (Experiment) Updated with loaded data
        %
        % Note: This method is called implicitly by the constructor.
        %       Calling it again will reload all files.
        %
        % Example:
        %   experiment = experiment.LoadExperiment();
        %
        % See also: LoadTargets, LoadConditions, LoadData
        
        function obj = LoadExperiment(obj)
            obj.targets = [];
            obj.conditions = [];
            obj.obstacles = [];
            obj.taskSettings = [];
            obj.rawData = [];
            obj.fixPts = [];
            obj.segmentedData = [];
            obj.metaData = [];
            obj.epochDuration = [];

            [obj.targets, ~] = obj.loadIfFileExists(@() obj.LoadTargets());
            [obj.conditions, ~] = obj.loadIfFileExists(@() obj.LoadConditions());
            [obj.obstacles, ~] = obj.loadIfFileExists(@() obj.LoadObstacles());
            [obj.taskSettings, ~] = obj.loadIfFileExists(@() obj.LoadTaskSettings());
            [obj.rawData, ~] = obj.loadIfFileExists(@() obj.LoadData());
            [obj.fixPts, ~] = obj.loadIfFileExists(@() obj.loadFixationPointsData());
        end

        %% LoadTargets - Load target objects from JSON file
        % Loads target definitions (positions, IDs, sizes) from {ID}_targets.json
        %
        % Output:
        %   targets - (1xN struct array) Each struct contains fields:
        %             - tgID: Target ID
        %             - x, y: Target coordinates
        %             - TargetSize: Size of the target
        %
        % Example:
        %   targets = obj.LoadTargets();
        %
        % See also: LoadObstacles, LoadTaskSettings, readJson
        
        function targets = LoadTargets(obj)
            filePath = obj.buildFilePath("targets", ".json");
            data = obj.readJson(filePath);

            if ~isstruct(data) || ~isfield(data, 'targets')
                error('Experiment:MissingTargetsField', ...
                    'File "%s" does not contain a "targets" field.', filePath);
            end

            targets = data.targets;
        end

        %% LoadConditions - Load experimental conditions from CSV file
        % Loads the conditions table from {ID}_conditions.csv
        %
        % Output:
        %   conditions - (table) Experimental conditions with columns such as:
        %                - condID: Condition ID
        %                - Target*, FixPt*: Target and fixation point assignments
        %                - Other experiment parameters
        %
        % Example:
        %   conditions = obj.LoadConditions();
        %
        % See also: LoadData, LoadExperiment
        
        function conditions = LoadConditions(obj)
            filePath = obj.buildFilePath("conditions", ".csv");
            % conditions = obj.readCsv(filePath);
            conditions = readtable(filePath);
        end

        %% LoadData - Load raw sensor data from CSV file
        % Loads the raw experimental data from {ID}_data.csv
        %
        % Output:
        %   rawData - (table) Raw data with columns typically including:
        %             - currentTime: Timestamp of measurement
        %             - state: Current experimental state
        %             - condID: Associated condition ID
        %             - x, y: Position coordinates (if tracked)
        %             - Other sensor readings
        %
        % Example:
        %   rawData = obj.LoadData();
        %
        % See also: LoadConditions, LoadExperiment, SegmentData
        
        function rawData = LoadData(obj)
            filePath = obj.buildFilePath("data", ".csv");
            rawData = readtable(filePath);
            % rawData = obj.readCsv(filePath);
        end

        %% LoadObstacles - Load obstacle objects from JSON file
        % Loads obstacle definitions from {ID}_obstacles.json
        %
        % Output:
        %   obstacles - (1xN struct array) Each struct defines an obstacle
        %
        % Example:
        %   obstacles = obj.LoadObstacles();
        %
        % See also: LoadTargets, LoadTaskSettings
        
        function obstacles = LoadObstacles(obj)
            filePath = obj.buildFilePath("obstacles", ".json");
            data = obj.readJson(filePath);

            if ~isstruct(data) || ~isfield(data, 'obstacles')
                error('Experiment:MissingObstaclesField', ...
                    'File "%s" does not contain an "obstacles" field.', filePath);
            end

            obstacles = data.obstacles;
        end

        %% LoadTaskSettings - Load task settings from JSON file
        % Loads task configuration and parameters from {ID}_taskSettings.json
        %
        % Output:
        %   taskSettings - (struct) Task configuration parameters
        %
        % Example:
        %   taskSettings = obj.LoadTaskSettings();
        %
        % See also: LoadTargets, LoadObstacles
        
        function taskSettings = LoadTaskSettings(obj)
            filePath = obj.buildFilePath("taskSettings", ".json");
            taskSettings = obj.readJson(filePath);
        end

        %% LoadFixationPoints - Load or assign fixation points
        % Updates experiment with fixation point data. Can load from file
        % or use default values. Also removes related computed columns from conditions.
        % Supports method chaining.
        %
        % Syntax:
        %   experiment = experiment.LoadFixationPoints()
        %   experiment = experiment.LoadFixationPoints(defaultX, defaultY, defaultSize)
        %
        % Inputs (optional):
        %   defaultX    - (numeric) Default X coordinate
        %   defaultY    - (numeric) Default Y coordinate
        %   defaultSize - (numeric) Default target size
        %
        % Output:
        %   experiment - (Experiment) Updated experiment with fixPts field and refreshed metadata
        %
        % Example:
        %   experiment = experiment.LoadFixationPoints();
        %   experiment = experiment.LoadFixationPoints(0, 0, 100);
        %
        % See also: loadFixationPointsData, LoadTargets
        
        function obj = LoadFixationPoints(obj, defaultX, defaultY, defaultSize)
            if nargin == 1
                fixPts = obj.loadFixationPointsData();
            elseif nargin == 4
                fixPts = obj.loadFixationPointsData(defaultX, defaultY, defaultSize);
            else
                error('Experiment:InvalidFixationLoadInput', ...
                    'Use LoadFixationPoints() or LoadFixationPoints(defaultX, defaultY, defaultSize).');
            end

            obj.fixPts = fixPts;
            if istable(obj.conditions)
                obj.conditions = obj.removeConditionVars( ...
                    obj.conditions, ["FixPtTargetDistance", "FixPtTargetDist", "Difficulty"]);
            end
            obj = obj.refreshExperimentMetaData();
        end

        %% SegmentData - Segment raw data into trials based on state sequences
        % Divides raw data into individual trials based on recurring state patterns.
        % This creates segmentedData{condition}{trial} structure.
        % Supports method chaining.
        %
        % Syntax:
        %   experiment = experiment.SegmentData(stateSequence, strictOpt)
        %
        % Inputs:
        %   stateSequence - (numeric vector) State sequence marking trial boundaries
        %   strictOpt     - (char) "strict" or "first&last"
        %                   "strict": Matches exact complete state sequence
        %                   "first&last": Matches first and last states, allows intermediate variation
        %
        % Output:
        %   experiment - (Experiment) Updated with segmentedData and refreshed metaData
        %                segmentedData{cond}{trial} = table of trial data
        %
        % Example:
        %   experiment = experiment.SegmentData([1 2 5], 'strict');
        %   % Looks for complete sequences of states [1 -> 2 -> 5]
        %
        % See also: BuildMetaData, ComputeEpochDuration, LoadData
        
        function obj = SegmentData(obj, stateSequence, strictOpt)
            if isempty(obj.rawData)
                error('Experiment:MissingExperimentRawData', ...
                    'rawData is empty. Load raw data first.');
            end
            if isempty(obj.conditions)
                error('Experiment:MissingExperimentConditions', ...
                    'conditions is empty. Load conditions first.');
            end

            rawData = obj.rawData;
            conditions = obj.conditions;

            states = rawData.state(:);
            seq = stateSequence(:);
            if isempty(seq)
                error('Experiment:EmptyStateSequence', ...
                    'stateSequence cannot be empty.');
            end
            if isempty(states)
                obj.segmentedData = cell(1, size(conditions, 1));
                obj.epochDuration = [];
                obj = obj.refreshExperimentMetaData();
                return;
            end 

            % Compress consecutive equal states into state "runs".
            runStarts = [1; find(diff(states) ~= 0) + 1];
            runEnds = [runStarts(2:end) - 1; numel(states)];
            runStates = states(runStarts);

            opt = lower(string(strictOpt));
            startRunIdx = [];
            endRunIdx = [];
            StartsOfTrials = [];
            EndsOfTrials = [];

            switch opt
                case "strict"
                    startRunIdx = strfind(runStates.', seq.').';
                    endRunIdx = startRunIdx + numel(seq) - 1;
                    valid = endRunIdx <= numel(runEnds);
                    startRunIdx = startRunIdx(valid);
                    endRunIdx = endRunIdx(valid);

                    StartsOfTrials = runStarts(startRunIdx);
                    EndsOfTrials = runEnds(endRunIdx);

                case "first&last"
                    if isempty(stateSequence)
                        error('Experiment:EmptyStateSequence', ...
                            'stateSequence cannot be empty.');
                    end
                    if isempty(rawData.state)
                        obj.segmentedData = cell(1, size(conditions, 1));
                        return;
                    end

                    stateChanges = [true; (diff(states(1:end-1)) ~= 0); true];
                    rewards = find(stateChanges & rawData.state == seq(end));

                    for t = 1:numel(rewards)
                        % look for trial beginning
                        i = rewards(t); found1state = false;
                        while true
                            i = i-1;
                            if i == 0
                                ii = i+1;
                                break;
                            elseif (states(i) ~= seq(1) && found1state)
                                ii = i+1;
                                break;
                            end

                            if states(i) == seq(1) && ~found1state
                                found1state = true;
                            end
                        end
                        % look for trial end
                        j = rewards(t);
                        while true
                            j = j+1;
                            if j > numel(states)
                                jj = j-1;
                                break;
                            elseif j == numel(states) || states(j) ~= seq(end)
                                jj = j-1;
                                break;
                            end
                        end
                        if found1state
                            StartsOfTrials(end + 1) = ii; %#ok<AGROW>
                            EndsOfTrials(end + 1) = jj; %#ok<AGROW>
                        end
                    end

                otherwise
                    error('strictOpt must be "strict" or "first&last".');
            end

            segmentedData = cell(1, size(conditions,1));

            for t = 1:numel(StartsOfTrials)
                cond = rawData.condID(StartsOfTrials(t));
                if ~isempty(segmentedData{1,cond})
                    trial = size(segmentedData{1,cond},2)+1;
                else
                    trial = 1;
                end
                segmentedData{1,cond}{trial} = rawData(StartsOfTrials(t):EndsOfTrials(t),:);
            end

            obj.segmentedData = segmentedData;
            obj.epochDuration = [];
            obj = obj.refreshExperimentMetaData();
        end

        %% BuildMetaData - Build metadata structure from segmented data
        % Creates metaData{condition}{trial} structs containing trial information
        % including target details and condition parameters.
        % Supports method chaining.
        %
        % Syntax:
        %   experiment = experiment.BuildMetaData()
        %
        % Output:
        %   experiment - (Experiment) Updated with metaData field
        %                metaData{cond}{trial} = struct with:
        %                  - target1, target2, ...: Target objects referenced in condition
        %                  - condition: Condition parameters as struct
        %
        % Example:
        %   experiment = experiment.BuildMetaData();
        %   metadata = experiment.metaData{1}{1}  % First trial of first condition
        %   disp(metadata.target1.x)  % Access target position
        %
        % See also: SegmentData, buildMetaDataFromParts
        
        function obj = BuildMetaData(obj)
            if isempty(obj.segmentedData)
                error('Experiment:MissingExperimentSegmentedData', ...
                    'segmentedData is empty. Segment the data first.');
            end
            if isempty(obj.conditions)
                error('Experiment:MissingExperimentConditions', ...
                    'conditions is empty. Load conditions first.');
            end
            if isempty(obj.targets)
                error('Experiment:MissingExperimentTargets', ...
                    'targets is empty. Load targets first.');
            end

            obj.metaData = obj.buildMetaDataFromParts( ...
                obj.segmentedData, obj.conditions, obj.targets);
        end

        %% ComputeEpochDuration - Compute duration of specific states per trial
        % Calculates the total duration (in seconds or time units) that specified
        % states were active for each trial.
        % Supports method chaining.
        %
        % Syntax:
        %   experiment = experiment.ComputeEpochDuration(states2comp)
        %
        % Inputs:
        %   states2comp - (numeric vector) States to compute duration for
        %
        % Output:
        %   experiment - (Experiment) Updated with epochDuration field
        %                epochDuration{cond} = [dur_trial1, dur_trial2, ...]
        %                (empty if no trials for that condition)
        %
        % Algorithm:
        %   For each specified state, finds all indices where it occurs,
        %   sums duration between first and last occurrence (or catches errors).
        %
        % Example:
        %   experiment = experiment.ComputeEpochDuration([2 3]);
        %   % Sum durations when states 2 and 3 were active
        %
        % See also: SegmentData, LoadData
        
        function obj = ComputeEpochDuration(obj, states2comp)
            if isempty(obj.segmentedData)
                error('Experiment:MissingExperimentSegmentedData', ...
                    'segmentedData is empty. Segment the data first.');
            end

            segmentedData = obj.segmentedData;
            if ~iscell(segmentedData)
                error('Experiment:InvalidSegmentedData', ...
                    'segmentedData must be a cell array like segmentedData{1,cond}{trial}.');
            end
            if ~(isnumeric(states2comp) && isvector(states2comp) && ~isempty(states2comp))
                error('Experiment:InvalidStatesToCompute', ...
                    'states2comp must be a non-empty numeric vector of integer states.');
            end

            states2comp = double(states2comp(:).');
            if any(~isfinite(states2comp) | states2comp ~= fix(states2comp))
                error('Experiment:InvalidStatesToCompute', ...
                    'states2comp must contain finite integer values.');
            end
            states2comp = unique(states2comp);

            epochDuration = cell(size(segmentedData));

            for cond = 1:numel(segmentedData)
                trials = segmentedData{cond};
                if isempty(trials)
                    epochDuration{cond} = [];
                    continue;
                end
                if ~iscell(trials)
                    error('Experiment:InvalidSegmentedDataEntry', ...
                        'segmentedData{1,%d} must be a cell array of trial tables.', cond);
                end

                condDurations = zeros(1, numel(trials));
                for trial = 1:numel(trials)
                    trialData = trials{trial};
                    if isempty(trialData)
                        condDurations(trial) = 0;
                        continue;
                    end
                    if ~istable(trialData)
                        error('Experiment:InvalidTrialData', ...
                            'segmentedData{1,%d}{%d} must be a table.', cond, trial);
                    end

                    stateVar = obj.findVarNameCaseInsensitive(trialData, "state");
                    timeVar = obj.findVarNameCaseInsensitive(trialData, "currentTime");
                    if strlength(stateVar) == 0 || strlength(timeVar) == 0
                        error('Experiment:MissingStateOrTimeColumn', ...
                            'Each trial table must contain "state" and "currentTime" columns.');
                    end

                    trialStates = obj.toNumericColumn(trialData.(stateVar), stateVar);
                    trialTime = obj.toNumericColumn(trialData.(timeVar), timeVar);

                    totalDur = 0;
                    for s = states2comp
                        stateDur = 0;
                        runIdx = find(trialStates == s);
                        if isempty(runIdx)
                            continue;
                        end
                        try
                            stateDur = trialTime(max(runIdx)+1)-trialTime(min(runIdx));
                        catch
                            stateDur = trialTime(max(runIdx))-trialTime(min(runIdx));
                        end
                        totalDur = totalDur + stateDur;
                    end
                    condDurations(trial) = totalDur;
                end

                epochDuration{cond} = condDurations;
            end

            obj.epochDuration = epochDuration;
        end

        %% ComputeFixPtTargetDistance - Compute distance between fixation points and targets
        % Calculates Euclidean distance from each fixation point to its associated
        % target for each condition. Results stored in conditions.FixPtTargetDistance.
        % Supports method chaining.
        %
        % Syntax:
        %   experiment = experiment.ComputeFixPtTargetDistance()
        %   experiment = experiment.ComputeFixPtTargetDistance(defaultX, defaultY, defaultSize)
        %
        % Inputs (optional):
        %   defaultX, defaultY, defaultSize - (numeric) Fallback values if no fixPts file
        %
        % Output:
        %   experiment - (Experiment) Updated with:
        %                - fixPts: Loaded or created fixation points
        %                - conditions.FixPtTargetDistance: Computed distances
        %
        % Algorithm:
        %   distance = sqrt((targetX - fixX)^2 + (targetY - fixY)^2)
        %   Uses Target6 and FixPt6 columns from conditions to find pairs.
        %
        % Example:
        %   experiment = experiment.ComputeFixPtTargetDistance();
        %   experiment = experiment.ComputeFixPtTargetDistance(0, 0, 100);
        %
        % See also: ComputeDifficulty, LoadFixationPoints
        
        function obj = ComputeFixPtTargetDistance(obj, defaultX, defaultY, defaultSize)
            if isempty(obj.conditions)
                error('Experiment:MissingExperimentConditions', ...
                    'conditions is empty. Load conditions first.');
            end
            if isempty(obj.targets)
                error('Experiment:MissingExperimentTargets', ...
                    'targets is empty. Load targets first.');
            end

            conditions = obj.conditions;
            targets = obj.targets;
            if nargin == 1
                if ~isempty(obj.fixPts)
                    fixPts = obj.fixPts;
                else
                    fixPts = obj.loadFixationPointsData();
                end
            elseif nargin == 4
                fixPts = obj.loadFixationPointsData(defaultX, defaultY, defaultSize);
            else
                error('Experiment:InvalidFixPtTargetDistanceInput', ...
                    'Use ComputeFixPtTargetDistance() or ComputeFixPtTargetDistance(defaultX, defaultY, defaultSize).');
            end

            if ~istable(conditions)
                error('Experiment:InvalidConditionsInput', ...
                    'conditions must be a table (as returned by LoadConditions).');
            end
            if ~isstruct(targets)
                error('Experiment:InvalidTargetsInput', ...
                    'targets must be a struct array (as returned by LoadTargets).');
            end
            if ~isstruct(fixPts)
                error('Experiment:InvalidFixPtsInput', ...
                    'fixPts must be a struct array (as returned by LoadFixationPoints).');
            end

            target6Var = obj.findVarNameCaseInsensitive(conditions, "Target6");
            if strlength(target6Var) == 0
                error('Experiment:MissingTarget6Column', ...
                    'Conditions file must contain a "Target6" column.');
            end

            targetIDs = obj.toNumericColumn(conditions.(target6Var), target6Var);

            if isempty(targets)
                error('Experiment:EmptyTargets', ...
                    'Targets file does not contain any target entry.');
            end
            targetLookupIDs = arrayfun(@(s) double(s.tgID), targets(:));
            targetLookupX = arrayfun(@(s) double(s.x), targets(:));
            targetLookupY = arrayfun(@(s) double(s.y), targets(:));

            [isTargetFound, targetIdx] = ismember(targetIDs, targetLookupIDs);
            if any(~isTargetFound)
                missingTargetIDs = unique(targetIDs(~isTargetFound));
                error('Experiment:UnknownTargetID', ...
                    'Target6 references unknown target ID(s): %s.', strjoin(string(missingTargetIDs), ', '));
            end

            targetX = targetLookupX(targetIdx);
            targetY = targetLookupY(targetIdx);

            if isempty(fixPts)
                error('Experiment:EmptyFixationPoints', ...
                    'Fixation points file does not contain any fixation point entry.');
            end
            fixLookupIDs = arrayfun(@(s) double(s.fixPtID), fixPts(:));
            fixLookupX = arrayfun(@(s) double(s.x), fixPts(:));
            fixLookupY = arrayfun(@(s) double(s.y), fixPts(:));

            fixPt6Var = obj.findVarNameCaseInsensitive(conditions, "FixPt6");
            if strlength(fixPt6Var) > 0
                fixIDs = obj.toNumericColumn(conditions.(fixPt6Var), fixPt6Var);
            else
                if numel(fixLookupIDs) ~= 1
                    error('Experiment:MissingFixPt6Column', ...
                        ['Conditions file has no "FixPt6" column and fixationPoints has %d entries. ' ...
                        'Expected exactly one fixation point in this case.'], numel(fixLookupIDs));
                end
                fixIDs = repmat(fixLookupIDs(1), height(conditions), 1);
            end

            [isFixFound, fixIdx] = ismember(fixIDs, fixLookupIDs);
            if any(~isFixFound)
                missingFixIDs = unique(fixIDs(~isFixFound));
                error('Experiment:UnknownFixPointID', ...
                    'FixPt6 references unknown fixation point ID(s): %s.', strjoin(string(missingFixIDs), ', '));
            end

            fixX = fixLookupX(fixIdx);
            fixY = fixLookupY(fixIdx);

            conditions = obj.removeConditionVars(conditions, ["FixPtTargetDistance", "FixPtTargetDist", "Difficulty"]);
            conditions.FixPtTargetDistance = hypot(targetX - fixX, targetY - fixY);

            obj.fixPts = fixPts;
            obj.conditions = conditions;
            obj = obj.refreshExperimentMetaData();
        end

        %% ComputeDifficulty - Compute task difficulty using Fitts' Law
        % Calculates difficulty index based on distance and target size using either
        % Classic or Adjusted Fitts' Law formulation.
        % Supports method chaining.
        %
        % Syntax:
        %   experiment = experiment.ComputeDifficulty(method)
        %
        % Inputs:
        %   method     - (char) "ClassicFittsLaw" or "AdjustedFittsLaw"
        %                Classic: ID = log2(2D/W)
        %                Adjusted: ID = log2(D/W + 1)
        %                where D = distance, W = target width
        %
        % Output:
        %   experiment - (Experiment) Updated with conditions.Difficulty column
        %
        % Requirements:
        %   - Must call ComputeFixPtTargetDistance first
        %   - FixPtTargetDistance must be in conditions
        %   - All target sizes must be positive
        %
        % Example:
        %   experiment = experiment.ComputeFixPtTargetDistance();
        %   experiment = experiment.ComputeDifficulty('ClassicFittsLaw');
        %
        % References:
        %   - Fitts, P. M. (1954). "The information capacity of the human motor system
        %     in controlling the amplitude of movement."
        %
        % See also: ComputeFixPtTargetDistance
        
        function obj = ComputeDifficulty(obj, method)
            if isempty(obj.conditions)
                error('Experiment:MissingExperimentConditions', ...
                    'conditions is empty. Load conditions first.');
            end
            if isempty(obj.targets)
                error('Experiment:MissingExperimentTargets', ...
                    'targets is empty. Load targets first.');
            end

            conditions = obj.conditions;
            targets = obj.targets;

            if ~istable(conditions)
                error('Experiment:InvalidConditionsInput', ...
                    'conditions must be a table (as returned by LoadConditions).');
            end
            if ~isstruct(targets)
                error('Experiment:InvalidTargetsInput', ...
                    'targets must be a struct array (as returned by LoadTargets).');
            end

            method = string(method);
            if ~isscalar(method) || strlength(strtrim(method)) == 0
                error('Experiment:InvalidDifficultyMethod', ...
                    'method must be "ClassicFittsLaw" or "AdjustedFittsLaw".');
            end
            method = strtrim(method);

            distVar = obj.findVarNameCaseInsensitive(conditions, "FixPtTargetDistance");
            if strlength(distVar) == 0
                error('Experiment:MissingFixPtTargetDistanceColumn', ...
                    'Distance not computed yet. First call ComputeFixPtTargetDistance().');
            end
            distanceVals = obj.toNumericColumn(conditions.(distVar), distVar);
            if any(distanceVals < 0)
                error('Experiment:InvalidFixPtTargetDistance', ...
                    '"FixPtTargetDistance" must contain non-negative values.');
            end

            target6Var = obj.findVarNameCaseInsensitive(conditions, "Target6");
            if strlength(target6Var) == 0
                error('Experiment:MissingTarget6Column', ...
                    'conditions must contain a "Target6" column.');
            end
            targetIDs = obj.toNumericColumn(conditions.(target6Var), target6Var);

            if isempty(targets)
                error('Experiment:EmptyTargets', ...
                    'Targets input does not contain any target entry.');
            end
            targetLookupIDs = arrayfun(@(s) double(s.tgID), targets(:));
            targetLookupSize = arrayfun(@(s) double(s.TargetSize), targets(:));

            [isTargetFound, targetIdx] = ismember(targetIDs, targetLookupIDs);
            if any(~isTargetFound)
                missingTargetIDs = unique(targetIDs(~isTargetFound));
                error('Experiment:UnknownTargetID', ...
                    'Target6 references unknown target ID(s): %s.', strjoin(string(missingTargetIDs), ', '));
            end

            targetSizeVals = targetLookupSize(targetIdx);
            if any(targetSizeVals <= 0)
                error('Experiment:InvalidTargetSize', ...
                    'All referenced target sizes must be positive.');
            end

            switch lower(method)
                case "classicfittslaw"
                    difficultyVals = log2(2 * distanceVals ./ targetSizeVals);
                case "adjustedfittslaw"
                    difficultyVals = log2((distanceVals ./ targetSizeVals) + 1);
                otherwise
                    error('Experiment:InvalidDifficultyMethod', ...
                        'Unknown method "%s". Use "ClassicFittsLaw" or "AdjustedFittsLaw".', method);
            end

            conditions.Difficulty = difficultyVals;
            obj.conditions = conditions;
            obj = obj.refreshExperimentMetaData();
        end

        %% PlotFixPtsTargets - Visualize fixation points and targets
        % Creates a scatter plot showing fixation point and target locations.
        % Targets are shown as blue circles, fixation points as red diamonds.
        %
        % Syntax:
        %   fig = experiment.PlotFixPtsTargets(conds2plot)
        %
        % Inputs:
        %   conds2plot  - (numeric array or "all") Condition ID(s) to plot
        %
        % Output:
        %   fig - (figure handle) Figure containing the plot
        %
        % Plot Details:
        %   - Blue circles: Targets (labeled T{ID})
        %   - Red diamonds: Fixation points (labeled F{ID})
        %   - Marker size reflects TargetSize property
        %
        % Example:
        %   fig = experiment.PlotFixPtsTargets([1 2 3]);
        %   fig = experiment.PlotFixPtsTargets("all");
        %
        % See also: LoadFixationPoints, LoadTargets
        
        function fig = PlotFixPtsTargets(obj, conds2plot)
            if isempty(obj.conditions)
                error('Experiment:MissingExperimentConditions', ...
                    'conditions is empty. Load conditions first.');
            end
            if isempty(obj.fixPts)
                error('Experiment:MissingExperimentFixPts', ...
                    'fixPts is empty. Load fixation points first.');
            end
            if isempty(obj.targets)
                error('Experiment:MissingExperimentTargets', ...
                    'targets is empty. Load targets first.');
            end

            conditions = obj.conditions;
            fixPts = obj.fixPts;
            targets = obj.targets;

            if isempty(targets)
                error('Experiment:EmptyTargets', ...
                    'Targets file does not contain any target entry.');
            end
            if isempty(fixPts)
                error('Experiment:EmptyFixationPoints', ...
                    'Fixation points file does not contain any fixation point entry.');
            end

            allTargetIDs = arrayfun(@(s) double(s.tgID), targets(:));
            allFixIDs = arrayfun(@(s) double(s.fixPtID), fixPts(:));

            isAll = false;
            if ischar(conds2plot) || (isstring(conds2plot) && isscalar(conds2plot))
                isAll = strcmpi(strtrim(string(conds2plot)), "all");
            end
            if isAll
                targetIDs = allTargetIDs;
                fixIDs = allFixIDs;
            else
                condIDVar = obj.findVarNameCaseInsensitive(conditions, "condID");
                if strlength(condIDVar) == 0
                    error('Experiment:MissingCondIDColumn', ...
                        'Conditions file must contain a "condID" column.');
                end

                requestedConds = obj.toNumericColumn(conds2plot, "conds2plot");
                requestedConds = requestedConds(:);
                condIDs = obj.toNumericColumn(conditions.(condIDVar), condIDVar);
                isSelected = ismember(condIDs, requestedConds);

                if ~any(isSelected)
                    error('Experiment:NoMatchingConditions', ...
                        'None of the requested condition IDs are present in conditions.');
                end

                selectedConds = conditions(isSelected, :);
                names = string(selectedConds.Properties.VariableNames);

                targetVars = names(startsWith(lower(names), "target"));
                if isempty(targetVars)
                    error('Experiment:MissingTargetColumns', ...
                        'Conditions file has no columns starting with "Target".');
                end

                targetIDs = [];
                for iVar = 1:numel(targetVars)
                    vals = obj.toNumericColumn(selectedConds.(targetVars(iVar)), targetVars(iVar));
                    targetIDs = [targetIDs; vals(:)]; %#ok<AGROW>
                end
                targetIDs = unique(targetIDs).';

                lowerNames = lower(names);
                isFixIDVar = false(size(lowerNames));
                for iName = 1:numel(lowerNames)
                    nm = char(lowerNames(iName));
                    isFixIDVar(iName) = ~isempty(regexp(nm, ...
                        '^(fixpt|fixationpt|fixpoint|fixationpoint)\d*$', 'once'));
                end
                fixVars = names(isFixIDVar);
                if isempty(fixVars)
                    if numel(allFixIDs) ~= 1
                        error('Experiment:MissingFixPtColumns', ...
                            ['No fixation-point ID column found in conditions and fixation points has %d entries. ' ...
                            'Expected exactly one fixation point in this case.'], numel(allFixIDs));
                    end
                    fixIDs = allFixIDs(1);
                else
                    fixIDs = [];
                    for iVar = 1:numel(fixVars)
                        vals = obj.toNumericColumn(selectedConds.(fixVars(iVar)), fixVars(iVar));
                        fixIDs = [fixIDs; vals(:)]; %#ok<AGROW>
                    end
                    fixIDs = unique(fixIDs).';
                end
            end

            [isTargetFound, targetIdx] = ismember(targetIDs, allTargetIDs);
            if any(~isTargetFound)
                missingTargetIDs = unique(targetIDs(~isTargetFound));
                error('Experiment:UnknownTargetID', ...
                    'Requested conditions reference unknown target ID(s): %s.', strjoin(string(missingTargetIDs), ', '));
            end
            plotTargets = targets(targetIdx);

            [isFixFound, fixIdx] = ismember(fixIDs, allFixIDs);
            if any(~isFixFound)
                missingFixIDs = unique(fixIDs(~isFixFound));
                error('Experiment:UnknownFixPointID', ...
                    'Requested conditions reference unknown fixation point ID(s): %s.', strjoin(string(missingFixIDs), ', '));
            end
            plotFixPts = fixPts(fixIdx);

            fig = figure('Color', 'w', 'Name', 'Fixation Points and Targets');
            ax = axes(fig);
            hold(ax, 'on');
            grid(ax, 'on');
            axis(ax, 'equal');
            xlabel(ax, 'X');
            ylabel(ax, 'Y');
            title(ax, sprintf('ID %s - Fixation Points and Targets', obj.ID));

            targetX = arrayfun(@(s) double(s.x), plotTargets(:));
            targetY = arrayfun(@(s) double(s.y), plotTargets(:));
            targetSize = arrayfun(@(s) double(s.TargetSize), plotTargets(:));
            targetIDsTxt = arrayfun(@(s) sprintf('T%d', s.tgID), plotTargets(:), 'UniformOutput', false);
            targetMarker = max(36, targetSize * 4);
            scatter(ax, targetX, targetY, targetMarker, ...
                'o', 'MarkerFaceColor', [0.40, 0.70, 1.00], 'MarkerEdgeColor', [0.00, 0.35, 0.85], 'LineWidth', 1.0);
            text(ax, targetX, targetY, targetIDsTxt, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.00, 0.20, 0.50], 'FontWeight', 'bold');

            fixX = arrayfun(@(s) double(s.x), plotFixPts(:));
            fixY = arrayfun(@(s) double(s.y), plotFixPts(:));
            fixSize = arrayfun(@(s) double(s.TargetSize), plotFixPts(:));
            fixIDsTxt = arrayfun(@(s) sprintf('F%d', s.fixPtID), plotFixPts(:), 'UniformOutput', false);
            fixMarker = max(49, fixSize * 4);
            scatter(ax, fixX, fixY, fixMarker, ...
                'd', 'MarkerFaceColor', [1.00, 0.75, 0.75], 'MarkerEdgeColor', [0.75, 0.00, 0.00], 'LineWidth', 1.2);
            text(ax, fixX, fixY, fixIDsTxt, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', [0.60, 0.00, 0.00], 'FontWeight', 'bold');

            legend(ax, {'Targets', 'Fixation Points'}, 'Location', 'best');
        end
    end

    methods (Access = private)
        %% normalizeID - Convert and validate experiment ID
        % Converts numeric IDs to string with 'ID' prefix, validates format.
        %
        % Input:
        %   inputID - (numeric or string) ID value or string
        %
        % Output:
        %   id - (string) Normalized ID (always starts with 'ID')
        %
        % Examples:
        %   normalizeID(2002) -> "ID2002"
        %   normalizeID("2002") -> "ID2002"
        %   normalizeID("ID2002") -> "ID2002"
        %
        % Error cases:
        %   - Numeric ID is negative or non-integer
        %   - String ID is empty
        %
        % See also: Experiment (constructor)
        
        function id = normalizeID(~, inputID)
            if isnumeric(inputID) && isscalar(inputID)
                if ~isfinite(inputID) || inputID < 0 || inputID ~= fix(inputID)
                    error('Experiment:InvalidID', ...
                        'Numeric ID must be a non-negative integer.');
                end
                id = "ID" + string(inputID);
                return;
            end

            id = strtrim(string(inputID));
            if ~isscalar(id) || strlength(id) == 0
                error('Experiment:InvalidID', ...
                    'ID must be text or a non-negative integer.');
            end

            if ~startsWith(upper(id), "ID")
                id = "ID" + id;
            end
        end

        function filePath = buildFilePath(obj, stem, ext)
            fileName = obj.ID + "_" + stem + ext;
            filePath = fullfile(obj.Folder, fileName);

            if ~isfile(filePath)
                error('Experiment:FileNotFound', ...
                    'Expected file not found: "%s".', filePath);
            end
        end

        function [value, wasLoaded] = loadIfFileExists(~, loaderFcn)
            value = [];
            wasLoaded = false;

            try
                value = loaderFcn();
                wasLoaded = true;
            catch ME
                if strcmp(ME.identifier, 'Experiment:FileNotFound')
                    return;
                end
                rethrow(ME);
            end
        end

        function fixationPoints = loadFixationPointsData(obj, defaultX, defaultY, defaultSize)
            filePath = fullfile(obj.Folder, obj.ID + "_fixationPoints.json");

            if nargin == 4
                if ~(isnumeric(defaultX) && isscalar(defaultX) && isfinite(defaultX) && ...
                     isnumeric(defaultY) && isscalar(defaultY) && isfinite(defaultY) && ...
                     isnumeric(defaultSize) && isscalar(defaultSize) && isfinite(defaultSize) && defaultSize > 0)
                    error('Experiment:InvalidFixationFallback', ...
                        'defaultX, defaultY must be finite scalars and defaultSize must be a positive finite scalar.');
                end

                fixationPoints = struct( ...
                    'fixPtID', 1, ...
                    'x', double(defaultX), ...
                    'y', double(defaultY), ...
                    'TargetSize', double(defaultSize) ...
                );
                return;
            end

            if nargin ~= 1
                error('Experiment:InvalidFixationFallback', ...
                    'Use loadFixationPointsData() or loadFixationPointsData(defaultX, defaultY, defaultSize).');
            end

            if ~isfile(filePath)
                error('Experiment:FileNotFound', ...
                    'Expected file not found: "%s".', filePath);
            end

            data = obj.readJson(filePath);
            if ~isstruct(data) || ~isfield(data, 'fixationPoints')
                error('Experiment:MissingFixationPointsField', ...
                    'File "%s" does not contain a "fixationPoints" field.', filePath);
            end
            fixationPoints = data.fixationPoints;
        end

        function metaData = buildMetaDataFromParts(obj, segmentedData, conditions, targets)
            if ~iscell(segmentedData)
                error('Experiment:InvalidSegmentedData', ...
                    'segmentedData must be a cell array like segmentedData{1,cond}{trial}.');
            end
            if ~istable(conditions)
                error('Experiment:InvalidConditionsInput', ...
                    'conditions must be a table (as returned by LoadConditions).');
            end
            if ~isstruct(targets)
                error('Experiment:InvalidTargetsInput', ...
                    'targets must be a struct array (as returned by LoadTargets).');
            end

            targetVars = obj.findOrderedTargetVarNames(conditions);
            if isempty(targetVars)
                error('Experiment:MissingTargetColumns', ...
                    'Conditions file must contain at least one "TargetX" column.');
            end

            if isempty(targets)
                error('Experiment:EmptyTargets', ...
                    'Targets input does not contain any target entry.');
            end

            targetLookupIDs = arrayfun(@(s) double(s.tgID), targets(:));
            metaData = cell(size(segmentedData));

            for cond = 1:numel(segmentedData)
                trials = segmentedData{1,cond};
                if isempty(trials)
                    metaData{1,cond} = [];
                    continue;
                end
                if ~iscell(trials)
                    error('Experiment:InvalidSegmentedDataEntry', ...
                        'segmentedData{1,%d} must be a cell array of trial tables.', cond);
                end

                condRow = obj.getConditionRowForIndex(conditions, cond);
                metaTemplate = struct();

                for iVar = 1:numel(targetVars)
                    targetVar = targetVars(iVar);
                    targetID = obj.toOptionalNumericScalar(condRow.(targetVar), targetVar);
                    if isempty(targetID)
                        continue;
                    end

                    [isTargetFound, targetIdx] = ismember(targetID, targetLookupIDs);
                    if ~isTargetFound
                        error('Experiment:UnknownTargetID', ...
                            'Condition %d references unknown target ID %g in column "%s".', ...
                            cond, targetID, targetVar);
                    end

                    fieldName = matlab.lang.makeValidName(lower(char(targetVar)));
                    metaTemplate.(fieldName) = targets(targetIdx);
                end

                metaTemplate.condition = table2struct(condRow);

                metaTrials = cell(size(trials));
                for trial = 1:numel(trials)
                    metaTrials{trial} = metaTemplate;
                end
                metaData{1,cond} = metaTrials;
            end
        end

        function conditions = removeConditionVars(~, conditions, requestedNames)
            if ~istable(conditions) || isempty(conditions)
                return;
            end

            requestedNames = string(requestedNames(:));
            varNames = string(conditions.Properties.VariableNames);
            keepMask = true(size(varNames));

            for iName = 1:numel(requestedNames)
                keepMask = keepMask & ~strcmpi(varNames, requestedNames(iName));
            end

            conditions = conditions(:, keepMask);
        end

        function obj = refreshExperimentMetaData(obj)
            obj.metaData = [];

            if ~istable(obj.conditions)
                return;
            end
            if isempty(obj.targets)
                return;
            end
            if isempty(obj.segmentedData)
                return;
            end

            obj.metaData = obj.buildMetaDataFromParts( ...
                obj.segmentedData, obj.conditions, obj.targets);
        end

        function data = readJson(~, filePath)
            try
                data = jsondecode(fileread(filePath));
            catch ME
                error('Experiment:InvalidJson', ...
                    'Could not parse JSON file "%s": %s', filePath, ME.message);
            end
        end

        function tbl = readCsv(~, filePath)
            % Use MATLAB's optimized readtable with European decimal format support
            try
                opts = delimitedTextImportOptions("NumVariables", 1, ...
                    "Delimiter", ";", ...
                    "Encoding", "UTF-8")%, ...
                    %"DecimalSeparator", ",");  % Handle European comma decimal format
                tbl = readtable(filePath, opts);
            catch ME
                error('Experiment:MalformedCsv', ...
                    'Could not read CSV file "%s": %s', filePath, ME.message);
            end
        end

        function varName = findVarNameCaseInsensitive(~, tbl, requestedName)
            names = string(tbl.Properties.VariableNames);
            idx = find(strcmpi(names, requestedName), 1, 'first');
            if isempty(idx)
                varName = "";
            else
                varName = names(idx);
            end
        end

        function targetVars = findOrderedTargetVarNames(~, tbl)
            names = string(tbl.Properties.VariableNames);
            targetIdx = false(size(names));
            targetOrder = zeros(size(names));

            for iName = 1:numel(names)
                tokens = regexp(char(names(iName)), '^target(\d+)$', 'tokens', 'once', 'ignorecase');
                if isempty(tokens)
                    continue;
                end

                targetIdx(iName) = true;
                targetOrder(iName) = str2double(tokens{1});
            end

            targetVars = names(targetIdx);
            targetOrder = targetOrder(targetIdx);
            [~, order] = sort(targetOrder);
            targetVars = targetVars(order);
        end

        function condRow = getConditionRowForIndex(obj, conditions, condIndex)
            condIDVar = obj.findVarNameCaseInsensitive(conditions, "condID");
            if strlength(condIDVar) > 0
                condIDs = obj.toNumericColumn(conditions.(condIDVar), condIDVar);
                rowIdx = find(condIDs == condIndex, 1, 'first');
                if ~isempty(rowIdx)
                    condRow = conditions(rowIdx, :);
                    return;
                end
            end

            if condIndex < 1 || condIndex > height(conditions)
                error('Experiment:UnknownConditionIndex', ...
                    'Condition index %d is not present in the conditions table.', condIndex);
            end

            condRow = conditions(condIndex, :);
        end

        function values = toNumericColumn(~, data, colName)
            if isnumeric(data)
                values = double(data(:));
            elseif isstring(data)
                values = str2double(strrep(data, ",", "."));
            elseif iscell(data)
                values = str2double(strrep(string(data), ",", "."));
            elseif iscategorical(data)
                values = str2double(strrep(string(data), ",", "."));
            else
                error('Experiment:UnsupportedColumnType', ...
                    'Column "%s" has unsupported type: %s.', colName, class(data));
            end

            if any(~isfinite(values))
                error('Experiment:InvalidNumericColumn', ...
                    'Column "%s" contains non-numeric or missing values.', colName);
            end
        end

        function value = toOptionalNumericScalar(~, data, colName)
            if isnumeric(data)
                values = double(data(:));
            elseif isstring(data)
                values = str2double(strrep(data(:), ",", "."));
            elseif iscell(data)
                values = str2double(strrep(string(data(:)), ",", "."));
            elseif iscategorical(data)
                values = str2double(strrep(string(data(:)), ",", "."));
            else
                error('Experiment:UnsupportedColumnType', ...
                    'Column "%s" has unsupported type: %s.', colName, class(data));
            end

            if isempty(values)
                value = [];
                return;
            end
            if numel(values) ~= 1
                error('Experiment:InvalidNumericScalar', ...
                    'Column "%s" must contain a single scalar value.', colName);
            end
            if ~isfinite(values) || isnan(values)
                value = [];
                return;
            end

            value = values;
        end
    end
end
