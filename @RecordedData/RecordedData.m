%% RECORDEDDATA CLASS
% RecordedData - Load, process, and analyze neural recording data
%
% This class provides tools for loading neural recording data (spike times, behavioral markers,
% kinematic data) and performing comprehensive analysis including binning, visualization,
% decoding, and dimensionality reduction via dPCA.
%
% USAGE:
%   recordedData = RecordedData.loadData(filename);
%   recordedData = RecordedData.makeCSMS(filenameSpikes, filenameEvents, ...);
%
% For full documentation, see RecordedData_Documentation.md
%
% Author: [Specify if known]
% Date: [Specify if known]

classdef RecordedData
    %% ===== CLASS PROPERTIES =====
    properties (SetAccess = public)
        % Input data
        Filename (1,1) string           % Path to loaded MAT file
        CS = [];                         % Spike times (neurons × conditions × trials)
        MS = [];                         % State marker times (neurons × conditions × trials)
        KIN = [];                        % Kinematic data (optional)
        EY = [];                         % Eye tracking data (optional)
        TC = [];                         % Touchscreen coordinates (optional)
        TCTS = [];                       % Touchscreen timestamps (optional)
        TRIAL_OUTCOME = [];              % Trial outcome metadata (optional)

        % Variable-bin tensors
        TensorWithVariableBinInfo = struct();      % Metadata for variable-bin tensor
        TensorWithVariableBin = [];                % Trial-by-trial activity (neurons × conditions × trials × time)
        TensorWithVariableBinCondAvg = [];         % Condition-averaged activity (neurons × conditions × time)
        MarkerTensorForVariableBin = [];           % Marker times trial-by-trial (neurons × conditions × trials × marker)
        MarkerTensorForVariableBinCondAvg = [];    % Marker times condition averaged (neurons × conditions × marker)
        TrialNumWithVariableBin = [];              % Trial counts for variable-bin pipeline (neurons × conditions)

        % Fixed-bin tensors
        TensorWithFixBinInfo = struct();           % Metadata for fixed-bin tensor
        TensorWithFixBin = [];                     % Trial-by-trial activity (neurons × conditions × trials × time)
        TensorWithFixBinCondAvg = [];              % Condition-averaged activity (neurons × conditions × time)
        MarkerTensorForFixBin = [];                % Marker times trial-by-trial (neurons × conditions × trials × marker)
        MarkerTensorForFixBinCondAvg = [];         % Marker times condition averaged (neurons × conditions × marker)
        TrialNumWithFixBin = [];                   % Trial counts for fixed-bin pipeline (neurons × conditions)

        % Parameter-split tensor
        TensorWithParamsInfo = struct();         % Metadata for parameter-split tensor
        TensorWithParams = [];                   % Trial-by-trial activity (neurons × param1 x param2 x ... × trials × time)
        TensorWithParamsCondAvg = [];            % Condition-averaged tensor (neurons × param1 x param2 x ... × time)
        MarkerTensorForParams = [];              % State marker times in parameter-split tensor
        TrialNumWithParams = [];                 % Trial counts (neurons × param1 × param2 × ...)
     
        
    end

    %% ===== PUBLIC METHODS =====
    methods
        %% Constructor
        function obj = RecordedData()
            % RECORDEDDATA - Initialize an empty RecordedData object
            %
            % USAGE:
            %   recordedData = RecordedData();
            %
            % Note: Use static methods loadData() or makeCSMS() to populate data.
            %
            % See also: loadData, makeCSMS
        end

        %% UTILITY METHODS
        obj = mergeConditionsCSMS(obj, cond2merge)

        %% TENSOR PREPARATION METHODS
        obj = prepareTensorWithVariableBin(obj, window, expectedBinWidth, neus2consider, conds2consider)
        
        obj = prepareTensorWithFixBin(obj, window, binWidth, neus2consider, conds2consider)

        obj = splitTensorWithParams(obj, params, fixedBinFlag)

        %% DATA REFINEMENT METHODS
        obj = refineNeuralData(obj, refType, refOpt, fixedFlag)
        
        %% VISUALIZATION METHODS
        plotPSTH(obj, neu2consider, conds2consider, window, bin, vars2plot, varargin)
        
        plotSDF(obj, neus2consider, conds2consider, window, binWidth, order, displayErrorFactor, varargin)

        %% DECODING AND ANALYSIS METHODS
        decodingScores = decodeCond(obj, neus2consider, conds2consider, window)

        %% computeAnovaOnCond
        anovaResults = computeAnovaOnCond(obj, neus2consider, conds2consider, refWindow, anaWindow, binWidth, varargin)

        %% computeAnovaOnParam
        anovaResults = computeAnovaOnParam(obj, neus2consider, conds2consider, params, refWindow, anaWindow, binWidth, varargin)

        %% decodeClassification
        % [decodingScores, Xtrain, labelsTrain, Xtest, labelsTest, mdl, predictedLabelsTrain, predictedLabelsTest, trainCVResults, testResults] = decodeClassification(obj, neus2consider, binWidth, conds2consider4train, trainParams, trainWindow, conds2consider4test, testParams, testWindow, varargin)
        [neuralDecodingScores, neuralDecodingData, neuralDecodingModels] = neuralDecodingClassification(obj, neus2consider, binWidth, conds2consider4train, trainParams, trainWindow, conds2consider4test, testParams, testWindow, varargin)
       
        %% dPCA ANALYSIS METHODS
        [dPCAmodel, dPCAprojTrain, dPCAprojTest] = computedPCA(obj, dPCAOpt, processingOpt, trainOpt, testOpt, visOpt)

    end

    %% ===== STATIC METHODS (CLASS CONSTRUCTION) =====
    methods (Static)
        %% loadData
        obj = loadData(filename)

        %% makeCSMS
        obj = makeCSMS(filenameSpikes, filenameEvents, eventType, stateSeqs, matchOpt, otherMarkers)
    end

    %% ===== PRIVATE HELPER METHODS =====
    methods (Access = private, Static)

        %% validateFilename
        function filename = validateFilename(filename)
            % VALIDATEFILENAME - Validate and process filename input (private method)
            %
            % SYNTAX:
            %   filename = RecordedData.validateFilename(filename)
            %
            % INPUTS:
            %   filename    - Input filename (string | char)
            %
            % OUTPUT:
            %   filename    - Validated and processed filename (string)
            %
            % DESCRIPTION:
            %   Validates that filename is a non-empty string, file exists, and has .mat extension.
            %
            % Errors if:
            %   - filename is empty or scalar mismatch
            %   - file does not exist
            %   - file extension is not .mat
            filename = string(filename);

            if ~isscalar(filename) || strlength(strtrim(filename)) == 0
                error('RecordedData:InvalidFilename', ...
                    'filename must be a non-empty text value.');
            end

            filename = strtrim(filename);

            if ~isfile(filename)
                error('RecordedData:FileNotFound', ...
                    'MAT file not found: "%s".', filename);
            end

            [~, ~, ext] = fileparts(filename);
            if ~strcmpi(ext, '.mat')
                error('RecordedData:InvalidFileExtension', ...
                    'filename must point to a .mat file.');
            end
        end

        %% getLoadedVariableOrDefault
        function value = getLoadedVariableOrDefault(loadedData, variableName)
            % GETLOADEDVARIABLEORDEFAULT - Get loaded variable or empty if missing (private method)
            %
            % SYNTAX:
            %   value = RecordedData.getLoadedVariableOrDefault(loadedData, variableName)
            %
            % INPUTS:
            %   loadedData      - Struct from load() function (struct)
            %   variableName    - Variable name to retrieve (string)
            %
            % OUTPUT:
            %   value           - Variable value or empty array if not found
            if isfield(loadedData, variableName)
                value = loadedData.(variableName);
            else
                value = [];
            end
        end

    end
end
