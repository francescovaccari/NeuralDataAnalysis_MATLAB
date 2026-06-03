function [obj] = mergeConditionsCSMS(obj, cond2merge)
    % MERGECONDITIONSCSMS - Merge multiple conditions into collapsed conditions
    %
    % INPUTS:
    %   obj             - RecordedData object
    %   cond2merge      - Cell array specifying condition merging
    %
    % OUTPUT:
    %   obj             - Updated RecordedData object
    %
    % DESCRIPTION:
    %   Restructures neural data by merging multiple conditions into fewer conditions.
    %   All trial data from merged conditions concatenated and sorted by trial onset.

    % ===== INPUT VALIDATION =====
    if ~iscell(cond2merge)
        error('RecordedData:InvalidInput', ...
            'cond2merge must be a cell array. Example: {[1 2], [3 4 5], [6]}');
    end

    if isempty(cond2merge)
        error('RecordedData:InvalidInput', ...
            'cond2merge must not be empty.');
    end

    if isempty(obj.CS) || isempty(obj.MS)
        error('RecordedData:MissingData', ...
            'CS and MS must be populated before merging conditions.');
    end

    % Extract data to be merged
    CS = obj.CS;
    MS = obj.MS;
    TC = obj.TC;
    TCTS = obj.TCTS;
    TRIAL_OUTCOME = obj.TRIAL_OUTCOME;

    % Get original number of conditions
    [nNeu, nOrigCond] = size(CS);

    % ===== VALIDATE CONDITION INDICES =====
    allIndices = [];
    for i = 1:numel(cond2merge)
        if ~isvector(cond2merge{i}) || ~isnumeric(cond2merge{i})
            error('RecordedData:InvalidInput', ...
                'Each element of cond2merge must be a vector of positive integers. Found element %d: %s', ...
                i, class(cond2merge{i}));
        end

        if any(cond2merge{i} < 1) || any(cond2merge{i} > nOrigCond)
            error('RecordedData:InvalidConditionIndex', ...
                'Condition indices must be between 1 and %d. Found invalid index in element %d: %s', ...
                nOrigCond, i, mat2str(cond2merge{i}));
        end

        if any(cond2merge{i} ~= round(cond2merge{i}))
            error('RecordedData:InvalidInput', ...
                'Condition indices must be integers. Found non-integer values in element %d.', i);
        end

        allIndices = [allIndices, cond2merge{i}];
    end

    % Check for duplicate indices
    if numel(allIndices) ~= numel(unique(allIndices))
        duplicates = unique(allIndices(duplicates(allIndices)));
        error('RecordedData:DuplicateIndices', ...
            'Each original condition can only appear once in cond2merge. Duplicates found: %s', ...
            mat2str(duplicates));
    end

    % Check if all original conditions are covered
    nNewCond = numel(cond2merge);
    allUsedIndices = sort(allIndices);
    expectedIndices = 1:nOrigCond;
    unusedIndices = expectedIndices(~ismember(expectedIndices, allUsedIndices));

    if ~isempty(unusedIndices)
        warning('RecordedData:UnusedConditions', ...
            'Not all original conditions are included in cond2merge. Unused: %s', ...
            mat2str(unusedIndices));
    end

    % ===== MERGE OPERATION =====
    mergedCS = cell(nNeu, nNewCond);
    mergedMS = cell(nNeu, nNewCond);
    mergedTC = [];
    mergedTCTS = [];
    mergedTRIAL_OUTCOME = [];

    % Check if optional properties exist and are populated
    hasTouchScreen = ~isempty(TC);
    hasTrialOutcome = ~isempty(TRIAL_OUTCOME);

    if hasTouchScreen
        mergedTC = cell(nNeu, nNewCond);
        mergedTCTS = cell(nNeu, nNewCond);
    end

    if hasTrialOutcome
        mergedTRIAL_OUTCOME = cell(nNeu, nNewCond);
    end

    % For each new condition, merge the original conditions
    for newCond = 1:nNewCond
        condList = cond2merge{newCond};

        % For each neuron
        for neu = 1:nNeu
            % Collect all trials
            allTrials_CS = {};
            allTrials_MS = {};
            allTrials_TC = {};
            allTrials_TCTS = {};
            allTrials_OUTCOME = {};
            allOnsetTimes = [];

            % Concatenate trials from all conditions to merge
            for origCond = condList
                nTrialsInCond = numel(CS{neu, origCond});

                for trial = 1:nTrialsInCond
                    allTrials_CS{end+1} = CS{neu, origCond}{trial};
                    allTrials_MS{end+1} = MS{neu, origCond}{trial};

                    % Extract onset time (first marker)
                    onsetTime = MS{neu, origCond}{trial}(1);
                    allOnsetTimes(end+1) = onsetTime;

                    if hasTouchScreen
                        iTrial = size(allTrials_TC, 2) + 1;
                        allTrials_TC{iTrial} = TC{neu, origCond}{trial};
                        allTrials_TCTS{iTrial} = TCTS{neu, origCond}{trial};
                    end

                    if hasTrialOutcome
                        allTrials_OUTCOME{end+1} = TRIAL_OUTCOME{neu, origCond}{trial};
                    end
                end
            end

            % Sort trials by onset time
            [~, sortIdx] = sort(allOnsetTimes);

            % Reorder all trials according to onset time
            mergedCS{neu, newCond} = allTrials_CS(sortIdx);
            mergedMS{neu, newCond} = allTrials_MS(sortIdx);

            if hasTouchScreen
                orderedTC = cell(1, numel(sortIdx));
                for trialIdx = 1:numel(sortIdx)
                    orderedTC{trialIdx} = allTrials_TC{sortIdx(trialIdx)};
                end
                mergedTC{neu, newCond} = orderedTC;
                mergedTCTS{neu, newCond} = allTrials_TCTS(sortIdx);
            end

            if hasTrialOutcome
                mergedTRIAL_OUTCOME{neu, newCond} = allTrials_OUTCOME(sortIdx);
            end
        end
    end

    % Display summary
    disp(['Successfully merged ' num2str(nOrigCond) ' original conditions into ' num2str(nNewCond) ' new conditions.']);
    for i = 1:nNewCond
        disp(['  New Condition ' num2str(i) ': merged original conditions ' mat2str(cond2merge{i}) ' and sorted by trial onset.']);
    end

    obj.CS = mergedCS;
    obj.MS = mergedMS;
    obj.TC = mergedTC;
    obj.TCTS = mergedTCTS;
    obj.TRIAL_OUTCOME = mergedTRIAL_OUTCOME;

end
