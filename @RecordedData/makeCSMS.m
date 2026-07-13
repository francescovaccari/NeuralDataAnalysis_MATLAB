function obj = makeCSMS(filenameSpikes, filenameEvents, eventType, stateSeqs, matchOpt, otherMarkers)
% MAKECSMS - Create RecordedData from spike times and behavioral events
%
% SYNTAX:
%   recordedData = makeCSMS(filenameSpikes, filenameEvents, ...
%       eventType, stateSeqs, matchOpt, otherMarkers)
%
% INPUTS:
%   filenameSpikes  - Path to spike data file (string | char)
%   filenameEvents  - Path to behavioral events CSV (string | char)
%   eventType       - 'per_frame_events' | 'simple_events' (string)
%   stateSeqs       - State sequences to extract (cell array for 'strict', matrix for 'first&last')
%   matchOpt        - 'strict' | 'first&last' (string)
%   otherMarkers    - Additional state indices to track (numeric vector | [])
%
% OUTPUT:
%   obj             - RecordedData object with CS, MS, TC, TCTS, TRIAL_OUTCOME populated
%
% DESCRIPTION:
%   Creates RecordedData from raw spike data and behavioral event files.
%   Extracts trials matching specified state sequences and segments spike times and
%   behavioral variables accordingly.
%
% REMARKS:
%   - CSV expected columns: state, condID, currentTime, optionally xCursor, yCursor, errorType
%   - Supports two matching strategies:
%     - 'strict': match exact state sequences
%     - 'first&last': match first→last state transitions
%   - Automatically segments touchscreen data (TC) if available
%   - Tracks trial outcomes (error types, final states)
%
% EXAMPLE:
%   stateSeqs = {[1 2 3 4], [1 2 4]};
%   recordedData = makeCSMS(...
%       '/path/to/spikes.mat', '/path/to/events.csv', ...
%       'per_frame_events', stateSeqs, 'strict', [2 3]);
%
% See also: RecordedData

obj = RecordedData();

% Load spikes
load(filenameSpikes, 'results');  % Assuming 'results' variable contains the spike data
spikes = results; % why 'results' if they are spikes ???
clear results;

events = readtable(filenameEvents);  % Load behavioral events from CSV file
states = events.state;

if strcmp(matchOpt, "strict");
    totalSeq = numel(stateSeqs);
elseif strcmp(matchOpt, "first&last");
    totalSeq = size(stateSeqs,1);
else
    error(['matchOpt' matchOpt 'cannot be handled. Set "strict" or "first&last"'])
end

StartsOfTrials = [];
EndsOfTrials = [];
TrialOnsets = [];
allPenultimateStates = [];

for s = 1:totalSeq

    if strcmp(matchOpt, "strict");
        seq = stateSeqs{s}';
    elseif strcmp(matchOpt, "first&last");
        seq = stateSeqs(s,:);
    end

    switch eventType
        case 'per_frame_events'

            stateStarts = [1; find(diff(states) ~= 0) + 1]; % Compress consecutive equal states into single events.
            stateEnds = [stateStarts(2:end) - 1; numel(states)];
            stateStates = states(stateStarts);
            prevStates = [NaN; stateStates(1:end-1)];

            startTrialIdx = [];
            endTrialIdx = [];


            switch matchOpt
                case "strict"

                    if numel(stateStates) >= numel(seq)
                        startTrialIdx = find(arrayfun(@(i) isequal(stateStates(i:i+numel(seq)-1), seq), 1:numel(stateStates)-numel(seq)+1)).';
                    else
                        startTrialIdx = [];
                    end
                    endTrialIdx = startTrialIdx + numel(seq) - 1;
                    valid = endTrialIdx <= numel(stateEnds);
                    startTrialIdx = startTrialIdx(valid);
                    endTrialIdx = endTrialIdx(valid);

                case "first&last"
                    tmp_idx = find(stateStates == seq(1) | stateStates == seq(2))';

                    filteredStateStates = stateStates(tmp_idx);
                    filteredStateStarts = stateStarts(tmp_idx);
                    filteredStateEnds = stateEnds(tmp_idx);
                    seq = seq(:);

                    validTrialIdxIdx = find(arrayfun(@(i) isequal(filteredStateStates(i:i+numel(seq)-1), seq), 1:numel(filteredStateStates)-numel(seq)+1)).';

                    [~, startTrialIdx] = intersect(stateStarts,filteredStateStarts(validTrialIdxIdx));
                    [~, endTrialIdx] = intersect(stateEnds,filteredStateEnds(validTrialIdxIdx+1));

                otherwise
                    error('matchOpt must be "strict" or "first&last".');

            end

            StartsOfTrials = [StartsOfTrials; stateStarts(startTrialIdx)];
            EndsOfTrials = [EndsOfTrials; stateEnds(endTrialIdx)];
            TrialOnsets = [TrialOnsets; events.currentTime(stateStarts(startTrialIdx))];
            allPenultimateStates = [allPenultimateStates; prevStates(endTrialIdx)];

        case 'simple_events'
            error('eventType = simple_events not yet implemented')

    end
end

[~, trialOrder] = sort(TrialOnsets);
StartsOfTrials = StartsOfTrials(trialOrder);
EndsOfTrials = EndsOfTrials(trialOrder);
allPenultimateStates = allPenultimateStates(trialOrder);

switch eventType
    case 'per_frame_events'
        stateStartsTS = events.currentTime(stateStarts);
        stateEndsTS = events.currentTime(stateEnds);

        MS = cell(1, max(events.condID));
        CS = cell(size(spikes,2), max(events.condID));
        TC = MS; %prepare for touchscreen touch coordinates
        TCTS = MS; %prepare for touchscreen touch coordinates timestamps
        TRIAL_OUTCOME = MS;

        for t = 1:numel(StartsOfTrials)
            cond = events.condID(StartsOfTrials(t));
            if ismember('errorType', events.Properties.VariableNames)
                errorType = events.errorType{EndsOfTrials(t)};
            else
                errorType = '';
            end
            lastState = events.state(EndsOfTrials(t));
            penultimateState = allPenultimateStates(t);

            if ~isempty(MS{1,cond})
                trial = size(MS{1,cond},2)+1;
            else
                trial = 1;
            end

            if isempty(otherMarkers)
                warning('Leaving otherMarkers empty, could lead to an inconsistent number of markers across different trials')
            else
                MS{1,cond}{trial} = nan(1,numel(otherMarkers)+3);
            end

            trialStates_idx = (stateStarts > StartsOfTrials(t) & stateStarts < EndsOfTrials(t) & stateEnds ~= EndsOfTrials(t));
            trialStates = stateStates(trialStates_idx)';
            tmp = find(ismember(otherMarkers, trialStates));

            MS{1,cond}{trial}(1) = stateStartsTS(stateStarts == StartsOfTrials(t));
            MS{1,cond}{trial}(end-1:end) = [stateStartsTS(stateEnds == EndsOfTrials(t)) stateEndsTS(stateEnds == EndsOfTrials(t))];
            MS{1,cond}{trial}(tmp+1) = stateStartsTS(trialStates_idx)';

            % segment touchscreen trajectories
            if ismember('xCursor', events.Properties.VariableNames) && ismember('yCursor', events.Properties.VariableNames)
                TC{1,cond}{trial}.xCursor = events.xCursor(events.currentTime > MS{1,cond}{trial}(1) & events.currentTime < MS{1,cond}{trial}(end))';
                TC{1,cond}{trial}.yCursor = events.yCursor(events.currentTime > MS{1,cond}{trial}(1) & events.currentTime < MS{1,cond}{trial}(end))';
            else
                TC{1,cond}{trial}.xCursor = [];
                TC{1,cond}{trial}.yCursor = [];
                if t == 1  % Only warn once per sequence
                    warning('xCursor and/or yCursor columns not found in events table. TC will be empty.');
                end
            end
            TCTS{1,cond}{trial} = events.currentTime(events.currentTime > MS{1,cond}{trial}(1) & events.currentTime < MS{1,cond}{trial}(end))';

            % segment spikes
            for neu = 1:size(spikes,2)
                CS{neu, cond}{trial} = spikes{1, neu}(spikes{1, neu} > MS{1,cond}{trial}(1) & spikes{1, neu} < MS{1,cond}{trial}(end));
            end

            % keep information of the trial outcome (important for errors)
            TRIAL_OUTCOME{1,cond}{trial}.lastState = lastState;
            TRIAL_OUTCOME{1,cond}{trial}.penultimateState = penultimateState;
            TRIAL_OUTCOME{1,cond}{trial}.errorType = errorType;

        end

        MS = repmat(MS,size(spikes,2),1);
        TC = repmat(TC,size(spikes,2),1);
        TCTS = repmat(TCTS,size(spikes,2),1);
        TRIAL_OUTCOME = repmat(TRIAL_OUTCOME,size(spikes,2),1);

    case 'simple_events'
        error('eventType = simple_events not yet implemented')
end
obj.MS = MS;
obj.CS = CS;
obj.TC = TC;
obj.TCTS = TCTS;
obj.TRIAL_OUTCOME = TRIAL_OUTCOME;
end
