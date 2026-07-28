function obj = loadData(filename, simultaneousRecording)
% LOADDATA - Load neural recording data from MAT file (static method)
%
% SYNTAX:
%   recordedData = RecordedData.loadData(filename)
%   recordedData = RecordedData.loadData(filename, simultaneousRecording)
%
% INPUTS:
%   filename              - Path to MAT file (string | char)
%   simultaneousRecording - (optional) logical flag indicating whether units
%                           were recorded simultaneously; if true, MS, TC and
%                           KIN are shared across units, which speeds up
%                           computations (default: false)
%
% OUTPUT:
%   obj         - RecordedData object with loaded data
%
% DESCRIPTION:
%   Static method to create and populate RecordedData object from MAT file.
%   Expected variables: CS, MS, and optionally KIN, EY, TC, TCTS.
%   Missing variables are set to empty.
%
% REMARKS:
%   - File must be valid .mat file
%   - CS and MS are core variables (warning if missing)
%   - Additional variables (KIN, EY, TC, TCTS) are optional
%
% EXAMPLE:
%   recordedData = RecordedData.loadData('/path/to/mydata.mat');
%   recordedData = RecordedData.loadData('/path/to/mydata.mat', true);
%   disp(size(recordedData.CS));   % neurons × conditions
%
% See also: makeCSMS

if nargin < 2 || isempty(simultaneousRecording)
    simultaneousRecording = false;
end

filename = RecordedData.validateFilename(filename);
if ~simultaneousRecording 
    loadedData = load(filename, 'CS', 'MS', 'KIN', 'EY', 'TC', 'TCTS');
else
    mFile = matfile(filename);
    loadedData = struct();
    loadedData.CS = mFile.CS;
    loadedData.MS = repmat(mFile.MS(1,:), size(mFile.CS, 1), 1);
    if isprop(mFile, 'KIN')
        loadedData.KIN = repmat(mFile.KIN(1,:), size(mFile.CS, 1), 1);
    end
    if isprop(mFile, 'EY')
        loadedData.EY = repmat(mFile.EY(1,:), size(mFile.CS, 1), 1);
    end
    if isprop(mFile, 'TC')
        loadedData.TC = repmat(mFile.TC(1,:), size(mFile.CS, 1), 1);
    end
    if isprop(mFile, 'TCTS')
        loadedData.TCTS = repmat(mFile.TCTS(1,:), size(mFile.CS, 1), 1);
    end
end

requiredVariables = ["CS", "MS", "KIN", "EY"];
loadedVariables = string(fieldnames(loadedData));
missingVariables = requiredVariables(~ismember(requiredVariables, loadedVariables));

if ~isempty(missingVariables)
    warning('RecordedData:MissingVariables', ...
        'MAT file "%s" is missing variable(s): %s. Missing properties were set to [].', ...
        filename, strjoin(cellstr(missingVariables), ', '));
end

obj = RecordedData();
obj.Filename = filename;
obj.simultaneousRecording = simultaneousRecording;
obj.CS = RecordedData.getLoadedVariableOrDefault(loadedData, 'CS');
obj.MS = RecordedData.getLoadedVariableOrDefault(loadedData, 'MS');
obj.KIN = RecordedData.getLoadedVariableOrDefault(loadedData, 'KIN');
obj.EY = RecordedData.getLoadedVariableOrDefault(loadedData, 'EY');
obj.TC = RecordedData.getLoadedVariableOrDefault(loadedData, 'TC');
obj.TCTS = RecordedData.getLoadedVariableOrDefault(loadedData, 'TCTS');
end
