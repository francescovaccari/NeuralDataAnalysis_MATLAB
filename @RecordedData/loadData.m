function obj = loadData(filename)
% LOADDATA - Load neural recording data from MAT file (static method)
%
% SYNTAX:
%   recordedData = RecordedData.loadData(filename)
%
% INPUTS:
%   filename    - Path to MAT file (string | char)
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
%   disp(size(recordedData.CS));   % neurons × conditions
%
% See also: makeCSMS

filename = RecordedData.validateFilename(filename);
loadedData = load(filename, 'CS', 'MS', 'KIN', 'EY', 'TC', 'TCTS');

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
obj.CS = RecordedData.getLoadedVariableOrDefault(loadedData, 'CS');
obj.MS = RecordedData.getLoadedVariableOrDefault(loadedData, 'MS');
obj.KIN = RecordedData.getLoadedVariableOrDefault(loadedData, 'KIN');
obj.EY = RecordedData.getLoadedVariableOrDefault(loadedData, 'EY');
obj.TC = RecordedData.getLoadedVariableOrDefault(loadedData, 'TC');
obj.TCTS = RecordedData.getLoadedVariableOrDefault(loadedData, 'TCTS');
end
