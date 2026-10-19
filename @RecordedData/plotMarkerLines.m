function plotMarkerLines(obj, window, markerNames, varargin)
    % PLOTMARKERLINES - Draw labeled vertical marker lines on the current axes
    %
    % SYNTAX:
    %   obj.plotMarkerLines(window)
    %   obj.plotMarkerLines(window, markerNames)
    %   obj.plotMarkerLines(window, markerNames, 'Color', [0.6 0.6 0.6], ...)
    %
    % INPUTS:
    %   obj          - RecordedData object, after prepareTensorWithFixBin has
    %                  been called for this window (so AlignedMarkerTensorForFixBin
    %                  is populated).
    %   window       - [start(s), alignment_marker, end(s)] (numeric vector,
    %                  length 3), same convention as prepareTensorWithFixBin.
    %   markerNames  - (optional) cell array of marker labels. May either span
    %                  all experiment markers (indexed the same way as
    %                  AlignedMarkerTensorForFixBin) or have one entry per
    %                  marker actually plotted within window. Default: {}
    %                  (lines are drawn unlabeled). Pass [] or {} to skip.
    %   varargin     - extra Name,Value pairs forwarded to xline (e.g. 'Color').
    %
    % DESCRIPTION:
    %   Single shared implementation of the marker-line plotting logic
    %   previously duplicated (and inconsistently reimplemented) across
    %   plotSDF, computeAnovaOnParam, computeAnovaOnCond, and
    %   neuralDecodingClassification. Computes marker positions via
    %   computeAlignedMarkerTimes (exact, no bin-index interpolation) and
    %   draws them with xline on the current axes (gca).
    %
    % See also: computeAlignedMarkerTimes, prepareTensorWithFixBin

    if nargin < 3 || isempty(markerNames)
        markerNames = {};
    end
    if ischar(markerNames) || isstring(markerNames)
        markerNames = cellstr(markerNames(:))';
    elseif iscell(markerNames)
        markerNames = cellfun(@char, markerNames(:)', 'UniformOutput', false);
    end

    [markerTimes, markerIdx] = obj.computeAlignedMarkerTimes(window);
    if isempty(markerTimes)
        return
    end

    if isempty(markerNames)
        xline(markerTimes, varargin{:});
        return
    end

    if ~isempty(markerIdx) && numel(markerNames) >= max(markerIdx)
        % markerNames spans (at least) all experiment markers referenced by
        % markerIdx; pick out the names of the markers actually plotted.
        labels = markerNames(markerIdx);
    elseif numel(markerNames) == numel(markerTimes)
        % markerNames has one entry per plotted marker.
        labels = markerNames;
    else
        xline(markerTimes, varargin{:});
        return
    end

    xline(markerTimes, '-', labels, 'LabelOrientation', 'aligned', ...
        'LabelHorizontalAlignment', 'left', varargin{:});
end
