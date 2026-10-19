function [markerTimes, markerIdx] = computeAlignedMarkerTimes(obj, window)
    % COMPUTEALIGNEDMARKERTIMES - Marker times relative to the window's alignment marker
    %
    % SYNTAX:
    %   [markerTimes, markerIdx] = obj.computeAlignedMarkerTimes(window)
    %
    % INPUTS:
    %   obj         - RecordedData object, after prepareTensorWithFixBin has
    %                 been called for this window (so AlignedMarkerTensorForFixBin
    %                 is populated).
    %   window      - [start(s), alignment_marker, end(s)] (numeric vector,
    %                 length 3), same convention as prepareTensorWithFixBin.
    %
    % OUTPUTS:
    %   markerTimes - Row vector of marker times (s), relative to the
    %                 window's alignment marker, averaged across neurons,
    %                 conditions and trials, restricted to markers that
    %                 fall inside [window(1) window(3)].
    %   markerIdx   - Absolute experiment-marker index (into the nMarker
    %                 dimension of AlignedMarkerTensorForFixBin) for each
    %                 entry in markerTimes. Lets a markerNames list that
    %                 spans all experiment markers be indexed consistently
    %                 by callers (e.g. plotMarkerLines).
    %
    % DESCRIPTION:
    %   Replaces the old pattern of converting BinnedMarkerTensorForFixBin
    %   (bin indices) to time via interpolation onto a manually reconstructed
    %   time axis. AlignedMarkerTensorForFixBin already stores each marker's
    %   exact timing relative to the alignment marker, so no reconstruction
    %   or interpolation is needed, and values stay perfectly aligned with
    %   TensorWithFixBinInfo.relativeBinCenters (same reference frame).
    %
    % See also: prepareTensorWithFixBin, plotMarkerLines

    if isempty(obj.AlignedMarkerTensorForFixBin)
        markerTimes = [];
        markerIdx = [];
        return
    end

    allMarkerTimes = squeeze(mean(obj.AlignedMarkerTensorForFixBin, [1 2 3], 'omitnan'))';
    mask = isfinite(allMarkerTimes) & allMarkerTimes >= window(1) & allMarkerTimes <= window(3);
    markerTimes = allMarkerTimes(mask);
    markerIdx = find(mask);
end
