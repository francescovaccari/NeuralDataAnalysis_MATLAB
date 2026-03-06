function dam_index = computeDAMindex(A)
% COMPUTEDAMINDEX  Compute a data-availability/adjacency-missingness (DAM) index
%
%   dam_index = computeDAMindex(A) computes a normalized index in the range
%   [0, 1] that summarizes how "missing" a 2-D matrix is, taking into account
%   the number of missing elements and their local spatial adjacency.
%   More missing elements or greater clustering results in a lower DAM index.
%
%   Syntax
%   ------
%   dam_index = computeDAMindex(A)
%
%   Inputs
%   ------
%   A : 2-D numeric or logical matrix
%       Matrix containing data values and NaNs. NaNs are interpreted as
%       missing samples; non-NaN entries are interpreted as available data.
%
%   Outputs
%   -------
%   dam_index : scalar double in [0, 1]
%       Normalized data availability index:
%       0  -> No missing samples (matrix full of data)
%       1  -> All samples missing (matrix full of NaNs)
%
%   Notes
%   -----
%   The algorithm counts valid (non-NaN) neighbors in the 8-connected neighborhood
%   for each cell, computing a penalty (8 - nNeighbors)^2 per cell. This penalty
%   is averaged over the matrix and normalized using reference cases (all data vs.
%   no data). Uses zero-padding outside the matrix boundaries (neighbors beyond
%   borders treated as missing). For very small matrices where normalization bounds
%   coincide, falls back to the simple fraction of NaNs.
%
%   Examples
%   --------
%   A = rand(20); A(A < 0.2) = NaN;
%   dam = computeDAMindex(A);
%
%   Please cite
%   -----------
%     Hadjidimitrakis, K., Vaccari, F. E., De Vitis, M., Filippini, M., Diomedi, S., & Fattori, P. (2026).
%     Spontaneous oculomotor behavior sharpens eye position signals in parietal cortex.
%     [Manuscript under review].
%

    % -----------------------
    % Input validation
    % -----------------------
    if nargin ~= 1
        error('computeDAMindex:InvalidNumInputs', 'Exactly one input argument is required.');
    end

    validateattributes(A, {'numeric','logical'}, {'2d','nonempty','real'}, mfilename, 'A', 1);

    % Treat Inf as invalid for this metric.
    if any(isinf(A(:)))
        error('computeDAMindex:InvalidValues', 'Input matrix A contains Inf values. Use finite values or NaN for missing data.');
    end

    % Special case: single element -> DAM is 1 if missing, else 0.
    if numel(A) == 1
        dam_index = double(isnan(A));
        return;
    end

    % -----------------------
    % Main computation
    % -----------------------
    dataMask = ~isnan(A);                     % true where data are available
    sampleIndex = meanLocalPenalty(dataMask);

    % Upper bound: no data (all missing)
    UpperBoundIndex = meanLocalPenalty(false(size(A)));

    % Lower bound: all data
    LowerBoundIndex = meanLocalPenalty(true(size(A)));

    denom = (UpperBoundIndex - LowerBoundIndex);

    % If bounds coincide (can happen for very small sizes), fall back to NaN fraction.
    if denom == 0 || ~isfinite(denom)
        dam_index = sampleIndex;
        return;
    end

    dam_index = (sampleIndex - LowerBoundIndex) / denom;

end

% -------------------------------------------------------------------------
% Helper: compute mean penalty given a logical mask of available data.
% -------------------------------------------------------------------------
function m = meanLocalPenalty(dataMask)
    % dataMask: logical matrix, true where data exist (available samples)

    % Count neighbors with data (8-connected) using a 3x3 convolution kernel.
    % 'same' uses zero-padding outside the borders, matching the original intent.
    neighborCount = conv2(double(dataMask), ones(3), 'same') - double(dataMask);

    % Penalty per cell: max 8 neighbors. Larger penalty when fewer data neighbors exist.
    index = (8 - neighborCount).^2;

    % Mean penalty over all cells.
    m = mean(index, 'all');
end