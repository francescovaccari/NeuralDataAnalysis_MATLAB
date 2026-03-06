function pi = computePI(avgSC_FIX)
% COMPUTEPI  Compute Preference Index (PI) from matrix or vector
%
%   pi = computePI(avgSC_FIX)
%   Computes the Population Index (PI) as a measure of population favorability.
%
%   Syntax
%   ------
%   pi = computePI(avgSC_FIX)
%
%   Inputs
%   ------
%   avgSC_FIX : numeric matrix or vector
%       Activity values (e.g., average spike counts, firing rates)
%
%   Outputs
%   -------
%   pi : scalar double in [-1, +1] or NaN
%       Population Index = (n - sum/max) / (n - 1)
%       where n is the number of non-NaN values.
%
%   Notes
%   -----
%   Implements the formula from Moody et al. (1998):
%       n = length(find(~isnan(avgSC_FIX)));
%       PI = (n - nansum(avgSC_FIX,'all')/max(avgSC_FIX,[],'all')) / (n - 1);
%   Returns NaN when computation is not possible (all values NaN, n <= 1,
%   or maximum value is zero).
%
%   References
%   ----------
%   Moody, S. L., Wise, S. P., di Pellegrino, G., & Zipser, D. (1998). 
%   A model that accounts for activity in primate frontal cortex during a delayed matching-to-sample task. 
%   The Journal of neuroscience : the official journal of the Society for Neuroscience, 18(1), 399–410. 
%
%   Hadjidimitrakis, K., Selen, L. P., Bonaiuto, J. J., Breveglieri, R., Bosco, A., Galletti, C., & Fattori, P. (2011).
%   Variability of multisensory neurons in the anterior intraparietal area during a reaching task.
%   PLoS ONE, 6(12), e29619.
%
%   Please cite
%   -----------
%     Hadjidimitrakis, K., Vaccari, F. E., De Vitis, M., Filippini, M., Diomedi, S., & Fattori, P. (2026).
%     Spontaneous oculomotor behavior sharpens eye position signals in parietal cortex.
%     [Manuscript under review].

% number of non-NaN entries
n = length(find(~isnan(avgSC_FIX)));

% handle degenerate cases
if n <= 1
    pi = NaN;
    return
end

% sum over non-NaNs and global maximum 
sumVal = sum(avgSC_FIX, 'all','omitmissing');
maxVal = max(avgSC_FIX, [], 'all','omitmissing');

% if maxVal is NaN or zero, avoid division by zero / invalid result
if isnan(maxVal) || maxVal == 0
    pi = NaN;
    return
end

% compute PI
pi = (n - sumVal / maxVal) / (n - 1);
end