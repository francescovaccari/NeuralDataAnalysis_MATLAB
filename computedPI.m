function dpi = computedPI(avgSC_FIX)
% COMPUTEDPI  Compute depth Preference Index (dPI) from matrix or vector
%
%   dpi = COMPUTEDPI(avgSC_FIX)
%   Computes the directional Population Index (dPI) as a measure of directional
%   selectivity in neural population responses.
%
%   Syntax
%   ------
%   dpi = computedPI(avgSC_FIX)
%
%   Inputs
%   ------
%   avgSC_FIX : numeric matrix or vector
%       Activity values (e.g., average spike counts, firing rates for different directions)
%
%   Outputs
%   -------
%   dpi : scalar double in [0, 1] or NaN
%       Directional Population Index = range / (max + min)
%       where range = max - min.
%
%   Notes
%   -----
%   Implements the formula from Hadjidimitrakis et al. (2011):
%       dPI = range(avgSC_FIX,'all') / (max(avgSC_FIX,[],'all') + min(avgSC_FIX,[],'all'));
%   Returns NaN when computation is not possible (all values NaN or denominator = 0).
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

% get global max and min 
maxVal = max(avgSC_FIX, [], 'all');
minVal = min(avgSC_FIX, [], 'all');

% handle degenerate cases
if isnan(maxVal) || isnan(minVal)
    dpi = NaN;
    return
end

% range: difference between max and min (equivalent to range(...,'all'))
r = maxVal - minVal;

denom = maxVal + minVal;
if denom == 0
    dpi = NaN;
    return
end

% compute dPI
dpi = r / denom;
end