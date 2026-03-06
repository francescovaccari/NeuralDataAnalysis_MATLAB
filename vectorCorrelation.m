function [rho, beta, theta, p] = vectorCorrelation(x, y, u, v, varargin)
% VECCORRELATION  2D vector correlation between two vector sets
%
%   [rho, beta, theta, p] = VECCORRELATION(x, y, u, v)
%   Computes the correlation between two paired 2D vector sets [x,y] and
%   [u,v] using the approach from Hanson et al., 1992.
%
%   Inputs
%     x, y    - Numeric vectors (Nx1 or 1xN) containing the first set of
%               paired coordinates. x(i) and y(i) form the i-th 2D vector.
%     u, v    - Numeric vectors (same length as x and y) containing the
%               second set of paired coordinates. u(i) and v(i) form the
%               i-th 2D vector to compare to [x(i),y(i)].
%     Name/value optional inputs (use quoted names):
%       'pvalue'   - 'default' (compute p with t-transform) or 'shuffle'
%                    (compute empirical p by random shuffling). Default:
%                    'default'.
%       'nShuffle' - Positive integer number of permutations used when
%                    'pvalue' is 'shuffle'. Default: 1000.
%
%   Outputs
%     rho     - Scalar correlation coefficient in [-1, +1]. Positive values
%               indicate rotation-like alignment; negative values indicate
%               reflection-like alignment.
%     beta    - Scale factor relating the magnitudes of the two vector sets. 
%               Example: Beta = 2 indicates that [x,y] vector fields must
%               be scaled by a factor of 2 to match [u,v].
%     theta   - Rotation/reflection angle in degrees (computed via atand).
%     p       - Two-tailed p-value for rho. If 'pvalue' is 'default', p is
%               computed from a t-transform. If 'pvalue' is 'shuffle', p
%               is an empirical two-tailed p-value from the permutation
%               null distribution (add-one correction applied).
%
%   Notes
%     - All input vectors must have equal length N >= 3 for the t-transform.
%     - When using 'shuffle', the implementation applies the same random
%       permutation to x and y (preserving their pairing) and a separate
%       same random permutation to u and v. Each shuffle yields a rho; the
%       empirical p is the fraction of |rho_null| >= |rho_observed|.
%     - For reproducible permutations, set the RNG state before calling:
%           rng(seed)
%
%   Examples
%     % Default p via t-transform
%     [rho,beta,theta,p] = veccorrelation(x,y,u,v);
%
%     % Empirical p with 5000 permutations
%     [rho,beta,theta,p] = veccorrelation(x,y,u,v,'pvalue','shuffle','nShuffle',5000);
%
%   References
%   ----------
%     Hanson, B., Klink, K., Matsuura, K., Robeson, S. M., & Willmott, C. J. (1992). 
%     Vector correlation: Review, exposition, and geographic application. 
%     Annals of the Association of American Geographers, 82(1), 103–116. 
%     https://doi.org/10.1111/j.1467-8306.1992.tb01900.x
%	
%   Please cite
%     Hadjidimitrakis, K., Vaccari, F. E., De Vitis, M., Filippini, M., Diomedi, S., & Fattori, P. (2026). 
%     Spontaneous oculomotor behavior sharpens eye position signals in parietal cortex.
%     [Manuscript under review].

narginchk(4, inf);

% Default options
opts.pvalue = 'default';
opts.nShuffle = 1000;

% Parse optional name/value pairs
if ~isempty(varargin)
    if mod(length(varargin),2) ~= 0
        error('Optional parameters must be provided as name/value pairs.');
    end
    for k = 1:2:length(varargin)
        name = varargin{k};
        val  = varargin{k+1};
        validateattributes(name, {'char','string'}, {'scalartext'});
        switch lower(char(name))
            case 'pvalue'
                validateattributes(val, {'char','string'}, {'scalartext'});
                valstr = lower(char(val));
                if ~ismember(valstr, {'default','shuffle'})
                    error('''pvalue'' must be ''default'' or ''shuffle''.');
                end
                opts.pvalue = valstr;
            case 'nshuffle'
                validateattributes(val, {'numeric'}, {'scalar','integer','positive'});
                opts.nShuffle = double(val);
            otherwise
                error('Unknown option name ''%s''.', char(name));
        end
    end
end

n = length(x);
if length(y) ~= n || length(u) ~= n || length(v) ~= n
    error('All input vectors must have the same length.');
end

% Compute observed rho, beta, theta
[rho, beta, theta] = veccorrelation_core_compute(x, y, u, v);

switch opts.pvalue
    case 'default'
        % Default p-value by t-transform (two-tailed)
        tstat = rho * sqrt(length(x) - 2) / sqrt(1 - rho^2);
        p = 2 * (1 - tcdf(abs(tstat), length(x) - 2)); % two-tailed
    case 'shuffle'
        % Shuffle-based p-value


        nShuffle = opts.nShuffle;
        rho_null = zeros(nShuffle, 1);

        % For reproducibility a user can set rng before calling this function.
        for i = 1:nShuffle
            perm_uv = randperm(n);   % same permutation for u and v
            u_sh = u(perm_uv);
            v_sh = v(perm_uv);
            [rho_i, ~, ~] = veccorrelation_core_compute(x, y, u_sh, v_sh);
            rho_null(i) = rho_i;
        end

        % Two-tailed empirical p-value: proportion of |rho_null| >= |rho_observed|
        count_extreme = sum(abs(rho_null) >= abs(rho));
        p = (count_extreme + 1) / (nShuffle + 1); % add-one correction

end

end

function [rho, beta, theta] = veccorrelation_core_compute(x, y, u, v)
    % Core computations (returns rho, beta, theta)
    sigmax = std(x, 1);
    sigmay = std(y, 1);
    sigmau = std(u, 1);
    sigmav = std(v, 1);

    tmp = cov(x, u, 1);
    sigmaxu = tmp(1,2);
    tmp = cov(x, v, 1);
    sigmaxv = tmp(1,2);
    tmp = cov(y, u, 1);
    sigmayu = tmp(1,2);
    tmp = cov(y, v, 1);
    sigmayv = tmp(1,2);

    ksi = (sigmaxu * sigmayv) - (sigmaxv * sigmayu);
    % handle degenerate ksi == 0 (avoid division by zero)
    if ksi == 0
        s = 0;
        disp("WARNING: (sigmaxu * sigmayv) - (sigmaxv * sigmayu) == 0, computations might be altered")
    else
        s = ksi / abs(ksi);
    end

    a = sigmaxu^2 + sigmayv^2 + sigmaxv^2 + sigmayu^2 + (2 * s * ksi);
    b = (sigmax^2 + sigmay^2) * (sigmau^2 + sigmav^2);

    if b <= 0
        rho = NaN;
    else
        rho = s * sqrt(a / b);
    end

    if (sigmax^2 + sigmay^2) <= 0
        beta = NaN;
    else
        beta = s * rho * sqrt(sigmau^2 + sigmav^2) / sqrt(sigmax^2 + sigmay^2);
    end

    denom = (sigmaxu - s * sigmayv);
    if denom == 0
        theta = NaN;
    else
        theta = atand((sigmaxv - s * sigmayu) / denom);
    end
end
