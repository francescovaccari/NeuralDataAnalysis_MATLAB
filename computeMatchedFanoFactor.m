% COMPUTEMATCHEDFANOFACTOR  Compute Fano Factor for matched distributions
%
%   [FF, FF_SE] = computeMatchedFanoFactor(meanA, meanB, varA, varB, num_resampling)
%   Computes the Fano Factor for two datasets with matched mean distributions
%   through iterative resampling. The Fano Factor represents the slope of the
%   variance vs. mean relationship, serving as a measure of noise characteristics.
%
%   Syntax
%   ------
%   [FF, FF_SE] = computeMatchedFanoFactor(meanA, meanB, varA, varB, num_resampling)
%
%   Inputs
%   ------
%   meanA (N x 1)          - Mean values for dataset A
%   meanB (M x 1)          - Mean values for dataset B
%   varA (N x 1)           - Variance values for dataset A
%   varB (M x 1)           - Variance values for dataset B
%   num_resampling (scalar) - Number of resampling iterations
%
%   Outputs
%   -------
%   FF (num_resampling x 2)    - Fano Factor estimates [FF_A, FF_B] for each iteration
%   FF_SE (num_resampling x 2) - Standard errors of FF estimates [SE_A, SE_B]
%
%   Notes
%   -----
%   For each resampling iteration:
%   1. Uses matchDistributions() to obtain paired indices with matched distributions
%   2. Selects corresponding mean and variance pairs from both datasets
%   3. Fits linear model: variance ~ mean
%   4. Extracts slope (Fano Factor) and standard error for each dataset
%   5. Repeats num_resampling times to obtain distribution of FF estimates
%
%   Examples
%   --------
%   [FF, FF_SE] = computeMatchedFanoFactor(meanA, meanB, varA, varB, 100);
%   FF_mean = mean(FF);  % Average Fano Factor across iterations
%
%   References
%   ----------
%   Churchland, M., Yu, B., Cunningham, J. et al. (2010)
%   Stimulus onset quenches neural variability: a widespread cortical phenomenon. 
%   Nat Neurosci 13, 369–378. 
%   https://doi.org/10.1038/nn.2501
%
%   Please cite
%   -----------
%     Hadjidimitrakis, K., Vaccari, F. E., De Vitis, M., Filippini, M., Diomedi, S., & Fattori, P. (2026).
%     Spontaneous oculomotor behavior sharpens eye position signals in parietal cortex.
%     [Manuscript under review].

function [FF, FF_SE] = computeMatchedFanoFactor(meanA, meanB, varA, varB, num_resampling)

    % Initialize output arrays to store FF estimates and standard errors
    FF = nan(num_resampling, 1); 
    FF_SE = FF;

    % Loop over resampling iterations
    for i = 1:num_resampling
        % Match the distributions between dataset A and B
        [~, ~, idx2pickA, idx2pickB] = match_distributions(meanA, meanB);

        % Extract matched mean values from both datasets
        matched_meanA = meanA(idx2pickA);
        matched_meanB = meanB(idx2pickB);

        % Extract corresponding variance values for matched pairs
        matched_varA = varA(idx2pickA);
        matched_varB = varB(idx2pickB);

        % Fit linear regression models to variance vs. mean relationship
        % Slope corresponds to the Fano Factor
        mdl_A = fitlm(matched_meanA, matched_varA);
        mdl_B = fitlm(matched_meanB, matched_varB);

        % Extract Fano Factor (slope, 2nd coefficient) and standard error for dataset A
        FF(i, 1) = mdl_A.Coefficients.Estimate(2);
        FF_SE(i, 1) = mdl_A.Coefficients.SE(2);

        % Extract Fano Factor (slope, 2nd coefficient) and standard error for dataset B
        FF(i, 2) = mdl_B.Coefficients.Estimate(2);
        FF_SE(i, 2) = mdl_B.Coefficients.SE(2);

    end


end
