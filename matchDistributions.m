% MATCHDISTRIBUTIONS  Resample two distributions to match empirical probabilities
%
%   [subsample_A, subsample_B, idx2pickA, idx2pickB] = matchDistributions(A, B)
%   Resamples two distributions to match their empirical probability distributions
%   across bins, implementing the distribution matching procedure from Churchland et al.
%   (2010). Both datasets are downsampled to have equal sampling probability in
%   each bin, enabling fair statistical comparison.
%
%   Syntax
%   ------
%   [subsample_A, subsample_B, idx2pickA, idx2pickB] = matchDistributions(A, B)
%   [subsample_A, subsample_B, idx2pickA, idx2pickB] = matchDistributions(A, B, display)
%
%   Inputs
%   ------
%   A (N x 1)         - First distribution (column vector of numeric values)
%   B (M x 1)         - Second distribution (column vector of numeric values)
%   display (logical) - Optional. If true, displays histograms before and after
%                       matching. Default: false
%
%   Outputs
%   -------
%   subsample_A (K x 1) - Resampled values from A with matched distribution
%   subsample_B (K x 1) - Resampled values from B with matched distribution
%   idx2pickA (K x 1)   - Indices into original A: A(idx2pickA) = subsample_A
%   idx2pickB (K x 1)   - Indices into original B: B(idx2pickB) = subsample_B
%
%   Notes
%   -----
%   Algorithm:
%   1. Divides combined range [A; B] into 20 bins (5% quantile intervals)
%   2. Calculates probability in each bin for both distributions
%   3. Takes minimum probability across both distributions for each bin
%   4. Resamples proportionally to match the minimum bin probabilities
%   5. Returns matched subsamples and their corresponding indices
%
%   Examples
%   --------
%   A = randn(1000, 1);
%   B = randn(800, 1) * 1.5 + 0.5;
%   [sub_A, sub_B, idx_A, idx_B] = matchDistributions(A, B, true);
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
%
% =========================================================================

function [subsample_A, subsample_B, idx2pickA, idx2pickB] = matchDistributions(A, B, display)

    % Set default value for display parameter
    if nargin < 3
        display = false;
    end
    
    % Validate input: check that A and B are column vectors
    if ~iscolumn(A) || ~iscolumn(B)
        disp('Error: A and B must be column vectors')
        return
    end
    
    % Define bin edges based on the combined range of A and B
    % Using 5% quantile intervals results in 20 bins
    edges = quantile([A; B], 0:0.05:1);
    
    % Compute probability histograms for both distributions
    [hisA, ~, binsA] = histcounts(A, 'Normalization', 'probability', 'BinEdges', edges);
    [hisB, ~, binsB] = histcounts(B, 'Normalization', 'probability', 'BinEdges', edges);
    
    % Calculate matching probabilities: take minimum in each bin
    % This ensures both distributions are downsampled to the same level
    % Reference: Churchland et al., 2010
    new_prob = min([hisA; hisB], [], 1);
    
    % Resample distribution A to match the target probabilities
    subsample_A = nan(size(A));
    for ii = 1:length(hisA)
        % Calculate resampling probability for this bin
        if hisA(ii) == new_prob(ii)
            prob2pick = 1;  % Keep all samples if already at minimum
        else
            prob2pick = new_prob(ii) / hisA(ii);  % Downsample if above minimum
        end
        
        % Find all indices belonging to this bin and randomly resample
        idx = find(binsA == ii);
        num2pick = round(prob2pick * numel(idx));
        idx2pick = idx(randperm(numel(idx), num2pick));
        subsample_A(idx2pick) = A(idx2pick);
    end
    
    % Extract indices of retained samples and remove NaN padding
    idx2pickA = find(~isnan(subsample_A));
    subsample_A(isnan(subsample_A)) = [];
    
    % Resample distribution B to match the target probabilities
    subsample_B = nan(size(B));
    for ii = 1:length(hisB)
        % Calculate resampling probability for this bin
        if hisB(ii) == new_prob(ii)
            prob2pick = 1;  % Keep all samples if already at minimum
        else
            prob2pick = new_prob(ii) / hisB(ii);  % Downsample if above minimum
        end
        
        % Find all indices belonging to this bin and randomly resample
        idx = find(binsB == ii);
        num2pick = round(prob2pick * numel(idx));
        idx2pick = idx(randperm(numel(idx), num2pick));
        subsample_B(idx2pick) = B(idx2pick);
    end
    
    % Extract indices of retained samples and remove NaN padding
    idx2pickB = find(~isnan(subsample_B));
    subsample_B(isnan(subsample_B)) = [];
    
    % Display histograms if requested
    if display
        % Calculate percentage of data dropped for each distribution
        percent_dropped_A = (length(A) - length(subsample_A)) / length(A) * 100;
        percent_dropped_B = (length(B) - length(subsample_B)) / length(B) * 100;
        
        figure('Name', 'Distribution Matching Results', 'NumberTitle', 'off')
        
        % Plot original distributions
        subplot(2, 1, 1)
        histogram(A, 'BinEdges', edges, 'DisplayName', 'Distribution A', 'FaceAlpha', 0.6, 'Normalization', 'probability')
        hold on
        histogram(B, 'BinEdges', edges, 'DisplayName', 'Distribution B', 'FaceAlpha', 0.6, 'Normalization', 'probability')
        hold off
        xlabel('Value')
        ylabel('Probability')
        title('Original Distributions')
        legend
        grid on
        
        % Plot resampled distributions
        subplot(2, 1, 2)
        histogram(subsample_A, 'BinEdges', edges, 'DisplayName', 'Matched A', 'FaceAlpha', 0.6, 'Normalization', 'probability')
        hold on
        histogram(subsample_B, 'BinEdges', edges, 'DisplayName', 'Matched B', 'FaceAlpha', 0.6, 'Normalization', 'probability')
        xlabel('Value')
        ylabel('Probability')
        title(sprintf('Resampled Distributions (Matched) - A dropped: %.1f%%, B dropped: %.1f%%', ...
            percent_dropped_A, percent_dropped_B))
        legend
        grid on
    end

end