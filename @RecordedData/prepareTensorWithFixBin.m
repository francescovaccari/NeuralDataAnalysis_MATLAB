function obj = prepareTensorWithFixBin(obj, window, binWidth, smoothGaussianStdInBin)
    % PREPARETENSORWITHFIXBIN - Prepare tensor with fixed-width binning
    %
    % INPUTS:
    %   obj                     - RecordedData object
    %   window                  - [start(s), alignment_marker, end(s)] (numeric vector, length 3)
    %   binWidth                - [bin_width, overlap] in seconds (numeric vector, length 2)
    %   smoothGaussianStdInBin  - Gaussian smoothing std in bins (numeric scalar ≥ 0)
    %
    % OUTPUT:
    %   obj                     - Updated RecordedData object
    %
    % DESCRIPTION:
    %   Creates binned neural activity with uniform bin widths. Supports overlapping bins
    %   for intrinsic temporal smoothing. State markers tracked in MarkerTensorForFixBin.
    %
    % REMARKS:
    %   - Firing rates in spikes/second
    %   - Overlapping bins computed using moving sum
    %   - Uses progress bar for monitoring
    %   - binWidth(1): bin size; binWidth(2): overlap (0 for no overlap)
    
    % Extract data from object
    CS = obj.CS;
    MS = obj.MS;

    TensorWithFixBinInfo = obj.TensorWithFixBinInfo;
    TensorWithFixBinInfo.window = window;
    TensorWithFixBinInfo.binWidth = binWidth;
    TensorWithFixBinInfo.smoothGaussianStdInBin = smoothGaussianStdInBin;

    [nNeu, nCond] = size(CS);
    nMaxTrial = max(cellfun(@numel, CS(:)));
    nMarker = numel(MS{5,1}{1,1});
    if binWidth(2) > 0 %partially overlapping bins
        edges = window(1) : binWidth(2) : window(3);
    else
        edges = window(1) : binWidth(1) : window(3);
    end

    TensorWithFixBin = nan(nNeu, nCond, nMaxTrial, numel(edges)-1);
    MarkerTensorForFixBin = nan(nNeu, nCond, nMaxTrial, nMarker);

    % Progress bar
    progressBar = waitbar(0, 'Preparing fixed-bin tensor: Neuron 0 / ' + string(nNeu) + ', Condition 0 / ' + string(nCond));

    for neu = 1:nNeu
        for cond = 1:nCond
            % Update progress bar
            progress = (neu - 1) / nNeu + (1 / nNeu) * (cond / nCond);
            waitbar(progress, progressBar, sprintf('Preparing fixed-bin tensor: Neuron %d / %d, Condition %d / %d', neu, nNeu, cond, nCond));

            for trial = 1:size(CS{neu, cond},2)
                markers = MS{neu,cond}{1,trial};
                for m = 1:nMarker
                    [~,~,I] = histcounts(markers(m), edges+markers(window(2)));
                    if I ~= 0 %if the marker is present in the interval
                        MarkerTensorForFixBin(neu, cond, trial, m) = I;
                    end

                end
                tmp_trial = histcounts(CS{neu,cond}{1,trial}, edges+markers(window(2)));
                if binWidth(2) > 0 %partially overlapping bins
                    tmp_trial = movsum(tmp_trial,floor(binWidth(1)/binWidth(2)));
                end
                if smoothGaussianStdInBin > 0
                    if binWidth(2)==0
                        tmp_trial = smoothdata(tmp_trial,'gaussian',smoothGaussianStdInBin*5);
                    else
                        tmp_trial = smoothdata(tmp_trial,'gaussian',smoothGaussianStdInBin*5*floor(binWidth(1)/binWidth(2)));
                    end

                end
                TensorWithFixBin(neu,cond,trial,:) = tmp_trial;
            end
        end
    end

    TensorWithFixBinCondAvg = squeeze(mean(TensorWithFixBin,3,'omitmissing'));
    close(progressBar);
    
    % Assign results to object properties
    obj.TensorWithFixBin = TensorWithFixBin;
    obj.TensorWithFixBinCondAvg = TensorWithFixBinCondAvg;
    obj.MarkerTensorForFixBin = MarkerTensorForFixBin;
    obj.TensorWithFixBinInfo = TensorWithFixBinInfo;

end
