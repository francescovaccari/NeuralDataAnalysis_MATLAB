function obj = prepareTensorWithVariableBin(obj, window, expectedBinWidth, smoothGaussianStdInBin)
    % PREPARETENSORWITHVARIABLEBIN - Prepare tensor with epoch-adaptive binning
    %
    % INPUTS:
    %   obj                     - RecordedData object
    %   window                  - [start_marker, end_marker] indices (numeric vector, length 2)
    %   expectedBinWidth        - Target bin duration in seconds (numeric scalar)
    %   smoothGaussianStdInBin  - Gaussian smoothing std in bins (numeric scalar ≥ 0)
    %
    % OUTPUT:
    %   obj                     - Updated RecordedData object
    %
    % DESCRIPTION:
    %   Creates binned neural activity where bin widths adapt to epoch durations. Bins are
    %   automatically merged if they fall below 50% of expected width. Gaussian smoothing
    %   is applied if smoothGaussianStdInBin > 0.
    %
    % REMARKS:
    %   - Firing rates in spikes/second
    %   - Uses progress bar for monitoring
    %   - Stores real bin widths in TensorWithVariableBinInfo.realBinWidth

    % Initialize struct
    TensorWithVariableBinInfo = obj.TensorWithVariableBinInfo;
    TensorWithVariableBinInfo.window = window;
    TensorWithVariableBinInfo.expectedBinWidth = expectedBinWidth;
    TensorWithVariableBinInfo.smoothGaussianStdInBin = smoothGaussianStdInBin;
    CS = obj.CS;
    MS = obj.MS;

    [nNeu, nCond] = size(CS);
    nMaxTrial = max(cellfun(@numel, CS(:)));

    MarkerTensorForVariableBin = nan(nNeu, nCond, nMaxTrial, range(window)+1);
    TensorWithVariableBinInfo.realBinWidth = nan(nNeu, nCond, nMaxTrial, range(window));


    for neu = 1:nNeu
        for cond = 1:nCond
            for trial = 1:size(CS{neu, cond},2)
                MarkerTensorForVariableBin(neu,cond,trial,:) = MS{neu,cond}{1,trial}(window(1):window(2));
            end
        end
    end

    MarkerTensorForVariableBin = squeeze(mean(MarkerTensorForVariableBin, 1,'omitnan'));

    meanEpochDurations = squeeze(mean(MarkerTensorForVariableBin(:,:,2:end)-MarkerTensorForVariableBin(:,:,1:end-1),[1 2], 'omitnan'));

    nBins = round(meanEpochDurations./expectedBinWidth);
    nBins(nBins==0) = 1;
    TensorWithVariableBin = nan(nNeu, nCond, nMaxTrial, sum(nBins));

    for neu = 1:nNeu
        for cond = 1:nCond
            for trial = 1:size(CS{neu, cond},2)
                tmp_trial = [];
                skip_epo = false;
                for ep = 1:length(nBins)
                    currentMarkers = fillmissing(MS{neu,cond}{1,trial},'next');
                    ff = currentMarkers(ep+1+window(1));
                    if skip_epo
                        ii = currentMarkers(ep-1+window(1));
                        realBinWidth = range([ii ff])/ (nBins(ep-1)+nBins(ep));
                        skip_epo = false;
                        deviation = 0;
                    else
                        ii = currentMarkers(ep+window(1));
                        realBinWidth = range([ii ff])/nBins(ep);
                        TensorWithVariableBinInfo.realBinWidth(neu,cond,trial,ep) = realBinWidth;
                        deviation = (realBinWidth-expectedBinWidth)/expectedBinWidth*100;
                    end

                    if deviation < -50 && (ep < length(nBins))
                        skip_epo = true;
                        continue
                    end

                    time_edges = ii:realBinWidth:ff;
                    tmp_trial = [tmp_trial histcounts(CS{neu,cond}{1,trial},time_edges)./realBinWidth];
                    TensorWithVariableBinInfo.realBinWidth(neu, cond, trial, ep) = realBinWidth;
                end

                if smoothGaussianStdInBin > 0
                    tmp_trial = smoothdata(tmp_trial,'gaussian', smoothGaussianStdInBin*5);
                end

                TensorWithVariableBin(neu,cond,trial,:) = tmp_trial;
            end
        end
    end
    TensorWithVariableBinCondAvg = squeeze(mean(TensorWithVariableBin,3,'omitmissing'));
    TensorWithVariableBinInfo.nBins = nBins;

    % Assign results to object properties
    obj.TensorWithVariableBin = TensorWithVariableBin;
    obj.TensorWithVariableBinCondAvg = TensorWithVariableBinCondAvg;
    obj.MarkerTensorForVariableBin = MarkerTensorForVariableBin;
    obj.TensorWithVariableBinInfo = TensorWithVariableBinInfo;
end
