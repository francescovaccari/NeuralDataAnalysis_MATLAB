function obj = prepareTensorWithVariableBin(obj, window, expectedBinWidth, neus2consider, conds2consider)
    % PREPARETENSORWITHVARIABLEBIN - Prepare tensor with epoch-adaptive binning
    %
    % INPUTS:
    %   obj                     - RecordedData object
    %   window                  - [start_marker, end_marker] indices (numeric vector, length 2)
    %   expectedBinWidth        - Target bin duration in seconds (numeric scalar)
    %   neus2consider           - Neuron indices to include (numeric/logical vector)
    %   conds2consider          - Condition indices to include (numeric/logical vector)
    %
    % OUTPUT:
    %   obj                     - Updated RecordedData object
    %
    % DESCRIPTION:
    %   Creates binned neural activity where bin widths adapt to epoch durations. Bins are
    %   automatically merged if they fall below 50% of expected width.
    %
    % REMARKS:
    %   - Firing rates in spikes/second
    %   - Uses progress bar for monitoring
    %   - Stores real bin widths in TensorWithVariableBinInfo.realBinWidth

    % Initialize struct
    TensorWithVariableBinInfo = obj.TensorWithVariableBinInfo;
    TensorWithVariableBinInfo.window = window;
    TensorWithVariableBinInfo.expectedBinWidth = expectedBinWidth;
    CS = obj.CS;
    MS = obj.MS;

    if nargin < 5 || isempty(neus2consider) || isempty(conds2consider)
        error('neus2consider and conds2consider are required and cannot be empty.');
    end
    if islogical(neus2consider)
        neus2consider = find(neus2consider);
    end
    if islogical(conds2consider)
        conds2consider = find(conds2consider);
    end
    neus2consider = reshape(neus2consider, 1, []);
    conds2consider = reshape(conds2consider, 1, []);

    nSelectedNeu = numel(neus2consider);
    nSelectedCond = numel(conds2consider);
    selectedCS = CS(neus2consider, conds2consider);
    nMaxTrial = max(cellfun(@numel, selectedCS(:)));

    nMarker = range(window)+1;
    nEpoch = range(window);
    allMarkers = nan(nSelectedNeu, nSelectedCond, nMaxTrial, nMarker);
    TensorWithVariableBinInfo.realBinWidth = nan(nSelectedNeu, nSelectedCond, nMaxTrial, nEpoch);
    TensorWithVariableBinInfo.neus2consider = neus2consider;
    TensorWithVariableBinInfo.conds2consider = conds2consider;


    for neuIdx = 1:nSelectedNeu
        neu = neus2consider(neuIdx);
        for condIdx = 1:nSelectedCond
            cond = conds2consider(condIdx);
            for trial = 1:size(CS{neu, cond},2)
                allMarkers(neuIdx,condIdx,trial,:) = MS{neu,cond}{1,trial}(window(1):window(2));
            end
        end
    end

    epochDurations = allMarkers(:,:,:,2:end)-allMarkers(:,:,:,1:end-1);
    medianEpochDurations = reshape(median(epochDurations, [1 2 3], 'omitnan'), 1, []);

    nBins = round(medianEpochDurations./expectedBinWidth);
    nBins(nBins==0) = 1;

    markerBinPositions = [0 cumsum(nBins)];
    MarkerTensorForVariableBin = repmat(reshape(markerBinPositions, 1, 1, 1, nMarker), ...
        nSelectedNeu, nSelectedCond, nMaxTrial, 1);
    MarkerTensorForVariableBinCondAvg = repmat(reshape(markerBinPositions, 1, 1, nMarker), ...
        nSelectedNeu, nSelectedCond, 1);

    TensorWithVariableBin = nan(nSelectedNeu, nSelectedCond, nMaxTrial, sum(nBins));

    for neuIdx = 1:nSelectedNeu
        neu = neus2consider(neuIdx);
        for condIdx = 1:nSelectedCond
            cond = conds2consider(condIdx);
            for trial = 1:size(CS{neu, cond},2)
                tmp_trial = [];
                skip_epo = false;
                for ep = 1:length(nBins)
                    currentMarkers = fillmissing(MS{neu,cond}{1,trial},'next');
                    % Epoch ep spans from marker (window(1)+ep-1) to marker (window(1)+ep)
                    ff = currentMarkers(ep+window(1));
                    if skip_epo
                        ii = currentMarkers(ep-2+window(1));
                        realBinWidth = range([ii ff])/ (nBins(ep-1)+nBins(ep));
                        skip_epo = false;
                        deviation = 0;
                    else
                        ii = currentMarkers(ep-1+window(1));
                        realBinWidth = range([ii ff])/nBins(ep);
                        TensorWithVariableBinInfo.realBinWidth(neuIdx,condIdx,trial,ep) = realBinWidth;
                        deviation = (realBinWidth-expectedBinWidth)/expectedBinWidth*100;
                    end

                    if deviation < -50 && (ep < length(nBins))
                        skip_epo = true;
                        continue
                    end

                    time_edges = ii:realBinWidth:ff;
                    tmp_trial = [tmp_trial histcounts(CS{neu,cond}{1,trial},time_edges)./realBinWidth];
                    TensorWithVariableBinInfo.realBinWidth(neuIdx, condIdx, trial, ep) = realBinWidth;
                end

                TensorWithVariableBin(neuIdx,condIdx,trial,:) = tmp_trial;
            end
        end
    end
    TensorWithVariableBinCondAvg = reshape(mean(TensorWithVariableBin,3,'omitmissing'), ...
        nSelectedNeu, nSelectedCond, sum(nBins));
    TensorWithVariableBinInfo.nBins = nBins;

    % Compute trial counts: N x C (number of valid trials per neuron per condition)
    TrialNumWithVariableBin = cellfun(@(c) size(c, 2), CS(neus2consider, conds2consider));

    % Assign results to object properties
    obj.TensorWithVariableBin = TensorWithVariableBin;
    obj.TensorWithVariableBinCondAvg = TensorWithVariableBinCondAvg;
    obj.MarkerTensorForVariableBin = MarkerTensorForVariableBin;
    obj.MarkerTensorForVariableBinCondAvg = MarkerTensorForVariableBinCondAvg;
    obj.TensorWithVariableBinInfo = TensorWithVariableBinInfo;
    obj.TrialNumWithVariableBin = TrialNumWithVariableBin;
end
