function obj = prepareTensorWithFixBin(obj, window, binWidth, neus2consider, conds2consider)
    % PREPARETENSORWITHFIXBIN - Prepare tensor with fixed-width binning
    %
    % INPUTS:
    %   obj                     - RecordedData object
    %   window                  - [start(s), alignment_marker, end(s)] (numeric vector, length 3)
    %   binWidth                - [bin_width, overlap] in seconds (numeric vector, length 2)
    %   neus2consider           - Neuron indices to include (numeric/logical vector)
    %   conds2consider          - Condition indices to include (numeric/logical vector)
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
    %   - Overlapping bins computed using moving mean
    %   - Uses progress bar for monitoring
    %   - binWidth(1): bin size; binWidth(2): overlap (0 for no overlap)
    
    % Extract data from object
    CS = obj.CS;
    MS = obj.MS;

    TensorWithFixBinInfo = obj.TensorWithFixBinInfo;
    TensorWithFixBinInfo.window = window;
    TensorWithFixBinInfo.binWidth = binWidth;

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
    nMarker = numel(MS{neus2consider(1), conds2consider(1)}{1,1});
    if binWidth(2) > 0 %partially overlapping bins
        edges = window(1) : binWidth(2) : window(3);
    else
        edges = window(1) : binWidth(1) : window(3);
    end

    relativeBinCenters = (edges(1:end-1) + edges(2:end)) / 2;
    TensorWithFixBin = nan(nSelectedNeu, nSelectedCond, nMaxTrial, numel(edges)-1);
    BinTimestampsForFixBin = nan(nSelectedNeu, nSelectedCond, nMaxTrial, numel(edges)-1);
    MarkerTensorForFixBin = nan(nSelectedNeu, nSelectedCond, nMaxTrial, nMarker);
    TensorWithFixBinInfo.neus2consider = neus2consider;
    TensorWithFixBinInfo.conds2consider = conds2consider;

    % Progress bar
    progressBar = waitbar(0, 'Preparing fixed-bin tensor: Neuron 0 / ' + string(nSelectedNeu) + ', Condition 0 / ' + string(nSelectedCond));

    for neuIdx = 1:nSelectedNeu
        neu = neus2consider(neuIdx);
        for condIdx = 1:nSelectedCond
            cond = conds2consider(condIdx);
            % Update progress bar
            progress = (neuIdx - 1) / nSelectedNeu + (1 / nSelectedNeu) * (condIdx / nSelectedCond);
            waitbar(progress, progressBar, sprintf('Preparing fixed-bin tensor: Neuron %d / %d, Condition %d / %d', neuIdx, nSelectedNeu, condIdx, nSelectedCond));

            for trial = 1:size(CS{neu, cond},2)
                markers = MS{neu,cond}{1,trial};
                mksAligned2WinStart = markers - (markers(window(2)) + window(1)); %align markers with the beginning of the window of interest
                if binWidth(2) > 0
                    binWidth2use = binWidth(2);
                else
                    binWidth2use = binWidth(1);
                end
                tmpMks = (floor(abs(mksAligned2WinStart)/ binWidth2use)+1).*sign(mksAligned2WinStart); %vectorized implementation
                tmpMks(tmpMks==0) = 1; %bin number 0 have no sense
                MarkerTensorForFixBin(neuIdx, condIdx, trial, :) = tmpMks;
                %% similar but with a for loop
                % for m = 1:nMarker
                %     mkAligned2WinStart = markers(m) - (markers(window(2)) + window(1));
                %     if binWidth(2) > 0
                %         MarkerTensorForFixBin(neuIdx, condIdx, trial, m) = round( mkAligned2WinStart / binWidth(2)); %compute the bin of the marker even if it is outside the window
                %     else
                %         MarkerTensorForFixBin(neuIdx, condIdx, trial, m) = round( mkAligned2WinStart / binWidth(1)); %compute the bin of the marker even if it is outside the window
                %     end
                % 
                % end
                tmp_trial = histcounts(CS{neu,cond}{1,trial}, edges+markers(window(2)))./binWidth(1); %compute Firing Rate (sp/s)
                if binWidth(2) > 0 %partially overlapping bins
                    tmp_trial = movmean(tmp_trial,floor(binWidth(1)/binWidth(2)));
                end
                TensorWithFixBin(neuIdx,condIdx,trial,:) = tmp_trial;
                BinTimestampsForFixBin(neuIdx, condIdx, trial, :) = relativeBinCenters + markers(window(2));
            end
        end
    end

    TensorWithFixBinCondAvg = reshape(mean(TensorWithFixBin,3,'omitmissing'), ...
        nSelectedNeu, nSelectedCond, numel(edges)-1);
    MarkerTensorForFixBinCondAvg = reshape(mean(MarkerTensorForFixBin,3,'omitmissing'), ...
        nSelectedNeu, nSelectedCond, nMarker);
    close(progressBar);
    
    % Compute trial counts: N x C (number of valid trials per neuron per condition)
    TrialNumWithFixBin = cellfun(@(c) size(c, 2), CS(neus2consider, conds2consider));

    % Assign results to object properties
    obj.TensorWithFixBin = TensorWithFixBin;
    obj.TensorWithFixBinCondAvg = TensorWithFixBinCondAvg;
    obj.MarkerTensorForFixBin = MarkerTensorForFixBin;
    obj.MarkerTensorForFixBinCondAvg = MarkerTensorForFixBinCondAvg;
    obj.TensorWithFixBinInfo = TensorWithFixBinInfo;
    obj.TrialNumWithFixBin = TrialNumWithFixBin;
    obj.BinTimestampsForFixBin = BinTimestampsForFixBin;

end
