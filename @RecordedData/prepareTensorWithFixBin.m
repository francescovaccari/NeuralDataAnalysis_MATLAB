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
    %   for intrinsic temporal smoothing. State markers tracked in BinnedMarkerTensorForFixBin
    %   (the bin index each marker falls into) and AlignedMarkerTensorForFixBin (the timing
    %   of each marker relative to the window's alignment marker, i.e. markers(window(2))).
    %   The relative bin-center time axis is exposed as
    %   TensorWithFixBinInfo.relativeBinCenters, in the same reference frame as
    %   AlignedMarkerTensorForFixBin, for consistent plotting downstream.
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
    if mod(binWidth(1),binWidth(2)) ~= 0 && binWidth(2) > 0
        warning(['Warning: when binning the spikes, the bin width should be a multiple of the step. ' ...
            'This is not the case: binWidth = ' num2str(binWidth(1)) ', step = ' num2str(binWidth(2)) '. '...
            'Bin width will be automatically adapted to allow vectorized computations.'] );
    end

    neus2consider = reshape(neus2consider, 1, []);
    conds2consider = reshape(conds2consider, 1, []);

    nSelectedNeu = numel(neus2consider);
    nSelectedCond = numel(conds2consider);
    selectedCS = CS(neus2consider, conds2consider);
    nMaxTrial = max(cellfun(@numel, selectedCS(:)));
    nMarker = numel(MS{neus2consider(1), conds2consider(1)}{1,1});
    if binWidth(2) > 0 %partially overlapping bins
        delta = mod(window(3)-window(1),binWidth(2))/2; %to avoid misallignemnt of central marker
        edges = window(1)+delta : binWidth(2) : window(3);
    else
        delta = mod(window(3)-window(1),binWidth(1))/2; %to avoid misallignemnt of central marker
        edges = window(1)+delta : binWidth(1) : window(3);
    end
    
    if binWidth(2) > 0
        relativeBinCenters = ([window(1)+delta:binWidth(2):window(3)-binWidth(2)] + [window(1)+binWidth(2)+delta:binWidth(2):window(3)]) / 2;
    else
        relativeBinCenters = (edges(1:end-1) + edges(2:end)) / 2;
    end
    % Expose the bin-center time axis (relative to the window's alignment
    % marker) so downstream plotting code (plotSDF, computeAnovaOnParam/Cond,
    % computedPCA, neuralDecodingClassification) can reuse this exact,
    % delta-corrected axis instead of re-deriving (and potentially
    % misaligning) their own.
    TensorWithFixBinInfo.relativeBinCenters = relativeBinCenters;
    TensorWithFixBin = nan(nSelectedNeu, nSelectedCond, nMaxTrial, numel(edges)-1);
    BinTimestampsForFixBin = nan(nSelectedNeu, nSelectedCond, nMaxTrial, numel(edges)-1);
    BinnedMarkerTensorForFixBin = nan(nSelectedNeu, nSelectedCond, nMaxTrial, nMarker);
    AlignedMarkerTensorForFixBin = nan(nSelectedNeu, nSelectedCond, nMaxTrial, nMarker); 
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
                mksAligned2WinStart = markers - (markers(window(2)) + window(1)); %align markers with the beginning of the window of interest (mks with negative times will before the window of interest)
                if binWidth(2) > 0
                    binWidth2use = binWidth(2);
                else
                    binWidth2use = binWidth(1);
                end
                tmpMks = (floor(abs(mksAligned2WinStart)/ binWidth2use)+1).*sign(mksAligned2WinStart); %vectorized implementation for assigning a bin to every mk (even those outside interest window)
                tmpMks(tmpMks==0) = 1; %bin number 0 have no sense
                BinnedMarkerTensorForFixBin(neuIdx, condIdx, trial, :) = tmpMks;
                %% similar but with a for loop
                % for m = 1:nMarker
                %     mkAligned2WinStart = markers(m) - (markers(window(2)) + window(1));
                %     if binWidth(2) > 0
                %         BinnedMarkerTensorForFixBin(neuIdx, condIdx, trial, m) = round( mkAligned2WinStart / binWidth(2)); %compute the bin of the marker even if it is outside the window
                %     else
                %         BinnedMarkerTensorForFixBin(neuIdx, condIdx, trial, m) = round( mkAligned2WinStart / binWidth(1)); %compute the bin of the marker even if it is outside the window
                %     end
                % 
                % end
                
                if binWidth(2) > 0 %partially overlapping bins
                    tmp_trial = histcounts(CS{neu,cond}{1,trial}, edges+markers(window(2)))./binWidth(2); %compute Firing Rate (sp/s)
                    tmp_trial = movmean(tmp_trial,[0 floor(binWidth(1)/binWidth(2))]); %1st binWidth(1)-wide bin will include n smaller bin from to 1 to n, and so on
                else
                    tmp_trial = histcounts(CS{neu,cond}{1,trial}, edges+markers(window(2)))./binWidth(1); %compute Firing Rate (sp/s)
                end
                TensorWithFixBin(neuIdx,condIdx,trial,:) = tmp_trial;
                BinTimestampsForFixBin(neuIdx, condIdx, trial, :) = relativeBinCenters + markers(window(2));
                AlignedMarkerTensorForFixBin(neuIdx, condIdx, trial, :) = markers - markers(window(2));
            end
        end
    end

    TensorWithFixBinCondAvg = reshape(mean(TensorWithFixBin,3,'omitmissing'), ...
        nSelectedNeu, nSelectedCond, numel(edges)-1);
    BinnedMarkerTensorForFixBinCondAvg = reshape(mean(BinnedMarkerTensorForFixBin,3,'omitmissing'), ...
        nSelectedNeu, nSelectedCond, nMarker);
    AlignedMarkerTensorForFixBinCondAvg = reshape(mean(AlignedMarkerTensorForFixBin,3,'omitmissing'), ...
        nSelectedNeu, nSelectedCond, nMarker);

    close(progressBar);
    
    % Compute trial counts: N x C (number of valid trials per neuron per condition)
    TrialNumWithFixBin = cellfun(@(c) size(c, 2), CS(neus2consider, conds2consider));

    % Assign results to object properties
    obj.TensorWithFixBin = TensorWithFixBin;
    obj.TensorWithFixBinCondAvg = TensorWithFixBinCondAvg;
    obj.BinnedMarkerTensorForFixBin = BinnedMarkerTensorForFixBin;
    obj.BinnedMarkerTensorForFixBinCondAvg = BinnedMarkerTensorForFixBinCondAvg;
    obj.AlignedMarkerTensorForFixBin = AlignedMarkerTensorForFixBin;
    obj.AlignedMarkerTensorForFixBinCondAvg = AlignedMarkerTensorForFixBinCondAvg;
    obj.TensorWithFixBinInfo = TensorWithFixBinInfo;
    obj.TrialNumWithFixBin = TrialNumWithFixBin;
    obj.BinTimestampsForFixBin = BinTimestampsForFixBin;

end
