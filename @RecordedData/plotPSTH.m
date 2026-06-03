function plotPSTH(obj, neu2consider, conds2consider, window, bin, vars2plot)
    % PLOTPSTH - Plot peri-stimulus time histogram with raster and optional variables
    %
    % SYNTAX:
    %   plotPSTH(obj, neu2consider, conds2consider, window, bin, vars2plot)
    %
    % INPUTS:
    %   obj                 - RecordedData object                     
    %   neu2consider        - Single neuron index (numeric scalar, positive integer)
    %   conds2consider      - Condition indices (numeric vector)
    %   window              - [onset, marker, offset] in seconds (numeric vector, length 3)
    %   bin                 - Bin size in seconds (numeric scalar, positive)
    %   vars2plot           - Cell array of variable names {'EY', 'TC', ...} or empty {}
    %
    %
    % OUTPUT:
    %   Figure with histogram, raster plot, and optional behavioral variables
    %
    % DESCRIPTION:
    %   Comprehensive visualization of neural activity aligned to a specific marker.
    %   Displays mean firing rate histogram with raster plot showing all trials.

    % Validate inputs
    if ~isscalar(neu2consider) || neu2consider < 1
        error('neu2consider must be a single positive neuron index');
    end
    
    if ~isvector(conds2consider) || any(conds2consider < 1)
        error('conds2consider must be a vector of positive condition indices');
    end
    
    if numel(window) ~= 3
        error('window must be a 3-element vector [onset, marker, offset]');
    end
    
    if bin <= 0
        error('bin must be positive (in seconds)');
    end
    
    if nargin < 6 || isempty(vars2plot)
        vars2plot = {};
    end

    % Get number of conditions
    nConds = numel(conds2consider);
    nVars = numel(vars2plot);
    
    % Calculate optimal grid layout
    nCols = ceil(sqrt(nConds));
    nRows = ceil(nConds / nCols);
    
    % Determine number of time ticks
    if nConds <= 4
        nTimeTicks = 15;
    elseif nConds <= 9
        nTimeTicks = 10;
    else
        nTimeTicks = 5;
    end

    % Exctract CS and MS from obj
    CS = obj.CS;
    MS = obj.MS;

    % Create edges for histogram binning
    edges = window(1):bin:window(3);
    nBins = numel(edges) - 1;
    bin_centers = edges(1:end-1) + bin/2;
    
    % Preallocate firing rate data
    meanFiringRates = [];
    
    % First pass: compute mean firing rates for all conditions to get global max
    for condIdx = 1:nConds
        cond = conds2consider(condIdx);
        nTrials = numel(CS{neu2consider, cond});
        
        condFiringRates = nan(nTrials, nBins);
        
        for trial = 1:nTrials
            markers = MS{neu2consider, cond}{trial};
            spikes = CS{neu2consider, cond}{trial};
            
            % Align to marker
            alignmentTime = markers(window(2));
            alignedSpikes = spikes - alignmentTime;
            
            % Compute firing rate (spikes per second)
            counts = histcounts(alignedSpikes, edges);
            condFiringRates(trial, :) = counts / bin;
        end
        
        % Compute mean for this condition
        condMeanFR = mean(condFiringRates, 1, 'omitnan');
        meanFiringRates = [meanFiringRates; condMeanFR];
    end
    
    % Get global maximum firing rate
    globalMaxFR = max(meanFiringRates(:), [], 'omitnan');
    
    % Create figure
    fig = figure('Name', sprintf('PSTH - Neuron %d (bin = %.3f s)', neu2consider, bin), ...
                 'NumberTitle', 'off');
    sgtitle(sprintf('PSTH - Neuron %d (bin = %.3f s)', neu2consider, bin));

    % Total rows per condition
    if nVars == 0
        rowsPerCond = 2;
    else
        rowsPerCond = 2 + nVars;
    end
    totalGridRows = nRows * rowsPerCond;
    totalGridCols = nCols;
    
    % Plot each condition
    for condIdx = 1:nConds
        cond = conds2consider(condIdx);
        
        % Calculate grid position
        [gridRow, gridCol] = ind2sub([nRows, nCols], condIdx);
        
        % Base subplot index
        baseGridIdx = (gridRow - 1) * rowsPerCond * totalGridCols + (gridCol - 1) + 1;
        
        nTrials = numel(CS{neu2consider, cond});
        condFiringRates = nan(nTrials, nBins);
        
        for trial = 1:nTrials
            markers = MS{neu2consider, cond}{trial};
            spikes = CS{neu2consider, cond}{trial};
            
            alignmentTime = markers(window(2));
            alignedSpikes = spikes - alignmentTime;
            
            counts = histcounts(alignedSpikes, edges);
            condFiringRates(trial, :) = counts / bin;
        end
        
        % Mean firing rate
        meanFR = mean(condFiringRates, 1, 'omitnan');
        
        % ===== HISTOGRAM =====
        ax_hist = subplot(totalGridRows, totalGridCols, baseGridIdx);
        bar(ax_hist, bin_centers, meanFR, 'FaceColor', 'k', 'EdgeColor', 'none', 'BarWidth', 1);
        hold on
        xline(0,'--')
        ax_hist.YLim = [0, (globalMaxFR * 1.1)+1];
        ax_hist.XLim = [window(1), window(3)];
        ax_hist.YLabel.String = 'Firing Rate (sp/s)';
        ax_hist.Title.String = ['Condition ' num2str(cond)];
        
        % Add fine grid
        ax_hist.XGrid = 'on';
        ax_hist.YGrid = 'on';
        ax_hist.GridLineStyle = ':';
        ax_hist.GridAlpha = 0.2;
        
        % Set x-ticks
        timeTicks = linspace(window(1), window(3), nTimeTicks);
        ax_hist.XTick = timeTicks;
        ax_hist.XTickLabel = arrayfun(@(x) sprintf('%.2f', x), timeTicks, 'UniformOutput', false);
        
        % ===== RASTER PLOT =====
        ax_raster = subplot(totalGridRows, totalGridCols, baseGridIdx + totalGridCols);
        hold(ax_raster, 'on');
        
        markerSize = max(5, min(5, 17 / sqrt(nTrials * nConds)));
        
        % Generate HSV colormap for markers
        hsvColors = hsv(size(MS{neu2consider, conds2consider(1)}{1}, 2));
        
        for trial = 1:nTrials
            markers = MS{neu2consider, cond}{trial};
            spikes = CS{neu2consider, cond}{trial};
            
            alignmentTime = markers(window(2));
            alignedSpikes = spikes - alignmentTime;
            alignedMarkers = markers - alignmentTime;
            
            % Plot spikes
            validSpikes = alignedSpikes(alignedSpikes >= window(1) & alignedSpikes <= window(3));
            for s = 1:numel(validSpikes)
                plot(ax_raster, [validSpikes(s) validSpikes(s)], ...
                     [trial - 0.4, trial + 0.4], 'k', 'LineWidth', 1);
            end
            
            % Plot markers
            for m = 1:numel(alignedMarkers)
                markerTime = alignedMarkers(m);
                if markerTime >= window(1) && markerTime <= window(3)
                    colorIdx = mod(m - 1, size(hsvColors, 1)) + 1;
                    plot(ax_raster, markerTime, trial, 'o', ...
                         'MarkerFaceColor', hsvColors(colorIdx, :), ...
                         'MarkerEdgeColor', hsvColors(colorIdx, :), ...
                         'MarkerSize', markerSize);
                end
            end
        end
        
        ax_raster.YDir = 'reverse';
        ax_raster.XLim = [window(1), window(3)];
        ax_raster.YLim = [0, nTrials + 1];
        ax_raster.YLabel.String = 'Trial';
        ax_raster.XLabel.String = 'Time (s)';
        
        % Add grid
        ax_raster.XGrid = 'on';
        ax_raster.YGrid = 'on';
        ax_raster.GridLineStyle = ':';
        ax_raster.GridAlpha = 0.15;
        
        % Set x-ticks
        timeTicks = linspace(window(1), window(3), nTimeTicks);
        ax_raster.XTick = timeTicks;
        ax_raster.XTickLabel = arrayfun(@(x) sprintf('%.2f', x), timeTicks, 'UniformOutput', false);
        
        box(ax_raster, 'on');
        
        % ===== ADDITIONAL VARIABLES =====
        for varIdx = 1:nVars
            subplotIdx = baseGridIdx + (varIdx + 1) * totalGridCols;
            ax_var = subplot(totalGridRows, totalGridCols, subplotIdx);
            
            varName = vars2plot{varIdx};
            
            if (isprop(obj, varName) || isfield(obj, varName))
                varData = obj.(varName);
                
                if ~isempty(varData) && iscell(varData)
                    hold(ax_var, 'on');
                    
                    for trial = 1:nTrials
                        markers = MS{neu2consider, cond}{trial};
                        alignmentTime = markers(window(2));
                        
                        % Handle TC (touchscreen coordinates)
                        if strcmp(varName, 'TC') && isprop(obj, 'TCTS')
                            if size(varData,1) >= neu2consider && size(varData,2) >= cond && ...
                               ~isempty(varData{neu2consider, cond}) && numel(varData{neu2consider, cond}) >= trial && ...
                               isstruct(varData{neu2consider, cond}{trial}) && isfield(varData{neu2consider, cond}{trial}, 'xCursor')
                                
                                xCoords = varData{neu2consider, cond}{trial}.xCursor;
                                yCoords = varData{neu2consider, cond}{trial}.yCursor;
                                
                                tcts = obj.TCTS{neu2consider, cond}{trial};
                                if ~isempty(tcts) && numel(tcts) == numel(xCoords)
                                    alignedTS = tcts - alignmentTime;
                                    
                                    validIdx = alignedTS >= window(1) & alignedTS <= window(3);
                                    alignedTS_valid = alignedTS(validIdx);
                                    xCoords_valid = xCoords(validIdx);
                                    yCoords_valid = yCoords(validIdx);
                                    
                                    if ~isempty(xCoords_valid)
                                        plot(ax_var, alignedTS_valid, xCoords_valid, 'Color', [1 0 0 0.3]);
                                    end
                                    if ~isempty(yCoords_valid)
                                        plot(ax_var, alignedTS_valid, yCoords_valid, 'Color', [0 0 1 0.3]);
                                    end
                                    xline(0,'--')
                                    legend({'x','y'})
                                end
                            end
                        else
                            % General case for other variables
                            if numel(varData) >= neu2consider && numel(varData{neu2consider}) >= cond && ...
                               ~isempty(varData{neu2consider, cond}) && numel(varData{neu2consider, cond}) >= trial
                                
                                varTrace = varData{neu2consider, cond}{trial};
                                
                                if ~isempty(varTrace)
                                    timeAxis = linspace(window(1), window(3), numel(varTrace));
                                    plot(ax_var, timeAxis, varTrace, 'Color', [0.5 0.5 0.5 0.3]);
                                end
                            end
                        end
                    end
                    
                    ax_var.Title.String = varName;
                    ax_var.XLabel.String = 'Time (s)';
                    ax_var.YLabel.String = varName;
                    ax_var.XLim = [window(1), window(3)];
                    ax_var.XGrid = 'on';
                    ax_var.YGrid = 'on';
                    ax_var.GridLineStyle = ':';
                    ax_var.GridAlpha = 0.15;
                    
                    timeTicks = linspace(window(1), window(3), nTimeTicks);
                    ax_var.XTick = timeTicks;
                    ax_var.XTickLabel = arrayfun(@(x) sprintf('%.2f', x), timeTicks, 'UniformOutput', false);
                    
                    box(ax_var, 'on');
                end
            end
        end
    end
end
