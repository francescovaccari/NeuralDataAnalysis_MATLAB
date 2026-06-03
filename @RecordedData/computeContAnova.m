function anovaResults = computeContAnova(obj, neus2consider, conds2consider, refWindow, anaWindow, binWidth)
    % COMPUTECONTANOVA - Compute ANOVA comparing baseline and analysis windows
    %
    % INPUTS:
    %   obj                     - RecordedData object
    %   neus2consider           - Neuron indices (numeric vector)
    %   conds2consider          - Condition indices (numeric vector)
    %   refWindow               - Reference [start(s), marker, end(s)] (numeric vector, length 3)
    %   anaWindow               - Analysis [start(s), marker, end(s)] (numeric vector, length 3)
    %   binWidth                - [bin_size, step] in seconds (numeric vector, length 2)
    %
    % OUTPUT:
    %   anovaResults            - Struct containing:
    %                             - p2way: (neurons × bins × 3) p-values [epoch, condition, epoch×condition]
    %                             - p1way: (neurons × bins) p-values [condition]
    %                             - sign2way: (bins × 3) % neurons with p<0.05 (2-way)
    %                             - sign1way: (bins) % neurons with p<0.05 (1-way)
    %
    % DESCRIPTION:
    %   Two-way ANOVA (epoch × condition) comparing baseline vs. analysis window.
    %   One-way ANOVA (condition only) within analysis window.

    if (refWindow(3)-refWindow(1)) ~= binWidth(1)
        disp("NOTE: The binWidth chosen for the analysis window doesn't match the size of the reference window")
    end
    
    % Call prepareTensorWithFixBin to prepare reference data
    refObj = obj.prepareTensorWithFixBin(refWindow, [refWindow(3)-refWindow(1) 0], 0);
    refData = refObj.TensorWithFixBin(neus2consider, conds2consider, :, :);

     % Call prepareTensorWithFixBin to prepare data to be analyzed
    anaObj = obj.prepareTensorWithFixBin(anaWindow, binWidth, 0);
    anaData = anaObj.TensorWithFixBin(neus2consider, conds2consider, :, :);


    p2way = nan(numel(neus2consider),size(anaData,4),3);
    p1way = nan(numel(neus2consider),size(anaData,4),1);
    
    for neu = 1:numel(neus2consider)
        for t = 1:size(anaData,4)
            refY = refData(neu,:,:);
            anaY = anaData(neu,:,:,t);
            conds = nan(size(anaY));
            for cond = 1:numel(conds2consider)
                conds(1,cond,:) = cond;
            end
            % 2way ANOVA (epoch, cond, epoch*cond)
            y = [refY(:); anaY(:)];
            g1 = [zeros(numel(refY),1); ones(numel(anaY),1)];
            g2 = [conds(:); conds(:)];
            p2way(neu,t,:) = anovan(y,{g1, g2}, 'model','interaction','varnames',{'epoch','condition'}, 'display', 'off')';

            % 1way ANOVA (cond)
            y = [anaY(:)];
            g = [conds(:)];
            p1way(neu,t) = anovan(y, g,'varnames',{'condition'}, 'display', 'off')';
        end
    end

    anovaResults.p2way = p2way;
    anovaResults.p1way = p1way;
    anovaResults.sign2way = nan(size(anaData,4),3);
    anovaResults.sign1way = nan(size(anaData,4),1);
    
    for t = 1:size(anaData,4)
        for vvar = 1:3
            anovaResults.sign2way(t,vvar) = length(find(p2way(:,t,vvar)<0.05))/numel(p2way(:,t,vvar))*100;
        end
        anovaResults.sign1way(t) = length(find(p1way(:,t)<0.05))/numel(p1way(:,t))*100;
    end

    figure
    plot(anovaResults.sign2way, 'LineWidth', 3)
    if ~isempty(obj.MarkerTensorForFixBin)
        hold on
        xline(squeeze(mean(obj.MarkerTensorForFixBin(neus2consider,conds2consider,:,:),[1 2 3],'omitmissing')))
    end
    ylabel("% of units modulated")
    xlabel("bins")
    title(["2wayAnova: baseline interval " num2str(refWindow) " / analysis interval " num2str(anaWindow)...
        " divided in bins of " num2str(binWidth(1)) " with an overlap of " num2str(binWidth(2))])
    legend({"epoch","cond","epoch*cond"})

    figure
    plot(anovaResults.sign1way, 'LineWidth', 3)
    if ~isempty(obj.MarkerTensorForFixBin)
        hold on
        xline(squeeze(mean(obj.MarkerTensorForFixBin(neus2consider,conds2consider,:,:),[1 2 3],'omitmissing')))
    end
    ylabel("% of units modulated")
    xlabel("bins")
    title(["1wayAnova: analysis interval " num2str(anaWindow)...
        " divided in bins of " num2str(binWidth(1)) " with an overlap of " num2str(binWidth(2))])
    legend({"cond"})
end
