function dPCAResults = computedPCA(obj, combinedParams, numComps, options)
    % COMPUTEDPCA - Compute demixed principal component analysis (dPCA)
    %
    % INPUTS:
    %   obj                 - RecordedData object
    %   combinedParams      - Parameter combinations (cell array)
    %   numComps            - Number of components (numeric scalar)
    %   options             - Optional analysis options (struct, optional)
    %
    % OUTPUT:
    %   dPCAResults         - Struct containing:
    %                         - W: decoder matrix weights (time × components)
    %                         - V: encoder matrix weights (features × components)
    %                         - explVar: explained variance per component
    %                         - whichMarg: component marginalization assignments
    %
    % DESCRIPTION:
    %   Performs demixed PCA to decompose neural population activity into components
    %   that selectively represent each parameter.
    
    Tensor4dPCA = obj.Tensor4dPCA;
    Tensor4dPCAInfo = obj.Tensor4dPCAInfo;

    if nargin < 5 || isempty(options)
        options = struct();
    end

    if ~isfield(options, 'margColours') || isempty(options.margColours)
        options.margColours = jet(numel(combinedParams));
    end

    if ~isfield(options, 'margNames') || isempty(options.margNames)
        for i=1:numel(combinedParams)
            options.margNames{i} = ['MargNum' num2str(i)];
        end
    end

    if ~isfield(options, 'timeMarginalization') || isempty(options.timeMarginalization)
        options.timeMarginalization = ndims(Tensor4dPCA);
    end

    if ~isfield(options, 'time') || isempty(options.time)
        if Tensor4dPCAInfo.fixedBinFlag
            options.time = (0:size(Tensor4dPCA,ndims(Tensor4dPCA))-1).*(Tensor4dPCAInfo.binWidth(1)/2);
        else
            options.time = (0:size(Tensor4dPCA,ndims(Tensor4dPCA))-1).*Tensor4dPCAInfo.expectedBinWidth;
        end
    end

    if ~isfield(options, 'timeEvents') || isempty(options.timeEvents)
        if isfield(Tensor4dPCAInfo, 'MarkersFixedBin') && Tensor4dPCAInfo.fixedBinFlag
            options.timeEvents = interp1(1:numel(options.time), options.time, Tensor4dPCAInfo.MarkersFixedBin)';
        elseif ~Tensor4dPCAInfo.fixedBinFlag && isfield(Tensor4dPCAInfo, 'nBins')
            options.timeEvents = cumsum(Tensor4dPCAInfo.nBins)'.*Tensor4dPCAInfo.expectedBinWidth;
        else
            options.timeEvents = [];
        end
    end

    [W, V, options.whichMarg] = dpca(obj.Tensor4dPCA, numComps, ...
        'combinedParams', combinedParams);

    explVar = dpca_explainedVariance(obj.Tensor4dPCA, W, V, ...
        'combinedParams', combinedParams);

    dpca_plot(obj.Tensor4dPCA, W, V, @dpca_plot_default, ...
        'whichMarg', options.whichMarg, ...
        'explainedVar', explVar, ...
        'marginalizationNames', options.margNames, ...
        'marginalizationColours', options.margColours, ...
        'time', options.time, ...
        'timeEvents', options.timeEvents, ...
        'timeMarginalization', options.timeMarginalization,...
        'legendSubplot', 16);

    dPCAResults.W = W;
    dPCAResults.V = V;
    dPCAResults.whichMarg = options.whichMarg;
    dPCAResults.explVar = explVar;
end
