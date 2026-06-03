function obj = prepareTensor4dPCA(obj, paramNames, conditions, fixedBinFlag, neus2consider, conds2consider)
    % PREPARETENSOR4DPCA - Prepare multi-dimensional tensor for dPCA
    %
    % INPUTS:
    %   obj                 - RecordedData object
    %   paramNames          - Parameter names (cell array of strings)
    %   conditions          - Condition table/struct with fields matching paramNames
    %   fixedBinFlag        - true for fixed bins, false for variable bins (logical)
    %   neus2consider       - Neuron indices (numeric vector)
    %   conds2consider      - Condition indices (numeric vector)
    %
    % OUTPUT:
    %   obj                 - Updated RecordedData object with newly filled:
    %            .Tensor4dPCA         - Multi-dimensional tensor (neurons × param1 × param2 × ... × time)
    %            .Tensor4dPCAInfo     - Struct with metadata (paramNames, params, fixedBinFlag, etc.)
    %
    % DESCRIPTION:
    %   Restructures condition-averaged tensor into multi-dimensional format where each
    %   parameter gets its own dimension. Allows dPCA to decompose neural activity by parameter.


    % Initialize Tensor4dPCAInfo struct
    obj.Tensor4dPCAInfo.paramNames = paramNames;
    obj.Tensor4dPCAInfo.fixedBinFlag = fixedBinFlag;

    if fixedBinFlag && ~isempty(obj.MarkerTensorForFixBin)
        obj.Tensor4dPCAInfo.MarkersFixedBin = squeeze(mean(obj.MarkerTensorForFixBin(neus2consider,conds2consider,:,:),[1 2 3],'omitnan'));
    end

    if fixedBinFlag
        obj.Tensor4dPCAInfo.binWidth = obj.TensorWithFixBinInfo.binWidth;
    else
        obj.Tensor4dPCAInfo.expectedBin = obj.TensorWithVariableBinInfo.expectedBin;
    end


    % Extract parameters from conditions based on paramNames
    params = [];
    for p = 1:numel(paramNames)
        if isstruct(conditions)
            params = [params, conditions.(paramNames{p})];
        elseif istable(conditions)
            params = [params, table2array(conditions(:, paramNames{p}))];
        else
            error('conditions must be a struct or table');
        end
    end

    obj.Tensor4dPCAInfo.params = params;

    % Select appropriate tensor based on fixedBinFlag
    if fixedBinFlag
        sourceTensor = obj.TensorWithFixBinCondAvg(neus2consider, conds2consider, :);
    else
        sourceTensor = obj.TensorWithVariableBinCondAvg(neus2consider, conds2consider, :);
    end

    % Initialize tensor with base dimensions (neuron x time)
    Tensor4dPCA = nan([size(sourceTensor, 1), size(sourceTensor, 3)]);

    % Add dimensions for each parameter
    for p = 1:numel(paramNames)
        Tensor4dPCA = repmat(Tensor4dPCA, [ones(1, ndims(Tensor4dPCA)), numel(unique(params(:, p)))]);
        if p==1
            Tensor4dPCA = permute(Tensor4dPCA, [1,  ndims(Tensor4dPCA), ndims(Tensor4dPCA)-1]);
        elseif p==2
            Tensor4dPCA = permute(Tensor4dPCA, [1:2, ndims(Tensor4dPCA), ndims(Tensor4dPCA)-1]);
        elseif p>2
            Tensor4dPCA = permute(Tensor4dPCA, [1, 2:p, ndims(Tensor4dPCA), ndims(Tensor4dPCA)-1]);
        end
    end

    % Populate tensor using dynamic indexing
    for cond = 1:size(sourceTensor,2)
        % Build dynamic index cell array
        indices = cell(1, numel(paramNames) + 2);
        indices{1} = ':';  % First dimension (neurons)
        for p = 1:numel(paramNames)
            indices{p+1} = params(cond, p);  % Parameter conditions
        end
        indices{end} = ':';  % Last dimension (time)

        % Assign data using dynamic indexing
        Tensor4dPCA(indices{:}) = sourceTensor(:, cond, :);
    end

    obj.Tensor4dPCA = Tensor4dPCA;
end
