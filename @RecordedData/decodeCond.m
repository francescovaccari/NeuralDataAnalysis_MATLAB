function decodingScores = decodeCond(obj, neus2consider, conds2consider, window)
    % DECODECOND - Decode conditions using Naïve Bayes classifier
    %
    % INPUTS:
    %   obj             - RecordedData object
    %   neus2consider   - Neuron indices (numeric vector)
    %   conds2consider  - Condition indices (numeric vector)
    %   window          - [start(s), marker, end(s)] relative to marker (numeric vector, length 3)
    %
    % OUTPUT:
    %   decodingScores  - Struct with performance metrics (confusion matrix displayed as figure)
    %
    % DESCRIPTION:
    %   Performs 10-fold cross-validation classification using Naïve Bayes decoder.
    %   Spike counts within specified window are used as features.
    
    % Extract data from object
    CS = obj.CS;
    MS = obj.MS;

    preprocessedData = [];
    labels = [];

    nNeus2consider = numel(neus2consider);

    for cond = conds2consider
        i = 0;
        nTrials = numel(CS{1, cond});
        condData = nan(nTrials, nNeus2consider);
        for neu = neus2consider
            i = i+1;

            for trial = 1:nTrials
                markers = MS{neu, cond}{trial};
                spikes = CS{neu, cond}{trial};

                ii = markers(window(2)) + window(1);
                ff = markers(window(2)) + window(3);
                if ff <= ii
                    error('RecordedData:InvalidWindowRange', ...
                        'Window end must be larger than window start for neuron %d, condition %d, trial %d.', ...
                        neu, cond, trial);
                end
                condData(trial, i) = length(find(spikes > ii & spikes < ff));
            end
        end
        preprocessedData = [preprocessedData; condData];
        labels = [labels; repelem(cond, nTrials, 1)];
    end

    hpartition = cvpartition(labels,KFold=10);

    X = preprocessedData+rand(size(preprocessedData))*10^-8;
    tbl = array2table(X); 
    tbl.Y = labels;

    for cross = 1:10
        idxTrain = training(hpartition,cross);
        tblTrain = tbl(idxTrain,:);
        idxTest = test(hpartition,cross);
        tblTest = tbl(idxTest,:);

        mdl = fitcnb(tblTrain,"Y");
        trueLabels{cross} = tblTest.Y;
        predictedLabels{cross} = predict(mdl,tblTest);
    end
    trueLabels = cell2mat(trueLabels');
    predictedLabels = cell2mat(predictedLabels');

    figure
    C = confusionchart(trueLabels, predictedLabels);
    decodingScores = statsOfMeasure(confusionmat(trueLabels,predictedLabels), false);

end
