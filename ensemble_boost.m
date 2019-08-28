%% I/O
X = input;
Y = categorical(output);
classOrder = unique(Y);
rng (1); % For reproducibility
classOrder()

%% Train a classifier
% This code specifies all the classifier options and trains the classifier.

t = templateTree('MaxNumSplits', 20);
enMdl = fitcensemble(X,Y,'Method', 'RUSBoost', 'Learners',t, ...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus'))

%% crossval
crossv = crossval(enMdl, 'KFold', 10)
[validationPredictions, validationScores] = kfoldPredict(crossv);

%% loss
validationAccuracy = 1 - kfoldLoss(crossv, 'LossFun', 'ClassifError');      
validationAccuracy()

%Y = cellstr(Y)
predictedY = resubPredict(enMdl);

plotconfusion(validationPredictions, Y)
plotroc(validationPredictions, Y)
