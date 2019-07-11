X = input;
Y = output;
classOrder = unique(Y);
rng (1); % For reproducibility
classOrder()

%% optm/mdlgenratn
t = templateSVM(...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true );
opts = struct('Optimizer','bayesopt','ShowPlots',true,...
    'AcquisitionFunctionName','expected-improvement-plus');
PMdl = fitcecoc(X,Y, 'Learners', t, 'Classnames',classOrder,...
    'OptimizeHyperparameters', 'all', 'HyperparameterOptimizationOptions' , opts)

%% crossval
crossv = crossval(PMdl, 'KFold', 10)
[validationPredictions, validationScores] = kfoldPredict(crossv);

%% loss
validationAccuracy = 1 - kfoldLoss(crossv, 'LossFun', 'ClassifError');
validationAccuracy()
