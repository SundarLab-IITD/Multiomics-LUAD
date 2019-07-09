X = features;
Y = categorical(labels);
classOrder = unique(Y);
rng(1); 
classOrder()

%% optm
t = templateSVM('Standardize',true);
opts= struct('Optimizer', 'bayesopt', 'ShowPlots',true, 'AcquisitionFunctionName', 'expected-improvement-plus' );
PMdl = fitcecoc(X,Y,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'))

%% model generation
%t = templateSVM('Standardize',true);
%PMdl = fitcecoc(X,Y,'Holdout',0.30,'Learners',t,'ClassNames',classOrder,...
 %    'HyperparameterOptimizationOptions' ,opts);
Mdl = PMdl.Trained{1};           % Extract trained, compact classifier
testInds = test(PMdl.Partition);  % Extract the test indices
XTest = X(testInds,:);
YTest = Y(testInds,:);
idx = randsample(sum(testInds),20);
table(YTest(idx),labels(idx),...
    'VariableNames',{'TrueLabels','PredictedLabels'})

%% crossval

crossv = crossval(PMdl, 'KFold', 10)
[validationPredictions, validationScores] = kfoldPredict(crossv);

%% loss{
validationAccuracy = 1 - kfoldLoss(PMdl, 'LossFun', 'ClassifError');
validationAccuracy()
