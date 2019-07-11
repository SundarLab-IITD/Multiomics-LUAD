X = features;
%labels = cellstr(labels)
Y = categorical(labels);
classOrder = unique(Y);
rng(1); 
classOrder()

%% Train a classifier
% This code specifies all the classifier options and trains the classifier.
opts = struct('Optimizer','bayesopt','ShowPlots',true,...
    'AcquisitionFunctionName','expected-improvement-plus');

rng(1)
Mdl = fitcsvm(...
    X, ...
    Y, ...
    'KernelFunction', 'polynomial', ...
    'Optimizehyperparameters', 'auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'),...
    'PolynomialOrder', 2, ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true, ...
    'ClassNames', classOrder);


%% Create the result struct with predict function

rng (10);
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
svmPredictFcn = @(x) predict(Mdl, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationSVM = Mdl;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2019a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 373 columns because this model was trained using 373 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationSVM, 'KFold', 10);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
validationAccuracy

%confmat=confusionmat(X(:,1),validationPredictions);
%% figure
predictedY = resubPredict(Mdl);

plotconfusion(validationPredictions, predictedY)
plotroc(validationPredictions, predictedY)


