%% load data
cnv = xlsread('cnv_sel_features_144genes.xlsx').';
mRna = xlsread('rna_sel_features_126genes.xlsx').';
mirna = xlsread('mirna_sel_features_78genes.xlsx').';
%meth = xlsread('meth_sel_features_1genes.xlsx').';
mut = xlsread('mutation_sel_features_7genes.xlsx').';
cli = xlsread('clinical_processed.xlsx').';

resp = xlsread('response.xlsx');

x = {cnv; mRna; mirna; mut; cli};

%% nn for selected_data

rng(1);
trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation.
% Create a Pattern Recognition Network
hiddenLayerSize = [9 5];
net = patternnet(hiddenLayerSize, trainFcn);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 75/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 10/100;

%set crossentropy as performance eval
net.performFcn = 'crossentropy';    
net.performParam.regularization = 0.0822;
net.performParam.normalization = 'none';

%net = feedforwardnet(10,trainFcn);
net.numinputs = 5;
net.inputConnect = [1 1 1 1 1 ; 0 0 0 0 0; 0 0 0 0 0];
%net.inputConnect = [1 1 1 1 1; 0 0 0 0 0];
% Train the Network
r = resp.';
[net,tr] = train(net,x,r);
%view(net)
y = net(x)
plotconfusion(r, y)
plotroc(r, y)
 
%% Test the Network  
val_cnv = readtable('cnv_val_out23.csv','ReadRowNames',true);
val_mut = readtable('mutation_val_out23.csv','ReadRowNames',true);

val_CNV= table2array(val_cnv);
val_MUT= table2array(val_mut);

val_resp = xlsread('val_response.xlsx');

net.numinputs = 1;
X = val_CNV.';
Y = val_resp.';
y_val = sim(net, X);
figure()

view(net)
%plotconfusion(Y, y_val)
%y = net(x);
%e = gsubtract(r,y);
%performance = perform(net,r,y)
