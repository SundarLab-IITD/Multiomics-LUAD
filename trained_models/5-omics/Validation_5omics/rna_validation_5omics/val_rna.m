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
net.performParam.regularization = 0.082;
net.performParam.normalization = 'none';

%net = feedforwardnet(10,trainFcn);
net.numinputs = 5;
net.inputConnect = [1 1 1 1 1 ; 0 0 0 0 0; 0 0 0 0 0];
%net.inputConnect = [1 1 1 1 1; 0 0 0 0 0];
% Train the Network
r = resp.';
[net,tr] = train(net,x,r);
view(net)
 
%% Test the Network(RNA)
val_rna = readtable('rna_out_val_final.csv','ReadRowNames',true);
save('val_rna.mat','val_rna')

load('val_rna.mat')
[val_d,val_t1] = xlsread('clin_out_val_for_rna.xlsx');
val_resp= categorical(val_d(:,2));
val_rna.class=val_resp;
validation_table_rna= val_rna;

T_val= xlsread('valresponse_rna.xlsx');
t_val= T_val';
val_age = val_d(:,[1])';
gender = categorical(val_t1(2:end,9));
stage=categorical(val_t1(2:end,5));
tumor=categorical(val_t1(2:end,6));
node=categorical(val_t1(2:end,7));
metastasis=categorical(val_t1(2:end,8));
t1=grp2idx(node);
t2=grp2idx(stage);
t3=grp2idx(tumor);
t4=grp2idx(gender);
t5=grp2idx(metastasis);

val_clinical_table = table(val_age',t1,t2,t3,t4,t5);
%vars = {'t1', 't2', 't3', 't4', 't5', 'gender', 'stage', 'tumor','node','metastasis','val_age'};
%clear 'vars';;
%F=zeros(6,163);

val_rnaclinical=(table2array(val_clinical_table))';
%val_cnv=table2array(val_cnv);

net.numinputs = 5;
net.inputConnect = [0 1 0 0 1 ;0 0 0 0 0 ;0 0 0 0 0];
view(net)
%net.inputs=3;
%X = {val_cnv.',val_cnvclinical'}
A=zeros(144,67);
B=zeros(78,67);
D=zeros(7,67);
E=zeros(6,67);
%F=zeros(6,163);


X=(table2array(validation_table_rna(:,1:126))');
rnaX={A; X; B; D; val_rnaclinical};
%rnaX={A; X; B; D; E};
%Y=grp2idx(val_resp)';
%y_val = sim(net, rnaX);
y_val = net(rnaX);
figure()

view(net)
plotconfusion(t_val, y_val)
plotroc(t_val, y_val)

%% Test the Network(CNV)


