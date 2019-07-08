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
 


%% Test the Network(Mutation)


val_mut = readtable('mutation_out_val_final.csv','ReadRowNames',true);
save('val_mut.mat','val_mut')

load('val_mut.mat')
[val_d,val_t1] = xlsread('clin_out_val_for_mutation.xlsx');
val_resp= categorical(val_d(:,2));
val_mut.class=val_resp;
validation_table_cnv= val_mut;

T_val= xlsread('valresponse_mut.xlsx');
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
%clear 'vars';
val_mutclinical=(table2array(val_clinical_table))';
%val_cnv=table2array(val_cnv);

net.numinputs = 5;
net.inputConnect = [0 0 0 1 0 ;0 0 0 0 0 ;0 0 0 0 0];
view(net)
%net.inputs=3;
%X = {val_cnv.',val_cnvclinical'}
A=zeros(144,64);
B=zeros(126,64);
D=zeros(78,64);
E=zeros(7,64);
F=zeros(6,64);


X=(table2array(validation_table_cnv(:,1:7))');
%mutX={A; B; D; X; val_mutclinical};
mutX={A; B; D; X; F};
%Y=grp2idx(val_resp)';
%y_val = sim(net, rnaX);
y_val = net(mutX);
figure()

view(net)
plotconfusion(t_val, y_val)
set(0,'DefaultLineLineWidth',4);
set(0,'DefaultaxesLineWidth',4);
set(0,'DefaultaxesFontSize',18);
set(0,'DefaultTextFontSize',18)
plotroc(t_val, y_val)

%% Test the Network(mirna)
val_mirna = readtable('mirna_out_val_final.csv','ReadRowNames',true);
save('val_mirna.mat','val_mirna')

load('val_mirna.mat')
[val_d,val_t1] = xlsread('clin_out_val_for_mirna.xlsx');
val_resp= categorical(val_d(:,2));
val_mirna.class=val_resp;
validation_table_rna= val_mirna;

T_val= xlsread('valresponse_mirna.xlsx');
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
%clear 'vars';
val_mirnaclinical=(table2array(val_clinical_table))';
%val_cnv=table2array(val_cnv);

net.numinputs = 5;
net.inputConnect = [0 0 1 0 1 ;0 0 0 0 0 ;0 0 0 0 0];
view(net)
%net.inputs=3;
%X = {val_cnv.',val_cnvclinical'}
A=zeros(144,8);
B=zeros(126,8);
D=zeros(7,8);
E=zeros(6,8);
%F=zeros(6,163);


X=(table2array(validation_table_rna(:,1:78))');
mirnaX={A; B; X; D; val_mirnaclinical};
%mirnaX={A; B; X; D; E};
%Y=grp2idx(val_resp)';
%y_val = sim(net, rnaX);
y_val = net(mirnaX);
figure()

view(net)
plotconfusion(t_val, y_val)
plotroc(t_val, y_val)