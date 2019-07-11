load('val_rna.mat')
[val_d,val_t1] = xlsread('clin_out_val_for_cnv.xlsx');
val_resp= categorical(val_d(:,2));