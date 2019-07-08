import os
import requests
import numpy as np
import pandas as pd
import functools

from decimal import Decimal

from sklearn import predfessing

## Input Files
miRna_file = os.path.join('LUAD_miRNA_RPM_log2_Zscaled_imputed.csv')
mRna_file = os.path.join('LUAD_mRNA_log2_zscale_imputed.txt')
protein_file = os.path.join('LUAD_protein_imputed.txt')
CNV_file = os.path.join('LUAD_CNV_gistic_imputed.txt')
mutation_file = os.path.join('LUAD_mutation_output_minus_germline.csv')
methylation_file = os.path.join('LUAD_methylation_imputed.txt')
clinical_file = os.path.join('clinical_out_filtered4NA.csv')
#mutation_file_250 = os.path.join('mutation_out.tsv')
## Output Files
###dfessing 20 removal
miRna_out_file = os.path.join( 'miRna_out.tsv')
mRna_out_file = os.path.join( 'mRna_out.tsv')
protein_out_file = os.path.join( 'protein_out.tsv')
CNV_out_file = os.path.join( 'CNV_gistic_out.tsv')
mutation_out_file = os.path.join( 'mutation_out.tsv')
methylation_out_file = os.path.join( 'methylation_out.tsv')
clinical_out_file = os.path.join( 'clinical_out.tsv')
CNV_out_file_val = os.path.join('data', 'CNV_gistic_out_val.csv')
### Load Data

rna_df = pd.read_csv(mRna_file, index_col=0, sep='\t')
mirna_df = pd.read_csv(miRna_file, index_col=0)
protein_df = pd.read_csv(protein_file, index_col=0, sep='\t')
CNV_df = pd.read_csv(CNV_file, index_col=0, sep='\t')
mutation_df = pd.read_csv(mutation_file, index_col=0)
methylation_df = pd.read_csv(methylation_file, index_col=0, sep='\t')
clinical_df = pd.read_csv(clinical_file, index_col=0)
#mutation_df_250= pd.read_csv(mutation_file_250, index_col=0, sep='\t')
clinical_df_T = clinical_df.T
##commom patient id
#common = set.intersection(set(rna_df.columns), set(mirna_df.columns), set(CNV_df.columns), set(mutation_df.columns), set(methylation_df.columns), set(clinical_df_T.columns))
# 294 patients
common = set.intersection(set(rna_df.columns), set(mirna_df.columns), set(protein_df.columns), set(CNV_df.columns), set(mutation_df.columns), set(clinical_df_T.columns))
# 198 patients

##print (common)
##Then keep only those samples whose values are within the set of common values:
df1= rna_df.loc[:, ~rna_df.columns.isin(common)]
df2= mirna_df.loc[:, ~mirna_df.columns.isin(common)]
df3= protein_df.loc[:, ~protein_df.columns.isin(common)]
df4= CNV_df.loc[:, ~CNV_df.columns.isin(common)]
df5= mutation_df.loc[:, ~mutation_df.columns.isin(common)]
df6= methylation_df.loc[:, ~methylation_df.columns.isin(common)]
df7= clinical_df_T.loc[:, ~clinical_df_T.columns.isin(common)]

## common between clinical and CNV
common = set.intersection(set(clinical_df_T.columns), set(CNV_df.columns))
df8= df4.loc[:, df4.columns.isin(common)]
df8_T = df8.T
df8_T.sort_index()
CNV_out_file_val=os.path.join('data', 'cnv_out_val.csv')
df8.to_csv(CNV_out_file_val)
df9= df7.loc[:, df7.columns.isin(common)]
df9_T = df9.T
df9_T.sort_index()
clin_file_val=os.path.join('data', 'clin_out_val_for_CNV.csv')
df9.to_csv(clin_file_val)


## common between clinical and RNA
common = set.intersection(set(clinical_df_T.columns), set(rna_df.columns))
df10= df1.loc[:, df1.columns.isin(common)]
df10_T = df10.T
df10_T.sort_index()
RNA_out_file_val=os.path.join('data', 'rna_out_val.csv')
df10.to_csv(RNA_out_file_val)
df11= df7.loc[:, df7.columns.isin(common)]
df11_T = df11.T
df11_T.sort_index()
clin_file_val_for_rna=os.path.join('data', 'clin_out_val_for_rna.csv')
df11.to_csv(clin_file_val_for_rna)

## common between clinical and miRNA
common = set.intersection(set(clinical_df_T.columns), set(mirna_df.columns))
df12= df2.loc[:, df2.columns.isin(common)]
df12_T = df12.T
df12_T.sort_index()
miRNA_out_file_val=os.path.join('data', 'mirna_out_val.csv')
df12.to_csv(miRNA_out_file_val)
df13= df7.loc[:, df7.columns.isin(common)]
df13_T = df13.T
df13_T.sort_index()
clin_file_val_for_mirna=os.path.join('data', 'clin_out_val_for_mirna.csv')
df13.to_csv(clin_file_val_for_mirna)

## common between clinical and Protein
common = set.intersection(set(clinical_df_T.columns), set(protein_df.columns))
df14= df3.loc[:, df3.columns.isin(common)]
df14_T = df14.T
df14_T.sort_index()
protein_out_file_val=os.path.join('data', 'protein_out_val.csv')
df14.to_csv(protein_out_file_val)
df15= df7.loc[:, df7.columns.isin(common)]
df15_T = df15.T
df15_T.sort_index()
clin_file_val_for_protein=os.path.join('data', 'clin_out_val_for_protein.csv')
df15.to_csv(clin_file_val_for_protein)

## common between clinical and mutation
common = set.intersection(set(clinical_df_T.columns), set(mutation_df.columns))
df16= df5.loc[:, df5.columns.isin(common)]
df16_T = df16.T
df16_T.sort_index()
mutation_out_file_val=os.path.join('data', 'mutation_out_val.csv')
df16.to_csv(mutation_out_file_val)
df17= df7.loc[:, df7.columns.isin(common)]
df17_T = df17.T
df17_T.sort_index()
clin_file_val_for_mutation=os.path.join('data', 'clin_out_val_for_mutation.csv')
df17.to_csv(clin_file_val_for_mutation)

##write dataframes in output files
#df1.to_csv(mRna_out_file)
#df2.to_csv(miRna_out_file)
#df3.to_csv(protein_out_file)
#df4.to_csv(CNV_out_file)
#df5.to_csv(mutation_out_file)
#df6.to_csv(methylation_out_file)
#df7.to_csv(clinical_out_file)
#df8.to_csv(CNV_out_file_val)
