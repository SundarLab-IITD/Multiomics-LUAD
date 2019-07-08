import os
import requests
import numpy as np
import pandas as pd
import functools

from decimal import Decimal

from sklearn import predfessing


#Input file
miRna_file = os.path.join('mirna_out_val.csv')
mRna_file = os.path.join('rna_out_val.csv')
protein_file = os.path.join('protein_out_val.csv')

mutation_file = os.path.join('mutation_out_val.csv')

#output files 
#miRna_out_file = os.path.join('NCA_outputs', 'mirna_gene_names_78.csv')
#mRna_out_file = os.path.join('NCA_outputs', 'rna_gene_names_126.csv')
#protein_out_file = os.path.join('NCA_outputs', 'protein_out_val.csv')
#CNV_out_file = os.path.join( 'NCA_outputs','cnv_gene_names_144.csv')
#mutation_out_file = os.path.join('NCA_outputs', 'mutation_gene_names_7.csv')

#Load gene files
CNV_gene_file = os.path.join('NCA_outputs','cnv_gene_names_144.csv')
cnv_gene=pd.read_csv(CNV_gene_file)
mirna_gene_file = os.path.join('NCA_outputs','mirna_gene_names_78.csv')
mirna_gene=pd.read_csv(mirna_gene_file)
mutation_gene_file = os.path.join('NCA_outputs','mutation_gene_names_7.csv')
mutation_gene=pd.read_csv(mutation_gene_file)
#protein_gene_file = os.path.join('NCA_outputs','protein_genes.csv')
#protein_gene=pd.read_csv(protein_gene_file)
rna_gene_file = os.path.join('NCA_outputs','rna_gene_names_126.csv')
rna_gene=pd.read_csv(rna_gene_file)


#Load exp files
rna_df = pd.read_csv(mRna_file, index_col=0)
mirna_df = pd.read_csv(miRna_file, index_col=0)
protein_df = pd.read_csv(protein_file, index_col=0)
CNV_df = pd.read_csv(CNV_file, index_col=0)
mutation_df = pd.read_csv(mutation_file, index_col=0)


#Extract Common Gene Names for CNV
CNV_file = os.path.join('cnv_val_genenameedited.csv')
CNV_gene_file = os.path.join('NCA_outputs','cnv_gene_names_144.csv')
CNV_df = pd.read_csv(CNV_file, index_col=0)
cnv_gene=pd.read_csv(CNV_gene_file)
#CNV_df= CNV_df.drop_duplicates(keep='first')
CNV_df_T = CNV_df.T
cnv_df_T=CNV_df_T.sort_index();

common = set.intersection(set(cnv_df_T.columns), set(cnv_gene.columns))
df1= cnv_df_T.loc[:, cnv_df_T.columns.isin(common)]
df1.shape
sequence = list(cnv_gene.columns.values)
df2 = df1.reindex(columns=sequence)
CNV_out_file_val=os.path.join('../', 'Output_val_files', 'cnv_out_val_final.csv')
df2.to_csv(CNV_out_file_val)

#Extract Common Gene Names for RNA
RNA_file = os.path.join('rna_val_genenameedited.csv')
rna_gene_file = os.path.join('NCA_outputs','rna_gene_names_126.csv')
RNA_df = pd.read_csv(RNA_file, index_col=0)
rna_gene=pd.read_csv(rna_gene_file)
#CNV_df= CNV_df.drop_duplicates(keep='first')
RNA_df_T = RNA_df.T
rna_df_T=RNA_df_T.sort_index();

common = set.intersection(set(rna_df_T.columns), set(rna_gene.columns))
df3= rna_df_T.loc[:, rna_df_T.columns.isin(common)]
df3.shape
sequence = list(rna_gene.columns.values)
df4 = df3.reindex(columns=sequence)
RNA_out_file_val=os.path.join('../', 'Output_val_files', 'rna_out_val_final.csv')
df4.to_csv(RNA_out_file_val)


#Extract Common Gene Names for miRNA
miRNA_file = os.path.join('mirna_val_genenameedited.csv')
mirna_gene_file = os.path.join('NCA_outputs','mirna_gene_names_78.csv')
miRNA_df = pd.read_csv(miRNA_file, index_col=0)
mirna_gene=pd.read_csv(mirna_gene_file)
#CNV_df= CNV_df.drop_duplicates(keep='first')

miRNA_df_T = miRNA_df.T
mirna_df_T=miRNA_df_T.sort_index();

common = set.intersection(set(mirna_df_T.columns), set(mirna_gene.columns))
df5= mirna_df_T.loc[:, mirna_df_T.columns.isin(common)]
df5.shape
sequence = list(mirna_gene.columns.values)
df6 = df5.reindex(columns=sequence)
miRNA_out_file_val=os.path.join('../', 'Output_val_files', 'mirna_out_val_final.csv')
df6.to_csv(miRNA_out_file_val)

#Extract Common Gene Names for mutation
mut_file = os.path.join('mutation_val_genenameedited.csv')
mut_gene_file = os.path.join('NCA_outputs','mutation_gene_names_7.csv')
MUT_df = pd.read_csv(mut_file, index_col=0)
mut_gene=pd.read_csv(mut_gene_file)
#CNV_df= CNV_df.drop_duplicates(keep='first')
MUT_df_T = MUT_df.T
mut_df_T=MUT_df_T.sort_index();

common = set.intersection(set(mut_df_T.columns), set(mut_gene.columns))
df7= mut_df_T.loc[:, mut_df_T.columns.isin(common)]
df7.shape
sequence = list(mut_gene.columns.values)
df8 = df7.reindex(columns=sequence)
MUT_out_file_val=os.path.join('../', 'Output_val_files', 'mutation_out_val_final.csv')
df8.to_csv(MUT_out_file_val)


#Extract Common Gene Names for protein
prot_file = os.path.join('protein_val_genenameedited.csv')
prot_gene_file = os.path.join('../', 'Input','protein_genes.csv')
PROT_df = pd.read_csv(prot_file, index_col=0)
prot_gene=pd.read_csv(prot_gene_file)
#CNV_df= CNV_df.drop_duplicates(keep='first')
PROT_df_T = PROT_df.T
prot_df_T=PROT_df_T.sort_index();

common = set.intersection(set(prot_df_T.columns), set(prot_gene.columns))
df9= prot_df_T.loc[:, prot_df_T.columns.isin(common)]
df9.shape
sequence = list(prot_gene.columns.values)
df9 = df8.reindex(columns=sequence)
PROT_out_file_val=os.path.join('../', 'Output_val_files', 'protein_out_val_final.csv')
df9.to_csv(PROT_out_file_val)





