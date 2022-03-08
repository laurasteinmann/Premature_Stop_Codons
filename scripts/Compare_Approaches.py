# Comparison of the two Approaches 
# Libraries
import pandas as pd 
import numpy as np
#Functions

#Analysis for small list
stop_list_geneexp = pd.read_csv('../data/processed/Genexpression_Differences/Type2_Stop_List_Significant_Bonferroni.csv', index_col=0)
stop_list_length = pd.read_csv('../data/processed/Relative_Lengths/Stop_Codons_Equal_Number_Approach1.csv', index_col=0)

codons_geneexp = []
for entry in stop_list_geneexp.index:
    id = str(stop_list_geneexp.loc[entry, 'Chrom']) + "-" + str(stop_list_geneexp.loc[entry, 'Pos'])
    codons_geneexp.append(id)

codons_length = []
for entry in stop_list_length.index:
    id = str(stop_list_length.loc[entry, 'Chrom']) + "-" + str(stop_list_length.loc[entry, 'Pos'])
    codons_length.append(id)

codons_overlap = stop_list_length[np.isin(codons_length, codons_geneexp)]
number_overlap = len(codons_overlap)
number_geneexp = len(codons_geneexp)
number_length = len(codons_length)


# Analysis with 30 percent list length

stop_list_length_full = pd.read_csv('../data/processed/Relative_Lengths/Stop_Codons_30percent_Threshold.csv', index_col=0)

codons_length_full = []
for entry in stop_list_length_full.index:
    id = str(stop_list_length_full.loc[entry, 'Chrom']) + "-" + str(stop_list_length_full.loc[entry, 'Pos'])
    codons_length_full.append(id)

codons_overlap_big = stop_list_length_full[np.isin(codons_length_full, codons_geneexp)]
number_overlap_big = len(codons_overlap_big)
number_length_full = len(codons_length_full)
