# Approach 2: Calculation of the remaining mRNA length
import pandas as pd 
import numpy as np
#Functions 

#Analysis
stop_codons = pd.read_csv('../data/preprocessed/Stop_Table_Full.csv', index_col=0)
ara11 = pd.read_csv('../data/processed/Ara11/Ara11_genes_gff.csv', index_col=0)
ara11['Short_Name'] = np.nan
for row in ara11.index:
    name = ara11.loc[row, 'Name']
    short_name = name[:9]
    ara11.loc[row, 'Short_Name'] = short_name 

results_table = pd.DataFrame(columns=['Chrom', 'Pos', 'Gene', 'Start_Protein', 'Stop_Protein', 'Length', 'Relative_Length'])
results_table['Chrom'] = stop_codons['Chrom']
results_table['Pos'] = stop_codons['Pos']
results_table['Gene'] = stop_codons['Gene']

for entry in results_table.index:
    gene = results_table.loc[entry, 'Gene']
    gene_attributes = ara11[ara11['Short_Name']==gene]
    if gene_attributes.shape[0] >= 1:
        gene_attributes = gene_attributes.iloc[0,]
        results_table.loc[entry,'Start_Protein'] = gene_attributes['Start_Protein']
        results_table.loc[entry, 'Stop_Protein'] = gene_attributes['Stop_Protein']
    else:
        pass

results_table = results_table.dropna(subset=['Start_Protein', 'Stop_Protein'])

for entry in results_table.index:
    stop_attributes = results_table.loc[entry,]
    if stop_attributes.Start_Protein < stop_attributes.Stop_Protein:
        length = stop_attributes.Pos - stop_attributes.Start_Protein
        full_length = stop_attributes.Stop_Protein - stop_attributes.Start_Protein
        relative_length = length/full_length * 100
    else:
        length = stop_attributes.Pos - stop_attributes.Stop_Protein
        full_length = stop_attributes.Start_Protein - stop_attributes.Stop_Protein 
        relative_length = length/full_length * 100

    results_table.loc[entry, 'Length'] = length
    results_table.loc[entry, 'Relative_Length'] = relative_length 

results_table.to_csv('../data/processed/Relative_Lengths/Stop_Codons_Full_List_Lengths.csv')
