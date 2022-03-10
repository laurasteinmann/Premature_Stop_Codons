# Analysis of Premature Stop Codons
# Libraries
import pandas as pd 

#Function

#Analysis of Single-Stop -Gene based Stop Distribution
single_stop_codons_table = pd.read_csv('../data/processed/Genexpression_Differences/Type2_Stop_List_Significant_Bonferroni.csv', index_col=0)

gene_based_stop_codons_table = pd.read_csv('../data/processed/Genexpression_Differences/Gene_Based_Bonferroni_Significant_Lower_Expression.csv', index_col=0)


stop_codons_attributes = pd.DataFrame(columns=['Stops_per_Gene'], index = gene_based_stop_codons_table.Gene)


for stop in gene_based_stop_codons_table.index:
    entry = gene_based_stop_codons_table.loc[stop,]
    gene = entry['Gene']
    subset = single_stop_codons_table[single_stop_codons_table['Gene']==gene]
    number_stops = subset.shape[0]
    stop_codons_attributes.loc[gene,'Stops_per_Gene'] = number_stops

print('Maximum Premature Stop Codons per Gene:',  stop_codons_attributes.Stops_per_Gene.max())
print('Minimum Premature Stop Codons per Gene:', stop_codons_attributes.Stops_per_Gene.min())

stop_codons_attributes.to_csv('../data/processed/Premature_Stop_Codons_Attributes/Premature_Stop_Attributes.csv')