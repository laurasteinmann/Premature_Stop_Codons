# Generation of the GWAS Phenotypes 
# Libraries
import pandas as pd
# Functions

    
# Analysis
coocurrence_table = pd.read_csv('../data/processed/Coocurrence/Coocurrence_Table_Bonferroni_Selected_Interactions.csv', index_col=0)

for entry in coocurrence_table.index:
    gene1 = coocurrence_table.loc[entry, 'Gene1']
    gene2 = coocurrence_table.loc[entry, 'Gene2']
    gene1_accessions = gene_premature_stop_list[gene_premature_stop_list['Gene']==gene1].iloc[0,5]
    gene2_accessions = gene_premature_stop_list[gene_premature_stop_list['Gene']==gene2].iloc[0,5]
    gene2_accessions = gene2_accessions.astype('int')
    accessions = gene1_accessions.astype('int')
    phenotype = pd.DataFrame(index=accessions, columns=['Phenotype'])
    for row in phenotype.index:
        if np.isin(row, gene2_accessions):
            phenotype.loc[row, 'Phenotype'] = 1
        else:
            phenotype.loc[row, 'Phenotype'] = 0
    
    filename = '../data/processed/Coocurrence/Phenotypes/Over_Coocurrence/' + gene1 + '_' + gene2 + '.csv'
    phenotype.to_csv(filename)