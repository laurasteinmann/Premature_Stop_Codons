# Coocurrence of Premature Stop Codons

#Libraries
import pandas as pd 
import numpy as np
from scipy import stats
#Functions
def generate_cooccurrence_table(gene_list):
    cooccurrence_table = pd.DataFrame(columns=['Gene1', 'Gene2', 'Gene1_count',
                                              'Gene2_count', 'Cooccurrence',
                                              'Lower_border', 'Upper_border',
                                              'pvalue_lower', 'pvalue_upper', 'Distance'])
    gene1_list = []
    gene2_list = []
    genes = gene_list.Gene
    for gene in genes:
        for gene2 in genes:
            gene1_list.append(gene)
            gene2_list.append(gene2)

    cooccurrence_table['Gene1'] = gene1_list
    cooccurrence_table['Gene2'] = gene2_list
    return cooccurrence_table

def prework_ara11(ara11_list):
    gene_names = ara11_list['Name']
    new_genes = []
    for gene in gene_names:
        new_gene = gene[:9]
        new_genes.append(new_gene)
    ara11_list['New_Genes'] = new_genes
    return ara11_list

def statistics_in_coocurrence(cooccurrence_table, gene_list, pvaluethreshold, ara11_list):
    for row in range(0, cooccurrence_table.shape[0]):
        gene = cooccurrence_table.iloc[row,0]
        attributes_gene1 = gene_list[gene_list.Gene == gene]
        cooccurrence_table.iloc[row,2] = int(attributes_gene1.KO_num.values)
    for row in range(0, cooccurrence_table.shape[0]):
        gene = cooccurrence_table.iloc[row,1]
        attributes_gene2 = gene_list[gene_list.Gene == gene]
        cooccurrence_table.iloc[row,3] = int(attributes_gene2.KO_num.values)
    for row in range(0, cooccurrence_table.shape[0]):
        gene1 = cooccurrence_table.iloc[row,0]
        gene2 = cooccurrence_table.iloc[row,1]
        attributes_gene1 = gene_list[gene_list.Gene == gene1]
        attributes_gene2 = gene_list[gene_list.Gene == gene2]
        gene1_accessions = np.asarray(attributes_gene1.KO_acc.values[0])
        gene2_accessions = np.asarray(attributes_gene2.KO_acc.values[0])
        cooccurrence = gene1_accessions[np.isin(gene1_accessions, gene2_accessions)]
        cooccurrence_table.iloc[row, 4] = int(len(cooccurrence))
        cooccurrence_table['Cooccurrence'] = cooccurrence_table['Cooccurrence'].fillna(0)
#    for row in range(0, cooccurrence_table.shape[0]):
#        table_attributes = cooccurrence_table.iloc[row,:]
#        total_accessions = 665
#        gene2_count = table_attributes.Gene2_count
#        gene1_count = table_attributes.Gene1_count
#        coocurrence = table_attributes.Cooccurrence
#        pvalue_under = stats.hypergeom(M=total_accessions, n=gene2_count, N=gene1_count).sf(coocurrence)
#        pvalue_over = 1 - stats.hypergeom(M=total_accessions, n=gene2_count, N=gene1_count).sf(coocurrence)
#        cooccurrence_table.iloc[row, 7] = pvalue_under
#        cooccurrence_table.iloc[row, 8] = pvalue_over
#        confidence_interval = stats.hypergeom(M=total_accessions, n=gene2_count, N=gene1_count).interval(pvaluethreshold)
#        cooccurrence_table.iloc[row, 5] = confidence_interval[0]
#        cooccurrence_table.iloc[row, 6] = confidence_interval[1]
    for row in range(0, cooccurrence_table.shape[0]):
        gene1 = cooccurrence_table.iloc[row,0]
        gene2 = cooccurrence_table.iloc[row,1]
        if gene1 == gene2:
            distance = 0
        elif gene1[2] == gene2[2]:
            dim_test1 = ara11_list[ara11_list['New_Genes'] == gene1]
            dim_test2 = ara11_list[ara11_list['New_Genes'] == gene2]
            if dim_test1.shape[0] == 0 or dim_test2.shape[0] == 0:
                distance = np.nan
            else:
                ara11_gene1 = ara11_list[ara11_list['New_Genes'] == gene1].iloc[0,]
                ara11_gene2 = ara11_list[ara11_list['New_Genes'] == gene2].iloc[0,]
                if ara11_gene1.Start_Protein < ara11_gene2.Start_Protein:
                    distance = np.abs(ara11_gene2.Start_Protein - ara11_gene1.Stop_Protein)
                else:
                    distance = np.abs(ara11_gene1.Start_Protein - ara11_gene2.Stop_Protein)
        elif gene1[2] != gene2[2]:
            distance = 'diff_Chr'
        cooccurrence_table.iloc[row, 9] = distance
    return cooccurrence_table

#Analysis
ara11 = pd.read_csv('../data/processed/Ara11/Ara11_genes_gff.csv', index_col=0)
coocurrence_table_significant = generate_cooccurrence_table(gene_premature_stop_list)
ara11 = prework_ara11(ara11)
coocurrence_table_significant = statistics_in_coocurrence(coocurrence_table_significant, gene_premature_stop_list, 0.05, ara11)
coocurrence_table_significant.to_csv('../data/processed/Coocurrence/Coocurrence_Table_Significant_Full.csv')

coocurrence_table_bonferroni = generate_cooccurrence_table(gene_premature_stop_list)
pvalue = 0.05/coocurrence_table_bonferroni.shape[0]
coocurrence_table_bonferroni = statistics_in_coocurrence(coocurrence_table_bonferroni, gene_premature_stop_list, pvalue, ara11)
coocurrence_table_bonferroni.to_csv('../data/processed/Coocurrence/Coocurrence_Table_Bonferroni_Significant_Full.csv')