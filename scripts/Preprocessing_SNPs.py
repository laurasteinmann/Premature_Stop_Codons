# Control with Synonmymous and Nonsynonymous SNP mutations
# Library
import pandas as pd
import numpy as np 
# Function
def generate_summary_table():
    table = pd.DataFrame(columns=['Chrom', 'Pos', 'ID', 'Type', 'ALT_Frequency'])
    return table

def fill_out_summary_table(table, data):
    table['Chrom'] = data['Chromosome']
    table['Pos'] = data['Position']
    table['ID'] = data['SNP']
    table['Type'] = data['SYN_vs_NON_SYN']
    return table

def filter_accessions(gt_section, expression_data):
    gt_section_first_filter = gt_section.iloc[:, 1:]
    expression_accessions = expression_data.index.values
    gt_accessions = gt_section_first_filter.columns.values
    expression_accessions = expression_accessions.astype(int)
    gt_accessions = gt_accessions.astype(int)
    gt_section_first_filter.columns = gt_accessions
    in_both = gt_accessions[np.isin(gt_accessions,expression_accessions)]
    gt_section_first_filter = gt_section_first_filter[in_both]
    return gt_section_first_filter

def make_ids(chrom, position):
    ids = np.empty(shape=0)
    for i in np.arange(0, len(chrom)):
        id = chrom[i] + '- ' + position[i]
        ids = np.append(ids, id)
    return ids
# Analysis
snps_attributes = pd.read_csv('../data/preprocessed/snps_2029.csv', index_col=0)
synonymous = snps_attributes[snps_attributes['SYN_vs_NON_SYN']=='SYN']
summary_table = generate_summary_table()
summary_table = fill_out_summary_table(summary_table, synonymous)
print('Summary table done:', summary_table.shape)
summary_table.to_csv('../data/processed/Synonymous_Table.csv')

gt_section = pd.read_csv('../data/preprocessed/GT_Section_1001_Genomes.csv', index_col=0)
print('GT_Section geladen:', gt_section.iloc[0:6,0:6])
expression_data = pd.read_csv('../data/preprocessed/727_Expressiondata.csv', index_col=0)
gt_section_665 = filter_accessions(gt_section, expression_data)
gt_section_665.to_csv('../data/preprocessed/665_GT_Section.csv')

fixed_section = pd.read_csv('../data/preprocessed/Fixed_Section_1001_Genomes.csv', index_col=0)
subset = fixed_section.iloc[400000:500000,:]
chrom = subset.CHROM.values
chrom = chrom.astype(str)
position = subset.POS.values
position = position.astype(str)
ids_array = make_ids(chrom, position)

fixed_section.iloc[400000:500000,2] = ids_array
fixed_section.to_csv('../data/preprocessed/Fixed_Section_1001_Genomes.csv')

#ids_synonymous = summary_table.ID.values
#ids_synonymous = ids_synonymous.astype(str)

#in_both = ids_synonymous[np.isin(ids_synonymous, id_array)]
#genomics2 = genomics.loc[in_both]
#genomics2.shape
#genomics2.to_csv('../data/preprocessed/Second_Filter_Step.csv')