# Control with Synonmymous and Nonsynonymous SNP mutations
# Library
import pandas as pd
# Function
def generate_summary_table():
    table = pd.DataFrame(columns=['Chrom', 'Pos', 'ID'])
    return table

def fill_snps(table, data):
    for snp in data.index:
        table.loc[snp,'Chrom'] = data.loc[snp,'CHROM']
        table.loc[snp, 'Pos'] = data.loc[snp, 'POS']
        table.loc[snp,'ID'] = str(table.loc[snp,'Chrom']) + '-' + str(table.loc[snp, 'Pos'])
    return table
# Analysis
fixed_section = pd.read_csv('data/preprocessed/Fixed_Section_1001_Genomes.csv', index_col=0)
summary_table = generate_summary_table()
summary_table = fill_snps(summary_table, fixed_section)
