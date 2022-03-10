# Analysis for the Coocurrence of Premature Stop Codons
import pandas as pd
# Libraries
# Functions 
# Analysis

coocurrence_table = pd.read_csv('../data/processed/Coocurrence/Coocurrence_Table_Bonferroni_Significant_Full.csv', index_col=0)

pvaluethreshold = 0.05/coocurrence_table.shape[0]

coocurrence_significant = coocurrence_table[coocurrence_table['pvalue_upper']<= pvaluethreshold]
coocurrence_significant = coocurrence_significant[coocurrence_significant['Gene1']!=coocurrence_significant['Gene2']]
coocurrence_significant = coocurrence_significant.replace('diff_Chr', 1000000)
coocurrence_significant["Distance"] = pd.to_numeric(coocurrence_significant["Distance"])
coocurrence_significant = coocurrence_significant[(coocurrence_significant['Distance'] > 10000)]
coocurrence_significant.to_csv('../data/processed/Coocurrence/Coocurrence_Table_Bonferroni_Selected_Interactions.csv')