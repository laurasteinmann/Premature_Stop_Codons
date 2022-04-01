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

coocurrence_table_less = pd.read_csv('../data/processed/Coocurrence/Coocurrence_Table_Significant_Full.csv', index_col=0)
pvaluethreshold2 = 0.05
coocurrence_significant_less = coocurrence_table_less[coocurrence_table_less['pvalue_upper']<= pvaluethreshold2]
coocurrence_significant_less = coocurrence_significant_less[coocurrence_significant_less['Gene1']!=coocurrence_significant_less['Gene2']]
coocurrence_significant_less = coocurrence_significant_less.replace('diff_Chr', 1000000)
coocurrence_significant_less["Distance"] = pd.to_numeric(coocurrence_significant_less["Distance"])

coocurrence_significant_under = coocurrence_table_less[coocurrence_table_less['pvalue_lower']<= pvaluethreshold2]
coocurrence_significant_under = coocurrence_significant_under[coocurrence_significant_under['Gene1']!=coocurrence_significant_under['Gene2']]
coocurrence_significant_under = coocurrence_significant_under.replace('diff_Chr', 1000000)
coocurrence_significant_under["Distance"] = pd.to_numeric(coocurrence_significant_under["Distance"])
