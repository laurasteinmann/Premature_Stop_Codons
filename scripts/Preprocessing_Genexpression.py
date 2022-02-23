# Preprocessing of Genexpression Differences Approach
# Libraries
import pandas as pd 
import numpy as np
#Function
def simple_gt(expression_data, gt_section):
    expression_accessions = expression_data.index.values
    genomic_accessions = gt_section.columns.values
    expression_accessions = expression_accessions.astype(int)
    genomic_accessions = genomic_accessions.astype(int)
    gt_section.columns = genomic_accessions
    in_both = genomic_accessions[np.isin(genomic_accessions,
                                         expression_accessions)]

    expression = expression_data.loc[in_both, ]
    genomics = gt_section[in_both]

    gt_simple = pd.DataFrame(columns=genomics.columns)

    for column in range(genomics.shape[1]):
        one_column = genomics.iloc[:, column]
        changed_strings = []
        for string in one_column:
            changed = string[:1]
            changed_strings = np.append(changed_strings, changed)
        gt_simple.iloc[:, column] = changed_strings

    gt_numeric = pd.DataFrame(columns=gt_simple.columns)

    for column in range(gt_simple.shape[1]):
        one_column = gt_simple.iloc[:, column]
        one_column = one_column.replace(".", np.nan)
        one_column = one_column.astype("float")
        gt_numeric.iloc[:, column] = one_column

    gt = gt_numeric.transpose()
    return gt, expression
#Analysis
gt_section = pd.read_csv('data/preprocessed/GT_Section_Stop_Codons_full_vcf.csv', index_col=0)
expression_data = pd.read_csv('data/preprocessed/727_Expressiondata.csv', index_col=0)
gt, expression = simple_gt(expression_data, gt_section)
gt.to_csv('data/preprocessed/GT_Section_Numeric_Overlap.csv')
expression.to_csv('data/preprocessed/Expression_Overlap.csv')