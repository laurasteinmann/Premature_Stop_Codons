# Functional Dataset based on Geneexpression Differences
#Libraries
import pandas as pd
import numpy as np
import copy as cp
from scipy import stats
#Functions

def gene_info(fixed_section):
    stop_attributes = pd.DataFrame(columns=["Chrom", "Pos", "Gene", 'WT_num', 'WT_mean', 'WT_std', "WT_acc", 'KO_num',
                                            "KO_mean", "KO_std", 'KO_acc', "na_Num","na_Acc", 'pvalue'])
    stops = fixed_section.shape[0]
    for stop in range(stops):
        string = fixed_section.iloc[stop, 7]
        split = string.partition("stop_gained")[2]
        split_start = split.find('AT')
        split_stop = split_start + 9
        gene = split[split_start:split_stop]
        chrom = fixed_section.iloc[stop, 0]
        pos = fixed_section.iloc[stop, 1]
        new_row = {'Chrom': [chrom], 'Pos': [pos], 'Gene': [gene], 'WT_num': [np.nan], 'WT_mean': [np.nan], 'WT_std': [np.nan],
                   'WT_acc': [np.nan], 'KO_num': [np.nan], "KO_mean": np.nan, "KO_std": [np.nan], "KO_acc": [np.nan],
                   "na_Num": [np.nan], "na_Acc": [np.nan], "pvalue": [np.nan]}
        new = pd.DataFrame.from_dict(new_row)
        stop_attributes = pd.concat([stop_attributes, new], ignore_index=True)
    return stop_attributes

def filter_Common_Genes(stop_list, expression, gt):
    filter = np.isin(stop_list.Gene, expression.columns.values)
    stop_list = stop_list.loc[filter]
    gt_t = gt.transpose()
    gt_new = gt_t[filter]
    gt_new = gt_new.transpose()
    return stop_list, gt_new

def fill_num(stop_list, gt, type, type_number):
    list_num = []
    for stop in range(stop_list.shape[0]):
        type_df = gt[gt.iloc[:, stop] == type_number]
        type_num = type_df.shape[0]
        list_num.append(type_num)
    stop_list = cp.copy(stop_list)
    stop_list[type] = list_num
    return stop_list

def filter_num(stop_list, gt, type):
    subset = stop_list[stop_list[type] <= 1]
    new_stop_list = stop_list.drop(subset.index)
    gt_t = gt.transpose()
    gt_t.index = gt_t.index.astype('int64')
    gt_new = gt_t.drop(subset.index)
    gt_new = gt_new.transpose()
    return new_stop_list, gt_new

def fill_stat(stop_list, expression, gt, type_mean, type_std, type_num):
    list_mean = []
    list_std = []
    for stop in range(stop_list.shape[0]):
        type_df = gt[gt.iloc[:, stop] == type_num]
        stop_att = stop_list.iloc[stop, ]
        gene = stop_att.Gene
        expression_sub = expression.loc[type_df.index]
        expression_gene = expression_sub[gene]
        mean = np.mean(expression_gene)
        std = np.std(expression_gene)
        list_mean.append(mean)
        list_std.append(std)

    stop_list[type_mean] = list_mean
    stop_list[type_std] = list_std
    return stop_list

def fill_acc(stop_list, gt, type, type_num):
    list_acc = []
    for stop in range(stop_list.shape[0]):
        type_df = gt[gt.iloc[:, stop] == type_num]
        type_acc = type_df.index.values
        list_acc.append(type_acc)
    stop_list[type] = list_acc
    return stop_list

def na_att(stop_list, gt):
    na_num_list = []
    na_acc_list = []
    for stop in range(stop_list.shape[0]):
        na = gt[gt.iloc[:, stop].isnull()]
        na_num = na.shape[0]
        na_num_list.append(na_num)
        na_acc_list.append(na.index.values)

    stop_list["na_Num"] = na_num_list
    stop_list["na_Acc"] = na_acc_list
    return stop_list

def p_values(stop_list, expression):
    list_pvalues = []
    for stop in range(stop_list.shape[0]):
        stop_att = stop_list.iloc[stop, ]
        wildtype_index = stop_att.WT_acc
        wildtype = expression.loc[wildtype_index, stop_att.Gene]
        ko_index = stop_att.KO_acc
        ko = expression.loc[ko_index, stop_att.Gene]
        trash, pval = stats.ttest_ind(wildtype, ko)
        list_pvalues.append(pval)

    stop_list['pvalue'] = list_pvalues
    return stop_list

def functional_datasets(stop_list, pvalue):
    stop_list_1 = stop_list[stop_list['pvalue'] <= pvalue]
    return stop_list_1

def filter_stop_significant_higher_expression(stop_list):
    stop_list_all_significant = functional_datasets(stop_list, 0.05)
    stop_list_significant_high_expression = stop_list_all_significant[
        stop_list_all_significant['WT_mean'] < stop_list_all_significant['KO_mean']]
    print('Significant premature stop codons with increase of gene expression ',
          stop_list_significant_high_expression.shape[0])

    return stop_list_significant_high_expression, stop_list_all_significant

def filter_stop_significant_lower_expression(stop_list_all_significant, stop_list_significant_high_expression):
    stop_list_significant_low_expression = stop_list_all_significant.drop(stop_list_significant_high_expression.index)
    print('Significant premature stop codons with decrease of gene expression ',
          stop_list_significant_low_expression.shape[0])

    return stop_list_significant_low_expression

def filter_gt_matrix(gt_filtered, stop_list_significant_low_expression):
    gt_short = gt_filtered.transpose()
    gt_short = gt_short.loc[stop_list_significant_low_expression.index]
    gt_short = gt_short.transpose()

    return gt_short

def filter_stop_bonferroni_higher_expression(stop_list):
    pvalue = 0.05 / stop_list.shape[0]
    stop_list_bonferroni_significant = functional_datasets(stop_list, pvalue)
    stop_list_bonferroni_high_expression = stop_list_bonferroni_significant[
        stop_list_bonferroni_significant['WT_mean'] < stop_list_bonferroni_significant['KO_mean']]
    print('Bonferroni corrected significant premature stop codons with increase of gene expression ',
          stop_list_bonferroni_high_expression.shape[0])

    return stop_list_bonferroni_high_expression, stop_list_bonferroni_significant

def filter_stop_bonferroni_low_expression(stop_list_bonferroni_significant, stop_list_bonferroni_high_expression):
    stop_list_bonferroni_low_expression = stop_list_bonferroni_significant.drop(
        stop_list_bonferroni_high_expression.index)
    print('Bonferroni corrected significant premature stop codons with decreased of gene expression ',
          stop_list_bonferroni_low_expression.shape[0])

    return stop_list_bonferroni_low_expression

def filter_gt_matrix_bonferroni(gt_filtered, stop_list_bonferroni_low_expression):
    gt_short = gt_filtered.transpose()
    gt_short = gt_short.loc[stop_list_bonferroni_low_expression.index]
    gt_short = gt_short.transpose()

    return gt_short

def gene_based_list(single_stop_list, gt):
    gt = gt.transpose()
    gt.index = gt.index.astype('int64')
    stop_attributes = pd.DataFrame(columns=['Gene', 'Prem_Stop_Codons_Num', 'WT_num', 'WT_acc', 'KO_num', 'KO_acc',
                                            'na_Num', 'na_Acc'])
    unique_genes = np.unique(single_stop_list.Gene)
    stop_attributes['Gene'] = unique_genes
    list_num_premature_stops = []
    for gene in unique_genes:
        stops_in_gene = single_stop_list[single_stop_list['Gene'] == gene]
        list_num_premature_stops.append(stops_in_gene.shape[0])

    stop_attributes['Prem_Stop_Codons_Num'] = list_num_premature_stops

    list_ko_num = []
    list_ko_accessions = []
    list_wt_num = []
    list_wt_accessions = []
    list_na_num = []
    list_na_accessions  = []
    for gene in unique_genes:
        stops_in_gene = single_stop_list[single_stop_list['Gene'] == gene]
        gt_stops = gt.loc[stops_in_gene.index]
        accessions_ko = np.empty(shape=[0, gt_stops.shape[0]])
        accessions_wt = np.empty(shape=[0, gt_stops.shape[0]])
        accessions_na = np.empty(shape=[0, gt_stops.shape[0]])
        for row in range(0, gt_stops.shape[0]):
            stop_codon = gt_stops.iloc[row,:]
            ko_accessions = stop_codon[stop_codon==1]
            for name in ko_accessions.index:
                accessions_ko = np.append(accessions_ko, name)
            wt_accessions = stop_codon[stop_codon == 0]
            for name in wt_accessions.index:
                accessions_wt = np.append(accessions_wt, name)
            na_accessions = stop_codon[stop_codon.isnull()]
            for name in na_accessions.index:
                accessions_na = np.append(accessions_na, name)

        unique_accessions_ko = np.unique(accessions_ko)
        unique_accessions_wt = np.unique(accessions_wt)
        unique_accessions_na = np.unique(accessions_na)
        overlap_corrected_wt = unique_accessions_wt[~np.isin(unique_accessions_wt, unique_accessions_ko)]
        overlap_corrected_ko_na = unique_accessions_na[~np.isin(unique_accessions_na, unique_accessions_ko)]
        overlap_corrected_na = overlap_corrected_ko_na[~np.isin(overlap_corrected_ko_na, unique_accessions_wt)]
        list_ko_num.append(len(unique_accessions_ko))
        list_ko_accessions.append(unique_accessions_ko)
        list_wt_num.append(len(overlap_corrected_wt))
        list_wt_accessions.append(overlap_corrected_wt)
        list_na_num.append(len(overlap_corrected_na))
        list_na_accessions.append(overlap_corrected_na)

    stop_attributes['KO_num'] = list_ko_num
    stop_attributes['KO_acc'] = list_ko_accessions
    stop_attributes['WT_num'] = list_wt_num
    stop_attributes['WT_acc'] = list_wt_accessions
    stop_attributes['na_Num'] = list_na_num
    stop_attributes['na_Acc'] = list_na_accessions
    print('After merging the single stop premature stop codons on a gene based level: ',
          stop_attributes.shape[0])

    return stop_attributes
#Analysis
gt_overlap = pd.read_csv('../data/preprocessed/GT_Section_Numeric_Overlap.csv', index_col=0)
fixed_section = pd.read_csv('../data/preprocessed/Fixed_Section_Stop_Codons_full_vcf.csv', index_col=0) 
expression_overlap = pd.read_csv('../data/preprocessed/Expression_Overlap.csv', index_col=0)

stop_table = gene_info(fixed_section)
print('Full List of Premature Stop Codons Ara11 contain: ', stop_table.shape[0])
stop_table.to_csv('../data/preprocessed/Stop_Table_Full.csv')

stop_list_filtered, gt_filtered = filter_Common_Genes(stop_table, expression_overlap, gt_overlap)
print('Filtered List of Premature Stop Codons after removing unknown genes of Ara11: ', stop_list_filtered.shape[0])
stop_table = fill_num(stop_list_filtered, gt_filtered, "WT_num", 0)
stop_table, gt_filtered = filter_num(stop_table, gt_filtered, 'WT_num')
stop_table = fill_stat(stop_table, expression_overlap, gt_filtered, 'WT_mean', 'WT_std', 0)
stop_table = fill_acc(stop_table, gt_filtered, 'WT_acc', 0)
stop_table = fill_num(stop_table, gt_filtered, "KO_num", 1)
stop_table, gt_filtered = filter_num(stop_table, gt_filtered, 'KO_num')
stop_table = fill_stat(stop_table, expression_overlap, gt_filtered, 'KO_mean', 'KO_std', 1)
stop_table = fill_acc(stop_table, gt_filtered, 'KO_acc', 1) 
stop_table = na_att(stop_table, gt_filtered)
print('The final list of Premature Stop Codons before calculating the pvalues includes: ', stop_table.shape[0])
stop_table = p_values(stop_table, expression_overlap)
stop_table.to_csv('../data/processed/Genexpression_Differences/Prem_Stops_WT_ko_fulllist.csv')

stop_list_significant_high_expression, stop_list_all_significant = filter_stop_significant_higher_expression(stop_table)
stop_list_significant_high_expression.to_csv('../data/processed/Genexpression_Differences/Type3_Stop_List_Significant.csv')
stop_list_significant_low_expression = filter_stop_significant_lower_expression(stop_list_all_significant, stop_list_significant_high_expression)
stop_list_significant_low_expression.to_csv('../data/processed/Genexpression_Differences/Type2_Stop_List_Significant.csv')
gt_short = filter_gt_matrix(gt_filtered, stop_list_significant_low_expression)
gt_short.to_csv('../data/preprocessed/GT_Significant_StopList.csv')

stop_list_bonferroni_high_expression, stop_list_bonferroni_significant = filter_stop_bonferroni_higher_expression(stop_table)
stop_list_bonferroni_high_expression.to_csv('../data/processed/Genexpression_Differences/Type3_Stop_List_Significant_Bonferroni.csv')
stop_list_bonferroni_low_expression = filter_stop_bonferroni_low_expression(stop_list_bonferroni_significant, stop_list_bonferroni_high_expression)
stop_list_bonferroni_low_expression.to_csv('../data/processed/Genexpression_Differences/Type2_Stop_List_Significant_Bonferroni.csv')
gt_short = filter_gt_matrix_bonferroni(gt_filtered, stop_list_bonferroni_low_expression)
gt_short.to_csv('../data/preprocessed/GT_Bonferroni_Significant_StopList.csv')

# Combining significant premature stop codons to gene level list
gene_premature_stop_list = gene_based_list(stop_list_bonferroni_low_expression, gt_short)
gene_premature_stop_list.to_csv('../data/processed/Genexpression_Differences/Gene_Based_Bonferroni_Significant_Lower_Expression.csv')
