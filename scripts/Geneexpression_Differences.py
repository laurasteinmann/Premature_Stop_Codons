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

#Analysis
gt_overlap = pd.read_csv('../data/preprocessed/GT_Section_Numeric_Overlap.csv', index_col=0)
fixed_section = pd.read_csv('../data/preprocessed/Fixed_Section_Stop_Codons_full_vcf.csv', index_col=0) 
expression_overlap = pd.read_csv('../data/preprocessed/Expression_Overlap.csv', index_col=0)

stop_table = gene_info(fixed_section)
print('Full List of Premature Stop Codons Ara11 contain: ', stop_table.shape[0])
stop_table.to_csv('../data/processed/Genexpression_Differences/Stop_Table_Full.csv')

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
stop_list = p_values(stop_table, expression_overlap)
stop_list.to_csv('../data/processed/Genexpression_Differences/Prem_Stops_WT_ko_fulllist.csv')