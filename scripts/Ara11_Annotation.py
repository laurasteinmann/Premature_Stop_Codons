# Libraries
import pandas as pd 
import copy as cp
import numpy as np
# Functions
def fill_gene_names(gff, length):
    genes = np.empty((gff.shape[0]), dtype='<U10')
    for entry in range(0, gff.shape[0]):
        string = gff.iloc[entry, 8]
        split = string.partition('ID')[2]
        split_start = split.find('AT')
        split_stop = split_start + length
        gene_name = split[split_start:split_stop]
        genes[entry] = gene_name

    gff = cp.copy(gff)
    gff['Genes'] = genes
    return gff 

def filter_genes(gff, filter):
    mask = []
    for entry in range(0, gff.shape[0]):
        string = gff.iloc[entry, 9]
        if string[3] == filter:
            mask.append(True)
        else:
            mask.append(False)
    gff_without_new_rnas = gff[mask]
    return gff_without_new_rnas

def find_id_in_description(table, row):
    string = table.loc[row, 8]
    split = string.partition("ID")[2]
    split_start = split.find('AT')
    split_stop = split.find('-')
    gene = split[split_start:split_stop]
    return gene 

def find_parent_in_description(table, row):
    string = table.loc[row, 8]
    split = string.partition("Parent")[2]
    split_start = split.find('AT')
    split_stop = split.find(';')
    number_names = split_stop/12
    names = []
    for i in range(0, int(number_names)):
        factor = i + 1
        gene_stop = 12 * factor
        gene_start = gene_stop - 11
        gene = split[gene_start:gene_stop]
        names.append(gene)
    return names

def protein_names_new_column(table):
    proteins_names = []
    for protein in table.index:
        name = find_id_in_description(table, protein)
        proteins_names.append(name)
    table = cp.copy(table)
    table['ID'] = proteins_names
    return table

def utr_new_arrangement(table):
    utr5 = pd.DataFrame(columns=['Name', "Start", "Stop"])
    for utr in table.index:
        name = find_parent_in_description(table, utr)
        start = table.loc[utr, 3]
        stop = table.loc[utr, 4]
        for entry in name:
            new_row = {"Name": [entry], "Start": [start], "Stop": [stop]}
            new = pd.DataFrame.from_dict(new_row)
            utr5 = pd.concat([utr5, new])
    utr5 = utr5.drop_duplicates(subset=['Name'])
    return utr5

def generate_ara11_protein_table(gff, results_table):
    gene_list = np.unique(gff.Genes)
    strange_genes = []
    for gene in gene_list:
        all_classes = gff[gff['Genes'] == gene]
        proteins = all_classes[all_classes[2] == 'protein']
        if proteins.shape[0] == 1:
            five_prime_utr = all_classes[all_classes[2] == 'five_prime_UTR']
            if five_prime_utr.shape[0] == 0:
                five_utr_start = np.nan
                five_utr_stop = np.nan
            else:
                five_utr_start = int(five_prime_utr.iloc[0, 3])
                five_utr_stop = int(five_prime_utr.iloc[0, 4])
            three_prime_utr = all_classes[all_classes[2] == 'three_prime_UTR']
            if three_prime_utr.shape[0] == 0:
                three_utr_start = np.nan
                three_utr_stop = np.nan
            else:
                three_utr_start = int(three_prime_utr.iloc[0, 3])
                three_utr_stop = int(three_prime_utr.iloc[0, 4])
            gene = find_id_in_description(proteins, proteins.index[0])
            orientation = proteins.iloc[0, 6]
            start_protein = int(proteins[3])
            stop_protein = int(proteins[4])
            new_row = {"Name": [gene], 'Orientation': [orientation], "Start_5UTR": [five_utr_start],
                       "Stop_5UTR": [five_utr_stop], "Start_Protein": [start_protein], 'Stop_Protein': [stop_protein],
                       'Start_3UTR': [three_utr_start], 'Stop_3UTR': [three_utr_stop]}
            new = pd.DataFrame.from_dict(new_row)
            results_table = pd.concat([results_table, new], ignore_index=True)
        elif proteins.shape[0] > 1:
            proteins = protein_names_new_column(proteins)
            five_prime_utr = all_classes[all_classes[2] == 'five_prime_UTR']
            new_five_utrs = utr_new_arrangement(five_prime_utr)
            three_prime_utr = all_classes[all_classes[2] == 'three_prime_UTR']
            new_three_utrs = utr_new_arrangement(three_prime_utr)
            proteins = proteins.sort_values(by='ID')
            for entry in proteins.index:
                gene = proteins.loc[entry, 'ID']
                orientation = proteins.loc[entry, 6]
                five_utr = new_five_utrs[new_five_utrs['Name'] == gene]
                if five_utr.shape[0] == 0:
                    five_start = np.nan
                    five_stop = np.nan
                else:
                    five_start = int(five_utr.Start)
                    five_stop = int(five_utr.Stop)
                start_protein = int(proteins.loc[entry, 3])
                stop_protein = int(proteins.loc[entry, 4])
                three_utr = new_three_utrs[new_three_utrs['Name'] == gene]
                if three_utr.shape[0] == 0:
                    three_start = np.nan
                    three_stop = np.nan
                else:
                    three_start = int(three_utr.Start)
                    three_stop = int(three_utr.Stop)
                new_row = {"Name": [gene], 'Orientation': [orientation], "Start_5UTR": [five_start],
                           "Stop_5UTR": [five_stop], "Start_Protein": [start_protein], 'Stop_Protein': [stop_protein],
                           'Start_3UTR': [three_start], 'Stop_3UTR': [three_stop]}
                new = pd.DataFrame.from_dict(new_row)
                results_table = pd.concat([results_table, new], ignore_index=True)
        else:
            strange_genes.append(gene)
    return results_table, strange_genes

def generate_genes_without_proteins_table(gff, genes_list):
    genes_without_proteins = pd.DataFrame(columns=["Name", 'Orientation', "Start_Gene", 'Stop_Gene'])
    strange_attributes_2 = []
    for gene in genes_list:
        all_classes = gff[gff['Genes'] == gene]
        gene_attributes = all_classes[all_classes[2] == 'gene']
        gene_attributes_2 = all_classes[all_classes[2] == 'pseudogene']
        gene_attributes_3 = all_classes[all_classes[2] == 'transposable_element_gene']
        if gene_attributes.shape[0] == 1:
            orientation = gene_attributes.iloc[0, 6]
            start_gene = int(gene_attributes.iloc[0, 3])
            stop_gene = int(gene_attributes.iloc[0, 4])
            new_row = {"Name": [gene], 'Orientation': [orientation], "Start_Gene": [start_gene], 'Stop_Gene': [stop_gene]}
            new = pd.DataFrame.from_dict(new_row)
            genes_without_proteins = pd.concat([genes_without_proteins, new], ignore_index=True)
        elif gene_attributes_2.shape[0] == 1:
            orientation = gene_attributes_2.iloc[0, 6]
            start_gene = int(gene_attributes_2.iloc[0, 3])
            stop_gene = int(gene_attributes_2.iloc[0, 4])
            new_row = {"Name": [gene], 'Orientation': [orientation], "Start_Gene": [start_gene], 'Stop_Gene': [stop_gene]}
            new = pd.DataFrame.from_dict(new_row)
            genes_without_proteins = pd.concat([genes_without_proteins, new], ignore_index=True)
        elif gene_attributes_3.shape[0] == 1:
            orientation = gene_attributes_3.iloc[0, 6]
            start_gene = int(gene_attributes_3.iloc[0, 3])
            stop_gene = int(gene_attributes_3.iloc[0, 4])
            new_row = {"Name": [gene], 'Orientation': [orientation], "Start_Gene": [start_gene], 'Stop_Gene': [stop_gene]}
            new = pd.DataFrame.from_dict(new_row)
            genes_without_proteins = pd.concat([genes_without_proteins, new], ignore_index=True)
        else:
            strange_attributes_2.append(gene)
    return genes_without_proteins, strange_attributes_2

def generate_transposons_table(transposons, gff_transposons):
    gene_list = np.unique(gff_transposons.Genes)
    strange_transposons = []
    for gene in gene_list:
        all_classes = gff_transposons[gff_transposons['Genes'] == gene]
        gene_attributes = all_classes[all_classes[2] == 'transposable_element']
        if gene_attributes.shape[0] == 1:
            orientation = gene_attributes.iloc[0, 6]
            start_gene = int(gene_attributes.iloc[0, 3])
            stop_gene = int(gene_attributes.iloc[0, 4])
            new_row = {"Name": [gene], 'Orientation': [orientation], "Start_Gene": [start_gene], 'Stop_Gene': [stop_gene]}
            new = pd.DataFrame.from_dict(new_row)
            transposons = pd.concat([transposons, new], ignore_index=True)
        else:
            strange_transposons.append(gene)

    return transposons, strange_transposons

def generate_complete_ara11_transposons_start_stop_table(gff):
    transposons = pd.DataFrame(columns=["Name", 'Orientation', "Start_Gene", 'Stop_Gene'])
    gff_filtered = fill_gene_names(gff, 10)
    gff_transposons = filter_genes(gff_filtered, 'T')
    transposons, strange_transposons = generate_transposons_table(transposons, gff_transposons)
    return transposons, strange_transposons
# Analysis
gff = pd.read_csv('../data/raw/Araport11_GFF3_genes_transposons.201606.gff', sep='\t', header=None)
ara11 = pd.DataFrame(columns=['Name', 'Orientation', 'Start_5UTR', 'Stop_5UTR', 'Start_Protein', 'Stop_Protein', 'Start_3UTR', 'Stop_3UTR'])
gff_without_nas = gff[~gff.iloc[:,8].isnull()]
gff_preprocessed = fill_gene_names(gff_without_nas, 9)
gff_genes = filter_genes(gff_preprocessed, 'G')
ara11_genes, strange_genes = generate_ara11_protein_table(gff_genes, ara11)
ara11_genes.to_csv('../data/processed/Ara11/Ara11_genes_gff.csv')
table_genes_without_proteins, strange_attributes_2 = generate_genes_without_proteins_table(gff_genes, strange_genes)
table_genes_without_proteins.to_csv('../data/processed/Ara11/Ara11_genes_without_proteins_gff.csv')
table_transposons, strange_transposons = generate_complete_ara11_transposons_start_stop_table(gff_without_nas)
table_transposons.to_csv('../data/processed/Ara11/Ara11_transposons_gff.csv')
