# Libraries
import pandas as pd 

# Functions


# Analysis
gff = pd.read_csv('data/raw/Araport11_GFF3_genes_transposons.201606.gff', sep='\t', header=None)

ara11 = pd.DataFrame(columns=['Name', 'Orientation', 'Start_5UTR', 'Stop_5UTR', 'Start_Protein', 'Stop_Protein', 'Start_3UTR', 'Stop_3UTR'])


