# Expressiondata
load('../data/raw/727_Expressiondata.rda')
write.csv(RNA, '../data/preprocessed/727_Expressiondata.csv')
# vcf Stop Codons 
library('vcfR', verbose=FALSE)
stop_codons_full_vcf <- read.vcfR('data/raw/Stopp_Codons_full_Ara11.vcf', verbose = FALSE)
extracted_fixed_section <- stop_codons_full_vcf@fix
write.csv(extracted_fixed_section, 'data/preprocessed/Fixed_Section_Stop_Codons_full_vcf.csv')
extracted_gt_section <- stop_codons_full_vcf@gt
write.csv(extracted_gt_section, 'data/preprocessed/GT_Section_Stop_Codons_full_vcf.csv')

# vcf Full Stop Codons 1035 Genomes
library('vcfR', verbose=FALSE)
genomes <- read.vcfR('../data/raw/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff')
extracted_fixed_section <- genomes@fix
write.csv(extracted_fixed_section, 'data/preprocessed/Fixed_Section_1001_Genomes.csv')
extracted_gt_section <- genomes@gt
write.csv(extracted_gt_section, '../data/preprocessed/GT_Section_1001_Genomes.csv') 
# SNPs R-Data
load('../data/raw/snps_2029.rda')
write.csv(SNPs, '../data/preprocessed/snps_2029.csv')

