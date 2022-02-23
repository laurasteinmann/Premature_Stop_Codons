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