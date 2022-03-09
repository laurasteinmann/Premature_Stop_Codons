cooccurrence_table <- read.csv('../data/processed/Coocurrence/Coocurrence_Table_Significant_Full.csv', row.names=1)

pvaluethreshold <- 0.05
interactions <- c(1 : dim(cooccurrence_table)[1])
for (i in interactions){
    table_attributes <- cooccurrence_table[i,]
    total_accessions <- 665
    gene2_count <- as.numeric(table_attributes['Gene2_count'])
    gene1_count <- as.numeric(table_attributes['Gene1_count'])
    coocurrence <- as.numeric(table_attributes["Cooccurrence"])
    pvalue_under <- dhyper(coocurrence,gene2_count,total_accessions-gene2_count,gene1_count)
    pvalue_over <-  1 - dhyper(coocurrence,gene2_count,total_accessions-gene2_count,gene1_count)
    cooccurrence_table[i, 9] <- pvalue_under
    cooccurrence_table[i, 8] <- pvalue_over
    lower_border <- qhyper(pvaluethreshold, gene2_count, total_accessions-gene2_count, gene1_count)
    upper_border <- qhyper(pvaluethreshold, gene2_count, total_accessions-gene2_count, gene1_count, lower.tail = FALSE)
    cooccurrence_table[i, 6] <- lower_border
    cooccurrence_table[i, 7] <- upper_border
}

write.csv(cooccurrence_table, '../data/processed/Coocurrence/Coocurrence_Table_Significant_Full.csv', row.names=TRUE)

cooccurrence_table <- read.csv('../data/processed/Coocurrence/Coocurrence_Table_Bonferroni_Significant_Full.csv', row.names=1)
pvaluethreshold <- 0.05/dim(cooccurrence_table)[1]
interactions <- c(1 : dim(cooccurrence_table)[1])
for (i in interactions){
    table_attributes <- cooccurrence_table[i,]
    total_accessions <- 665
    gene2_count <- as.numeric(table_attributes['Gene2_count'])
    gene1_count <- as.numeric(table_attributes['Gene1_count'])
    coocurrence <- as.numeric(table_attributes["Cooccurrence"])
    pvalue_under <- dhyper(coocurrence,gene2_count,total_accessions-gene2_count,gene1_count)
    pvalue_over <-  1 - dhyper(coocurrence,gene2_count,total_accessions-gene2_count,gene1_count)
    cooccurrence_table[i, 9] <- pvalue_under
    cooccurrence_table[i, 8] <- pvalue_over
    lower_border <- qhyper(pvaluethreshold, gene2_count, total_accessions-gene2_count, gene1_count)
    upper_border <- qhyper(pvaluethreshold, gene2_count, total_accessions-gene2_count, gene1_count, lower.tail = FALSE)
    cooccurrence_table[i, 6] <- lower_border
    cooccurrence_table[i, 7] <- upper_border
}
write.csv(cooccurrence_table, '../data/processed/Coocurrence/Coocurrence_Table_Bonferroni_Significant_Full.csv', row.names=TRUE)