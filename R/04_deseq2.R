# 04 - Differential abundance of genera, IBD vs healthy (DESeq2)
# Input : genus_counts.csv (ASV table summed per genus with the taxonomy in Table S3; unassigned ASVs -> "Unclassified")
# Filter: "Unclassified" removed; genera with >= 10 reads in total (186 genera tested).
# Output: DESeq2_All_Genera.csv (S11) and DESeq2_Significant_Genera.csv (S10).
# Result (run in R with DESeq2, sfType='poscounts'): 186 tested, 44 significant (12 higher / 32 lower in IBD).
# Most significant genera are absent from one group, so their fold changes are very large (presence/absence).
source("R/00_setup.R")
library(DESeq2)
counts <- read.csv("genus_counts.csv", row.names = 1, check.names = FALSE)
counts <- counts[rownames(counts) != "Unclassified", ]
counts <- counts[rowSums(counts) >= 10, ]
coldata <- data.frame(row.names = colnames(counts),
                      Group = factor(group_from_sample(colnames(counts)), levels = c("Healthy", "IBD")))
dds <- DESeqDataSetFromMatrix(round(counts), coldata, ~ Group)
dds <- DESeq(dds, sfType = "poscounts")
res <- as.data.frame(results(dds, contrast = c("Group", "IBD", "Healthy")))
res$Genus <- rownames(res)
write.csv(res, "DESeq2_All_Genera.csv", row.names = FALSE)
sig <- subset(res, !is.na(padj) & padj < 0.05)
write.csv(sig, "DESeq2_Significant_Genera.csv", row.names = FALSE)
cat("tested:", nrow(res), " significant:", nrow(sig), " higher in IBD:", sum(sig$log2FoldChange > 0),
    " lower:", sum(sig$log2FoldChange < 0), "\n")
