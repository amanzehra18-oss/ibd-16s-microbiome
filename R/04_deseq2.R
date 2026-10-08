# 04 - Differential abundance of genera, IBD vs healthy (DESeq2)  [RECONSTRUCTED]
# Input : genus count table (rows = genera, columns = samples) built from the ASV
#         table and DADA2 taxonomy by summing ASV counts per genus.
# Output: DESeq2_All_Genera.csv (S11) and DESeq2_Significant_Genera.csv (S10)
# Reference results: 175 genera tested, 49 with padj < 0.05 (45 lower, 4 higher in IBD).
library(DESeq2)
counts <- read.csv("genus_counts.csv", row.names = 1, check.names = FALSE)
coldata <- data.frame(row.names = colnames(counts),
                      Group = factor(group_from_sample(colnames(counts)), levels = c("Healthy", "IBD")))
dds <- DESeqDataSetFromMatrix(round(counts), coldata, ~ Group)
dds <- DESeq(dds, sfType = "poscounts")   # [CHECK] original size-factor setting not recorded
res <- as.data.frame(results(dds, contrast = c("Group", "IBD", "Healthy")))
res$Genus <- rownames(res)
write.csv(res, "DESeq2_All_Genera.csv", row.names = FALSE)
sig <- subset(res, !is.na(padj) & padj < 0.05)
write.csv(sig, "DESeq2_Significant_Genera.csv", row.names = FALSE)
cat("tested:", nrow(res), " significant:", nrow(sig), "\n")
