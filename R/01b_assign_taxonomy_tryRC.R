# 01b - Genus taxonomy for the ASVs (SILVA 138.2), trying both strands
# Why tryRC = TRUE: the ASVs of nine IBD samples (2-5, 7-11) are reverse-complemented
# relative to the SILVA reference. Without tryRC they receive no genus (0% of ASVs).
# Input : ASV_sequences_fixed.fasta (2,266 ASVs, 120 bp)
# Output: taxonomy_tryRC.csv (= Supplementary Table S3, 61% of ASVs with a genus)
library(dada2)
seqs <- getSequences("ASV_sequences_fixed.fasta")
tax  <- assignTaxonomy(seqs, "silva_nr99_v138.2_toGenus_trainset.fa.gz",
                       multithread = TRUE, tryRC = TRUE)
write.csv(data.frame(ASV = rownames(tax), tax), "taxonomy_tryRC.csv", row.names = FALSE)
cat("ASVs with genus (%):", mean(!is.na(tax[, "Genus"])) * 100, "\n")
# genus_counts.csv = ASV table summed by Genus (ASVs without a genus -> "Unclassified"), samples in columns.
