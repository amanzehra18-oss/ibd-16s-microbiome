# 02 - reverse-complement ASVs of the nine reverse-oriented IBD samples
# Input : ASV_table.csv  (rows = ASV sequence, columns = samples)
# Output: ASV_table_rc_fixed.tsv, ASV_sequences_rc_fixed.fasta, ASV_id_map.csv
library(Biostrings)
asv <- read.csv("ASV_table.csv", row.names = 1, check.names = FALSE)
seqs <- rownames(asv)
# ASVs starting with GGACTAC are reverse-complemented (806R-side reads).
rc_flag <- startsWith(seqs, "GGACTAC")
new_seq <- seqs
new_seq[rc_flag] <- as.character(reverseComplement(DNAStringSet(seqs[rc_flag])))
ids <- sprintf("ASV_%04d", seq_along(seqs))
map <- data.frame(ASV_ID = ids, original = seqs, rc_applied = rc_flag, sequence = new_seq)
write.csv(map, "ASV_id_map.csv", row.names = FALSE)
# merge ASVs that became identical after reverse-complementing
tab <- rowsum(as.matrix(asv), group = new_seq)
final_ids <- sprintf("ASV_%04d", seq_len(nrow(tab)))
out <- data.frame(`#OTU ID` = final_ids, tab, check.names = FALSE)
write.table(out, "ASV_table_rc_fixed.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeXStringSet(DNAStringSet(setNames(rownames(tab), final_ids)), "ASV_sequences_rc_fixed.fasta")
# Then run "PICRUSt2 full pipeline" on Galaxy Europe with these two files
# (the FASTA and the table converted to BIOM) and download the pathway abundances.
# NOTE: the numbers in the manuscript used the IDs/merging done in the original run;
# the ID map in the supplementary data (if provided) is authoritative.
