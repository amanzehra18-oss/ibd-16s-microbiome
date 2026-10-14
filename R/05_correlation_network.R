# 05 - Genus x pathway Spearman correlations (Table S12) and network tables
# Reproduces Table S12 (v8: 84 pairs; 47 LPSSYN-PWY, 37 PEPTIDOGLYCANSYN-PWY). The same settings reproduced the previous S12 (95 pairs) exactly.
# Settings that reproduce S12:
#  * ALL genera present in >= 2 samples (not only DESeq2-significant ones)
#  * relative abundance = genus reads / total reads of the sample, INCLUDING Unclassified
#  * pathway abundance = per-million table S8b
#  * Spearman; Benjamini-Hochberg FDR applied WITHIN each pathway; FDR < 0.05
#  * the "Unclassified" row is dropped from the final table
source("R/00_setup.R")
genus <- read.csv("genus_counts.csv", row.names = 1, check.names = FALSE)
genus_rel <- sweep(genus, 2, colSums(genus), "/")
genus_rel <- genus_rel[rowSums(genus_rel > 0) >= 2, ]
path <- read.csv(supp("Supplementary_Table_S8b_PICRUSt2_Pathway_Abundance_per_million.csv"),
                 row.names = 1, check.names = FALSE)
# match samples by their number (1..22), never by column position
snum <- function(x) as.integer(sub("^X?([0-9]+).*$", "\\1", x))
common <- intersect(snum(colnames(genus_rel)), snum(colnames(path)))
genus_rel <- genus_rel[, match(common, snum(colnames(genus_rel)))]
path      <- path[,      match(common, snum(colnames(path)))]
paths <- c("LPSSYN-PWY", "PEPTIDOGLYCANSYN-PWY")
rows <- list()
for (g in rownames(genus_rel)) for (p in paths) {
  ct <- suppressWarnings(cor.test(as.numeric(genus_rel[g, ]), as.numeric(path[p, ]), method = "spearman"))
  rows[[length(rows) + 1]] <- data.frame(Genus = g, Pathway = p, Rho = unname(ct$estimate), Pvalue = ct$p.value)
}
cors <- do.call(rbind, rows)
cors$FDR <- ave(cors$Pvalue, cors$Pathway, FUN = function(x) p.adjust(x, method = "BH"))
sigc <- subset(cors, FDR < 0.05 & Genus != "Unclassified")
write.csv(sigc, "Significant_Correlations_FDR.csv", row.names = FALSE)
cat("significant correlations:", nrow(sigc), "\n")
edges <- data.frame(Source = sigc$Genus, Target = sigc$Pathway, Interaction = "correlates_with",
                    Rho = sigc$Rho, Correlation_Sign = ifelse(sigc$Rho > 0, "Positive", "Negative"),
                    Weight = abs(sigc$Rho))
write.csv(edges, "Edge_Table.csv", row.names = FALSE)
