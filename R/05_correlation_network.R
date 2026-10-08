# 05 - Genus x pathway Spearman correlations and network tables  [RECONSTRUCTED]
# Inputs : genus relative abundance (ASV counts summed by genus, DADA2 taxonomy,
#          genera present in >= 2 samples), per-million pathway abundance (S8b),
#          DESeq2 significant genera (S10).
# Outputs: significant correlations (S12) and the edge/node tables for Cytoscape (S25/S26).
# Reference results: 95 significant correlations (59 LPSSYN-PWY, 36 PEPTIDOGLYCANSYN-PWY).
source("R/00_setup.R")
genus <- read.csv("genus_counts.csv", row.names = 1, check.names = FALSE)
genus_rel <- sweep(genus, 2, colSums(genus), "/")
genus_rel <- genus_rel[rowSums(genus_rel > 0) >= 2, ]
path <- read.csv(supp("Supplementary_Table_S8b_PICRUSt2_Pathway_Abundance_per_million.csv"),
                 row.names = 1, check.names = FALSE)
paths <- c("LPSSYN-PWY", "PEPTIDOGLYCANSYN-PWY")
common <- intersect(colnames(genus_rel), colnames(path))
sig <- read.csv(supp("Supplementary_Table_S10_DESeq2_Significant_Genera.csv"))$Genus
rows <- list()
for (g in intersect(sig, rownames(genus_rel))) for (p in paths) {
  ct <- suppressWarnings(cor.test(as.numeric(genus_rel[g, common]),
                                  as.numeric(path[p, common]), method = "spearman"))
  rows[[length(rows) + 1]] <- data.frame(Genus = g, Pathway = p, Rho = unname(ct$estimate), P = ct$p.value)
}
cors <- do.call(rbind, rows)
cors$FDR <- p.adjust(cors$P, method = "BH")
sigc <- subset(cors, FDR < 0.05)
write.csv(sigc, "Significant_Correlations_FDR.csv", row.names = FALSE)
cat("significant correlations:", nrow(sigc), "\n")
# Edge table for Cytoscape (import: Source = Source, Target = Target)
edges <- data.frame(Source = sigc$Genus, Target = sigc$Pathway,
                    Interaction = "correlates_with", Rho = sigc$Rho,
                    Correlation_Sign = ifelse(sigc$Rho > 0, "Positive", "Negative"),
                    Weight = abs(sigc$Rho))
write.csv(edges, "Edge_Table.csv", row.names = FALSE)
# Hub analysis (degree, betweenness, closeness): Cytoscape > Tools > NetworkAnalyzer
# (undirected), as reported in Tables S27/S27b.
