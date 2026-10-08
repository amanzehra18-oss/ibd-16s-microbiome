# 06 - Integrate DESeq2, taxonomy and correlations (Table S14)
# Cleaned version of the recovered original commands (recovered_original_R/integration.Rhistory).
source("R/00_setup.R")
deseq <- read_csv("DESeq2_Significant_Genera.csv", show_col_types = FALSE)
network <- read_csv("Significant_Correlations_FDR.csv", show_col_types = FALSE)
taxonomy <- read_csv("taxonomy.csv", show_col_types = FALSE) %>% distinct(Genus, .keep_all = TRUE)
final_table <- deseq %>%
  left_join(taxonomy, by = "Genus") %>%
  left_join(network, by = "Genus") %>%
  mutate(Pathway = recode(Pathway,
    "LPSSYN-PWY" = "Lipopolysaccharide biosynthesis",
    "PEPTIDOGLYCANSYN-PWY" = "Peptidoglycan biosynthesis")) %>%
  select(Genus, Phylum, Family, log2FoldChange, padj, Pathway, Rho, FDR)
write.csv(final_table, "Final_Integrated_Taxonomy_Function.csv", row.names = FALSE)
