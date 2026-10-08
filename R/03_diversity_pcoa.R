# 03 - PCoA / PERMANOVA plot from the Bray-Curtis results (Tables S5-S7)
source("R/00_setup.R")
library(vegan)
pcoa <- read_csv(supp("Supplementary_Table_S6_PCoA_Coordinates.csv"), show_col_types = FALSE)
colnames(pcoa)[1] <- "Sample"
pcoa$Group <- group_from_sample(pcoa$Sample)     # fixes the earlier inverted labels
ggplot(pcoa, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 4) + stat_ellipse(level = 0.95, linetype = 2) +
  scale_colour_manual(values = c(IBD = "#D55E00", Healthy = "#0072B2")) +
  theme_classic() + labs(title = "PCoA (Bray-Curtis)")
ggsave("Figure4_beta_diversity.png", width = 7, height = 5, dpi = 300)

# PERMANOVA on the distance matrix (S5)
d <- as.matrix(read.csv(supp("Supplementary_Table_S5_Bray_Curtis_Distance_Matrix.csv"), row.names = 1, check.names = FALSE))
grp <- group_from_sample(rownames(d))
print(adonis2(as.dist(d) ~ grp, permutations = 9999))
