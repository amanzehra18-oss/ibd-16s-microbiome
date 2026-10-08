# Packages used across the R scripts. Install once.
# install.packages(c("tidyverse","vegan","ggplot2","readr","dplyr","BiocManager"))
# BiocManager::install(c("dada2","phyloseq","DESeq2"))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

# Set this to the folder holding the supplementary tables (S1-S32).
supp_dir <- "supplementary"
supp <- function(f) file.path(supp_dir, f)

# Samples 1-11 are IBD ("Disease"), 12-22 are healthy controls (see Table S1).
# Always derive the group from the sample NUMBER, never from row order.
group_from_sample <- function(ids) {
  n <- as.integer(sub("^X?([0-9]+).*$", "\\1", ids))
  ifelse(n <= 11, "IBD", "Healthy")
}
