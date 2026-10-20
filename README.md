[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23359812.svg)](https://doi.org/10.5281/zenodo.23359812)
# Integrated 16S rRNA bioinformatics analysis of gut microbiota, predicted pathways and candidate biomarkers in IBD

Code accompanying the preprint "Integrated 16S rRNA Bioinformatics Analysis of Gut Microbiota, Predicted Pathways, and Candidate Biomarkers in Inflammatory Bowel Disease".

## Preprint
This work is available as a preprint on bioRxiv (not yet peer-reviewed):

**Integrated 16S rRNA Bioinformatics Analysis of Gut Microbiota, Predicted Pathways, and Candidate Biomarkers in Inflammatory Bowel Disease**

Aman Zehra et al. bioRxiv 2026. doi: [10.64898/2026.10.16.759573](https://doi.org/10.64898/2026.10.16.759573)
> **Exploratory study.** 22 public stool 16S samples (11 IBD, Shanghai; 11 healthy, Singapore). The group is confounded with source cohort, and pathways are PICRUSt2 predictions. Results are hypothesis-generating.

## Pipeline overview
1. **Preprocessing:** SRA Toolkit → FastQC → Trimmomatic → KneadData/Bowtie2 (hg38 host removal)
2. **Taxonomy and diversity:** Kraken2 + Bracken; alpha (Shannon, Simpson) and beta (Bray-Curtis, PCoA, PERMANOVA)
3. **ASVs and function:** DADA2 (R) → PICRUSt2 (Galaxy Europe), pathway abundances normalised per million
4. **Statistics:** DESeq2 differential abundance; Spearman genus-pathway correlations (BH-FDR)
5. **Networks:** Cytoscape and NetworkAnalyzer
6. **Machine learning (Python):** Logistic Regression, Random Forest, SVM, XGBoost; ROC-AUC; biomarker prioritisation; external validation; SCFA-related pathways

## Data
- Raw reads (NCBI SRA): IBD SRR33757848–SRR33757858; healthy controls SRR37842383–SRR37842393
- Supplementary Tables S1–S32 are in `supplementary/`
- External validation: 43 stool samples (23 Crohn's disease, 20 controls) from the MicrobiomeAnalyst example dataset [ibd_data.zip](https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/data/ibd_data.zip)

## Repository layout
| Path | Content |
|---|---|
| `notebooks/` | Google Colab notebook (ML, ROC, network integration, external validation, SCFA). Paths point to Google Drive; edit before running. |
| `python/` | Diversity script and Python taxonomy classifier |
| `R/` | R scripts (see reproducibility notes) |
| `cytoscape/` | Cytoscape sessions and network files |
| `supplementary/` | Supplementary Tables S1–S32 |
| `data/`, `genus_counts.csv` | Genus count table and ASV taxonomy (SILVA 138.2) |
| `recovered_original_R/` | Command histories from the original RStudio sessions |

## Software
Python 3 (see `requirements.txt`); R with dada2, phyloseq, DESeq2, vegan, ggplot2, dplyr, readr; Cytoscape 3.10.3; FastQC 0.11.9, Trimmomatic 0.39, Bowtie2 2.5.1, Kraken2 2.1.3, Bracken, PICRUSt2.

## Reproducibility
See [REPRODUCIBILITY_NOTES.md](REPRODUCIBILITY_NOTES.md) for the status of the R scripts and the taxonomy correction in v1.1.

## AI assistance
AI-assisted tools were used to help write analysis code; the authors verified the outputs against the supplementary tables.

## License
MIT (see `LICENSE`). Cite the manuscript if you use this code.

