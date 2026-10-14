[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23359812.svg)](https://doi.org/10.5281/zenodo.23359812)
# Integrated 16S rRNA bioinformatics analysis of gut microbiota, predicted pathways and candidate biomarkers in IBD

Code accompanying the manuscript *"Integrated 16S rRNA Bioinformatics Analysis of Gut Microbiota,
Predicted Pathways, and Candidate Biomarkers in Inflammatory Bowel Disease"* (preprint: [link to be added]).

> **Exploratory study.** 22 public stool 16S samples (11 IBD, Shanghai; 11 healthy, Singapore).
> Group is confounded with source cohort, and pathways are PICRUSt2 *predictions*. Results are hypothesis-generating.

## Data
- Raw reads (NCBI SRA): IBD `SRR33757848-SRR33757858`, healthy controls `SRR37842383-SRR37842393`.
- Supplementary Tables S1-S32 (incl. S8b, S27b) are provided with the manuscript (Supplementary_Data_1_v8, already included in `supplementary/`).
- External validation: 43 stool samples (23 Crohn's disease, 20 controls) from the MicrobiomeAnalyst example
  dataset `ibd_data.zip` (https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/data/ibd_data.zip).

## Pipeline overview
1. SRA Toolkit -> FastQC -> Trimmomatic -> KneadData/Bowtie2 (hg38 host removal).
2. Kraken2 + Bracken taxonomy and diversity (alpha: Shannon/Simpson; beta: Bray-Curtis, PCoA, PERMANOVA).
3. DADA2 ASVs and taxonomy (R); reverse-oriented ASVs of nine IBD samples reverse-complemented;
   PICRUSt2 on Galaxy Europe on all 22 samples; pathway abundances normalised per million.
4. DESeq2 differential abundance; Spearman genus-pathway correlations (BH-FDR); Cytoscape network and NetworkAnalyzer hubs.
5. Machine learning (LR, RF, SVM, XGBoost), ROC-AUC, biomarker prioritisation, external validation, SCFA-related pathways (Python).

## Repository layout
| Path | Content |
|---|---|
| `notebooks/` | Original Google Colab notebook (ML, ROC, network integration, external validation, SCFA). Paths point to Google Drive; edit them before running. Cell order reflects the order of execution in Colab, not step numbers. |
| `python/diversity_bracken.py` | Alpha/beta diversity and PERMANOVA for all 22 samples. |
| `R/` | R scripts, see status below. |
| `cytoscape/` | Cytoscape sessions (.cys) and network files (.xgmml) for Figures 5C and 7B. |
| `supplementary/` | Supplementary Tables S1-S32 (v8, corrected genus taxonomy). |
| `genus_counts.csv`, `data/` | Genus count table and ASV taxonomy (SILVA 138.2, tryRC) used by `R/04` and `R/05`. |
| `python/` | Diversity script plus the Python taxonomy classifier (RDP naive-Bayes re-implementation) used to cross-check the tryRC taxonomy. |
| `recovered_original_R/` | Command histories recovered from the original RStudio sessions (kept for transparency). |

## Status of the R scripts - please read
The original R working files were lost. Scripts `R/01`, `R/04` and `R/05` are **reconstructions**. In v1.1 the genus taxonomy was
re-assigned with `tryRC = TRUE` (`R/01b`) because the ASVs of nine IBD samples are reverse-oriented and had received no genus; the
genus-level tables (S3, S10-S14, S20, S24-S28, S27b) were regenerated from the corrected taxonomy. DESeq2 (`R/04`) was run in R on the corrected genus table (output in `data/DESeq2_All_Genera_R_output.csv`): 186 genera tested, 44 significant (12 higher, 32 lower in IBD). `R/05_correlation_network.R` reference: 84 significant correlations (47 LPSSYN-PWY, 37 PEPTIDOGLYCANSYN-PWY). Scripts `R/03` and `R/06` are cleaned versions of recovered original commands.
Parameters marked `[CHECK]` must be confirmed.

## Software
Python 3 (see `requirements.txt`); R with dada2, phyloseq, DESeq2, vegan, ggplot2, dplyr, readr; Cytoscape 3.10.3;
FastQC 0.11.9, Trimmomatic 0.39, Bowtie2 2.5.1, Kraken2 2.1.3, Bracken, PICRUSt2 (Galaxy Europe).
Record exact versions of R and each package with `sessionInfo()`.

## AI assistance
AI-assisted tools were used to help write analysis code; the authors verified the outputs against the supplementary tables.

## License
MIT (see `LICENSE`). Cite the manuscript if you use this code.
