# Reproducibility notes

## Status of the R scripts
The original R working files were not preserved. Scripts `R/01`, `R/04` and `R/05` are reconstructions; `R/03` and `R/06` are cleaned versions of commands recovered from the original sessions (see `recovered_original_R/`). Outputs were checked against the supplementary tables. Parameters marked [CHECK] in the scripts still need to be confirmed.

## Taxonomy correction in v1.1
ASVs of nine IBD samples were reverse-oriented and initially received no genus. In v1.1 the genus taxonomy was re-assigned with `tryRC = TRUE` (`R/01b`), and all genus-level tables were regenerated from the corrected taxonomy. DESeq2 (`R/04`) output is in `data/DESeq2_All_Genera_R_output.csv`.

## Software versions
Please record exact R and package versions with `sessionInfo()` when re-running.
