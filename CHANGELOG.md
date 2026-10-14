# Changelog

## v1.1 (October 2026)
- Fixed genus taxonomy: ASVs of nine IBD samples (2-5, 7-11) are reverse-oriented and had no genus assigned (0%); re-assigned with SILVA 138.2 and `tryRC = TRUE` (genus-level assignment 39% -> 61% of ASVs; 874 of 892 previously assigned ASVs keep the same genus, none changed).
- Regenerated S3, S10-S14, S20, S24-S28, S27b (see `supplementary/Supplementary_Data_Index.txt`, version note v8).
- DESeq2: 186 genera tested, 44 significant (12 higher, 32 lower in IBD; previously 175 tested and 49 significant), computed in R. Correlations: 84 significant (previously 95). Paludicola is no longer DESeq2-significant.
- Added `R/01b_assign_taxonomy_tryRC.R`, `genus_counts.csv`, `data/taxonomy_tryRC.csv`, `python/` taxonomy scripts, DESeq2 R output in `data/`.

## v1.0.1 / v1.0
- Initial release (deposited at Zenodo, concept DOI 10.5281/zenodo.23359812).
