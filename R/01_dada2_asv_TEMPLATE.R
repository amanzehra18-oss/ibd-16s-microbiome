# 01 - ASV inference with DADA2  (TEMPLATE - not the original script)
# The original script was lost. Parameters marked [CHECK] must be set to the
# values used for the published results before this is run.
# NOTE: the ASVs analysed in the manuscript were 120 bp, and the nine IBD samples
# 2-5 and 7-11 were in reverse-complement orientation. Inspect read orientation
# and primers (e.g. with cutadapt) on the raw reads before filtering.
library(dada2)
fq_dir <- "fastq_trimmed"                     # [CHECK] host-filtered reads
fnFs <- sort(list.files(fq_dir, pattern = "\\.fastq(\\.gz)?$", full.names = TRUE))
sample.names <- sub("\\.fastq.*$", "", basename(fnFs))
filtFs <- file.path("filtered", paste0(sample.names, "_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, truncLen = 120,   # [CHECK]
                     maxN = 0, maxEE = 2, truncQ = 2, compress = TRUE)
errF <- learnErrors(filtFs, multithread = TRUE)
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
seqtab <- makeSequenceTable(dadaFs)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138_train_set.fa.gz",  # [CHECK] database
                       multithread = TRUE)
write.csv(t(seqtab.nochim), "ASV_table.csv")
write.csv(taxa, "taxonomy.csv")
# ASV_table.csv + a FASTA of the ASV sequences were then run through PICRUSt2 on
# Galaxy Europe, after reverse-complementing the reverse-oriented ASVs (see 02).
