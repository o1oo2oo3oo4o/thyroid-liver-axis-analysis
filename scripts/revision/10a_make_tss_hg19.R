#!/usr/bin/env Rscript
# =============================================================================
# 10a_make_tss_hg19.R      (run ONCE, before 10_TRb_ChIP_annotation.R)
#
# Purpose
#   Creates data/tss_hg19.csv  (columns: gene, chr, tss, strand) that 10_TRb_ChIP_annotation.R reads.
#   The Ayers 2014 TR-beta peaks are hg19 (GRCh37), so the TSS table must be hg19 too.
#
# What this script does NOT do
#   It does not download or build a genome. It reads GENE COORDINATES (UCSC knownGene gene models, hg19)
#   from the Bioconductor annotation package TxDb.Hsapiens.UCSC.hg19.knownGene.
#   If that package is not installed yet, the script installs it from Bioconductor (internet needed, one time).
#
# Working directory: scripts/revision (the folder that contains data/).  Output: data/tss_hg19.csv
# Packages: BiocManager (CRAN); TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db, GenomicFeatures (Bioconductor)
# =============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
need <- c("GenomicFeatures", "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  message("installing from Bioconductor: ", paste(miss, collapse = ", "))
  BiocManager::install(miss, ask = FALSE, update = FALSE)
}
still <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(still)) stop("could not install: ", paste(still, collapse = ", "),
                        "\n  check the internet connection and that your R version matches the Bioconductor release.")
suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db)
})

out_file <- file.path("data", "tss_hg19.csv")
dir.create("data", showWarnings = FALSE)

# ---- 1. gene ranges (one range per Entrez gene; genes on several strands/chromosomes are dropped by genes()) ----
g <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
cat("genes in TxDb:", length(g), "\n")
# keep chr1-22, X, Y only: alternative haplotype / unplaced contigs would give a second TSS for the same symbol
# (plain filtering on the chromosome name; no GenomeInfoDb function is needed, so it works in every Bioconductor release)
std_chr <- paste0("chr", c(1:22, "X", "Y"))
g <- g[as.character(seqnames(g)) %in% std_chr]
cat("genes on standard chromosomes:", length(g), "\n")

# ---- 2. Entrez ID -> current gene symbol ----
sym <- suppressMessages(mapIds(org.Hs.eg.db, keys = g$gene_id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
st  <- as.character(strand(g))
tss <- data.frame(gene   = unname(sym[g$gene_id]),
                  chr    = as.character(seqnames(g)),
                  tss    = ifelse(st == "+", start(g), end(g)),    # TSS = gene start on +, gene end on -
                  strand = st,
                  stringsAsFactors = FALSE)
tss <- tss[!is.na(tss$gene) & tss$gene != "" & tss$strand %in% c("+", "-"), ]
n_dup <- sum(duplicated(tss$gene))
if (n_dup) {
  message(n_dup, " rows share a symbol with another Entrez gene; keeping the first of each (10_TRb_ChIP_annotation.R does the same)")
  tss <- tss[!duplicated(tss$gene), ]
}
tss <- tss[order(tss$chr, tss$tss), ]
write.csv(tss, out_file, row.names = FALSE)
cat("\nwritten:", out_file, "|", nrow(tss), "genes | build hg19\n")
print(table(tss$strand))

# ---- 3. checks: are the genes used in 10_TRb_ChIP_annotation.R present? ----
panel <- c("CUX2","BDH2","ADH4","LILRB3","SBK1","ABHD1","INCENP","CDC42EP5","CD19","GPCPD1",
           "BTC","CFHR1","CLIC5","COL11A2","COL15A1","DEF6","DNAJB4","F13B","LGALS4","PKLR",
           "POR","SPARCL1","TIMD4","GSTM1","HSPA1A","SLC39A8")                 # same 26 genes as panel_genes in config.R
want <- panel
set_file <- file.path("data", "redox_gene_sets_human.csv")
if (file.exists(set_file)) {
  rs <- read.csv(set_file, stringsAsFactors = FALSE)
  if ("human_symbol" %in% names(rs)) want <- unique(c(panel, trimws(rs$human_symbol)))
}
miss_g <- setdiff(want, tss$gene)
cat("\ngenes checked:", length(want), "| missing from the TSS table:", length(miss_g), "\n")
if (length(miss_g)) cat("  missing:", paste(miss_g, collapse = ", "),
                        "\n  (usually an older/newer symbol or a gene not in UCSC knownGene; 10_TRb_ChIP_annotation.R reports them again and skips them)\n")
print(tss[tss$gene %in% c("GPCPD1", "ADH4", "F13B"), ], row.names = FALSE)
