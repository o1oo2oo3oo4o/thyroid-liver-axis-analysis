################################################################################
#  03_deseq2_mouse_liver_T3.R
#  DESeq2 analysis of the hepatic T3 response in mouse liver (public dataset GEO GSE85793; vehicle C1-C4 vs T3 T1-T4, 4 mice per group)
#
#  Re-written for the repository from the original script 2_run_DESeq2.R (which asks for its settings at the console) and the R
#  console history of the run (250910_code_test/.Rhistory). Reproduction check against the original result file
#  (250912__DEG_all_genes.csv), run with R 4.3 / DESeq2 1.42.0:
#    - 12,135 genes after the low-count filter, identical baseMean and pvalue (relative difference < 1e-6),
#      the same genes with a missing padj (11,819 genes with a padj), padj identical within 3e-7;
#    - the Wald statistic of the file 261020__DEG_all_genes_with_stat.csv equals the statistic recomputed here (difference 0);
#    - the log2 fold changes of both files are the apeglm-shrunken estimates (apeglm is required, see below). apeglm could not be
#      installed in the checking environment, so the shrunken values themselves were not recomputed.
#
#  METHODS (as in the manuscript)
#    - genes with count >= 10 in at least half of the samples are kept (MIN_COUNT, half of the samples)
#    - DESeq2 default Wald test, design ~ condition, reference level Control
#    - log2 fold changes are shrunken with apeglm (lfcShrink(type = "apeglm")). lfcShrink calls results() with its default
#      alpha (0.1) for the p values; the padj of the output file is therefore the padj of results(dds) at alpha = 0.1, which also
#      fixes which genes are removed by independent filtering (padj = NA).
#    - responsive gene = padj < 0.05 and |log2FC| >= 1 (shrunken log2FC)
#
#  INPUTS (DATA_DIR)
#    1_cleaned_count_matrix.csv : genes x samples (C1-C4, T1-T4) raw counts      (made by 02_build_count_matrix_GSE85793.R)
#    1_cleaned_coldata.csv      : sample, condition (Control / T3)
#  OUTPUTS (OUT_DIR, default = DATA_DIR because later steps read them from there)
#    250912__DEG_all_genes.csv                        : gene, baseMean, log2FoldChange (apeglm), lfcSE, pvalue, padj; all 12,135 genes, ordered by padj
#    GSE85793_mouse_liver_T3_DESeq2_all_genes.csv     : identical copy under the name read by figures/Fig2_*.R
#    261020__DEG_all_genes_with_stat.csv              : the same table with the Wald statistic (stat; from the non-shrunken results) added; read by 07_redox_cameraPR.R
#    250912__DEG_sig_genes_lfc_1.csv                  : responsive genes (padj < 0.05, |log2FC| >= 1)
#    03_deseq2_mouse_liver_T3_sessionInfo.txt
#  (The file names keep the date prefixes of the files used for the manuscript.)
#
#  USAGE  DATA_DIR=/path/to/data Rscript 03_deseq2_mouse_liver_T3.R
#         SHRINK_TYPE=none Rscript ...   # no shrinkage (apeglm not installed): padj and stat are unchanged, but log2FC, lfcSE and the
#                                        # number of responsive genes differ from the manuscript (779 instead of 541 in the reproduction check)
################################################################################
suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr); library(tibble); library(stringr) })

DATA_DIR    <- Sys.getenv("DATA_DIR", "data")
OUT_DIR     <- Sys.getenv("OUT_DIR", DATA_DIR)
SHRINK_TYPE <- Sys.getenv("SHRINK_TYPE", "apeglm")   # "apeglm" (manuscript) or "none"
stopifnot(SHRINK_TYPE %in% c("apeglm", "none"))
if (SHRINK_TYPE == "apeglm" && !requireNamespace("apeglm", quietly = TRUE))
  stop("The apeglm package is required (BiocManager::install('apeglm')); set SHRINK_TYPE=none to run without shrinkage.")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

count_path <- file.path(DATA_DIR, "1_cleaned_count_matrix.csv")   # rows = genes, columns = samples
col_path   <- file.path(DATA_DIR, "1_cleaned_coldata.csv")        # columns: sample, condition
alpha_fdr  <- 0.05      # padj cut
lfc_cut    <- 1         # |log2FC| cut
min_count  <- 10L       # low-count filter: count >= min_count in at least half of the samples

## ===== 1) Load and tidy =====
cnt <- read_csv(count_path, show_col_types = FALSE) |> as.data.frame()
rownames(cnt) <- cnt[[1]]; cnt[[1]] <- NULL
cnt[] <- lapply(cnt, function(x) as.integer(round(as.numeric(x)))); cnt[is.na(cnt)] <- 0L
colnames(cnt) <- str_trim(gsub("\r", "", colnames(cnt)))
meta <- read_csv(col_path, show_col_types = FALSE) |>
  mutate(sample    = str_trim(gsub("\r", "", as.character(sample))),
         condition = str_trim(gsub("\r", "", as.character(condition))))

## ===== 2) Use the common samples only, in the same order =====
common <- intersect(colnames(cnt), meta$sample)
stopifnot(length(common) > 0)
cnt  <- cnt[, common, drop = FALSE]
meta <- meta |> filter(sample %in% common) |>
  arrange(match(sample, colnames(cnt))) |>
  column_to_rownames("sample")

## ===== 3) Condition as a factor (Control is the reference) =====
lvl <- unique(meta$condition)
if ("Control" %in% lvl) {
  meta$condition <- factor(meta$condition, levels = c("Control", setdiff(lvl, "Control")))
} else {
  meta$condition <- factor(meta$condition, levels = lvl)   # without a "Control" label the first level is the reference
}

## ===== 4) Low-count filter =====
keep <- rowSums(cnt >= min_count) >= ceiling(ncol(cnt) / 2)
cnt  <- cnt[keep, , drop = FALSE]

## ===== 5) DESeq2 =====
dds <- DESeqDataSetFromMatrix(countData = cnt, colData = meta, design = ~ condition)
dds <- DESeq(dds)

## ===== 6) Coefficient (explicit) + LFC shrinkage =====
rn <- resultsNames(dds)                                   # e.g. "Intercept", "condition_T3_vs_Control"
coef_name <- rn[grepl("^condition_.*_vs_.*$", rn, ignore.case = TRUE)]
stopifnot(length(coef_name) == 1)
res_wald <- results(dds, name = coef_name)                # Wald test, default alpha (0.1) as used inside lfcShrink
if (SHRINK_TYPE == "apeglm") {
  suppressPackageStartupMessages(library(apeglm))
  res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
} else res <- res_wald
res <- res[order(res$padj), ]

## ===== 7) Save =====
res_df <- as.data.frame(res) |> rownames_to_column("gene")
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
readr::write_csv(res_df, file.path(OUT_DIR, "250912__DEG_all_genes.csv"))
readr::write_csv(res_df, file.path(OUT_DIR, "GSE85793_mouse_liver_T3_DESeq2_all_genes.csv"))
stat_df <- as.data.frame(res_wald) |> rownames_to_column("gene") |> select(gene, stat)
with_stat <- res_df |> left_join(stat_df, by = "gene") |> select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
readr::write_csv(with_stat, file.path(OUT_DIR, "261020__DEG_all_genes_with_stat.csv"))
sig_df <- res_df |> filter(!is.na(padj), padj < alpha_fdr, abs(log2FoldChange) >= lfc_cut)
readr::write_csv(sig_df, file.path(OUT_DIR, "250912__DEG_sig_genes_lfc_1.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "03_deseq2_mouse_liver_T3_sessionInfo.txt"))
cat(sprintf("[OK] coef=%s | shrinkage=%s | genes after filter=%d | genes with padj=%d | responsive (padj < %.2f, |log2FC| >= %g)=%d (up %d, down %d)\n",
            coef_name, SHRINK_TYPE, nrow(res_df), sum(!is.na(res_df$padj)), alpha_fdr, lfc_cut, nrow(sig_df),
            sum(sig_df$log2FoldChange > 0), sum(sig_df$log2FoldChange < 0)))
cat("[OK] Manuscript values: 11,819 genes tested, 541 responsive (266 up, 275 down)\n")
cat("[OK] Written to", OUT_DIR, "\n")
