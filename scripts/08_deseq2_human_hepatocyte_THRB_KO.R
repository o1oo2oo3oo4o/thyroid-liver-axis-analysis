################################################################################
#  08_deseq2_human_hepatocyte_THRB_KO.R
#  2 x 2 DESeq2 analysis of the T3 response in human iPSC-derived hepatocytes, wild type (WT) versus THRB knockout (KO)
#  (public dataset GEO GSE270233; vehicle or T3, n = 3 per group, 12 samples). Output = Supplementary Table S9a.
#
#  METHODS (as in the manuscript)
#    - genes with count >= 10 in at least 3 samples are kept (the size of the smallest group)
#    - DESeq2 default Wald test, design ~ genotype + treatment + genotype:treatment (reference levels WT and Veh)
#    - four contrasts
#        WT_T3_effect    T3 versus vehicle in WT cells                       coefficient treatment_T3_vs_Veh
#        KO_T3_effect    T3 versus vehicle in THRB-KO cells                  treatment_T3_vs_Veh + genotypeKO.treatmentT3
#        Interaction     KO T3 effect minus WT T3 effect                     coefficient genotypeKO.treatmentT3
#        KO_vs_WT_basal  KO versus WT in vehicle-treated cells               coefficient genotype_KO_vs_WT
#    - the reported log2FoldChange, lfcSE, stat, pvalue and padj are the unshrunken DESeq2 values (results(dds, alpha = 0.05));
#      an additional column lfc_shrunk holds the ashr-shrunken log2 fold change (not used in any test of the manuscript)
#
#  INPUTS (DATA_DIR)
#    GSE270233_counts.csv          genes x 12 samples (columns S1-S12), integer counts. The column order is the order of the GEO
#                                  sample list: S1-S3 WT vehicle (GSM8337685-87), S4-S6 WT T3 (GSM8337688-90),
#                                  S7-S9 KO vehicle (GSM8337691-93), S10-S12 KO T3 (GSM8337694-96).
#                                  Gene names have hyphens replaced by dots (HLA-A -> HLA.A).
#    GSE270233_series_matrix.txt   (optional) GEO series matrix; if present, the sample table below is checked against it
#    candidate_genes.txt           (optional, CANDIDATE_FILE) one gene symbol per line; adds candidate_genes_response.csv
#
#  OUTPUTS (OUT_DIR, default DATA_DIR/GSE270233_DESeq2)
#    coldata_used.csv, DEG_WT_T3_effect.csv, DEG_KO_T3_effect.csv, DEG_Interaction.csv, DEG_KO_vs_WT_basal.csv,
#    DEG_summary_counts.csv, T3_effect_WT_KO_interaction_wide.csv (read by 11_THRB_candidate_enrichment.R and the figure scripts),
#    Volcano_<contrast>.png / .pdf, QC_PCA.pdf, 08_deseq2_human_hepatocyte_THRB_KO_sessionInfo.txt
#
#  REPRODUCTION CHECK  The original run used R 4.5.2 and DESeq2 1.50.2. Re-running this script in R 4.3.3 / DESeq2 1.42.0 on
#  GSE270233_counts.csv gave the same genes, baseMean, log2FoldChange, lfcSE, stat, pvalue and padj as the original files
#  DEG_*.csv and T3_effect_WT_KO_interaction_wide.csv (see docs/PROVENANCE.md). The lfc_shrunk column needs the ashr package,
#  which was not available in the checking environment, so it was not recomputed.
#
#  USAGE  DATA_DIR=/path/to/data Rscript 08_deseq2_human_hepatocyte_THRB_KO.R
################################################################################
suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(readr); library(tibble); library(ggplot2)
})

## ---- settings ----------------------------------------------------------------
DATA_DIR       <- Sys.getenv("DATA_DIR", "data")
OUT_DIR        <- Sys.getenv("OUT_DIR", file.path(DATA_DIR, "GSE270233_DESeq2"))
count_file     <- file.path(DATA_DIR, "GSE270233_counts.csv")
matrix_file    <- file.path(DATA_DIR, "GSE270233_series_matrix.txt")
CANDIDATE_FILE <- Sys.getenv("CANDIDATE_FILE", "")   # optional: one gene symbol per line, no header
min_count      <- 10     # a gene needs at least this count ...
min_samples    <- 3      # ... in at least this many samples (size of the smallest group)
alpha          <- 0.05
lfc_cut        <- 1      # colour threshold |log2FC| of the volcano plots (not a filter of the analysis)
n_label        <- 10     # number of top genes (by padj) labelled in each volcano plot; 0 = no labels
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
HAVE_ASHR <- requireNamespace("ashr", quietly = TRUE)
if (!HAVE_ASHR) message("Package 'ashr' not installed: the column lfc_shrunk will be NA (no other result depends on it).")

## ---- 1. counts and sample table ------------------------------------------------
cnt <- read.csv(count_file, row.names = 1, check.names = FALSE)
stopifnot(ncol(cnt) == 12, all(cnt == round(cnt)))

coldata <- data.frame(
  sample    = paste0("S", 1:12),
  genotype  = factor(rep(c("WT", "KO"), each = 6), levels = c("WT", "KO")),
  treatment = factor(rep(rep(c("Veh", "T3"), each = 3), 2), levels = c("Veh", "T3")),
  row.names = paste0("S", 1:12)
)
stopifnot(identical(colnames(cnt), rownames(coldata)))

# Check the assumed sample order against the GEO series matrix (samples are listed in the order of the count columns)
if (file.exists(matrix_file)) {
  ml  <- readLines(matrix_file, warn = FALSE)
  chr  <- ml[startsWith(ml, "!Sample_characteristics_ch1")]
  geno <- chr[grepl("genotype:", chr, fixed = TRUE)]
  trt  <- chr[grepl("treatment:", chr, fixed = TRUE)]
  split_q <- function(l) gsub('^"|"$', "", strsplit(l, "\t", fixed = TRUE)[[1]][-1])
  g_geo <- ifelse(grepl("knockout", split_q(geno[1]), fixed = TRUE), "KO", "WT")
  t_geo <- ifelse(grepl("T3", split_q(trt[1]), fixed = TRUE), "T3", "Veh")
  stopifnot(identical(g_geo, as.character(coldata$genotype)), identical(t_geo, as.character(coldata$treatment)))
  cat("Sample table agrees with the GEO series matrix (genotype and treatment of all 12 samples).\n")
} else {
  cat("No GSE270233_series_matrix.txt found; the sample order S1-S12 = GSM8337685-GSM8337696 was not checked.\n")
}
print(table(coldata$genotype, coldata$treatment))
write_csv(coldata, file.path(OUT_DIR, "coldata_used.csv"))

## ---- 2. DESeq2 -----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(as.matrix(cnt), coldata, design = ~ genotype + treatment + genotype:treatment)
keep <- rowSums(counts(dds) >= min_count) >= min_samples
cat(sprintf("Filter: %d -> %d genes (count >= %d in >= %d samples)\n", nrow(dds), sum(keep), min_count, min_samples))
dds <- dds[keep, ]
dds <- DESeq(dds)
rn <- resultsNames(dds); print(rn)
stopifnot(all(c("genotype_KO_vs_WT", "treatment_T3_vs_Veh", "genotypeKO.treatmentT3") %in% rn))

## ---- 3. contrasts --------------------------------------------------------------
get_res <- function(name = NULL, contrast = NULL) {
  r  <- if (!is.null(name)) results(dds, name = name, alpha = alpha) else results(dds, contrast = contrast, alpha = alpha)
  df <- as.data.frame(r) %>% rownames_to_column("gene")
  if (HAVE_ASHR) {
    s <- if (!is.null(name)) lfcShrink(dds, coef = name, res = r, type = "ashr")
         else lfcShrink(dds, contrast = contrast, res = r, type = "ashr")
    df$lfc_shrunk <- as.data.frame(s)[df$gene, "log2FoldChange"]
  } else {
    df$lfc_shrunk <- NA_real_
  }
  df %>% arrange(padj)
}

res_list <- list(
  WT_T3_effect   = get_res(name = "treatment_T3_vs_Veh"),
  KO_T3_effect   = get_res(contrast = list(c("treatment_T3_vs_Veh", "genotypeKO.treatmentT3"))),
  Interaction    = get_res(name = "genotypeKO.treatmentT3"),   # KO T3 effect minus WT T3 effect
  KO_vs_WT_basal = get_res(name = "genotype_KO_vs_WT")         # KO versus WT, vehicle
)
for (nm in names(res_list)) write_csv(res_list[[nm]], file.path(OUT_DIR, paste0("DEG_", nm, ".csv")))

summ <- bind_rows(lapply(names(res_list), function(nm) {
  d <- res_list[[nm]]
  tibble(contrast  = nm,
         tested    = sum(!is.na(d$padj)),
         padj_lt   = sum(d$padj < alpha, na.rm = TRUE),
         up        = sum(d$padj < alpha & d$log2FoldChange > 0, na.rm = TRUE),
         down      = sum(d$padj < alpha & d$log2FoldChange < 0, na.rm = TRUE),
         padj_lfc1 = sum(d$padj < alpha & abs(d$log2FoldChange) >= 1, na.rm = TRUE))
}))
print(summ); write_csv(summ, file.path(OUT_DIR, "DEG_summary_counts.csv"))

## ---- 4. the three T3 contrasts in one table -----------------------------------------
pick <- function(d, tag) d %>% select(gene, log2FoldChange, lfcSE, pvalue, padj) %>%
  rename_with(~ paste0(.x, "_", tag), -gene)
wide <- pick(res_list$WT_T3_effect, "WT") %>%
  full_join(pick(res_list$KO_T3_effect, "KO"), by = "gene") %>%
  full_join(pick(res_list$Interaction, "int"), by = "gene")
write_csv(wide, file.path(OUT_DIR, "T3_effect_WT_KO_interaction_wide.csv"))

cand <- NULL
if (nzchar(CANDIDATE_FILE) && file.exists(CANDIDATE_FILE)) {
  # gene names of the count file have hyphens replaced by dots (HLA-A -> HLA.A); bring both sides to the same form
  cand <- trimws(readLines(CANDIDATE_FILE)); cand <- cand[cand != ""]
  key  <- function(x) gsub("-", ".", x, fixed = TRUE)
  out  <- tibble(candidate = cand, gene = key(cand)) %>%
    left_join(wide, by = "gene") %>%
    mutate(in_tested_set = gene %in% wide$gene, in_raw_counts = gene %in% rownames(cnt))
  write_csv(out, file.path(OUT_DIR, "candidate_genes_response.csv"))
  cat(sprintf("%d candidates: %d tested, %d only in the raw counts, %d not in the raw counts\n",
              nrow(out), sum(out$in_tested_set), sum(!out$in_tested_set & out$in_raw_counts), sum(!out$in_raw_counts)))
}

## ---- 5. volcano plots (PNG and PDF) ------------------------------------------------
# x axis = unshrunken log2 fold change (the same quantity as the up/down counts of DEG_summary_counts.csv)
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
cand_keys <- if (!is.null(cand)) gsub("-", ".", cand, fixed = TRUE) else character(0)

volcano_titles <- c(
  WT_T3_effect   = "WT: T3 vs vehicle",
  KO_T3_effect   = "THRB KO: T3 vs vehicle",
  Interaction    = "Interaction (KO T3 effect - WT T3 effect)",
  KO_vs_WT_basal = "Vehicle: THRB KO vs WT"
)

make_volcano <- function(d, nm) {
  d <- d %>%
    filter(!is.na(padj)) %>%
    mutate(neg_log10_padj = -log10(pmax(padj, 1e-300)),
           status = case_when(padj < alpha & log2FoldChange >=  lfc_cut ~ "Up",
                              padj < alpha & log2FoldChange <= -lfc_cut ~ "Down",
                              TRUE ~ "NS"))
  n_up <- sum(d$status == "Up"); n_down <- sum(d$status == "Down")
  lab <- d %>% filter(status != "NS") %>% arrange(padj) %>% slice_head(n = n_label)
  if (length(cand_keys) > 0) lab <- bind_rows(lab, d %>% filter(gene %in% cand_keys)) %>% distinct(gene, .keep_all = TRUE)

  p <- ggplot(d, aes(x = log2FoldChange, y = neg_log10_padj, color = status)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c(Up = "#D55E00", Down = "#0072B2", NS = "grey70"),
                       breaks = c("Up", "Down"),
                       labels = c(sprintf("Up (%d)", n_up), sprintf("Down (%d)", n_down))) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", linewidth = 0.3) +
    labs(title = volcano_titles[[nm]],
         subtitle = sprintf("padj < %g, |log2FC| >= %g; %d genes tested", alpha, lfc_cut, nrow(d)),
         x = expression(log[2]~"fold change"), y = expression(-log[10]~"(adjusted P)"), color = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")

  if (nrow(lab) > 0 && n_label > 0) {
    p <- p + if (has_repel) {
      ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3, color = "black",
                               max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2)
    } else {
      geom_text(data = lab, aes(label = gene), size = 2.5, color = "black", vjust = -0.5, check_overlap = TRUE)
    }
  }
  base <- file.path(OUT_DIR, paste0("Volcano_", nm))
  ggsave(paste0(base, ".png"), p, width = 6, height = 5.5, dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), p, width = 6, height = 5.5)
  invisible(p)
}
for (nm in names(res_list)) make_volcano(res_list[[nm]], nm)

## ---- 6. QC ---------------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)
pdf(file.path(OUT_DIR, "QC_PCA.pdf"), width = 6, height = 5)
print(plotPCA(vsd, intgroup = c("genotype", "treatment")))
dev.off()

chk <- intersect(c("THRB", "THRA", "ANGPTL4", "DIO1"), rownames(vsd))
print(round(assay(vsd)[chk, ], 2))   # sanity check of the sample order (THRB should be low in the KO samples)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "08_deseq2_human_hepatocyte_THRB_KO_sessionInfo.txt"))
cat("\nDone:", normalizePath(OUT_DIR), "\n")
