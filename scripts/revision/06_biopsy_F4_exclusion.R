#!/usr/bin/env Rscript
# =============================================================================
# 06_biopsy_F4_exclusion.R        (revision analysis; Reviewer 1 comment 6, Reviewer 2 major comment 4)
#
# Reviewer 1 comment 6 and Reviewer 2 major comment 4:
#   * exact F0-F4 sample counts per cohort (GSE130970 contains only two F4 biopsies)
#   * all models of Table S12 recomputed from raw counts, with cohort-specific 95% CI,
#     fixed-effect meta-analysis, Cochran's Q and I2
#   * NEW: hepatocyte marker score KEPT, stage F4 EXCLUDED
#       M6  stage + hepatocyte score                  (F4 excluded)
#       M7  stage + hepatocyte score + NAS            (F4 excluded; only where NAS is available)
#     to answer whether the four hepatocyte-adjusted associations (ABHD1, F13B, ADH4, LGALS4)
#     depend on cirrhotic F4 samples.
#
# Models (cohort-specific; GSE135251 has no age/sex, so its covariates are reduced automatically)
#   M1 primary (stage [+ sex + age])          M2 + hepatocyte score        M3 primary + NAS
#   M4 primary, F4 excluded                   M5 F3-F4 versus F0-F1 (F2 left out)
#   M6 hepatocyte score, F4 excluded (NEW)    M7 hepatocyte score + NAS, F4 excluded (NEW)
# Input : config.R (cohorts, meta_cols, panel/candidate genes), count matrices and sample tables
# Output: results/biopsy_F4/*.csv, *.pdf
# =============================================================================
suppressPackageStartupMessages({ library(DESeq2); library(data.table); library(ggplot2) })
source(Sys.getenv("STAGE_CFG", "config.R"))
out_dir <- file.path(out_root, "06_biopsy_F4"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
panel <- if (exists("panel_genes")) panel_genes else candidate_genes
focus <- candidate_genes

read_cohort <- function(cfg) {
  cnt <- fread(cfg$counts); setnames(cnt, 1, "gene"); cnt[, gene := trimws(gene)]; cnt <- cnt[!duplicated(gene)]
  m <- as.matrix(cnt[, -1]); rownames(m) <- cnt$gene; storage.mode(m) <- "integer"
  meta <- fread(cfg$meta); mc <- meta_cols
  meta[, sample := as.character(get(mc$sample))]; meta[, stage := as.integer(get(mc$stage))]
  if (!is.na(mc$group) && mc$group %in% names(meta)) meta <- meta[get(mc$group) == mc$group_keep]
  for (nm in c("nas", "age", "sex")) { col <- mc[[nm]]; meta[[nm]] <- if (!is.na(col) && col %in% names(meta)) meta[[col]] else NA }
  meta <- as.data.frame(meta[!is.na(stage) & sample %in% colnames(m)]); m <- m[, meta$sample, drop = FALSE]
  list(counts = m, meta = meta)
}
add_covariates <- function(meta, vsd) {
  meta$stage_num <- as.numeric(meta$stage)
  mk <- intersect(hepatocyte_markers, rownames(vsd)); z <- t(scale(t(vsd[mk, , drop = FALSE])))
  meta$hep <- as.numeric(scale(colMeans(z)))
  if (!all(is.na(meta$age))) meta$age_s <- as.numeric(scale(as.numeric(meta$age)))
  if (!all(is.na(meta$sex))) meta$sex <- factor(meta$sex)
  if (!all(is.na(meta$nas))) meta$nas_s <- as.numeric(scale(as.numeric(meta$nas)))
  meta
}
has <- function(meta, v) v %in% names(meta) && !all(is.na(meta[[v]]))

run_model <- function(name, cohort, counts, meta, terms, coef, genes) {
  keep <- rowSums(counts >= min_count) >= ncol(counts) / 2
  dds <- tryCatch(DESeq(DESeqDataSetFromMatrix(counts[keep, , drop = FALSE], meta,
                        as.formula(paste("~", paste(terms, collapse = " + ")))), quiet = TRUE), error = function(e) NULL)
  if (is.null(dds)) return(NULL)
  r <- as.data.table(as.data.frame(results(dds, name = coef)), keep.rownames = "gene")[gene %in% genes]
  r[, .(model = name, cohort = cohort, n = ncol(counts), gene, log2FC = log2FoldChange, lfcSE, pvalue,
        ci_low = log2FoldChange - 1.96 * lfcSE, ci_high = log2FoldChange + 1.96 * lfcSE)]
}

all_res <- list(); stage_tab <- list()
for (cn in names(cohorts)) {
  d <- read_cohort(cohorts[[cn]]); counts <- d$counts
  keep0 <- rowSums(counts >= min_count) >= ncol(counts) / 2
  vsd <- assay(varianceStabilizingTransformation(DESeqDataSetFromMatrix(counts[keep0, ], d$meta, ~ 1), blind = TRUE))
  meta <- add_covariates(d$meta, vsd)
  stage_tab[[cn]] <- data.table(cohort = cn, stage = paste0("F", 0:4), n = as.integer(table(factor(meta$stage, levels = 0:4))))
  bt <- c(if (has(meta, "sex")) "sex", if (has(meta, "age_s")) "age_s"); genes <- unique(c(panel, focus))
  sub <- function(sel) list(c = counts[, sel, drop = FALSE], m = meta[sel, , drop = FALSE])
  noF4 <- meta$stage <= 3; M5s <- meta$stage %in% c(0, 1, 3, 4)
  run <- function(nm, sel, terms, coef) { s <- sub(sel); run_model(nm, cn, s$c, s$m, terms, coef, genes) }
  all_res[[paste(cn, 1)]] <- run("M1_primary", rep(TRUE, nrow(meta)), c(bt, "stage_num"), "stage_num")
  all_res[[paste(cn, 2)]] <- run("M2_hepatocyte_adjusted", rep(TRUE, nrow(meta)), c(bt, "hep", "stage_num"), "stage_num")
  if (has(meta, "nas_s")) all_res[[paste(cn, 3)]] <- run("M3_NAS_adjusted", rep(TRUE, nrow(meta)), c(bt, "nas_s", "stage_num"), "stage_num")
  all_res[[paste(cn, 4)]] <- run("M4_primary_F4_excluded", noF4, c(bt, "stage_num"), "stage_num")
  mm <- meta; mm$adv <- as.numeric(mm$stage >= 3)
  { s <- list(c = counts[, M5s, drop = FALSE], m = mm[M5s, , drop = FALSE])
    all_res[[paste(cn, 5)]] <- run_model("M5_F3F4_vs_F0F1", cn, s$c, s$m, c(bt, "adv"), "adv", genes) }
  all_res[[paste(cn, 6)]] <- run("M6_hepatocyte_F4_excluded_NEW", noF4, c(bt, "hep", "stage_num"), "stage_num")
  if (has(meta, "nas_s")) all_res[[paste(cn, 7)]] <- run("M7_hepatocyte+NAS_F4_excluded_NEW", noF4, c(bt, "hep", "nas_s", "stage_num"), "stage_num")
}
stage_tab <- rbindlist(stage_tab); fwrite(stage_tab, file.path(out_dir, "1_stage_distribution.csv"))
cat("[1] samples per fibrosis stage:\n"); print(dcast(stage_tab, cohort ~ stage, value.var = "n"))

res <- rbindlist(all_res, fill = TRUE); fwrite(res, file.path(out_dir, "2_cohort_specific_estimates.csv"))

## ---- meta-analysis, Q and I2, BH across candidates ----------------------------------------------
both <- res[, .(n_cohorts = uniqueN(cohort)), by = .(model, gene)][n_cohorts >= 2]
meta_tab <- res[both, on = .(model, gene), nomatch = 0][, {
  b <- log2FC; se <- lfcSE; w <- 1 / se^2; bm <- sum(w * b) / sum(w); sem <- sqrt(1 / sum(w)); k <- .N
  Q <- sum(w * (b - bm)^2)
  list(meta_log2FC = bm, meta_ci_low = bm - 1.96 * sem, meta_ci_high = bm + 1.96 * sem, meta_p = 2 * pnorm(-abs(bm / sem)),
       same_direction = length(unique(sign(b))) == 1, all_cohort_p_below_0.05 = all(pvalue < 0.05),
       Q = Q, Q_p = pchisq(Q, k - 1, lower.tail = FALSE), I2_percent = max(0, (Q - (k - 1)) / Q) * 100,
       cohort_estimates = paste(sprintf("%s %.3f [%.3f, %.3f] P=%.2g", cohort, b, ci_low, ci_high, pvalue), collapse = " | "))
}, by = .(model, gene)]
meta_tab[gene %in% panel, meta_BH_across_candidates := p.adjust(meta_p, "BH"), by = model]
fwrite(meta_tab, file.path(out_dir, "3_meta_analysis_Q_I2.csv"))

key <- meta_tab[gene %in% focus, .(model, gene, meta_log2FC = round(meta_log2FC, 3), meta_p = signif(meta_p, 3),
        meta_BH = signif(meta_BH_across_candidates, 3), same_direction, all_cohort_p_below_0.05,
        Q_p = signif(Q_p, 3), I2_percent = round(I2_percent, 0))]
setorder(key, gene, model); fwrite(key, file.path(out_dir, "4_four_candidates_summary.csv"))
cat("\n[4] four candidates, all models (the new F4-excluded models are M6 and M7):\n"); print(key)
cat("\nHeterogeneity tests have little power with two cohorts; report Q and I2 but do not over-interpret them.\n")

pdf(file.path(out_dir, "4_four_candidates_forest.pdf"), width = 10, height = 6)
fp <- res[gene %in% focus]
fp[, model := factor(model, levels = rev(unique(model)))]
print(ggplot(fp, aes(log2FC, model, colour = cohort)) + geom_vline(xintercept = 0, lty = 2) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high), position = position_dodge(0.6), size = 0.3, orientation = "y") +
  facet_wrap(~ gene, scales = "free_x") + labs(x = "log2 fold change per stage (95% CI)", y = NULL) + theme_bw(base_size = 9))
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDone. Results are in", out_dir, "\n")
