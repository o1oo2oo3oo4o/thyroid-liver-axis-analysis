#!/usr/bin/env Rscript
# =============================================================================
# 02_biopsy_stage_analysis.R
#
# Question (reviewer comment 6): are the fibrosis-stage associations of the
# candidate genes robust to fibrosis-stage heterogeneity?
#
# What it does for every cohort (GSE135251, GSE130970) from RAW COUNTS
#   M1  primary model            (numeric stage [+ sex + age where available])
#   M2  + hepatocyte score       (reproduces the hepatocyte-adjusted model, Table S12)
#   M3  + NAS
#   M4  + steatosis, lobular inflammation, ballooning (instead of the summed NAS)
#   M5  stage as a factor        (F2, F3, F4 versus F0-1; each step separately)
#   M6  likelihood-ratio test    (is a linear stage effect sufficient?  factor vs numeric)
#   M7  stage-stratified slopes  (early F0-F2 and advanced F3-F4)
# and then
#   * between-cohort heterogeneity (Cochran Q, I2) and fixed-effect meta-analysis
#   * figures: expression by stage, hepatocyte score by stage, step effects
#   * a reproduction check against the values reported in Table S12
#
# Input  : config.R  (paths, column names) + count matrices + metadata tables
# Output : results/biopsy_stage/*.csv, *.pdf
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(ggplot2)
})
source(Sys.getenv("STAGE_CFG", "config.R"))
out_dir <- file.path(out_root, "02_biopsy_stage"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
genes_of_interest <- unique(c(candidate_genes, control_genes))

## ---- helpers ------------------------------------------------------------------
read_cohort <- function(cfg) {
  cnt <- fread(cfg$counts); setnames(cnt, 1, "gene"); cnt[, gene := trimws(gene)]
  cnt <- cnt[!duplicated(gene)]
  m <- as.matrix(cnt[, -1]); rownames(m) <- cnt$gene; storage.mode(m) <- "integer"
  meta <- fread(cfg$meta); mc <- meta_cols
  meta[, sample := as.character(get(mc$sample))]
  meta[, stage := as.integer(get(mc$stage))]
  if (!is.na(mc$group) && mc$group %in% names(meta)) meta <- meta[get(mc$group) == mc$group_keep]
  for (nm in c("nas", "steatosis", "inflammation", "ballooning", "age", "sex")) {
    col <- mc[[nm]]
    meta[[nm]] <- if (!is.na(col) && col %in% names(meta)) meta[[col]] else NA
  }
  meta <- meta[!is.na(stage) & sample %in% colnames(m)]
  m <- m[, meta$sample, drop = FALSE]
  keep <- rowSums(m >= min_count) >= ncol(m) / 2
  list(counts = m[keep, ], meta = as.data.frame(meta))
}

prep_covariates <- function(meta, vsd_mat) {
  meta$stage_num <- as.numeric(meta$stage)
  meta$stage_grp <- factor(ifelse(meta$stage <= 1, "F0-1", paste0("F", meta$stage)),
                           levels = c("F0-1", "F2", "F3", "F4"))
  mk <- intersect(hepatocyte_markers, rownames(vsd_mat))
  if (length(mk) < 3) stop("fewer than 3 hepatocyte marker genes found")
  z <- t(scale(t(vsd_mat[mk, , drop = FALSE])))
  meta$hep <- as.numeric(scale(colMeans(z)))                 # hepatocyte score (z-scaled)
  if (!all(is.na(meta$age))) meta$age_s <- as.numeric(scale(as.numeric(meta$age)))
  if (!all(is.na(meta$sex))) meta$sex <- factor(meta$sex)
  for (v in c("nas", "steatosis", "inflammation", "ballooning"))
    if (!all(is.na(meta[[v]]))) meta[[paste0(v, "_s")]] <- as.numeric(scale(as.numeric(meta[[v]])))
  meta
}

has <- function(meta, v) v %in% names(meta) && !all(is.na(meta[[v]]))
base_terms <- function(meta) c(if (has(meta, "sex")) "sex", if (has(meta, "age_s")) "age_s")

# DESeq2 prints two harmless notes for these designs; they are muted on purpose:
#  (1) "numeric variables with integer values": stage_num (0-4) is numeric ON PURPOSE -> M1-M4, M6, M7
#      estimate one log2FC per fibrosis stage (linear trend). The factor version is M5.
#  (2) "levels of factors ... other characters": the stage_grp level "F0-1" contains a hyphen (label only).
# Neither note depends on how the input files are delimited, and neither changes any result.
# Only these two notes are muted; other DESeq2 messages (e.g. "rows did not converge") still show.
mute_design_notes <- function(expr) withCallingHandlers(expr, message = function(m) {
  if (grepl("numeric variables with integer values|levels of factors in the design", conditionMessage(m)))
    invokeRestart("muffleMessage")
})

fit_dds <- function(counts, meta, terms) {
  f <- as.formula(paste("~", paste(terms, collapse = " + ")))
  mute_design_notes(DESeq(DESeqDataSetFromMatrix(counts, meta, f), quiet = TRUE))
}
pick <- function(res, genes, model, cohort, coef_name) {
  r <- as.data.table(as.data.frame(res), keep.rownames = "gene")[gene %in% genes]
  r[, `:=`(model = model, cohort = cohort, coefficient = coef_name)]
  r[, .(cohort, model, coefficient, gene, log2FC = log2FoldChange, lfcSE, pvalue, padj, baseMean)]
}

## ---- run all models for one cohort ------------------------------------------------
run_cohort <- function(name, cfg) {
  cat("\n=== ", name, " ===\n")
  d <- read_cohort(cfg); counts <- d$counts
  cat(sprintf("samples: %d | genes kept: %d | stage table: %s\n", ncol(counts), nrow(counts),
              paste(names(table(d$meta$stage)), table(d$meta$stage), sep = "=", collapse = ", ")))
  vsd <- assay(varianceStabilizingTransformation(
    DESeqDataSetFromMatrix(counts, d$meta, ~ 1), blind = TRUE))
  meta <- prep_covariates(d$meta, vsd)
  bt <- base_terms(meta); out <- list()
  extract_num <- function(dds, model) pick(results(dds, name = "stage_num"), genes_of_interest, model, name, "stage_num")

  # M1, M2, M3, M4: numeric stage
  m1 <- fit_dds(counts, meta, c(bt, "stage_num"));                       out$M1 <- extract_num(m1, "M1_primary")
  m2 <- fit_dds(counts, meta, c(bt, "hep", "stage_num"));                out$M2 <- extract_num(m2, "M2_hepatocyte_adjusted")
  if (has(meta, "nas_s")) { m3 <- fit_dds(counts, meta, c(bt, "nas_s", "stage_num")); out$M3 <- extract_num(m3, "M3_NAS_adjusted") }
  comp <- paste0(c("steatosis", "inflammation", "ballooning"), "_s")
  if (all(sapply(comp, function(v) has(meta, v)))) {
    m4 <- fit_dds(counts, meta, c(bt, "hep", comp, "stage_num")); out$M4 <- extract_num(m4, "M4_hep+NAS_components")
  }

  # M5: stage as factor (F0-1 reference)  +  M6: LRT factor versus numeric stage
  grp_ok <- table(meta$stage_grp) >= 3
  if (all(grp_ok)) {
    m5 <- fit_dds(counts, meta, c(bt, "hep", "stage_grp"))
    for (cn in setdiff(resultsNames(m5), c("Intercept", grep("^(sex|age_s|hep)", resultsNames(m5), value = TRUE))))
      out[[cn]] <- pick(results(m5, name = cn), genes_of_interest, "M5_stage_factor", name, cn)
    # LRT: does the factor model fit better than the linear stage effect?
    m6 <- tryCatch({
      dds6 <- mute_design_notes(DESeqDataSetFromMatrix(counts, meta, as.formula(paste("~", paste(c(bt, "hep", "stage_num", "stage_grp"), collapse = " + ")))))
      dds6 <- mute_design_notes(DESeq(dds6, test = "LRT", reduced = as.formula(paste("~", paste(c(bt, "hep", "stage_num"), collapse = " + "))), quiet = TRUE))
      as.data.table(as.data.frame(results(dds6)), keep.rownames = "gene")[gene %in% genes_of_interest,
        .(cohort = name, model = "M6_LRT_nonlinearity", coefficient = "stage_grp beyond linear",
          gene, log2FC = NA_real_, lfcSE = NA_real_, pvalue, padj, baseMean)]
    }, error = function(e) { message("LRT skipped: ", conditionMessage(e)); NULL })
    out$M6 <- m6
  } else message("M5/M6 skipped: a stage group has fewer than 3 samples")

  # M7: stage-stratified slopes
  for (st in list(early = 0:2, advanced = 3:4)) {
    nm <- if (identical(st, 0:2)) "early_F0-F2" else "advanced_F3-F4"
    sel <- meta$stage %in% st
    if (sum(sel) >= min_n_stratum && length(unique(meta$stage[sel])) >= 2) {
      mm <- meta[sel, ]; cc <- counts[, sel, drop = FALSE]
      cc <- cc[rowSums(cc >= min_count) >= ncol(cc) / 2, ]
      m7 <- tryCatch(fit_dds(cc, mm, c(bt, "hep", "stage_num")), error = function(e) NULL)
      if (!is.null(m7)) out[[paste0("M7_", nm)]] <- pick(results(m7, name = "stage_num"), genes_of_interest,
                                                         paste0("M7_stratum_", nm), name, "stage_num")
    } else message("M7 stratum ", nm, " skipped (n = ", sum(sel), ")")
  }

  # data for plots
  expr <- data.table(sample = colnames(vsd), t(vsd[intersect(genes_of_interest, rownames(vsd)), , drop = FALSE]))
  pl <- merge(expr, as.data.table(meta[, c("sample", "stage", "stage_grp", "hep")]), by = "sample")
  pl[, cohort := name]
  list(res = rbindlist(out, fill = TRUE), plot = pl)
}

all_runs <- lapply(names(cohorts), function(n) run_cohort(n, cohorts[[n]]))
names(all_runs) <- names(cohorts)
res <- rbindlist(lapply(all_runs, `[[`, "res"), fill = TRUE)
fwrite(res, file.path(out_dir, "all_models_candidates_and_controls.csv"))

## ---- cross-cohort meta-analysis with heterogeneity ------------------------------
num_models <- res[coefficient == "stage_num" | grepl("^M7", model)]
meta_tab <- num_models[, {
  ok <- is.finite(log2FC) & is.finite(lfcSE) & lfcSE > 0
  if (sum(ok) < 2) NULL else {
    b <- log2FC[ok]; se <- lfcSE[ok]; w <- 1 / se^2
    bm <- sum(w * b) / sum(w); sem <- sqrt(1 / sum(w)); Q <- sum(w * (b - bm)^2); k <- sum(ok)
    list(n_cohorts = k, meta_log2FC = bm, meta_se = sem,
         meta_p = 2 * pnorm(-abs(bm / sem)), same_direction = length(unique(sign(b))) == 1,
         all_cohort_p_below_0.05 = all(pvalue[ok] < 0.05),
         Q = Q, het_p = pchisq(Q, k - 1, lower.tail = FALSE), I2 = max(0, (Q - (k - 1)) / Q) * 100,
         cohort_log2FC = paste(sprintf("%s %.3f", cohort[ok], b), collapse = " | "),
         cohort_P = paste(sprintf("%s %.2g", cohort[ok], pvalue[ok]), collapse = " | "))
  }
}, by = .(model, gene)]
fwrite(meta_tab, file.path(out_dir, "meta_analysis_with_heterogeneity.csv"))
cat("\n=== robustness summary for the candidates (same direction & P < 0.05 in all cohorts) ===\n")
print(dcast(meta_tab[gene %in% candidate_genes & !grepl("^M7", model)], gene ~ model,
            value.var = "all_cohort_p_below_0.05"))

## ---- reproduction check against Table S12 ---------------------------------------------
expected <- data.table(
  gene = rep(c("ABHD1", "F13B", "ADH4", "LGALS4"), each = 2),
  cohort = rep(c("GSE130970", "GSE135251"), 4),
  M1_expected = c(-0.129, -0.061, -0.055, -0.090, -0.042, -0.092, 0.245, 0.248),
  M2_expected = c(-0.139, -0.118, -0.075, -0.074, -0.081, -0.077, 0.217, 0.302))
chk <- merge(expected,
  res[model == "M1_primary", .(gene, cohort, M1_obs = log2FC)], by = c("gene", "cohort"), all.x = TRUE)
chk <- merge(chk, res[model == "M2_hepatocyte_adjusted", .(gene, cohort, M2_obs = log2FC)], by = c("gene", "cohort"), all.x = TRUE)
chk[, `:=`(M1_diff = M1_obs - M1_expected, M2_diff = M2_obs - M2_expected)]
fwrite(chk, file.path(out_dir, "reproduction_check_vs_TableS12.csv"))
cat("\n=== reproduction check (differences should be ~0 when the same samples and filters are used) ===\n")
print(chk[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

## ---- figures --------------------------------------------------------------------------
plot_dt <- rbindlist(lapply(all_runs, `[[`, "plot"), fill = TRUE)
pdf(file.path(out_dir, "stage_figures.pdf"), width = 10, height = 7)
# (a) hepatocyte score by stage: the composition shift behind the adjustment
print(ggplot(unique(plot_dt[, .(sample, cohort, stage_grp, hep)]), aes(stage_grp, hep)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") + geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
  facet_wrap(~ cohort) + labs(x = "fibrosis stage group", y = "hepatocyte marker score (z)",
  title = "Tissue composition changes with fibrosis stage") + theme_bw())
# (b) expression of the candidates and controls by stage group
long <- melt(plot_dt, id.vars = c("sample", "cohort", "stage", "stage_grp", "hep"),
             measure.vars = intersect(candidate_genes, names(plot_dt)), variable.name = "gene", value.name = "vst")
print(ggplot(long, aes(stage_grp, vst)) + geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.6) + facet_grid(gene ~ cohort, scales = "free_y") +
  labs(x = "fibrosis stage group", y = "VST expression", title = "Candidate genes by fibrosis stage") + theme_bw())
# (c) step effects: F2, F3, F4 versus F0-1
st <- res[model == "M5_stage_factor" & gene %in% candidate_genes]
if (nrow(st) > 0) {
  st[, coefficient := sub("stage_grp_", "", coefficient)]
  print(ggplot(st, aes(log2FC, coefficient, colour = cohort)) + geom_vline(xintercept = 0, lty = 2) +
    geom_pointrange(aes(xmin = log2FC - 1.96 * lfcSE, xmax = log2FC + 1.96 * lfcSE),
                    position = position_dodge(0.5), size = 0.35, orientation = "y") +
    facet_wrap(~ gene, scales = "free_x") + labs(x = "log2 fold change vs F0-1 (95% CI)", y = NULL,
    title = "Stage-resolved effects (hepatocyte-adjusted)") + theme_bw())
}
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDone. Results are in", out_dir, "\n")
