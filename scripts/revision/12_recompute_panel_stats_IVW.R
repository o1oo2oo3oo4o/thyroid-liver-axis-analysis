#!/usr/bin/env Rscript
# =============================================================================
# 12_recompute_panel_stats_IVW.R
#
# Purpose
#   Recompute, from RESULT TABLES THAT ALREADY EXIST (scripts 04, 06 and 10), the three statistics that depend on
#   how many genes are in the prioritized panel. Nothing is re-estimated (no DESeq2, no Mendelian randomization):
#   only the panel-dependent summary statistics are recomputed.
#
#     [1] Fisher enrichment of the panel among THRB-sensitive genes      (Figure 4D, Table S9b)
#     [2] Benjamini-Hochberg adjustment of the biopsy meta-analysis P     (Tables S12b, S18)
#     [3] Proportion of panel genes with a TRb ChIP peak near the TSS     (Table S19)
#
# Panels
#   Baseline (all genes of config.R, 26 genes)    STAGE_DROP="" (or "none")
#       must reproduce the numbers of the submitted manuscript.
#   IVW-centred definition (default, 22 genes)    no STAGE_DROP, or STAGE_DROP="BTC,COL11A2,HSPA1A,LGALS4"
#       A protein-genetic candidate has to be supported by an IVW estimate (FDR < 0.05). The five genes whose only
#       support came from the weighted-median method are ADH4, BTC, COL11A2, HSPA1A and LGALS4. ADH4 stays in the
#       panel because it is also a positional (colocalization) candidate; the other four are removed.
#   Any other list can be given in STAGE_DROP (comma-separated); the output folder is then called custom_drop.
#
# Run from scripts/revision (the folder that contains config.R, data/ and results/):
#   Rscript 12_recompute_panel_stats_IVW.R                          # IVW-centred panel (22 genes)
#   STAGE_DROP="" Rscript 12_recompute_panel_stats_IVW.R            # baseline panel (26 genes)
#   (in R: Sys.setenv(STAGE_DROP = ""); source("12_recompute_panel_stats_IVW.R"))
#
# Inputs (written by scripts 04, 06 and 10 of this folder, or supplied; see data/README.md)
#   config.R                                          panel_genes (the 26 genes)
#   data/DEG_hepT3_raw.csv                            DESeq2 output of GSE270233 (WT, KO, interaction)
#   results/06_biopsy_F4/3_meta_analysis_Q_I2.csv     per-gene meta-analysis P for all 7 biopsy models
#   results/10_TRb_ChIP/1_gene_level_peak_annotation.csv, 2_group_summary.csv
# Output: results/12_panel_recompute_<tag>/*.csv   (tag: baseline_26, IVW_drop4 or custom_drop; a new folder,
#         existing results are never overwritten)
# =============================================================================
suppressPackageStartupMessages(library(data.table))

## ---- settings ----------------------------------------------------------------------------------------
# Default: the IVW-centred definition (four genes removed). STAGE_DROP overrides the list:
#   STAGE_DROP="BTC,COL11A2"   remove these genes      STAGE_DROP="" or "none"   keep all genes of config.R (baseline)
ivw_drop   <- c("BTC", "COL11A2", "HSPA1A", "LGALS4")
drop_env   <- Sys.getenv("STAGE_DROP", unset = NA)
drop_genes <- if (is.na(drop_env)) ivw_drop else if (tolower(trimws(drop_env)) == "none") character(0) else trimws(strsplit(drop_env, ",")[[1]])
drop_genes <- drop_genes[nzchar(drop_genes)]
tag <- if (identical(sort(drop_genes), sort(ivw_drop))) "IVW_drop4" else if (!length(drop_genes)) "baseline_26" else "custom_drop"

cfg <- new.env(); sys.source("config.R", envir = cfg)
panel_all <- cfg$panel_genes
panel     <- setdiff(panel_all, drop_genes)
stopifnot(all(drop_genes %in% panel_all))
out_dir <- file.path("results", paste0("12_panel_recompute_", tag)); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("panel: %d genes in config.R -> %d after dropping [%s]\n", length(panel_all), length(panel),
            paste(drop_genes, collapse = ", ")))

## ---- [1] Fisher enrichment (hepatocyte THRB-KO experiment) ----------------------------------------------
hep <- fread("data/DEG_hepT3_raw.csv"); hep[, gene := trimws(gene)]; hep <- hep[!duplicated(gene)]
hep[, tested_both := !is.na(padj_WT) & !is.na(padj_int)]            # genes tested in the WT and the interaction analysis
hep[, wt_resp  := tested_both & padj_WT < 0.05 & abs(log2FoldChange_WT) >= 1]
hep[, altered  := wt_resp & padj_int < 0.05 & abs(log2FoldChange_int) >= 1 & !is.na(log2FoldChange_int)]
N_tested <- sum(hep$tested_both); N_wt <- sum(hep$wt_resp); N_alt <- sum(hep$altered)
cat(sprintf("\n[1] tested in both analyses: %d | WT-responsive: %d | altered response in KO: %d (%.1f%% of WT-responsive)\n",
            N_tested, N_wt, N_alt, 100 * N_alt / N_wt))
cand <- hep[tested_both == TRUE & gene %in% panel]                  # panel genes expressed/tested in hepatocytes
cat(sprintf("    panel genes expressed (tested): %d of %d; not detected: %s\n", nrow(cand), length(panel),
            paste(setdiff(panel, cand$gene), collapse = ", ")))
fisher_one <- function(a, b, c, d) {                                # a,b = candidates altered/not ; c,d = others altered/not
  m <- matrix(c(a, b, c, d), 2, byrow = TRUE); f <- fisher.test(m, alternative = "greater")
  c(p = f$p.value, odds_ratio = unname(f$estimate))
}
n_c <- nrow(cand); a1 <- sum(cand$altered); c1 <- N_alt - a1
t1 <- fisher_one(a1, n_c - a1, c1, (N_tested - n_c) - c1)
n_c2 <- sum(cand$wt_resp); a2 <- sum(cand$altered)                 # test 2: only among WT-responsive genes
c2 <- N_alt - a2
t2 <- fisher_one(a2, n_c2 - a2, c2, (N_wt - n_c2) - c2)
fisher_tab <- data.table(
  test = c("Test 1: candidates vs all tested genes", "Test 2: candidates vs wild-type responsive genes"),
  n_candidates = c(n_c, n_c2), candidates_altered = c(a1, a2),
  others_altered = c(c1, c2), n_others = c(N_tested - n_c, N_wt - n_c2),
  expected_by_chance = c(n_c * N_alt / N_tested, n_c2 * N_alt / N_wt),
  one_sided_Fisher_P = c(t1["p"], t2["p"]), odds_ratio = c(t1["odds_ratio"], t2["odds_ratio"]))
show <- copy(fisher_tab); show[, `:=`(expected_by_chance = round(expected_by_chance, 3),
                                   one_sided_Fisher_P = signif(one_sided_Fisher_P, 3), odds_ratio = round(odds_ratio, 2))]
print(show)
fwrite(fisher_tab, file.path(out_dir, "1_Fisher_enrichment.csv"))
fwrite(cand[altered | wt_resp, .(gene, log2FoldChange_WT, padj_WT, log2FoldChange_int, padj_int, wt_resp, altered)],
       file.path(out_dir, "1_panel_genes_WT_responsive_or_altered.csv"))

## ---- [2] BH of the biopsy meta-analysis P across the tested candidates ----------------------------------
mt <- fread("results/06_biopsy_F4/3_meta_analysis_Q_I2.csv")
mt <- mt[gene %in% panel_all]                                       # same candidate set the script used (config panel)
mt[, BH_baseline_stored := meta_BH_across_candidates]
mt_new <- mt[gene %in% panel]
mt_new[, BH_recomputed := p.adjust(meta_p, "BH"), by = model]       # BH within each model across the tested candidates
mt_base <- copy(mt)[, BH_check := p.adjust(meta_p, "BH"), by = model]
chk <- max(abs(mt_base$BH_check - mt_base$BH_baseline_stored), na.rm = TRUE)
cat(sprintf("\n[2] check: BH recomputed on the full set vs the value stored by script 06: max |difference| = %.2e\n", chk))
stopifnot(chk < 1e-10)                                              # the stored values are reproduced exactly
cat(sprintf("    candidates tested in both cohorts per model: %s\n",
            paste(unique(mt_new[, .N, by = model]$N), collapse = "/")))
heads <- c("ABHD1", "F13B", "ADH4", "LGALS4")
bh_tab <- mt_new[, .(model, gene, meta_log2FC, meta_p, BH_stored_26panel = BH_baseline_stored,
                     BH_recomputed, same_direction, all_cohort_p_below_0.05)]
setorder(bh_tab, model, meta_p)
fwrite(bh_tab, file.path(out_dir, "2_biopsy_BH_all_models.csv"))
cat("    headline genes (stored BH with the 26-gene panel | recomputed with this panel):\n")
print(bh_tab[gene %in% heads & model %in% c("M2_hepatocyte_adjusted", "M6_hepatocyte_F4_excluded_NEW", "M1_primary"),
             .(model, gene, meta_p = signif(meta_p, 3), BH_26panel = signif(BH_stored_26panel, 3),
               BH_this_panel = signif(BH_recomputed, 3))][order(model, meta_p)])
cat("    number of candidates with BH < 0.05 per model (this panel):\n")
print(bh_tab[, .(n_tested = .N, n_BH_lt_0.05 = sum(BH_recomputed < 0.05), genes = paste(gene[BH_recomputed < 0.05], collapse = ",")), by = model])
dropped_in_biopsy <- intersect(drop_genes, mt$gene)
cat(sprintf("    dropped genes that were part of the 20 tested: %s\n", paste(dropped_in_biopsy, collapse = ", ")))

## ---- [3] ChIP proportions ------------------------------------------------------------------------------
pk <- fread("results/10_TRb_ChIP/1_gene_level_peak_annotation.csv")
gs <- fread("results/10_TRb_ChIP/2_group_summary.csv")
bgrow <- gs[grepl("background", group)]
pp <- pk[gene %in% panel]
stopifnot(nrow(pp) == length(panel))                                # every panel gene must have a TSS row
chip_tab <- data.table(group = c(paste0("panel (", length(panel), " genes)"), "genome-wide background (all genes in TSS table)"),
                       genes_with_TSS = c(nrow(pp), bgrow$genes_with_TSS),
                       with_promoter_peak = c(sum(pp$promoter_m2kb_p1kb > 0), bgrow$with_promoter_peak),
                       with_peak_10kb = c(sum(pp$proximal_10kb > 0), bgrow$with_peak_10kb),
                       with_peak_50kb = c(sum(pp$proximal_50kb > 0), bgrow$with_peak_50kb))
chip_tab[, `:=`(pct_promoter = round(100 * with_promoter_peak / genes_with_TSS, 1), pct_10kb = round(100 * with_peak_10kb / genes_with_TSS, 1),
                pct_50kb = round(100 * with_peak_50kb / genes_with_TSS, 1))]
cat("\n[3] TRb BioChIP peaks (Ayers 2014, HepG2) around the panel:\n"); print(chip_tab)
three <- pk[gene %in% c("GPCPD1", "ADH4", "F13B"), .(gene, promoter_m2kb_p1kb, proximal_10kb, proximal_50kb, nearest_peak_distance_bp)]
cat("    the three THRB-sensitive candidates (not panel-size dependent):\n"); print(three)
cat("    panel genes with any peak within 50 kb:\n"); print(pp[proximal_50kb > 0, .(gene, promoter_m2kb_p1kb, proximal_10kb, proximal_50kb)])
fwrite(chip_tab, file.path(out_dir, "3_ChIP_panel_vs_background.csv"))
fwrite(pp[order(gene)], file.path(out_dir, "3_ChIP_panel_gene_level.csv"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDone. Results:", out_dir, "\n")
