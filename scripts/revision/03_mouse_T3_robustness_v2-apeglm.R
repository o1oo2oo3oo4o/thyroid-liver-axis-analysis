#!/usr/bin/env Rscript
# =============================================================================
# 03_mouse_T3_robustness_v2-apeglm.R
#
# v2-apeglm: identical to 03_mouse_T3_robustness.R EXCEPT that a gene is called "responsive"
# (adj P < 0.05 and |log2FC| >= 1) using the apeglm-SHRUNKEN log2FC, exactly as in the manuscript's DEG
# table (1,199 genes with adj P < 0.05; 541 of them with |shrunken log2FC| >= 1). Version 1 used the
# unshrunken log2FC and called 779 genes, so its Fisher/GO results were not comparable to the manuscript.
# P values and adj P are the DESeq2 Wald values (identical for both log2FC versions); the cameraPR score
# keeps the sign of the unshrunken log2FC, so cameraPR results are the same as in version 1.
# Needs the Bioconductor package apeglm.
#
# Question (reviewer comment 4): with 4 vs 4 mice (GSE85793), how robust are the
# pathway-level conclusions (glutathione set decreased, microsomal P450 set increased)?
#
# What it does from RAW COUNTS
#   1  re-runs the DESeq2 T3-vs-vehicle analysis and the set tests of the manuscript
#      (reproduction check: 541 responsive genes; glutathione Fisher P = 2.07e-4,
#       cameraPR P = 1.8e-6 at rho = 0.01)
#   2  ESTIMATES the inter-gene correlation from the data (instead of assuming
#      0.01 / 0.05 / 0.10) and repeats cameraPR; also limma::camera on voom values
#      with the correlation estimated per set (inter.gene.cor = NA)
#   3  leave-one-out resampling: every analysis is repeated after removing
#      1 mouse (8 subsets) or 1 T3 + 1 vehicle mouse (16 subsets of 3 vs 3)
#   4  an alternative differential-expression engine (limma-voom) to show the
#      conclusions do not depend on DESeq2
# Input  : config.R, count matrix (gene symbol x sample), metadata (sample, group)
# Output : results/mouse_robustness/*.csv, *.pdf
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(limma); library(ggplot2)
})
if (!requireNamespace("apeglm", quietly = TRUE))
  stop("Package 'apeglm' is required: BiocManager::install('apeglm'). Not falling back to unshrunken log2FC on purpose.")
source(Sys.getenv("STAGE_CFG", "config.R"))
out_dir <- file.path(out_root, "03_mouse_robustness_v2-apeglm"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rhos <- c(0.01, 0.05, 0.10)

## ---- data ------------------------------------------------------------------------------
cnt <- fread(mouse_counts_file); setnames(cnt, 1, "gene"); cnt[, gene := trimws(gene)]
cnt <- cnt[!duplicated(gene)]
counts <- as.matrix(cnt[, -1]); rownames(counts) <- cnt$gene; storage.mode(counts) <- "integer"
meta <- as.data.frame(fread(mouse_meta_file)); rownames(meta) <- meta$sample
counts <- counts[, rownames(meta), drop = FALSE]
meta$group <- factor(meta[[mouse_group_col]], levels = c(mouse_ctrl_label, mouse_treat_label))
stopifnot(!any(is.na(meta$group)))
keep <- rowSums(counts >= mouse_filter_min_count) >= ncol(counts) / 2
counts <- counts[keep, ]
cat(sprintf("samples: %d (%s) | genes kept: %d\n", ncol(counts),
            paste(names(table(meta$group)), table(meta$group), sep = "=", collapse = ", "), nrow(counts)))

sets <- fread(mouse_sets_file); sets[, mouse_symbol := trimws(mouse_symbol)]
sets_list <- split(sets$mouse_symbol, sets$set)

## ---- core function: DESeq2 + set tests for any subset of samples ---------------------------
analyse <- function(samples, label) {
  cc <- counts[, samples, drop = FALSE]; mm <- meta[samples, , drop = FALSE]
  dds <- DESeq(DESeqDataSetFromMatrix(cc, mm, ~ group), quiet = TRUE)
  coef_name <- paste0("group_", mouse_treat_label, "_vs_", mouse_ctrl_label)
  stopifnot(coef_name %in% resultsNames(dds))
  r_raw <- results(dds, name = coef_name)                                            # Wald test, unshrunken log2FC
  r_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)           # apeglm-shrunken log2FC (as in the manuscript)
  stopifnot(identical(rownames(r_raw), rownames(r_shr)))
  r <- as.data.table(as.data.frame(r_shr), keep.rownames = "gene")                  # log2FoldChange = SHRUNKEN; pvalue/padj = Wald
  r[, lfc_unshrunk := r_raw$log2FoldChange]
  r <- r[!is.na(padj)]
  r[, resp := padj < mouse_padj_cut & abs(log2FoldChange) >= mouse_lfc_cut]          # shrunken log2FC, as in the manuscript
  r[, score := -log10(pvalue) * sign(lfc_unshrunk)]            # signed score used in Methods 2.7 (sign of the unshrunken estimate)
  sc <- setNames(r$score, r$gene)
  rows <- lapply(names(sets_list), function(s) {
    mem <- sets_list[[s]]; ins <- r$gene %in% mem
    a <- sum(r$resp & ins); b <- sum(!r$resp & ins); c <- sum(r$resp & !ins); d <- sum(!r$resp & !ins)
    o <- list(subset = label, set = s, n_tested = sum(ins), n_responsive_in_set = a,
              n_responsive_total = sum(r$resp),
              n_set_down = sum(r$resp & ins & r$log2FoldChange < 0), n_set_up = sum(r$resp & ins & r$log2FoldChange > 0),
              fisher_p = fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE), alternative = "greater")$p.value)
    idx <- which(names(sc) %in% mem)
    for (rho in rhos) {
      cm <- if (length(idx) >= 3) cameraPR(sc, list(set = idx), inter.gene.cor = rho, use.ranks = FALSE) else NULL
      o[[sprintf("cameraPR_p_rho%.2f", rho)]] <- if (is.null(cm)) NA_real_ else cm$PValue
      o[[sprintf("cameraPR_dir_rho%.2f", rho)]] <- if (is.null(cm)) NA_character_ else as.character(cm$Direction)
    }
    as.data.table(o)
  })
  list(table = rbindlist(rows), res = r, dds = dds)
}

## ---- 1. full-data reproduction -------------------------------------------------------------
all_samples <- rownames(meta)
full <- analyse(all_samples, "full")
cat(sprintf("\n[1] responsive genes (adj P < %.2f and apeglm |log2FC| >= %g): %d  (manuscript: 541 of 11,819 tested; here %d tested)\n",
            mouse_padj_cut, mouse_lfc_cut, sum(full$res$resp), nrow(full$res)))
cat(sprintf("    (same data with the unshrunken log2FC would give %d genes; this is what version 1 reported)\n",
            sum(full$res$padj < mouse_padj_cut & abs(full$res$lfc_unshrunk) >= mouse_lfc_cut)))
if (file.exists(mouse_deg_file)) {                      # compare with the manuscript's own DEG table
  ms <- fread(mouse_deg_file); ms[, gene := trimws(gene)]
  ms_resp <- ms[responsive == TRUE, gene]; my_resp <- full$res[resp == TRUE, gene]
  cat(sprintf("    manuscript table: %d responsive genes | overlap with this run: %d | only here: %d | only in manuscript: %d\n",
              length(ms_resp), length(intersect(ms_resp, my_resp)), length(setdiff(my_resp, ms_resp)), length(setdiff(ms_resp, my_resp))))
}
print(full$table[, .(set, n_tested, n_responsive_in_set, n_set_down, n_set_up,
                     fisher_p = signif(fisher_p, 3), cameraPR_p_rho0.01 = signif(cameraPR_p_rho0.01, 3),
                     cameraPR_p_rho0.10 = signif(cameraPR_p_rho0.10, 3))])
fwrite(full$table, file.path(out_dir, "1_full_data_set_tests.csv"))

## ---- 2. estimated inter-gene correlation ------------------------------------------------------
vst_mat <- assay(varianceStabilizingTransformation(full$dds, blind = FALSE))
resid <- vst_mat - t(apply(vst_mat, 1, function(x) ave(x, meta$group)))   # remove group means
mean_cor <- function(genes) {
  g <- intersect(genes, rownames(resid)); if (length(g) < 3) return(NA_real_)
  cm <- cor(t(resid[g, , drop = FALSE])); mean(cm[upper.tri(cm)], na.rm = TRUE)
}
set.seed(1)
rho_tab <- rbindlist(lapply(names(sets_list), function(s) {
  g <- intersect(sets_list[[s]], rownames(resid))
  rnd <- replicate(1000, mean_cor(sample(rownames(resid), length(g))))
  data.table(set = s, n = length(g), rho_hat = mean_cor(g),
             random_set_rho_mean = mean(rnd, na.rm = TRUE),
             percentile_among_random_sets = mean(rnd <= mean_cor(g), na.rm = TRUE) * 100)
}))
fwrite(rho_tab, file.path(out_dir, "2a_estimated_intergene_correlation.csv"))
cat("\n[2a] mean pairwise correlation of residuals inside each set (6 residual df; noisy):\n"); print(rho_tab)

sc_full <- setNames(full$res$score, full$res$gene)
cam_est <- rbindlist(lapply(names(sets_list), function(s) {
  idx <- which(names(sc_full) %in% sets_list[[s]])
  rho_s <- max(0, min(0.5, rho_tab[set == s, rho_hat]))                  # clamp to a sensible range
  cm <- cameraPR(sc_full, list(set = idx), inter.gene.cor = rho_s, use.ranks = FALSE)
  data.table(set = s, rho_used = rho_s, cameraPR_dir = as.character(cm$Direction), cameraPR_p = cm$PValue)
}))
# limma::camera on voom values: the correlation is estimated per set (inter.gene.cor = NA)
y <- edgeR::calcNormFactors(edgeR::DGEList(counts))
design <- model.matrix(~ group, meta)
v <- voom(y, design)
cam_voom <- rbindlist(lapply(names(sets_list), function(s) {
  idx <- which(rownames(v) %in% sets_list[[s]])
  if (length(idx) < 3) return(NULL)
  cm <- camera(v, list(set = idx), design, contrast = 2, inter.gene.cor = NA)
  data.table(set = s, camera_voom_dir = as.character(cm$Direction), camera_voom_p_estimated_cor = cm$PValue)
}))
est <- merge(cam_est, cam_voom, by = "set", all = TRUE)
fwrite(est, file.path(out_dir, "2b_cameraPR_with_estimated_correlation.csv"))
cat("\n[2b] set tests with data-estimated inter-gene correlation:\n"); print(est)

## ---- 3. leave-one-out resampling --------------------------------------------------------------
t3  <- all_samples[meta$group == mouse_treat_label]; ctl <- all_samples[meta$group == mouse_ctrl_label]
subsets <- c(
  lapply(all_samples, function(x) list(label = paste0("drop_", x), s = setdiff(all_samples, x))),
  unlist(lapply(t3, function(a) lapply(ctl, function(b)
    list(label = paste0("drop_", a, "+", b), s = setdiff(all_samples, c(a, b))))), recursive = FALSE))
cat(sprintf("\n[3] leave-one-out: %d subsets (%d single-sample, %d one-per-group)\n", length(subsets),
            length(all_samples), length(t3) * length(ctl)))
loo <- rbindlist(lapply(subsets, function(z) {
  a <- analyse(z$s, z$label)
  t <- a$table
  t[, `:=`(n_t3 = sum(meta[z$s, "group"] == mouse_treat_label), n_ctl = sum(meta[z$s, "group"] == mouse_ctrl_label),
           jaccard_with_full = length(intersect(a$res[resp == TRUE, gene], full$res[resp == TRUE, gene])) /
                               length(union(a$res[resp == TRUE, gene], full$res[resp == TRUE, gene])))]
  g <- a$res[gene %in% c("Gclc", "Por"), .(gene, log2FoldChange, padj)]
  for (i in seq_len(nrow(g))) { set(t, j = paste0(g$gene[i], "_log2FC"), value = g$log2FoldChange[i])
                                set(t, j = paste0(g$gene[i], "_padj"),   value = g$padj[i]) }
  t
}), fill = TRUE)
fwrite(loo, file.path(out_dir, "3_leave_one_out_all_subsets.csv"))
loo_sum <- loo[, .(n_subsets = .N,
                   frac_fisher_p_below_0.05 = mean(fisher_p < 0.05),
                   frac_fisher_bonferroni_below_0.05 = mean(fisher_p * length(sets_list) < 0.05),
                   median_fisher_p = median(fisher_p), max_fisher_p = max(fisher_p),
                   frac_cameraPR_rho0.01_below_0.05 = mean(cameraPR_p_rho0.01 < 0.05, na.rm = TRUE),
                   frac_cameraPR_rho0.10_below_0.05 = mean(cameraPR_p_rho0.10 < 0.05, na.rm = TRUE),
                   median_jaccard_with_full_gene_list = median(jaccard_with_full),
                   frac_set_direction_unchanged = mean(sign(n_set_up - n_set_down) == sign(full$table[set == .BY$set, n_set_up - n_set_down]))),
               by = set]
fwrite(loo_sum, file.path(out_dir, "3_leave_one_out_summary.csv"))
cat("\n[3] summary over the leave-one-out subsets:\n"); print(loo_sum)

## ---- 4. alternative engine: limma-voom ---------------------------------------------------------
fit <- eBayes(lmFit(v, design))
tt <- as.data.table(topTable(fit, coef = 2, number = Inf, sort.by = "none"), keep.rownames = "gene")
tt[, resp := adj.P.Val < mouse_padj_cut & abs(logFC) >= mouse_lfc_cut]
cmp <- merge(tt[, .(gene, logFC, resp_voom = resp)], full$res[, .(gene, log2FoldChange = lfc_unshrunk, resp_deseq2 = resp)], by = "gene")  # log2FC compared on the unshrunken scale
alt <- rbindlist(lapply(names(sets_list), function(s) {
  ins <- tt$gene %in% sets_list[[s]]
  a <- sum(tt$resp & ins); b <- sum(!tt$resp & ins); c <- sum(tt$resp & !ins); d <- sum(!tt$resp & !ins)
  data.table(set = s, n_responsive_voom = a, n_responsive_deseq2 = full$table[set == s, n_responsive_in_set],
             fisher_p_voom = fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE), alternative = "greater")$p.value)
}))
# note: voom responsive genes use the limma log2FC (unshrunken); DESeq2 responsive genes use the apeglm log2FC
fwrite(alt, file.path(out_dir, "4_limma_voom_set_tests.csv"))
cat(sprintf("\n[4] limma-voom: %d responsive genes (DESeq2 %d); log2FC correlation with DESeq2 = %.3f; overlap = %d genes\n",
            sum(tt$resp), sum(full$res$resp), cor(cmp$logFC, cmp$log2FoldChange), sum(cmp$resp_voom & cmp$resp_deseq2)))
print(alt)

## ---- figures ---------------------------------------------------------------------------------------
pdf(file.path(out_dir, "mouse_robustness_figures.pdf"), width = 9, height = 6)
pl <- melt(loo[, .(set, subset, fisher = -log10(fisher_p), cameraPR_rho0.01 = -log10(cameraPR_p_rho0.01),
                   cameraPR_rho0.10 = -log10(cameraPR_p_rho0.10))],
           id.vars = c("set", "subset"), variable.name = "test", value.name = "neglog10p")
print(ggplot(pl, aes(set, neglog10p)) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6, colour = "#3973A5") +
  geom_hline(yintercept = -log10(0.05), lty = 2) + facet_wrap(~ test, ncol = 1) +
  labs(x = NULL, y = "-log10 P across leave-one-out subsets", title = "Stability of set-level results (dashed: P = 0.05)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1)))
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDone. Results are in", out_dir, "\n")
