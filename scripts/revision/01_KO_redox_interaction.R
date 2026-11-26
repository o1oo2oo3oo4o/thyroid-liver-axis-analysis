#!/usr/bin/env Rscript
# =============================================================================
# 01_KO_redox_interaction.R
#
# Question (reviewer comment 5): is the T3 response of redox / glutathione genes
# attenuated by THRB loss in human iPSC-derived hepatocytes (GSE270233)?
#
# Input  : DEG_hepT3_raw.csv            gene-level DESeq2 output of the
#                                       genotype x treatment model (WT, KO,
#                                       interaction, baseline contrast)
#          redox_gene_sets_human.csv    columns: set, human_symbol
# Output : results/KO_redox/*.csv, *.pdf
#
# Sign convention used throughout
#   z_*  = log2FoldChange / lfcSE (unshrunken Wald statistic)
#   att_z = -sign(log2FC_WT) * z_int
#           > 0  : the KO response is smaller than (or opposite to) the WT
#                  response, i.e. the response is ATTENUATED by THRB loss
#           < 0  : the KO response is larger than the WT response
#
# IMPORTANT - selection bias (winner's curse)
# Attenuation is only interpretable where the WT cells respond, but choosing
# genes on the WT result makes the interaction look attenuated for almost every
# gene (WT is over-estimated, KO - WT is therefore negative on average).
# In GSE270233 (n = 3 per group) ~98% of the WT-nominal genes look attenuated.
# Consequently
#   * section B compares set members with ALL OTHER genes in the SAME selected
#     universe (the bias cancels) - never read the raw proportion of attenuated genes;
#   * section B4 avoids the bias altogether by using the T3 direction of the
#     INDEPENDENT mouse experiment (GSE85793) instead of the human WT result.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(ggplot2)
})

## ---- configuration ----------------------------------------------------------
in_dir   <- "data"                                # folder with the input files
hep_file <- file.path(in_dir, "DEG_hepT3_raw.csv")
set_file <- file.path(in_dir, "redox_gene_sets_human.csv")
mouse_file <- file.path(in_dir, "mouse_T3_results.csv")   # optional (section B4)
out_dir  <- file.path("results", "01_KO_redox")
rhos     <- c(0.01, 0.05, 0.10)   # assumed inter-gene correlations (as in the manuscript)
padj_cut <- 0.05                  # response criteria used in the manuscript
lfc_cut  <- 1
wt_nominal_p <- 0.05              # "WT responds" for the attenuation analysis
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- read data ---------------------------------------------------------------
hep <- fread(hep_file)
hep[, gene := trimws(gene)]
hep <- hep[!duplicated(gene)]
sets <- fread(set_file)
sets[, human_symbol := trimws(human_symbol)]
sets_list <- split(sets$human_symbol, sets$set)

hep[, `:=`(
  z_wt  = log2FoldChange_WT  / lfcSE_WT,
  z_ko  = log2FoldChange_KO  / lfcSE_KO,
  z_int = log2FoldChange_int / lfcSE_int
)]
# signed score used for the set-level WT tests in the manuscript (-log10 P x sign of log2FC)
hep[, ss_wt := -log10(pvalue_WT) * sign(log2FoldChange_WT)]
hep[, att_z := -sign(log2FoldChange_WT) * z_int]
hep[, tested_wt  := !is.na(padj_WT)]
hep[, tested_int := !is.na(padj_int)]
hep[, resp_wt  := tested_wt  & padj_WT  < padj_cut & abs(log2FoldChange_WT)  >= lfc_cut]
hep[, int_crit := tested_int & padj_int < padj_cut & abs(log2FoldChange_int) >= lfc_cut]
hep[, altered  := resp_wt & int_crit]
hep[, wt_nominal := tested_wt & tested_int & pvalue_WT < wt_nominal_p]

cat(sprintf("Genes tested in WT: %d | WT-responsive: %d | altered (interaction criteria): %d\n",
            sum(hep$tested_wt), sum(hep$resp_wt), sum(hep$altered)))
cat(sprintf("WT-nominal universe for the attenuation analysis: %d genes\n", sum(hep$wt_nominal)))

## ---- helpers --------------------------------------------------------------------
camera_set <- function(stat_vec, members, rho) {
  idx <- which(names(stat_vec) %in% members)
  if (length(idx) < 3) return(list(n = length(idx), dir = NA_character_, p = NA_real_))
  r <- cameraPR(stat_vec, list(set = idx), inter.gene.cor = rho, use.ranks = FALSE)
  list(n = r$NGenes, dir = as.character(r$Direction), p = r$PValue)
}
named <- function(dt, col) { v <- dt[[col]]; names(v) <- dt$gene; v[is.finite(v)] }

## =============================================================================
## A. Reproduce the WT result for every set (manuscript: glutathione Fisher
##    P = 0.75, cameraPR P = 0.85 in human hepatocytes)
## =============================================================================
uni_wt <- hep[tested_wt == TRUE]
stat_wt <- named(uni_wt, "ss_wt")     # signed score, as in Methods 2.7 / 2.10
resA <- rbindlist(lapply(names(sets_list), function(s) {
  mem <- sets_list[[s]]
  in_set  <- uni_wt$gene %in% mem
  a <- sum(uni_wt$resp_wt &  in_set); b <- sum(!uni_wt$resp_wt &  in_set)
  c <- sum(uni_wt$resp_wt & !in_set); d <- sum(!uni_wt$resp_wt & !in_set)
  fp <- fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE), alternative = "greater")$p.value
  out <- list(set = s, n_tested = sum(in_set), n_responsive = a, fisher_p = fp)
  for (rho in rhos) {
    cm <- camera_set(stat_wt, mem, rho)
    out[[sprintf("cameraPR_dir_rho%.2f", rho)]] <- cm$dir
    out[[sprintf("cameraPR_p_rho%.2f",   rho)]] <- cm$p
  }
  as.data.table(out)
}))
fwrite(resA, file.path(out_dir, "A_WT_set_level_reproduction.csv"))
cat("\n[A] WT set-level results (compare with Table S11 / Table 6):\n"); print(resA)

## =============================================================================
## B. THRB dependence of the set response (the new analysis)
## =============================================================================
uni_att <- hep[wt_nominal == TRUE & is.finite(att_z)]
stat_att <- named(uni_att, "att_z")
resB <- rbindlist(lapply(names(sets_list), function(s) {
  mem <- sets_list[[s]]
  in_set <- uni_att$gene %in% mem
  n_set <- sum(in_set)
  if (n_set < 3) return(data.table(set = s, n_wt_nominal = n_set))
  out <- list(set = s, n_wt_nominal = n_set,
              median_att_z_set = median(uni_att$att_z[in_set]),
              median_att_z_other = median(uni_att$att_z[!in_set]))
  # B1: cameraPR on the attenuation score (Direction "Up" = more attenuated than other genes)
  for (rho in rhos) {
    cm <- camera_set(stat_att, mem, rho)
    out[[sprintf("B1_cameraPR_dir_rho%.2f", rho)]] <- cm$dir
    out[[sprintf("B1_cameraPR_p_rho%.2f",   rho)]] <- cm$p
  }
  # B3: Wilcoxon (one-sided) on att_z, set versus all other genes of the same universe
  out$B3_wilcox_p <- wilcox.test(uni_att$att_z[in_set], uni_att$att_z[!in_set],
                                 alternative = "greater")$p.value
  as.data.table(out)
}), fill = TRUE)
fwrite(resB, file.path(out_dir, "B_attenuation_set_level.csv"))
cat("\n[B] Is the set MORE attenuated in KO than other WT-responsive genes? (small P = yes)\n")
cat("    (the overall level of attenuation is inflated by selection bias; only the comparison counts)\n")
print(resB[, .(set, n_wt_nominal, median_att_z_set, median_att_z_other,
               B1_p_rho0.01 = B1_cameraPR_p_rho0.01, B1_p_rho0.10 = B1_cameraPR_p_rho0.10,
               B3_wilcox_p)])


## ---- B4. selection-free test using the independent mouse experiment --------------
## Hypothesis: if the conserved T3 response of a set is mediated by THRB, genes that
## change in mouse liver should have an interaction of OPPOSITE sign in human cells
## (e.g. mouse down -> positive interaction). The human WT result is not used for
## selection, so there is no winner's-curse bias.
if (file.exists(mouse_file)) {
  mo <- fread(mouse_file)
  mo <- mo[homology_types == "ortholog_one2one" & n_human == 1 & is.finite(log2FC)]
  mo[, gene := trimws(human_genes)]
  mo <- mo[!duplicated(gene), .(gene, mouse_log2FC = log2FC, mouse_responsive = responsive)]
  b4 <- merge(hep[is.finite(z_int) & tested_int == TRUE, .(gene, z_int)], mo, by = "gene")
  b4[, aligned_z := -sign(mouse_log2FC) * z_int]          # > 0 : attenuated in KO, given the mouse direction
  stat_b4 <- named(b4, "aligned_z")
  cat(sprintf("\n[B4] genes with a one-to-one mouse ortholog tested in human KO: %d\n", nrow(b4)))
  resB4 <- rbindlist(lapply(names(sets_list), function(s) {
    mem <- sets_list[[s]]; in_set <- b4$gene %in% mem
    if (sum(in_set) < 3) return(data.table(set = s, n = sum(in_set)))
    out <- list(set = s, n_all = sum(in_set), mean_aligned_z_set = mean(b4$aligned_z[in_set]),
                mean_aligned_z_other = mean(b4$aligned_z[!in_set]))
    for (rho in rhos) { cm <- camera_set(stat_b4, mem, rho)
      out[[sprintf("cameraPR_dir_rho%.2f", rho)]] <- cm$dir
      out[[sprintf("cameraPR_p_rho%.2f",   rho)]] <- cm$p }
    # restricted to genes responsive in mouse (independent selection)
    r <- b4[mouse_responsive == TRUE]; ins <- r$gene %in% mem
    out$n_mouse_responsive_in_set <- sum(ins)
    out$mouse_resp_wilcox_p <- if (sum(ins) >= 3)
      wilcox.test(r$aligned_z[ins], r$aligned_z[!ins], alternative = "greater")$p.value else NA_real_
    as.data.table(out)
  }), fill = TRUE)
  fwrite(resB4, file.path(out_dir, "B4_mouse_informed_attenuation.csv"))
  cat("[B4] selection-free attenuation test (Up / small P = response attenuated in THRB-KO):\n"); print(resB4)
} else message("mouse_T3_results.csv not found: skipping section B4")

## ---- C. strict criteria: altered classification among WT-responsive genes -------
resC <- rbindlist(lapply(names(sets_list), function(s) {
  mem <- sets_list[[s]]
  r <- hep[resp_wt == TRUE]
  in_set <- r$gene %in% mem
  a <- sum(r$altered & in_set); b <- sum(!r$altered & in_set)
  c <- sum(r$altered & !in_set); d <- sum(!r$altered & !in_set)
  data.table(set = s, wt_responsive_in_set = sum(in_set), altered_in_set = a,
             altered_rate_other = c / (c + d),
             fisher_p = if (sum(in_set) > 0)
               fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE), alternative = "greater")$p.value else NA)
}))
fwrite(resC, file.path(out_dir, "C_altered_among_WT_responsive.csv"))
cat("\n[C] Strict criteria (adj P < 0.05 and |log2FC| >= 1 for WT and interaction):\n"); print(resC)

## =============================================================================
## D. Gene-level table for every set member
## =============================================================================
gene_tab <- hep[gene %in% sets$human_symbol & tested_wt == TRUE,
  .(gene, WT_log2FC = log2FoldChange_WT, WT_se = lfcSE_WT, WT_p = pvalue_WT, WT_padj = padj_WT,
    KO_log2FC = log2FoldChange_KO, KO_se = lfcSE_KO, KO_p = pvalue_KO, KO_padj = padj_KO,
    INT_log2FC = log2FoldChange_int, INT_p = pvalue_int, INT_padj = padj_int,
    basal_KOvsWT_log2FC = log2FoldChange_KOvsWT_basal, basal_padj = padj_KOvsWT_basal,
    att_z, wt_nominal, resp_wt, altered)]
gene_tab[, KO_to_WT_ratio := KO_log2FC / WT_log2FC]
gene_tab <- merge(sets[, .(set, gene = human_symbol)], gene_tab, by = "gene")
setorder(gene_tab, set, WT_p)
fwrite(gene_tab, file.path(out_dir, "D_gene_level_redox_members.csv"))
cat("\n[D] WT-nominal members of the glutathione and P450 sets:\n")
print(gene_tab[wt_nominal == TRUE & set %in% c("Glutathione-related", "Microsomal P450 oxidoreduction"),
               .(set, gene, WT_log2FC = round(WT_log2FC, 2), WT_padj = signif(WT_padj, 2),
                 KO_log2FC = round(KO_log2FC, 2), INT_log2FC = round(INT_log2FC, 2),
                 INT_p = signif(INT_p, 2), ratio = round(KO_to_WT_ratio, 2))])


## =============================================================================
## G. Prespecified set-level tests in each contrast (WT, KO, interaction)
##    Background = all genes with a valid statistic in that contrast.
##    cameraPR (two-sided), Wald statistic and signed -log10 P score, rho = 0.01 is the
##    primary setting, 0.05 and 0.10 are sensitivity settings; Bonferroni over the 6 sets.
##    Reproduces the glutathione interaction result of the response letter
##    (rho = 0.01: P = 0.994 with Wald statistics, P = 0.725 with signed scores).
## =============================================================================
contr <- list(WT  = c("log2FoldChange_WT",  "lfcSE_WT",  "pvalue_WT"),
              KO  = c("log2FoldChange_KO",  "lfcSE_KO",  "pvalue_KO"),
              INT = c("log2FoldChange_int", "lfcSE_int", "pvalue_int"))
resG <- rbindlist(lapply(names(contr), function(cn) {
  cc <- contr[[cn]]
  u <- hep[is.finite(get(cc[1])) & is.finite(get(cc[2])) & is.finite(get(cc[3])) & get(cc[2]) > 0]
  stats <- list(wald = setNames(u[[cc[1]]] / u[[cc[2]]], u$gene),
                signed_score = setNames(-log10(u[[cc[3]]]) * sign(u[[cc[1]]]), u$gene))
  rbindlist(lapply(names(stats), function(st) rbindlist(lapply(names(sets_list), function(sn) {
    rbindlist(lapply(rhos, function(rho) {
      cm <- camera_set(stats[[st]], sets_list[[sn]], rho)
      data.table(contrast = cn, statistic = st, set = sn, rho = rho, n_universe = nrow(u),
                 n_set_genes = cm$n, direction = cm$dir, cameraPR_p = cm$p)
    })) }))))
}))
resG[, bonferroni_6sets := pmin(1, cameraPR_p * length(sets_list)), by = .(contrast, statistic, rho)]
fwrite(resG, file.path(out_dir, "G_prespecified_set_tests_WT_KO_INT.csv"))
cat("\n[G] Glutathione set, all contrasts (rho = 0.01):\n")
print(resG[set == "Glutathione-related" & rho == 0.01,
           .(contrast, statistic, n_set_genes, direction, cameraPR_p = signif(cameraPR_p, 3), bonferroni_6sets = signif(bonferroni_6sets, 3))])

## =============================================================================
## E. Figures
## =============================================================================
pdf(file.path(out_dir, "E_KO_redox_figures.pdf"), width = 9, height = 7)

# E1: attenuation score of set members versus all other WT-nominal genes
plot_dt <- rbindlist(lapply(names(sets_list), function(s) {
  d <- copy(uni_att[, .(gene, att_z)]); d[, in_set := gene %in% sets_list[[s]]]; d[, set := s]; d }))
plot_dt[, att_z_c := pmax(pmin(att_z, 8), -8)]
print(ggplot(plot_dt, aes(x = in_set, y = att_z_c)) +
  geom_violin(data = plot_dt[in_set == FALSE], fill = "grey85", colour = NA) +
  geom_boxplot(data = plot_dt[in_set == FALSE], width = 0.15, outlier.shape = NA) +
  geom_jitter(data = plot_dt[in_set == TRUE], width = 0.1, colour = "#B65F55", size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~ set, ncol = 3) +
  scale_x_discrete(labels = c("FALSE" = "other genes", "TRUE" = "set members")) +
  labs(x = NULL, y = "attenuation score (>0: response attenuated in THRB-KO)",
       title = "THRB loss and redox sets (WT-nominal genes)") + theme_bw(base_size = 10))

# E2: WT versus KO log2 fold change (95% CI) for WT-nominal members of the key sets
key <- gene_tab[wt_nominal == TRUE & set %in% c("Glutathione-related", "Microsomal P450 oxidoreduction")]
if (nrow(key) > 0) {
  lng <- rbind(
    key[, .(set, gene, genotype = "WT",      lfc = WT_log2FC, se = WT_se)],
    key[, .(set, gene, genotype = "THRB-KO", lfc = KO_log2FC, se = KO_se)])
  lng[, gene := factor(gene, levels = rev(unique(key[order(WT_log2FC)]$gene)))]
  print(ggplot(lng, aes(x = lfc, y = gene, colour = genotype)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_pointrange(aes(xmin = lfc - 1.96 * se, xmax = lfc + 1.96 * se),
                    position = position_dodge(0.5), size = 0.35, orientation = "y") +
    scale_colour_manual(values = c(WT = "#3973A5", "THRB-KO" = "#B65F55")) +
    facet_grid(set ~ ., scales = "free_y", space = "free_y") +
    labs(x = "T3 effect (log2 fold change, 95% CI)", y = NULL,
         title = "T3 response of WT-nominal redox genes in WT and THRB-KO hepatocytes") +
    theme_bw(base_size = 10))
}
dev.off()

## =============================================================================
## F. Which functional themes characterise the THRB-altered genes? (GO BP ORA)
## =============================================================================
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  bg <- hep[tested_wt == TRUE & tested_int == TRUE, gene]
  groups <- list(
    altered_all   = hep[altered == TRUE, gene],
    altered_up    = hep[altered == TRUE & log2FoldChange_WT > 0, gene],   # THRB-dependent induction
    altered_down  = hep[altered == TRUE & log2FoldChange_WT < 0, gene])   # THRB-dependent repression
  for (g in names(groups)) {
    ego <- tryCatch(clusterProfiler::enrichGO(groups[[g]], OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
            ont = "BP", universe = bg, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
            pvalueCutoff = 1, qvalueCutoff = 1), error = function(e) NULL)
    if (!is.null(ego)) {
      fwrite(as.data.table(as.data.frame(ego)), file.path(out_dir, sprintf("F_GO_BP_%s.csv", g)))
      cat(sprintf("\n[F] %s (n = %d): top terms\n", g, length(groups[[g]])))
      print(head(as.data.frame(ego)[, c("Description", "Count", "p.adjust")], 8))
    }
  }
} else {
  message("clusterProfiler / org.Hs.eg.db not installed: skipping section F")
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDone. Results are in", out_dir, "\n")
