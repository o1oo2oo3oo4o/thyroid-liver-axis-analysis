#!/usr/bin/env Rscript
# =============================================================================
# 05_mouse_QC_leave_one_animal_out_v3-mainmap.R        (revision analysis; Reviewer 1, comment 4)
#
# v3-mainmap = v2-apeglm + ONE change: the mouse-to-human ortholog table used for the GO term.
#   v2 read data/mouse_T3_results.csv, whose human_genes column is the synonym-resolved sensitivity mapping
#   (498 foreground / 11,290 background genes; xenobiotic adj P 4.63e-18 for all animals).
#   v3 reads data/ortholog_mouse_status_all.tsv, the mapping of the manuscript's main analysis (Table S10:
#   448 foreground / 10,085 background genes; xenobiotic adj P 3.44e-18), so the eight leave-one-animal-out
#   values are directly comparable with the manuscript. The gene counts used for each GO test are now written to
#   the output (GO_n_foreground, GO_n_universe). Nothing else differs from v2.
#
# v2-apeglm: identical to 05_mouse_QC_leave_one_animal_out.R EXCEPT that a gene is called "responsive"
# (adj P < 0.05 and |log2FC| >= 1) using the apeglm-SHRUNKEN log2FC, as in the manuscript's DEG table
# (541 responsive genes). Version 1 used the unshrunken log2FC (779 genes), so its Fisher P and GO results
# were not on the manuscript's gene pool. adj P / P are the DESeq2 Wald values; the cameraPR score keeps the
# sign of the unshrunken log2FC, so cameraPR results equal version 1. median_log2FC_set and frac_set_down
# now use the shrunken log2FC (the unshrunken median is kept in median_log2FC_set_unshrunk).
# Needs the Bioconductor package apeglm.
#
# Reviewer 1, comment 4: does any single animal dominate the pathway-level results?
#   1  sample QC: PCA, library sizes, Cook's distance
#   2  eight leave-one-animal-out (LOAO) DESeq2 fits with the identical thresholds and fixed gene sets
#      reporting, for each fit: number of responsive genes, glutathione set size / direction / median
#      effect / Fisher P / cameraPR P, and the adjusted P of "xenobiotic metabolic process"
#   3  a summary of the ranges
# Interpretation rule: the analysis tests whether ONE animal drives the direction;
# it does NOT require all eight fits to be significant and does not increase the effective sample size.
#
# Input : config.R, GSE85793 counts and metadata, mouse gene sets,
#         data/ortholog_mouse_status_all.tsv (main-analysis ortholog mapping, needed for the GO term; the file written
#         by scripts/04_reortholog_synonym.R, copied into data/)
# GO term (optional): requires clusterProfiler + org.Hs.eg.db (human orthologs, as in Methods 2.6)
# =============================================================================
suppressPackageStartupMessages({ library(DESeq2); library(data.table); library(limma); library(ggplot2) })
if (!requireNamespace("apeglm", quietly = TRUE))
  stop("Package 'apeglm' is required: BiocManager::install('apeglm'). Not falling back to unshrunken log2FC on purpose.")
source(Sys.getenv("STAGE_CFG", "config.R"))
out_dir <- file.path(out_root, "05_mouse_QC_LOAO_v3-mainmap"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rhos <- c(0.01, 0.10)
go_term <- "xenobiotic metabolic process"
ortholog_main_file <- file.path(data_dir, "ortholog_mouse_status_all.tsv")   # main-analysis mapping (manuscript Table S10)

## ---- data ------------------------------------------------------------------------------------
cnt <- fread(mouse_counts_file); setnames(cnt, 1, "gene"); cnt[, gene := trimws(gene)]; cnt <- cnt[!duplicated(gene)]
counts_all <- as.matrix(cnt[, -1]); rownames(counts_all) <- cnt$gene; storage.mode(counts_all) <- "integer"
meta <- as.data.frame(fread(mouse_meta_file)); rownames(meta) <- meta$sample
counts_all <- counts_all[, rownames(meta)]
meta$group <- factor(meta[[mouse_group_col]], levels = c(mouse_ctrl_label, mouse_treat_label))
counts <- counts_all[rowSums(counts_all >= mouse_filter_min_count) >= ncol(counts_all) / 2, ]
msets <- fread(mouse_sets_file); msets[, mouse_symbol := trimws(mouse_symbol)]
sets_list <- split(msets$mouse_symbol, msets$set)
if (!file.exists(ortholog_main_file))
  stop("Missing ", ortholog_main_file, ": copy ortholog_mouse_status_all.tsv (written by scripts/04_reortholog_synonym.R) into data/. ",
       "v3 does not skip the GO term silently, because the comparison with the manuscript depends on it.")
orth <- fread(ortholog_main_file, select = c("gene", "human_genes"), na.strings = c("", "NA"))
orth[, gene := trimws(gene)]

## ---- 1. sample QC -------------------------------------------------------------------------------
dds_full <- DESeq(DESeqDataSetFromMatrix(counts, meta, ~ group), quiet = TRUE)
vst_mat <- assay(varianceStabilizingTransformation(dds_full, blind = FALSE))
pca <- prcomp(t(vst_mat), scale. = FALSE)
pc <- data.table(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = meta$group)
var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
cook <- assays(dds_full)[["cooks"]]
qc <- data.table(sample = colnames(counts), group = meta$group,
                 library_size = colSums(counts_all), library_size_filtered = colSums(counts),
                 size_factor = sizeFactors(dds_full),
                 max_cooks_distance = apply(cook, 2, max, na.rm = TRUE),
                 n_genes_with_cooks_above_F_cutoff = colSums(cook > qf(0.99, 2, ncol(counts) - 2), na.rm = TRUE),
                 PC1 = pc$PC1, PC2 = pc$PC2)
fwrite(qc, file.path(out_dir, "1_sample_QC.csv"))
cat("[1] sample QC (PC1 explains", var_expl[1], "%, PC2", var_expl[2], "%):\n"); print(qc[, .(sample, group, library_size, size_factor = round(size_factor, 2),
      max_cooks = round(max_cooks_distance, 2), PC1 = round(PC1, 1), PC2 = round(PC2, 1))])
pdf(file.path(out_dir, "1_sample_QC.pdf"), width = 9, height = 4)
print(ggplot(pc, aes(PC1, PC2, colour = group, label = sample)) + geom_point(size = 3) + geom_text(vjust = -1, size = 3) +
        labs(x = paste0("PC1 (", var_expl[1], "%)"), y = paste0("PC2 (", var_expl[2], "%)"), title = "GSE85793 PCA (VST)") + theme_bw())
print(ggplot(qc, aes(sample, library_size / 1e6, fill = group)) + geom_col() + labs(y = "library size (millions)", x = NULL) + theme_bw())
dev.off()

## ---- 2. leave-one-animal-out ---------------------------------------------------------------------
go_available <- requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)
if (!go_available) message("GO adjusted P skipped (needs clusterProfiler and org.Hs.eg.db); the gene counts are still reported")

go_info <- function(res) {
  map <- orth[gene %in% res$gene & !is.na(human_genes) & human_genes != "", .(gene, human_genes)]
  hum <- function(genes) unique(trimws(unlist(strsplit(map[gene %in% genes, human_genes], ";"))))
  resp <- hum(res[resp == TRUE, gene]); univ <- hum(res$gene)
  out <- list(p = NA_real_, n_fg = length(resp), n_univ = length(univ))
  if (!go_available || length(resp) < 5) return(out)
  e <- tryCatch(clusterProfiler::enrichGO(resp, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", universe = univ,
          pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1), error = function(e) NULL)
  if (is.null(e)) return(out)
  d <- as.data.frame(e); v <- d$p.adjust[d$Description == go_term]
  if (length(v)) out$p <- v
  out
}

fit_one <- function(samples, label) {
  cc <- counts[, samples, drop = FALSE]; mm <- meta[samples, , drop = FALSE]
  cc <- cc[rowSums(cc >= mouse_filter_min_count) >= ncol(cc) / 2, ]
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
  r[, score := -log10(pvalue) * sign(lfc_unshrunk)]            # sign of the unshrunken estimate (same as version 1)
  sc <- setNames(r$score, r$gene)
  rows <- lapply(names(sets_list), function(s) {
    mem <- sets_list[[s]]; ins <- r$gene %in% mem
    a <- sum(r$resp & ins); b <- sum(!r$resp & ins); c <- sum(r$resp & !ins); d <- sum(!r$resp & !ins)
    idx <- which(names(sc) %in% mem)
    o <- list(fit = label, set = s, n_tested = sum(ins), n_responsive_total = sum(r$resp),
              n_set_responsive = a, n_set_down = sum(r$resp & ins & r$log2FoldChange < 0),
              n_set_up = sum(r$resp & ins & r$log2FoldChange > 0),
              median_log2FC_set = median(r$log2FoldChange[ins]),
              median_log2FC_set_unshrunk = median(r$lfc_unshrunk[ins]),
              frac_set_down = mean(r$log2FoldChange[ins] < 0),
              fisher_p = fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE), alternative = "greater")$p.value)
    for (rho in rhos) o[[sprintf("cameraPR_p_rho%.2f", rho)]] <- cameraPR(sc, list(set = idx), inter.gene.cor = rho, use.ranks = FALSE)$PValue
    as.data.table(o)
  })
  t <- rbindlist(rows)
  g <- go_info(r)
  t[, `:=`(xenobiotic_GO_adjP = g$p, GO_n_foreground = g$n_fg, GO_n_universe = g$n_univ)]
  t
}
animals <- rownames(meta)
fits <- c(list(list(label = "all_animals", s = animals)),
          lapply(animals, function(a) list(label = paste0("without_", a), s = setdiff(animals, a))))
res <- rbindlist(lapply(fits, function(z) fit_one(z$s, z$label)))
fwrite(res, file.path(out_dir, "2_leave_one_animal_out_all_fits.csv"))
ra <- res[fit == "all_animals" & set == "Glutathione-related"]
cat(sprintf("\n[2] reproduction check, all animals: %d responsive genes (manuscript 541); glutathione Fisher P = %.3g (manuscript 2.07e-4); xenobiotic GO adj P = %s (manuscript 3.44e-18); GO foreground / background genes = %d / %d (manuscript 448 / 10,085)\n",
            ra$n_responsive_total, ra$fisher_p, format(ra$xenobiotic_GO_adjP, digits = 3), ra$GO_n_foreground, ra$GO_n_universe))
cat("\n[2] glutathione set across the fits:\n")
print(res[set == "Glutathione-related", .(fit, n_responsive_total, n_set_responsive, n_set_down, n_set_up,
          median_log2FC = round(median_log2FC_set, 2), fisher_p = signif(fisher_p, 3),
          cameraPR_rho0.01 = signif(cameraPR_p_rho0.01, 3), xenobiotic_GO_adjP = signif(xenobiotic_GO_adjP, 3),
          GO_fg = GO_n_foreground, GO_bg = GO_n_universe)])

## ---- 3. summary of ranges for the manuscript ------------------------------------------------------
loao <- res[fit != "all_animals"]
summ <- loao[, .(n_fits = .N,
                 fisher_p_min = min(fisher_p), fisher_p_max = max(fisher_p),
                 cameraPR_rho0.01_min = min(cameraPR_p_rho0.01), cameraPR_rho0.01_max = max(cameraPR_p_rho0.01),
                 cameraPR_rho0.10_max = max(cameraPR_p_rho0.10),
                 median_log2FC_min = min(median_log2FC_set), median_log2FC_max = max(median_log2FC_set),
                 fits_with_majority_down = sum(frac_set_down > 0.5),
                 fits_fisher_p_below_0.05 = sum(fisher_p < 0.05),
                 fits_cameraPR_rho0.01_below_0.05 = sum(cameraPR_p_rho0.01 < 0.05),
                 xenobiotic_GO_adjP_min = if (all(is.na(xenobiotic_GO_adjP))) NA_real_ else min(xenobiotic_GO_adjP, na.rm = TRUE),
                 xenobiotic_GO_adjP_max = if (all(is.na(xenobiotic_GO_adjP))) NA_real_ else max(xenobiotic_GO_adjP, na.rm = TRUE)),
             by = set]
fwrite(summ, file.path(out_dir, "3_leave_one_animal_out_summary.csv"))
cat("\n[3] ranges over the eight leave-one-animal-out fits (report them as they are):\n"); print(summ[set == "Glutathione-related"])

pdf(file.path(out_dir, "2_leave_one_animal_out.pdf"), width = 8, height = 5)
g <- res[set == "Glutathione-related"]; g[, fit := factor(fit, levels = unique(fit))]
print(ggplot(g, aes(fit, -log10(fisher_p))) + geom_col(fill = "#3973A5") + geom_hline(yintercept = -log10(0.05), lty = 2) +
        coord_flip() + labs(x = NULL, y = "-log10 Fisher P (glutathione set)", title = "Leave-one-animal-out") + theme_bw())
print(ggplot(g, aes(fit, median_log2FC_set)) + geom_col(fill = "#B65F55") + coord_flip() +
        labs(x = NULL, y = "median log2 fold change of the glutathione set", title = "Direction and size of the set effect") + theme_bw())
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDone. Results are in", out_dir, "\n")
