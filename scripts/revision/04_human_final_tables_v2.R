#!/usr/bin/env Rscript
# =============================================================================
# 04_human_final_tables_v2.R     (revision analysis; Reviewer 1, comments 2, 3, 5)
#   v2: T3 baseMean now uses the size factors of the submitted DESeq2 run (see T3 comment below)
#
# Produces the tables that extend Table S9 / S11 for the response to
# Reviewer 1 (comments 2, 3, 5):
#   T1  26-gene panel: WT, KO and interaction effects with 95% CI in human hepatocytes
#       (23 genes are expressed; the others are listed as "not detected")
#   T2  annotation of all 36 curated glutathione genes: 28 testable and 8 untested in mouse,
#       functional module, mouse effect, human ortholog effect
#   T3  baseline maturity markers (AFP, DLK1, ALB, ASGR1, CYP3A4, CYP2E1, CYP1A2, UGT1A1, GSTM1)
#       from the raw counts (optional; needs GSE270233 counts); baseMean equals the submitted Table S9a values
# The six-set interaction table (WT, KO, interaction) is produced by 01_KO_redox_interaction.R (section G).
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
source(Sys.getenv("STAGE_CFG", "config.R"))
out_dir <- file.path(out_root, "04_human_final_tables_v2"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hep <- fread(hep_deg_file); hep[, gene := trimws(gene)]; hep <- hep[!duplicated(gene)]
mo  <- fread(mouse_deg_file); mo[, gene := trimws(gene)]
msets <- fread(mouse_sets_file); msets[, mouse_symbol := trimws(mouse_symbol)]

## ---- T1: panel table -----------------------------------------------------------------
ci <- function(l, se) sprintf("%.2f (%.2f, %.2f)", l, l - 1.96 * se, l + 1.96 * se)
t1 <- data.table(gene = panel_genes)
t1 <- merge(t1, hep, by = "gene", all.x = TRUE, sort = FALSE)
t1[, expressed_in_hepatocytes := !is.na(padj_WT)]
t1[, `:=`(
  WT_log2FC_CI  = ifelse(expressed_in_hepatocytes, ci(log2FoldChange_WT,  lfcSE_WT),  "not detected"),
  KO_log2FC_CI  = ifelse(expressed_in_hepatocytes, ci(log2FoldChange_KO,  lfcSE_KO),  "not detected"),
  INT_log2FC_CI = ifelse(expressed_in_hepatocytes, ci(log2FoldChange_int, lfcSE_int), "not detected"),
  WT_responsive = expressed_in_hepatocytes & padj_WT < 0.05 & abs(log2FoldChange_WT) >= 1,
  altered_response = expressed_in_hepatocytes & padj_WT < 0.05 & abs(log2FoldChange_WT) >= 1 &
                     padj_int < 0.05 & abs(log2FoldChange_int) >= 1)]
t1 <- t1[match(panel_genes, gene),
         .(gene, expressed_in_hepatocytes, WT_log2FC_CI, WT_padj = padj_WT, KO_log2FC_CI, KO_padj = padj_KO,
           INT_log2FC_CI, INT_padj = padj_int, baseline_KOvsWT_log2FC = log2FoldChange_KOvsWT_basal,
           baseline_padj = padj_KOvsWT_basal, WT_responsive, altered_response)]
fwrite(t1, file.path(out_dir, "T1_panel_WT_KO_interaction_with_CI.csv"))
cat(sprintf("[T1] panel genes: %d | expressed: %d | WT-responsive: %d | altered: %d\n",
            nrow(t1), sum(t1$expressed_in_hepatocytes), sum(t1$WT_responsive, na.rm = TRUE),
            sum(t1$altered_response, na.rm = TRUE)))
print(t1[expressed_in_hepatocytes == FALSE, gene])

## ---- T2: 36-gene glutathione annotation --------------------------------------------------
gsh <- msets[set == "Glutathione-related", .(mouse_symbol)]
module_of <- function(g) {
  g1 <- tolower(g)
  if (g1 %in% c("gclc", "gclm", "gss"))                      return("de novo synthesis")
  if (grepl("^(gst|mgst)", g1))                               return("Phase II conjugation")
  if (grepl("^gpx", g1))                                      return("hydroperoxide reduction")
  if (g1 %in% c("gsr", "idh1"))                               return("disulfide reduction / NADPH supply")
  if (grepl("^ggt", g1))                                      return("extracellular turnover")
  "UNASSIGNED (classify manually)"
}
gsh[, module := vapply(mouse_symbol, module_of, "")]
t2 <- merge(gsh, mo[, .(mouse_symbol = gene, human_genes, homology_types, mouse_log2FC = log2FC,
                        mouse_padj = padj, mouse_responsive = responsive)], by = "mouse_symbol", all.x = TRUE)
t2[, testable_in_mouse := !is.na(mouse_padj)]
t2[, mouse_direction := fifelse(!testable_in_mouse, NA_character_,
                                fifelse(mouse_log2FC < 0, "down", "up"))]
t2[, met_2fold_criterion := fifelse(!testable_in_mouse, NA, mouse_padj < 0.05 & abs(mouse_log2FC) >= 1)]
# first human ortholog (one-to-one mapped genes have exactly one symbol)
t2[, human_symbol := trimws(sub(";.*", "", human_genes))]
t2 <- merge(t2, hep[, .(human_symbol = gene, human_WT_log2FC = log2FoldChange_WT, human_WT_padj = padj_WT,
                        human_INT_log2FC = log2FoldChange_int, human_INT_padj = padj_int)],
            by = "human_symbol", all.x = TRUE)
t2[, reason_untested := fifelse(testable_in_mouse, "", "not retained in the mouse analysis (low expression or not mapped)")]
mod_order <- c("de novo synthesis", "Phase II conjugation", "hydroperoxide reduction",
               "disulfide reduction / NADPH supply", "extracellular turnover")
setorder(t2, module, -testable_in_mouse, mouse_padj)
t2 <- t2[, .(module, mouse_symbol, human_symbol, homology_types, testable_in_mouse, mouse_log2FC, mouse_padj,
             mouse_direction, met_2fold_criterion, human_WT_log2FC, human_WT_padj, human_INT_log2FC,
             human_INT_padj, reason_untested)]
fwrite(t2, file.path(out_dir, "T2_glutathione_36_gene_annotation.csv"))
cat(sprintf("\n[T2] glutathione genes: %d | testable in mouse: %d | untested: %d\n", nrow(t2),
            sum(t2$testable_in_mouse), sum(!t2$testable_in_mouse)))
cat("[T2] module counts (testable / untested):\n")
print(t2[, .(testable = sum(testable_in_mouse), untested = sum(!testable_in_mouse),
             met_2fold = sum(met_2fold_criterion, na.rm = TRUE)), by = module])
if (any(grepl("UNASSIGNED", t2$module))) message("Classify manually: ", paste(t2[grepl("UNASSIGNED", module), mouse_symbol], collapse = ", "))

## ---- T3: maturity markers (optional, raw counts) ---------------------------------------------
# baseMean = mean of size-factor-normalized counts over all 12 samples, as in DESeq2.
# The submitted GSE270233 run (261019_GSE270233_DEG.R) first kept genes with count >= 10 in >= 3 samples
# (33,551 -> 15,915 genes) and DESeq2 estimated the size factors on that gene set only. v1 of this script
# estimated them on every gene of the count file, which made baseMean 0.04-0.06% lower (size factors differ by
# up to 0.13%). Here the size factors come from the same pre-filtered set and are then applied to all genes,
# so genes removed by the pre-filter (e.g. CYP3A4) still get a baseMean on the same scale.
hep_prefilter_min_samples <- 3      # smallest group size of the 2 x 2 design, as in the submitted run
if (file.exists(hep_counts_file)) {
  suppressPackageStartupMessages(library(DESeq2))
  cnt <- fread(hep_counts_file); setnames(cnt, 1, "gene"); cnt[, gene := trimws(gene)]; cnt <- cnt[!duplicated(gene)]
  m <- as.matrix(cnt[, -1]); rownames(m) <- cnt$gene; storage.mode(m) <- "integer"
  keep_deseq <- rowSums(m >= min_count) >= hep_prefilter_min_samples
  cat(sprintf("\n[T3] pre-filter (count >= %d in >= %d samples): %d -> %d genes\n", min_count,
              hep_prefilter_min_samples, nrow(m), sum(keep_deseq)))
  dds <- DESeqDataSetFromMatrix(m[keep_deseq, ], data.frame(row.names = colnames(m), x = factor(rep("a", ncol(m)))), ~ 1)
  sf  <- sizeFactors(estimateSizeFactors(dds))
  nc  <- sweep(m, 2, sf, "/")
  cpm <- t(t(m) / colSums(m) * 1e6)
  g <- intersect(maturity_genes, rownames(m))
  n10 <- rowSums(m[g, , drop = FALSE] >= min_count)
  # a gene is "tested" in the submitted run when it has a non-missing adjusted P (DESeq2 also sets padj to NA
  # for genes below its independent-filtering threshold, so passing the pre-filter is not enough)
  padj_sub <- hep$padj_WT[match(g, hep$gene)]
  t3 <- data.table(gene = g,
                   baseMean_all_samples = rowMeans(nc[g, , drop = FALSE]),
                   mean_CPM = rowMeans(cpm[g, , drop = FALSE]),
                   n_samples_count_ge10 = as.integer(n10),
                   passes_DESeq2_prefilter_ge3_samples = keep_deseq[g],
                   passes_half_of_samples_filter = n10 >= ncol(m) / 2,
                   tested_in_submitted_DESeq2 = !is.na(padj_sub))
  t3 <- t3[match(maturity_genes, gene)][!is.na(gene)]
  fwrite(t3, file.path(out_dir, "T3_maturity_markers.csv"))
  cat("\n[T3] baseline maturity markers (baseMean = mean of normalized counts over all samples):\n"); print(t3)
} else message("\n[T3] skipped: ", hep_counts_file, " not found (baseMean for these genes is also in Table S9a)")
