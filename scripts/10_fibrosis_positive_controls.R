################################################################################
#  10_fibrosis_positive_controls.R
#  Technical check that the fibrosis stages and the samples are paired correctly in the two human liver cohorts (not a hypothesis test):
#  twenty established fibrosis-associated genes should increase with fibrosis stage in the primary model of both cohorts.
#  Output = sheet "stage_check" of the results workbook of the manuscript (Supplementary Table S13d).
#
#  NOTE  The stage_check sheet was added to the results workbook of 09_human_liver_fibrosis.R after that script had been run; the code that
#        produced it was not kept. This script was re-written for the repository and gives the same values as the stage_check sheet of
#        the delivered workbook (261021_human_liver_fibrosis_results_S12.xlsx): same 20 genes, same log2FC and P values (checked to
#        machine precision), same flags.
#
#  INPUTS (OUT_DIR of 09_human_liver_fibrosis.R)
#    261021_human_liver_fibrosis_gene_results_GSE130970_primary.csv
#    261021_human_liver_fibrosis_gene_results_GSE135251_primary.csv     (columns: gene, baseMean, log2FC, lfcSE, stat, pvalue, padj)
#  OUTPUT
#    261021_human_liver_fibrosis_stage_check.csv   (columns as in the sheet stage_check)
#
#  USAGE  FIBROSIS_DIR=results/09_human_liver_fibrosis Rscript 10_fibrosis_positive_controls.R
################################################################################
FIBROSIS_DIR <- Sys.getenv("FIBROSIS_DIR", "results/09_human_liver_fibrosis")
PREFIX <- "261021_human_liver_fibrosis"
POSITIVE_CONTROLS <- c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A1", "TIMP1", "TIMP2", "ACTA2", "LUM", "DCN",
                       "SPP1", "THY1", "PDGFRB", "LOXL2", "AKR1B10", "GDF15", "CD44", "VIM", "CCL2", "TGFB1")
read_primary <- function(cohort) {
  f <- file.path(FIBROSIS_DIR, sprintf("%s_gene_results_%s_primary.csv", PREFIX, cohort))
  if (!file.exists(f)) stop("Not found: ", f)
  read.csv(f, check.names = FALSE)
}
a <- read_primary("GSE130970"); b <- read_primary("GSE135251")
ia <- match(POSITIVE_CONTROLS, a$gene); ib <- match(POSITIVE_CONTROLS, b$gene)
if (anyNA(ia) || anyNA(ib)) stop("Genes missing in the results: ", paste(POSITIVE_CONTROLS[is.na(ia) | is.na(ib)], collapse = ", "))
out <- data.frame(gene = POSITIVE_CONTROLS,
                  log2FC_per_stage_GSE130970 = a$log2FC[ia], pvalue_GSE130970 = a$pvalue[ia],
                  log2FC_per_stage_GSE135251 = b$log2FC[ib], pvalue_GSE135251 = b$pvalue[ib])
out$both_positive <- out$log2FC_per_stage_GSE130970 > 0 & out$log2FC_per_stage_GSE135251 > 0
out[["both_P_lt_0.05"]] <- out$pvalue_GSE130970 < 0.05 & out$pvalue_GSE135251 < 0.05
write.csv(out, file.path(FIBROSIS_DIR, paste0(PREFIX, "_stage_check.csv")), row.names = FALSE)
cat(sprintf("%d of %d genes increased with stage in both cohorts; %d had P < 0.05 in both cohorts\n",
            sum(out$both_positive), nrow(out), sum(out[["both_P_lt_0.05"]])))
