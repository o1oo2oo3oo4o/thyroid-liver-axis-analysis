################################################################################
#  11_THRB_candidate_enrichment.R
#  Two one-sided Fisher exact tests asking whether the prioritised candidate genes more often show a THRB-dependent T3 response
#  than other genes (human hepatocytes, wild-type vs THRB-knockout, GSE270233). Output = Supplementary Table S9b.
#
#  DEFINITIONS
#    responsive in wild-type cells = adjusted P < 0.05 and |log2FC| >= 1 for the T3 effect in wild-type cells (799 genes)
#    significant interaction       = adjusted P < 0.05 and |log2FC| >= 1 for the genotype-by-treatment interaction term
#    altered response (in KO)      = responsive in wild-type cells AND significant interaction (194 genes)
#    panel                         = 26 prioritised candidates (10 locus-based + 16 protein-based, ADH4 in both, + SLC39A8);
#                                    the same definition as in figures/Fig4_THRB_dependency_supp.R; 23 are testable
#                                    (CUX2, ABHD1 and CD19 are not detected)
#    Test 1: candidates vs all genes tested in both analyses (15,297 genes) for "altered response"
#    Test 2: among the 799 wild-type responsive genes, candidates vs the other wild-type responsive genes
#
#  NOTE  The code that produced Table S9b was not kept. This script was re-written for the repository and reproduces
#        Table S9b exactly (2 x 2 counts and both P values; checked against the submitted workbook).
#
#  INPUT (DATA_DIR)
#    T3_effect_WT_KO_interaction_wide.csv  (output of the GSE270233 DESeq2 analysis; columns gene, log2FoldChange_WT, padj_WT,
#                                           log2FoldChange_int, padj_int, ... one row per gene)
#  OUTPUT (OUT_DIR)
#    THRB_candidate_enrichment_2x2.csv
#
#  USAGE  DATA_DIR=data Rscript 11_THRB_candidate_enrichment.R
################################################################################
DATA_DIR <- Sys.getenv("DATA_DIR", "data")
OUT_DIR  <- Sys.getenv("OUT_DIR", "results/11_THRB_candidate_enrichment")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
PADJ_CUT <- 0.05; LFC_CUT <- 1

POS10  <- c("CUX2","BDH2","ADH4","LILRB3","SBK1","ABHD1","INCENP","CDC42EP5","CD19","GPCPD1")
PROT16 <- c("ADH4","BTC","CFHR1","CLIC5","COL11A2","COL15A1","DEF6","DNAJB4","F13B","LGALS4","PKLR","POR","SPARCL1","TIMD4","GSTM1","HSPA1A")
PANEL26 <- sort(unique(c(POS10, PROT16, "SLC39A8")))

w <- read.csv(file.path(DATA_DIR, "T3_effect_WT_KO_interaction_wide.csv"), check.names = FALSE)
sig <- function(p, l) !is.na(p) & p < PADJ_CUT & !is.na(l) & abs(l) >= LFC_CUT
w$wt_resp  <- sig(w$padj_WT,  w$log2FoldChange_WT)
w$int_sig  <- sig(w$padj_int, w$log2FoldChange_int)
w$altered  <- w$wt_resp & w$int_sig
tested <- w[!is.na(w$padj_WT) & !is.na(w$padj_int), ]            # tested in both analyses (15,297 genes)
cand   <- tested[tested$gene %in% PANEL26, ]
cat(sprintf("tested in both: %d | WT responsive: %d | altered: %d | interaction significant (among genes tested in both): %d | candidates testable: %d (not testable: %s)\n",
            nrow(tested), sum(tested$wt_resp), sum(tested$altered), sum(tested$int_sig), nrow(cand), paste(setdiff(PANEL26, cand$gene), collapse = ", ")))

two_by_two <- function(label, cand_alt, cand_n, other_alt, other_n, ncand, nother) {
  m <- matrix(c(cand_alt, cand_n - cand_alt, other_alt, other_n - other_alt), 2, byrow = TRUE)
  p <- fisher.test(m, alternative = "greater")$p.value
  data.frame(test = label,
             group = c(sprintf("candidate (n=%d)", cand_n), sprintf("other (n=%d)", other_n)),
             altered_response_in_KO = c(cand_alt, other_alt), not_altered = c(cand_n - cand_alt, other_n - other_alt),
             one_sided_Fisher_P = p)
}
t1 <- two_by_two("Test 1: candidates vs all tested genes", sum(cand$altered), nrow(cand), sum(tested$altered) - sum(cand$altered), nrow(tested) - nrow(cand))
wr <- tested[tested$wt_resp, ]; cw <- wr[wr$gene %in% PANEL26, ]
t2 <- two_by_two("Test 2: candidates vs wild-type responsive genes", sum(cw$altered), nrow(cw), sum(wr$altered) - sum(cw$altered), nrow(wr) - nrow(cw))
out <- rbind(t1, t2)
print(out, row.names = FALSE)
write.csv(out, file.path(OUT_DIR, "THRB_candidate_enrichment_2x2.csv"), row.names = FALSE)
