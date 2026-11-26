################################################################################
#  11_glutathione_name_matched_concordance.R
#
#  Purpose: from the panel of 36 glutathione genes, select the genes that were tested in BOTH mouse and human by
#           matching gene symbols ("name matching"; this replaces the earlier set of 18 genes, common_18), and
#           recompute the mouse-human direction agreement and the human set tests.
#           The output is the input of Figure 6 (figures/Fig6_glutathione_name_matched.R), Table 7, Table S9c and
#           Table S11.
#
#  Stand-alone script: no other R file of this repository is needed.
#
#  Inputs (all relative to scripts/revision; see data/README.md)
#    MOUSE_FILE     : data/250912__DEG_all_genes.csv
#                     gene, baseMean, log2FoldChange, lfcSE, pvalue, padj  (mouse T3, GSE85793; output of scripts/03_deseq2_mouse_liver_T3.R)
#    HUMAN_WT_FILE  : data/DEG_WT_T3_effect.csv
#                     gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, lfc_shrunk  (GSE270233 wild type, T3 vs vehicle;
#                     output of scripts/08_deseq2_human_hepatocyte_THRB_KO.R)
#    HUMAN_INT_FILE : data/T3_effect_WT_KO_interaction_wide.csv   (optional; without it the interaction check is skipped)
#    ORTHOLOG_FILE  : data/ortholog_mouse_status_all.tsv           (optional; written by scripts/04_reortholog_synonym.R)
#
#  Name-matching rule (fixed before the results were inspected)
#    1) the upper-case form of the mouse symbol (KEY) of the 36 genes
#    2) is identical to a gene symbol of the human wild-type result, and
#    3) has a non-missing padj in both the mouse and the human result (= passed the expression filter and was tested).
#    Only genes that meet all three conditions enter the direction analysis.
#    Note: an identical symbol does not mean an Ensembl ortholog (GSTM2, GSTM3 and GSTA4 are such cases, see table 11g).
#
#  Direction agreement (same method as scripts/07_redox_cameraPR.R)
#    - agreement: sign(mouse log2FC) == sign(human log2FC)
#    - one-sided binomial test against 0.5 and against the agreement rate of all name-matched genes of the transcriptome (background)
#    - shift of the human set: human log2FC of the set genes against all other tested genes, one-sided Mann-Whitney (less)
#    - human set tests: Fisher (over-representation of responsive genes) and cameraPR; main statistic = signed -log10 P
#      (signed = -log10(pvalue) x sign(log2FC), rho 0.01); sensitivity: Wald stat and rho 0 / 0.05 / 0.10
#    - reproduction check: the same code is run for the earlier common_18 set and compared with the earlier values
#      (14/18, P 0.0154, ...)
#
#  Output (OUT_DIR)
#    11a_membership_36genes.csv              for each of the 36 genes: included or not, and the reason
#    11b_concordance_genes_name_matched.csv  read by figures/Fig6_glutathione_name_matched.R
#    11c_Table7_name_matched.csv             Table 7 layout
#    11d_concordance_summary.csv             name_matched and legacy_common_18 summary (agreement, Mann-Whitney)
#    11e_human_set_tests.csv                 human set Fisher and cameraPR tests (name_matched and legacy_common_18)
#    11f_interaction_check.csv               WT vs THRB-KO interaction (when HUMAN_INT_FILE exists)
#    11g_name_vs_Ensembl_ortholog.csv        (when ORTHOLOG_FILE exists) whether the same-name pairs are also Ensembl orthologs
#    11_sessionInfo.txt
#
#  Usage: Rscript 11_glutathione_name_matched_concordance.R   (run from scripts/revision)
#  R packages: limma, dplyr, readr
################################################################################

suppressPackageStartupMessages({
  library(limma)    # cameraPR
  library(dplyr)
  library(readr)
})

# ==============================================================================
# [settings]  every path can also be set with an environment variable (otherwise the defaults below are used)
# ==============================================================================
MOUSE_FILE     <- Sys.getenv("MOUSE_FILE",     "data/250912__DEG_all_genes.csv")
HUMAN_WT_FILE  <- Sys.getenv("HUMAN_WT_FILE",  "data/DEG_WT_T3_effect.csv")
HUMAN_INT_FILE <- Sys.getenv("HUMAN_INT_FILE", "data/T3_effect_WT_KO_interaction_wide.csv")
ORTHOLOG_FILE  <- Sys.getenv("ORTHOLOG_FILE",  "data/ortholog_mouse_status_all.tsv")   # optional; creates 11g when present
OUT_DIR        <- Sys.getenv("OUT_DIR", "results/11_glutathione_name_matched")

PADJ_CUT <- 0.05
LFC_CUT  <- 1
RHO_MAIN <- 0.01
RHO_SENS <- c(0, 0.05, 0.10)
USE_RANKS <- FALSE      # cameraPR; same as the earlier analysis
STAT_MAIN <- "signed"   # main statistic: the values of Table S9 / the manuscript (0.337, 0.85) were computed with "signed"; the other one is written as a sensitivity analysis

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_path <- function(f) file.path(OUT_DIR, f)

# ==============================================================================
# The 36 glutathione genes (mouse symbols; the same list as redox_sets_mouse in scripts/07_redox_cameraPR.R)
# ==============================================================================
glutathione_36 <- c(
  "Gclc", "Gclm", "Ggt1", "Ggt5", "Gpx1", "Gpx2", "Gpx3", "Gpx4", "Gpx6", "Gpx7", "Gsr",
  "Gss", "Gsta1", "Gsta2", "Gsta3", "Gsta4", "Gstk1", "Gstm1", "Gstm2", "Gstm3", "Gstm4",
  "Gstm5", "Gstm6", "Gstm7", "Gsto1", "Gsto2", "Gstp1", "Gstp2", "Gstt1", "Gstt2", "Gstt3",
  "Gstz1", "Idh1", "Mgst1", "Mgst2", "Mgst3"
)
stopifnot(length(glutathione_36) == 36, !anyDuplicated(glutathione_36))

# The 18 genes of the first submission (used only for comparison and the reproduction check; the selection rule was not recorded)
legacy_common_18 <- c("GSTA4","GSTM1","GSTM2","GSTM3","GSTM4","GSTP1","GSTO1","GSTZ1",
                      "GSTK1","MGST1","MGST3","GPX3","GPX4","GSR","GCLC","GCLM","GSS","IDH1")

# ==============================================================================
# Read the inputs
# ==============================================================================
read_deg <- function(path, need) {
  if (!file.exists(path)) stop("Input file not found: ", path, " (see data/README.md)")
  d <- read_csv(path, show_col_types = FALSE)
  miss <- setdiff(need, names(d))
  if (length(miss)) stop(path, ": required columns are missing: ", paste(miss, collapse = ", "))
  d
}
mouse_raw <- read_deg(MOUSE_FILE,    c("gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"))
human_raw <- read_deg(HUMAN_WT_FILE, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))

tested <- function(d) { d <- d[!is.na(d$padj), ]; d[!duplicated(d$gene), ] }   # genes that passed the expression filter and were tested
mouse <- tested(mouse_raw)
human <- tested(human_raw)
if (anyNA(human$stat)) stop("The human result has rows with a padj but a missing stat")
human$signed <- -log10(human$pvalue) * sign(human$log2FoldChange)    # signed -log10 P score
human$responsive <- human$padj < PADJ_CUT & abs(human$log2FoldChange) >= LFC_CUT
cat(sprintf("mouse: %d rows, %d tested | human WT: %d rows, %d tested\n",
            nrow(mouse_raw), nrow(mouse), nrow(human_raw), nrow(human)))

# ==============================================================================
# 11a. does each of the 36 genes enter the name matching?
# ==============================================================================
status_of <- function(g, raw, tst) {
  ifelse(!(g %in% raw$gene), "no row",
         ifelse(g %in% tst$gene, "tested", "padj NA (removed by the expression filter)"))
}
membership <- tibble(mouse_symbol = glutathione_36, KEY = toupper(glutathione_36)) %>%
  mutate(mouse_status = status_of(mouse_symbol, mouse_raw, mouse),
         human_status = status_of(KEY, human_raw, human),
         matched = mouse_status == "tested" & human_status == "tested",
         reason = case_when(
           matched ~ "included (name matched)",
           mouse_status != "tested" & human_status != "tested" ~ "not tested in mouse and not tested in human",
           mouse_status != "tested" ~ "not tested in mouse",
           TRUE ~ "the human gene of the same name was not tested"))
write_csv(membership, out_path("11a_membership_36genes.csv"))
matched_keys <- membership$KEY[membership$matched]
cat(sprintf("\n[name matching] %d of 36 genes included\n", length(matched_keys)))
print(table(membership$reason))

# ==============================================================================
# Direction agreement / human set tests (same calculation as scripts/07_redox_cameraPR.R)
# ==============================================================================
mo <- mouse %>% transmute(KEY = toupper(gene), mouse_symbol = gene, mouse_lfc = log2FoldChange,
                          mouse_padj = padj, mouse_base = baseMean)
hu <- human %>% transmute(KEY = gene, human_lfc = log2FoldChange, human_lfcSE = lfcSE,
                          human_padj = padj, human_base = baseMean)
both <- inner_join(mo, hu, by = "KEY")                                   # all name-matched mouse-human pairs
bg   <- mean(sign(both$mouse_lfc) == sign(both$human_lfc))               # background agreement rate

concordance_for <- function(keys, label) {
  conc <- both %>% filter(KEY %in% keys) %>%
    mutate(concordant = sign(mouse_lfc) == sign(human_lfc), direction = ifelse(concordant, "Concordant", "Discordant"))
  k <- sum(conc$concordant); n <- nrow(conc)
  in_set <- human$gene %in% keys
  summ <- tibble(
    set = label, n_genes = n, n_concordant = k,
    p_vs_half = binom.test(k, n, 0.5, alternative = "greater")$p.value,
    background_rate = round(bg, 3), n_background_pairs = nrow(both),
    p_vs_background = binom.test(k, n, bg, alternative = "greater")$p.value,
    median_mouse_lfc = median(conc$mouse_lfc), median_human_lfc = median(conc$human_lfc),
    median_abs_mouse = median(abs(conc$mouse_lfc)), median_abs_human = median(abs(conc$human_lfc)),
    n_human_abs_ge_1 = sum(abs(conc$human_lfc) >= 1),
    n_mouse_abs_lt_0.05 = sum(abs(conc$mouse_lfc) < 0.05),     # genes whose mouse change is close to 0, so the direction carries little meaning (descriptive)
    n_other_genes = sum(!in_set), median_human_lfc_others = median(human$log2FoldChange[!in_set]),
    wilcox_one_sided_p = wilcox.test(human$log2FoldChange[in_set], human$log2FoldChange[!in_set],
                                     alternative = "less")$p.value)
  list(genes = conc, summary = summ)
}

camera_for <- function(keys, label) {
  idx <- list(set = which(human$gene %in% keys))
  a <- sum(human$gene %in% keys & human$responsive); b <- sum(human$gene %in% keys & !human$responsive)
  cc <- sum(!(human$gene %in% keys) & human$responsive); dd <- sum(!(human$gene %in% keys) & !human$responsive)
  fp <- fisher.test(matrix(c(a, b, cc, dd), nrow = 2, byrow = TRUE), alternative = "greater")$p.value
  stats <- c(STAT_MAIN, setdiff(c("signed", "stat"), STAT_MAIN))          # main statistic first
  bind_rows(lapply(stats, function(st) bind_rows(lapply(c(RHO_MAIN, RHO_SENS), function(r) {
    cm <- cameraPR(human[[st]], idx, use.ranks = USE_RANKS, inter.gene.cor = r, sort = FALSE)
    tibble(set = label, n_tested = a + b, n_responsive = a, fisher_p = fp, statistic = st,
           main_analysis = (st == STAT_MAIN & r == RHO_MAIN),
           rho = r, camera_direction = cm$Direction, camera_p = cm$PValue)
  }))))
}

res_new <- concordance_for(matched_keys,      "name_matched")
res_old <- concordance_for(legacy_common_18,  "legacy_common_18")
summary_tab <- bind_rows(res_new$summary, res_old$summary)
set_tests   <- bind_rows(camera_for(matched_keys, "name_matched"), camera_for(legacy_common_18, "legacy_common_18"))

# ==============================================================================
# Save
# ==============================================================================
genes_out <- res_new$genes %>% arrange(human_lfc) %>%
  select(KEY, mouse_symbol, mouse_lfc, mouse_padj, mouse_base, human_lfc, human_lfcSE, human_padj, human_base,
         concordant, direction)
write_csv(genes_out, out_path("11b_concordance_genes_name_matched.csv"))

table7 <- genes_out %>%
  transmute(Gene = KEY, `Mouse gene` = mouse_symbol, `Mouse log2FC` = round(mouse_lfc, 3),
            `Human log2FC` = round(human_lfc, 3), `Human adj. P` = signif(human_padj, 3),
            `Human baseMean` = round(human_base, 1), Direction = direction)
write_csv(table7, out_path("11c_Table7_name_matched.csv"))
write_csv(summary_tab, out_path("11d_concordance_summary.csv"))
write_csv(set_tests,   out_path("11e_human_set_tests.csv"))

# ---- reproduction check: the values of the earlier common_18 analysis (S9_* files of the earlier version of the tables) ----
expected <- c(n_concordant = 14, p_vs_half = 0.015441894531249997, p_vs_background = 0.1424202901805731,
              median_mouse_lfc = -0.694477330435022, median_human_lfc = -0.2655578144474382,
              wilcox_one_sided_p = 0.010769380095965906)
got <- unlist(res_old$summary[1, names(expected)])
cam_exp <- c(`0.01` = 0.337244, `0` = 0.2993159, `0.05` = 0.4452887, `0.1` = 0.5274601)
cam_got <- set_tests %>% filter(set == "legacy_common_18", statistic == "signed") %>% { setNames(.$camera_p, as.character(.$rho)) }
ok <- all(abs(got - expected) < 1e-6) && all(abs(cam_got[names(cam_exp)] - cam_exp) < 1e-6)
cat("\n[reproduction check] are the earlier common_18 values reproduced?", if (ok) "yes (identical)" else "no (compare the values below)", "\n")
if (!ok) { print(rbind(expected = expected, got = got)); print(rbind(expected = cam_exp, got = cam_got[names(cam_exp)])) }

# ---- WT vs THRB-KO interaction (optional) ----
if (file.exists(HUMAN_INT_FILE)) {
  w <- read_csv(HUMAN_INT_FILE, show_col_types = FALSE)
  need <- c("gene", "log2FoldChange_WT", "padj_WT", "log2FoldChange_KO", "padj_KO", "log2FoldChange_int", "padj_int")
  if (all(need %in% names(w))) {
    inter <- w %>% filter(gene %in% matched_keys) %>%
      transmute(KEY = gene, WT_lfc = log2FoldChange_WT, WT_padj = padj_WT, KO_lfc = log2FoldChange_KO,
                KO_padj = padj_KO, INT_lfc = log2FoldChange_int, INT_padj = padj_int,
                meets_interaction_criterion = !is.na(padj_int) & padj_int < PADJ_CUT & abs(log2FoldChange_int) >= LFC_CUT) %>%
      arrange(INT_padj)
    write_csv(inter, out_path("11f_interaction_check.csv"))
    cat(sprintf("[interaction] of %d name-matched genes, %d meet the interaction criterion (adj P<%.2f, |log2FC|>=%g); smallest adj P = %.3g (%s)\n",
                length(matched_keys), sum(inter$meets_interaction_criterion), PADJ_CUT, LFC_CUT,
                min(inter$INT_padj, na.rm = TRUE), inter$KEY[which.min(inter$INT_padj)]))
  } else cat("[note] the column names of HUMAN_INT_FILE differ from the expected ones; the interaction check is skipped.\n")
} else cat("[note] HUMAN_INT_FILE not found; the interaction check is skipped.\n")

# ---- are the same-name pairs also Ensembl orthologs? (optional) ----
if (nzchar(ORTHOLOG_FILE) && file.exists(ORTHOLOG_FILE)) {
  o <- read_tsv(ORTHOLOG_FILE, show_col_types = FALSE, col_types = cols(.default = "c"))
  chk <- membership %>% filter(matched) %>%
    left_join(o %>% select(gene, human_genes, homology_types), by = c("mouse_symbol" = "gene")) %>%
    mutate(same_name_is_ortholog = mapply(function(k, h) !is.na(h) && k %in% strsplit(h, ";")[[1]], KEY, human_genes)) %>%
    select(mouse_symbol, KEY, ensembl_human_genes = human_genes, homology_types, same_name_is_ortholog)
  write_csv(chk, out_path("11g_name_vs_Ensembl_ortholog.csv"))
  cat(sprintf("[ortholog] of %d same-name pairs, %d are also Ensembl orthologs of that human gene\n",
              nrow(chk), sum(chk$same_name_is_ortholog)))
}

writeLines(capture.output(sessionInfo()), out_path("11_sessionInfo.txt"))

# ==============================================================================
# Console summary (numbers for the text and the legend)
# ==============================================================================
cat("\n================ direction agreement ================\n")
print(as.data.frame(summary_tab), digits = 4)
cat("\n================ human set tests ================\n")
print(as.data.frame(set_tests), digits = 4)
cat("\nOutput folder:", normalizePath(OUT_DIR), "\n")
