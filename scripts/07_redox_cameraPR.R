################################################################################
#  07_redox_cameraPR.R   (redox gene-set analysis in one file)
#  (Original: 261020_redox_cameraPR_1.R, derived from the co-author's 261019_redox_gene_sets_and_cameraPR_1.R)
#
#  STEP 1  Symbol check : checks that the mouse and human symbols of the redox gene sets match Ensembl (no P values are produced)
#  STEP 2  Set analysis : Fisher over-representation and cameraPR for the six sets, glutathione sensitivity analyses,
#                         mouse-human direction concordance
#
#  INPUTS (DATA_DIR)
#    MOUSE_FILE    : 261020__DEG_all_genes_with_stat.csv  (gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
#                    produced by 03_deseq2_mouse_liver_T3.R
#    HUMAN_WT_FILE : DEG_WT_T3_effect.csv  - T3 effect in wild-type human hepatocytes (GSE270233; same columns; gene = upper-case human symbol)
#    HUMAN_KO_FILE : DEG_KO_T3_effect.csv  - T3 effect in THRB-knockout hepatocytes (same columns). If either file is absent only the mouse analysis is run.
#    ENSEMBL_FILE  : Ensembl BioMart table of all mouse genes (GRCm39); needed for STEP 1 only
#    RAW_COUNT     : raw count matrix (only to label the test status in STEP 1; optional)
#
#  DECIDED BEFORE THE RUN (fixed before the results were seen and stated in the Methods)
#    STAT_INPUT          : "stat" (DESeq2 Wald statistic) or "signed" (-log10 p x sign(LFC)); the main analysis of the manuscript uses "signed",
#                          the other statistic is run automatically as a sensitivity analysis (S9_camera_stat_sensitivity.csv)
#    GSTM_REPRESENTATIVE : prespecified rule used when the Gstm family is collapsed to one gene (do not pick the strongest gene)
#    human sets (redox_sets_human) : symbol rule agreed with the co-author after inspecting the STEP 1 output
#    use.ranks, inter.gene.cor and package versions are recorded with S9_sessionInfo.txt
#
#  OUTPUTS (OUT_DIR)
#    STEP 1: redox_symbol_check.xlsx, redox_human_sets_strict.csv, redox_human_sets_all.csv
#    STEP 2: S9_gene_set_results.csv, S9_camera_rho_sensitivity.csv, S9_camera_stat_sensitivity.csv,
#            S9_mouse_human_concordance_genes.csv, S9_mouse_human_concordance_summary.csv,
#            S9_dataset_summary.csv, S9_gene_level_all_sets.csv, S9_glutathione_conditions_membership.csv,
#            S9_glutathione_condition_rules.csv, S9_fisher_2x2_tables.rds, S9_sessionInfo.txt
#    (The "S9" file names are the internal numbering of the analysis files; they correspond to Table S11 of the manuscript.)
#
#  CHANGES relative to the co-author's original script
#    1) the mouse input is a file that has a stat column (the original file had none, so STAT_INPUT = "stat" stopped)
#    2) checks for the required columns and for missing stat values were added
#    3) Bonferroni is applied only to the table in which the six sets are tested together (the original multiplied by the number of sets at every call)
#    4) the mouse analysis is repeated with the other statistic and saved as S9_camera_stat_sensitivity.csv
#    5) if the human input files are absent only the human part is skipped
#    6) STEP 1 (symbol check) was added so that the set definitions are written once, in this file
#
#  USAGE   Rscript 07_redox_cameraPR.R                    # main analysis, signed statistic
#          STAT_INPUT=stat Rscript 07_redox_cameraPR.R    # Wald-statistic run (the signed statistic becomes the sensitivity analysis)
################################################################################

suppressPackageStartupMessages({
  library(limma)       # cameraPR
  library(dplyr)
  library(readr)
  library(data.table)
  library(openxlsx)
})

# ==============================================================================
# [Settings]
# ==============================================================================
RUN_STEP1 <- TRUE      # symbol check (when the Ensembl file is present)
RUN_STEP2 <- TRUE      # Fisher, cameraPR

STAT_INPUT          <- Sys.getenv("STAT_INPUT", "signed")      # "stat" or "signed"
USE_RANKS           <- FALSE       # record this value
RHO_MAIN            <- 0.01
RHO_SENS            <- c(0, 0.05, 0.10)
PADJ_CUT            <- 0.05
LFC_CUT             <- 1
GSTM_REPRESENTATIVE <- "Gstm1"     # prespecified rule, e.g. first alphabetically

DATA_DIR      <- Sys.getenv("DATA_DIR", "data")
MOUSE_FILE    <- file.path(DATA_DIR, "261020__DEG_all_genes_with_stat.csv")
HUMAN_WT_FILE <- file.path(DATA_DIR, "DEG_WT_T3_effect.csv")
HUMAN_KO_FILE <- file.path(DATA_DIR, "DEG_KO_T3_effect.csv")
ENSEMBL_FILE  <- file.path(DATA_DIR, "261019_ensembl_GRCm39.txt")
RAW_COUNT     <- file.path(DATA_DIR, "1_raw_count_matrix_liver_T3.csv")      # NA if absent
OUT_DIR       <- Sys.getenv("OUT_DIR", file.path("results/07_redox_cameraPR", STAT_INPUT))

stopifnot(STAT_INPUT %in% c("stat", "signed"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_path <- function(f) file.path(OUT_DIR, f)

# ---------------------------------------------------------------- gene sets --
# Assembled for this analysis, fixed before it was run, not used in candidate
# selection. Mouse symbols.
redox_sets_mouse <- list(
  "Glutathione-related" = c(
    "Gclc", "Gclm", "Ggt1", "Ggt5", "Gpx1", "Gpx2", "Gpx3", "Gpx4", "Gpx6", "Gpx7", "Gsr",
    "Gss", "Gsta1", "Gsta2", "Gsta3", "Gsta4", "Gstk1", "Gstm1", "Gstm2", "Gstm3", "Gstm4",
    "Gstm5", "Gstm6", "Gstm7", "Gsto1", "Gsto2", "Gstp1", "Gstp2", "Gstt1", "Gstt2", "Gstt3",
    "Gstz1", "Idh1", "Mgst1", "Mgst2", "Mgst3"
  ),
  "NRF2 target genes" = c(
    "Abcc4", "Akr1b8", "Akr1c18", "Cat", "Cbr1", "Cbr3", "Ephx1", "Gclc", "Gclm", "Gsta3",
    "Hmox1", "Keap1", "Mt1", "Mt2", "Nfe2l2", "Nqo1", "Prdx1", "Prdx2", "Prdx3", "Prdx4",
    "Prdx5", "Prdx6", "Sod1", "Sod2", "Sod3", "Srxn1", "Txn1", "Txn2", "Txnrd1"
  ),
  "Iron homeostasis" = c(
    "Aco1", "Bdh2", "Cp", "Fech", "Fth1", "Ftl1", "Hamp", "Hamp2", "Heph", "Hfe", "Hmox1",
    "Ireb2", "Ncoa4", "Slc11a2", "Slc39a14", "Slc39a8", "Slc40a1", "Steap3", "Tf", "Tfr2",
    "Tfrc"
  ),
  "Lipid peroxidation / ferroptosis" = c(
    "Acsl4", "Aifm2", "Alox15", "Alox5", "Ch25h", "Chac1", "Cyb5r3", "Cybb", "Elovl5", "Fads1",
    "Fads2", "Gpx4", "Lpcat3", "Ncoa4", "Nox1", "Nox4", "Pebp1", "Por", "Sat1", "Slc7a11",
    "Tfrc"
  ),
  "Microsomal P450 oxidoreduction" = c(
    "Cyb5a", "Cyb5r3", "Cyp1a1", "Cyp1a2", "Cyp27a1", "Cyp2a4", "Cyp2a5", "Cyp2b10", "Cyp2c29",
    "Cyp2e1", "Cyp3a11", "Cyp3a16", "Cyp3a25", "Cyp4a10", "Cyp4a14", "Cyp4a31", "Cyp7a1",
    "Cyp8b1", "Fmo1", "Fmo3", "Fmo5", "Por"
  ),
  "Fatty acid beta-oxidation (mitochondrial and peroxisomal)" = c(
    "Acadl", "Acadm", "Acads", "Acadvl", "Acox1", "Acox2", "Bdh1", "Cpt1a", "Cpt2", "Echs1",
    "Etfa", "Etfb", "Etfdh", "Hadha", "Hadhb", "Hmgcl", "Hmgcs2", "Oxct1", "Slc25a20"
  )
)

# Human orthologue symbols for the GSE270233 analysis.
redox_sets_human <- list(
  "Glutathione-related" = c(
    "GCLC", "GCLM", "GGT1", "GGT5", "GPX1", "GPX2", "GPX3", "GPX4", "GPX6", "GPX7", "GSR",
    "GSS", "GSTA1", "GSTA2", "GSTA3", "GSTA4", "GSTK1", "GSTM1", "GSTM2", "GSTM3", "GSTM4",
    "GSTM5", "GSTO1", "GSTO2", "GSTP1", "GSTP2", "GSTT1", "GSTT2", "GSTZ1", "IDH1", "MGST1",
    "MGST2", "MGST3"
  ),
  "NRF2 target genes" = c(
    "ABCC4", "AKR1B10", "AKR1C1", "CAT", "CBR1", "CBR3", "EPHX1", "GCLC", "GCLM", "GSTA3",
    "HMOX1", "KEAP1", "MT1X", "MT2A", "NFE2L2", "NQO1", "PRDX1", "PRDX2", "PRDX3", "PRDX4",
    "PRDX5", "PRDX6", "SOD1", "SOD2", "SOD3", "SRXN1", "TXN", "TXN2", "TXNRD1"
  ),
  "Iron homeostasis" = c(
    "ACO1", "BDH2", "CP", "FECH", "FTH1", "FTL1", "HAMP", "HEPH", "HFE", "HMOX1", "IREB2",
    "NCOA4", "SLC11A2", "SLC39A14", "SLC39A8", "SLC40A1", "STEAP3", "TF", "TFR2", "TFRC"
  ),
  "Lipid peroxidation / ferroptosis" = c(
    "ACSL4", "AIFM2", "ALOX15", "ALOX5", "CH25H", "CHAC1", "CYB5R3", "CYBB", "ELOVL5", "FADS1",
    "FADS2", "GPX4", "LPCAT3", "NCOA4", "NOX1", "NOX4", "PEBP1", "POR", "SAT1", "SLC7A11",
    "TFRC"
  ),
  "Microsomal P450 oxidoreduction" = c(
    "CYB5A", "CYB5R3", "CYP1A1", "CYP1A2", "CYP27A1", "CYP2A6", "CYP2A7", "CYP2B6", "CYP2C9",
    "CYP2E1", "CYP3A4", "CYP3A5", "CYP3A7", "CYP4A11", "CYP4A22", "CYP4F11", "CYP7A1", "CYP8B1",
    "FMO1", "FMO3", "FMO5", "POR"
  ),
  "Fatty acid beta-oxidation (mitochondrial and peroxisomal)" = c(
    "ACADL", "ACADM", "ACADS", "ACADVL", "ACOX1", "ACOX2", "BDH1", "CPT1A", "CPT2", "ECHS1",
    "ETFA", "ETFB", "ETFDH", "HADHA", "HADHB", "HMGCL", "HMGCS2", "OXCT1", "SLC25A20"
  )
)

# Genes testable in both species, used for the concordance analysis (Table 6).
common_18 <- c("GSTA4","GSTM1","GSTM2","GSTM3","GSTM4","GSTP1","GSTO1","GSTZ1",
               "GSTK1","MGST1","MGST3","GPX3","GPX4","GSR","GCLC","GCLM","GSS","IDH1")


# Glutathione sensitivity conditions. The fourth removes every cytosolic AND
# microsomal transferase; MGST genes are transferases and must be removed too.
gsh        <- redox_sets_mouse[["Glutathione-related"]]
gstm_fam   <- grep("^Gstm", gsh, value = TRUE)
transferase<- grep("^Gst|^Mgst", gsh, value = TRUE)
gsh_conditions <- list(
  "full"                  = gsh,
  "Gstm collapsed"        = c(setdiff(gsh, gstm_fam), GSTM_REPRESENTATIVE),
  "Gstm removed"          = setdiff(gsh, gstm_fam),
  "all transferases removed" = setdiff(gsh, transferase)
)


# ==============================================================================
# STEP 1  Symbol check of the sets (against Ensembl)
# ==============================================================================
if (RUN_STEP1 && !is.na(ENSEMBL_FILE) && file.exists(ENSEMBL_FILE)) {
  cat("\n[STEP 1] Symbol check of the redox sets\n")
  e <- data.table::fread(ENSEMBL_FILE, sep = "\t", quote = "", colClasses = "character", na.strings = c("", "NA"))
  pk <- function(pat) { i <- grep(pat, names(e), ignore.case = TRUE)[1]; if (is.na(i)) stop("Column not found in the Ensembl file: ", pat); names(e)[i] }
  setnames(e, c(pk("^Gene stable ID$"), pk("^Gene name$"), pk("^Human gene name$"), pk("^Human homology type$"), pk("^Human orthology confidence")),
           c("mid", "mname", "hname", "htype", "hconf"))
  ep    <- e[!is.na(mname)]
  pairs <- unique(ep[!is.na(hname), .(mname, hname, htype, hconf)])
  
  raw_genes <- if (!is.na(RAW_COUNT) && file.exists(RAW_COUNT)) unique(as.character(data.table::fread(RAW_COUNT, select = 1)[[1]])) else NULL
  deg1 <- if (file.exists(MOUSE_FILE)) data.table::fread(MOUSE_FILE) else NULL
  
  rows <- list()
  for (sn in names(redox_sets_mouse)) for (g in redox_sets_mouse[[sn]]) {
    pr <- pairs[mname == g]
    in_ens <- g %in% ep$mname
    deg_status <- if (is.null(deg1)) NA_character_ else if (!(g %in% deg1$gene)) {
      if (!is.null(raw_genes) && !(g %in% raw_genes)) "absent from count matrix" else "excluded by pre-filter"
    } else if (is.na(deg1$padj[match(g, deg1$gene)])) "padj NA" else "tested"
    resp <- if (is.null(deg1) || !(g %in% deg1$gene)) NA else { d <- deg1[match(g, gene)]; !is.na(d$padj) && d$padj < PADJ_CUT && abs(d$log2FoldChange) >= LFC_CUT }
    code <- if (!in_ens) "M1 name not in Ensembl" else if (!nrow(pr)) "M2 no human orthologue" else if (uniqueN(pr$hname) == 1L && all(pr$htype == "ortholog_one2one")) "M0 one-to-one" else "M3 not one-to-one"
    rows[[length(rows) + 1]] <- data.table(set = sn, mouse_gene = g, deg_status = deg_status, responsive = resp, ensembl_code = code,
                                           human_ortholog = if (nrow(pr)) paste(sort(unique(pr$hname)), collapse = ";") else NA_character_,
                                           homology_type = if (nrow(pr)) paste(sort(unique(pr$htype)), collapse = ";") else NA_character_,
                                           low_confidence = if (nrow(pr)) any(pr$hconf == "0") else NA)
  }
  gene_tab <- rbindlist(rows)
  
  set_rows <- list(); rule_rows <- list()
  for (sn in names(redox_sets_mouse)) {
    pm <- pairs[mname %in% redox_sets_mouse[[sn]]]
    strict <- pm[htype == "ortholog_one2one", .(mouse_genes = paste(sort(unique(mname)), collapse = ";")), by = .(human_symbol = hname)]
    all_h  <- pm[, .(mouse_genes = paste(sort(unique(mname)), collapse = ";"), homology_types = paste(sort(unique(htype)), collapse = ";")), by = .(human_symbol = hname)]
    rule_rows[[length(rule_rows) + 1]] <- rbind(cbind(set = sn, rule = "strict", strict), cbind(set = sn, rule = "all", all_h), fill = TRUE)
    set_rows[[length(set_rows) + 1]] <- data.table(
      set = sn, n_mouse_set = length(redox_sets_mouse[[sn]]), n_prof_human_set = length(redox_sets_human[[sn]]),
      n_ensembl_strict = nrow(strict), n_ensembl_all = nrow(all_h),
      prof_only_vs_all = paste(sort(setdiff(redox_sets_human[[sn]], all_h$human_symbol)), collapse = ";"),
      ensembl_all_only = paste(sort(setdiff(all_h$human_symbol, redox_sets_human[[sn]])), collapse = ";"),
      prof_not_in_strict = paste(sort(setdiff(redox_sets_human[[sn]], strict$human_symbol)), collapse = ";"))
  }
  set_tab <- rbindlist(set_rows); rule_tab <- rbindlist(rule_rows, fill = TRUE)
  
  summ <- gene_tab[, .N, by = "ensembl_code"]; setorderv(summ, "ensembl_code")
  cat("  Ensembl reason codes of the mouse set genes:\n"); print(summ)
  cat("  Set-size comparison (co-author's human set / Ensembl strict / Ensembl all):\n")
  print(set_tab[, c("set", "n_prof_human_set", "n_ensembl_strict", "n_ensembl_all"), with = FALSE])
  
  setnames(gene_tab, c("set", "mouse_gene", "deg_status", "responsive", "ensembl_code", "human_ortholog", "homology_type", "low_confidence"),
           c("set", "mouse_gene", "test_status", "responsive", "Ensembl_reason_code", "human_counterpart", "homology_type", "low_confidence"))
  setnames(set_tab, c("set", "n_mouse_set", "n_prof_human_set", "n_ensembl_strict", "n_ensembl_all", "prof_only_vs_all", "ensembl_all_only", "prof_not_in_strict"),
           c("set", "n_mouse_set", "n_coauthor_human_set", "n_Ensembl_strict", "n_Ensembl_all", "coauthor_only_vs_all", "Ensembl_all_only", "coauthor_symbols_not_in_strict"))
  rule_strict <- rule_tab[rule_tab$rule == "strict"]; rule_all <- rule_tab[rule_tab$rule == "all"]
  setnames(rule_tab, c("set", "rule", "human_symbol", "mouse_genes", "homology_types"), c("set", "rule", "human_symbol", "mouse_genes", "homology_types"))
  
  hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", halign = "center", border = "Bottom")
  add_sheet <- function(wb, name, df) {
    addWorksheet(wb, name); df <- as.data.frame(df)
    writeData(wb, name, df, headerStyle = hs); freezePane(wb, name, firstRow = TRUE)
    if (ncol(df)) setColWidths(wb, name, cols = seq_len(ncol(df)),
                               widths = vapply(seq_len(ncol(df)), function(j) min(max(nchar(c(names(df)[j], as.character(head(df[[j]], 200)))), na.rm = TRUE) + 2, 45), numeric(1)))
  }
  wb <- createWorkbook()
  add_sheet(wb, "README", data.table(item = c("1_per_gene", "2_set_comparison", "3_rule_based_human_sets", "note"),
                                     value = c("per mouse set gene: test status, Ensembl reason code (M0-M3), human counterpart, homology type",
                                           "difference between the co-author's human sets and the Ensembl-rule sets (strict, all)",
                                           "human sets built with the Ensembl rules (strict: one2one only, all: every orthologue)",
                                           "the rule for the human sets was agreed with the co-author and then written into redox_sets_human at the top of this file")))
  add_sheet(wb, "1_per_gene", gene_tab); add_sheet(wb, "2_set_comparison", set_tab); add_sheet(wb, "3_rule_based_human_sets", rule_tab)
  saveWorkbook(wb, out_path("redox_symbol_check.xlsx"), overwrite = TRUE)
  data.table::fwrite(rule_strict, out_path("redox_human_sets_strict.csv"))
  data.table::fwrite(rule_all,    out_path("redox_human_sets_all.csv"))
  cat("  Saved:", out_path("redox_symbol_check.xlsx"), "\n")
} else if (RUN_STEP1) cat("\n[STEP 1] Ensembl file not found; symbol check skipped:", ENSEMBL_FILE, "\n")


# ==============================================================================
# STEP 2  Set analysis (Fisher, cameraPR)
# ==============================================================================
if (RUN_STEP2) {
  cat("\n[STEP 2] Redox set analysis\n")
  
  # ---------------------------------------------------------------- helpers --
  prep <- function(path) {
    d <- read_csv(path, show_col_types = FALSE)
    need <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", if (STAT_INPUT == "stat") "stat")
    miss <- setdiff(need, names(d))
    if (length(miss)) stop(path, ": required columns missing: ", paste(miss, collapse = ", "))
    d <- d %>% filter(!is.na(padj))
    if (STAT_INPUT == "stat" && anyNA(d$stat)) stop(path, ": some genes with a padj have a missing stat value")
    d$signed <- -log10(d$pvalue) * sign(d$log2FoldChange)
    d$responsive <- d$padj < PADJ_CUT & abs(d$log2FoldChange) >= LFC_CUT
    d[!duplicated(d$gene), ]
  }
  
  stat_vector <- function(d) if (STAT_INPUT == "stat") d$stat else d$signed
  
  # One-sided Fisher over-representation, with the 2x2 table saved.
  ora <- function(d, set) {
    inset <- d$gene %in% set
    a <- sum(inset &  d$responsive); b <- sum(inset & !d$responsive)
    cc<- sum(!inset &  d$responsive); dd<- sum(!inset & !d$responsive)
    tab <- matrix(c(a, b, cc, dd), nrow = 2, byrow = TRUE,
                  dimnames = list(c("in_set","out_set"), c("responsive","not")))
    ft <- fisher.test(tab, alternative = "greater")
    list(table = tab, n_tested = a + b, n_responsive = a,
         expected = (a + b) * sum(d$responsive) / nrow(d), p = ft$p.value)
  }
  
  run_sets <- function(d, sets, label, bonf_n = NULL) {
    idx <- lapply(sets, function(s) which(d$gene %in% s))
    cam <- cameraPR(statistic = stat_vector(d), index = idx,
                    use.ranks = USE_RANKS, inter.gene.cor = RHO_MAIN, sort = FALSE)
    o <- lapply(sets, function(s) ora(d, s))
    res <- tibble(
      dataset     = label,
      set         = names(sets),
      n_tested    = vapply(o, function(x) x$n_tested, 1L),
      n_responsive= vapply(o, function(x) x$n_responsive, 1L),
      expected    = round(vapply(o, function(x) x$expected, 1), 2),
      fold        = round(vapply(o, function(x) x$n_responsive / x$expected, 1), 2),
      n_resp_up   = vapply(sets, function(s) sum(d$gene %in% s & d$responsive & d$log2FoldChange > 0), 1L),
      n_resp_down = vapply(sets, function(s) sum(d$gene %in% s & d$responsive & d$log2FoldChange < 0), 1L),
      fisher_p    = vapply(o, function(x) x$p, 1),
      camera_dir  = cam$Direction,
      camera_p    = cam$PValue
    ) %>%
      mutate(fisher_bonf = if (is.null(bonf_n)) NA_real_ else pmin(fisher_p * bonf_n, 1),
             camera_bonf = if (is.null(bonf_n)) NA_real_ else pmin(camera_p * bonf_n, 1))
    attr(res, "tables") <- lapply(o, function(x) x$table)
    res
  }
  
  camera_rho_sweep <- function(d, set, label) {
    idx <- list(set = which(d$gene %in% set))
    bind_rows(lapply(c(RHO_MAIN, RHO_SENS), function(r) {
      cm <- cameraPR(stat_vector(d), idx, use.ranks = USE_RANKS,
                     inter.gene.cor = r, sort = FALSE)
      tibble(dataset = label, rho = r, direction = cm$Direction, p = cm$PValue)
    }))
  }
  
  # Values per set gene (gene-level output of the supplementary analysis). Untested genes are kept as rows
  gene_level <- function(d, sets, label) {
    bind_rows(lapply(names(sets), function(sn) {
      g <- sets[[sn]]; m <- d[match(g, d$gene), ]
      tibble(dataset = label, set = sn, gene = g, tested = g %in% d$gene,
             log2FC = m$log2FoldChange, padj = m$padj, stat = if ("stat" %in% names(d)) m$stat else NA_real_,
             responsive = ifelse(is.na(m$padj), NA, m$responsive),
             direction = ifelse(is.na(m$padj), "not tested",
                                ifelse(m$responsive, ifelse(m$log2FoldChange > 0, "responsive up", "responsive down"),
                                       ifelse(m$log2FoldChange > 0, "up (below criteria)", "down (below criteria)"))))
    }))
  }
  
  # ------------------------------------------------------------------- run ---
  mouse <- prep(MOUSE_FILE)
  cat("cameraPR input statistic (STAT_INPUT):", STAT_INPUT, "| use.ranks:", USE_RANKS, "| inter.gene.cor:", RHO_MAIN, "\n")
  cat("mouse: tested", nrow(mouse), " responsive", sum(mouse$responsive), "\n")
  have_human <- file.exists(HUMAN_WT_FILE) && file.exists(HUMAN_KO_FILE)
  if (have_human) {
    hu_wt <- prep(HUMAN_WT_FILE); hu_ko <- prep(HUMAN_KO_FILE)
    cat("human WT: tested", nrow(hu_wt), " responsive", sum(hu_wt$responsive), "\n")
    cat("human KO: tested", nrow(hu_ko), " responsive", sum(hu_ko$responsive), "\n")
  } else cat("[note] Human input files not found; the human analysis (WT, KO, common 18, direction concordance) is skipped:", HUMAN_WT_FILE, "/", HUMAN_KO_FILE, "\n")
  
  # Bonferroni is applied only to the table in which the six sets are tested together (for the sensitivity and common-18 tables only the raw P is interpreted)
  res_mouse    <- run_sets(mouse, redox_sets_mouse, "mouse", bonf_n = length(redox_sets_mouse))
  res_gsh_sens <- run_sets(mouse, gsh_conditions, "mouse_glutathione_sensitivity")
  rho_mouse    <- camera_rho_sweep(mouse, gsh, "mouse")
  
  # Repeat the same mouse analysis with the other statistic (to check how the choice for the main analysis affects the conclusions)
  STAT_MAIN <- STAT_INPUT
  STAT_ALT  <- setdiff(c("stat", "signed"), STAT_MAIN)
  STAT_INPUT <- STAT_ALT
  alt_six <- run_sets(mouse, redox_sets_mouse, paste0("mouse_alt_", STAT_ALT), bonf_n = length(redox_sets_mouse))
  alt_gsh <- run_sets(mouse, gsh_conditions, paste0("mouse_glutathione_sensitivity_alt_", STAT_ALT))
  STAT_INPUT <- STAT_MAIN
  stat_sens <- bind_rows(res_mouse, alt_six, res_gsh_sens, alt_gsh) %>%
    select(dataset, set, n_tested, n_responsive, fisher_p, camera_dir, camera_p, camera_bonf)
  
  out_list <- list(res_mouse, res_gsh_sens)
  tab_list <- list(mouse = attr(res_mouse, "tables"), sensitivity = attr(res_gsh_sens, "tables"))
  rho_list <- list(rho_mouse)
  concordance <- NULL; shift_p <- NA_real_; conc <- NULL
  
  if (have_human) {
    res_huWT <- run_sets(hu_wt, redox_sets_human, "human_WT", bonf_n = length(redox_sets_human))
    res_huKO <- run_sets(hu_ko, redox_sets_human, "human_KO", bonf_n = length(redox_sets_human))
    # Human over-representation restricted to the 18 genes testable in both species.
    res_hu18 <- run_sets(hu_wt, list("Glutathione common 18" = common_18), "human_WT_18")
    rho_huWT <- camera_rho_sweep(hu_wt, redox_sets_human[["Glutathione-related"]], "human_WT")
    rho_hu18 <- camera_rho_sweep(hu_wt, common_18, "human_WT_18")
    out_list <- list(res_mouse, res_huWT, res_huKO, res_hu18, res_gsh_sens)
    tab_list <- c(tab_list, list(human_WT = attr(res_huWT, "tables"), human_KO = attr(res_huKO, "tables"), human_18 = attr(res_hu18, "tables")))
    rho_list <- list(rho_mouse, rho_huWT, rho_hu18)
    
    # ----------------------------------- mouse to human direction concordance --
    mo <- mouse %>% mutate(KEY = toupper(gene)) %>% select(KEY, mouse_lfc = log2FoldChange)
    hu <- hu_wt  %>% mutate(KEY = gene) %>% select(KEY, human_lfc = log2FoldChange,
                                                   human_padj = padj, human_base = baseMean)
    conc <- inner_join(mo, hu, by = "KEY") %>% filter(KEY %in% common_18) %>%
      mutate(concordant = sign(mouse_lfc) == sign(human_lfc))
    
    k <- sum(conc$concordant); n <- nrow(conc)
    bg <- inner_join(mo, hu, by = "KEY") %>%
      summarise(rate = mean(sign(mouse_lfc) == sign(human_lfc))) %>% pull(rate)
    
    concordance <- tibble(
      n_genes          = n,
      n_concordant     = k,
      p_vs_half        = binom.test(k, n, 0.5, alternative = "greater")$p.value,
      background_rate  = round(bg, 3),
      p_vs_background  = binom.test(k, n, bg, alternative = "greater")$p.value,
      median_mouse_lfc = median(conc$mouse_lfc),
      median_human_lfc = median(conc$human_lfc),
      median_abs_mouse = median(abs(conc$mouse_lfc)),
      median_abs_human = median(abs(conc$human_lfc)),
      n_human_abs_ge_1 = sum(abs(conc$human_lfc) >= 1)
    )
    
    # Human set median shift against the rest of the transcriptome, one sided.
    shift_p <- wilcox.test(hu_wt$log2FoldChange[hu_wt$gene %in% common_18],
                           hu_wt$log2FoldChange[!hu_wt$gene %in% common_18],
                           alternative = "less")$p.value
  }
  
  # ---------------------------------------------------------------- save ---
  ds_list <- list(tibble(dataset = "mouse", n_tested = nrow(mouse), n_responsive = sum(mouse$responsive),
                         n_up = sum(mouse$responsive & mouse$log2FoldChange > 0), n_down = sum(mouse$responsive & mouse$log2FoldChange < 0)))
  gl_list <- list(gene_level(mouse, redox_sets_mouse, "mouse"))
  if (have_human) {
    ds_list <- list(ds_list[[1]],
                    tibble(dataset = "human_WT", n_tested = nrow(hu_wt), n_responsive = sum(hu_wt$responsive), n_up = sum(hu_wt$responsive & hu_wt$log2FoldChange > 0), n_down = sum(hu_wt$responsive & hu_wt$log2FoldChange < 0)),
                    tibble(dataset = "human_KO", n_tested = nrow(hu_ko), n_responsive = sum(hu_ko$responsive), n_up = sum(hu_ko$responsive & hu_ko$log2FoldChange > 0), n_down = sum(hu_ko$responsive & hu_ko$log2FoldChange < 0)))
    gl_list <- c(gl_list, list(gene_level(hu_wt, redox_sets_human, "human_WT"), gene_level(hu_ko, redox_sets_human, "human_KO")))
  }
  write_csv(bind_rows(ds_list), out_path("S9_dataset_summary.csv"))
  write_csv(bind_rows(gl_list), out_path("S9_gene_level_all_sets.csv"))
  write_csv(bind_rows(lapply(names(gsh_conditions), function(cn) tibble(condition = cn, gene = gsh, included = gsh %in% gsh_conditions[[cn]], tested = gsh %in% mouse$gene))),
            out_path("S9_glutathione_conditions_membership.csv"))
  write_csv(tibble(item = c("Gstm representative (prespecified rule)", "Gstm family members in the set", "Genes removed as transferases (^Gst|^Mgst)"),
                   value = c(GSTM_REPRESENTATIVE, paste(gstm_fam, collapse = ";"), paste(transferase, collapse = ";"))),
            out_path("S9_glutathione_condition_rules.csv"))
  out <- bind_rows(out_list)
  write_csv(out,        out_path("S9_gene_set_results.csv"))
  write_csv(bind_rows(rho_list), out_path("S9_camera_rho_sensitivity.csv"))
  write_csv(stat_sens,  out_path("S9_camera_stat_sensitivity.csv"))
  if (have_human) {
    write_csv(conc,     out_path("S9_mouse_human_concordance_genes.csv"))
    write_csv(concordance %>% mutate(wilcox_one_sided_p = shift_p), out_path("S9_mouse_human_concordance_summary.csv"))
  }
  saveRDS(tab_list, out_path("S9_fisher_2x2_tables.rds"))
  writeLines(capture.output(sessionInfo()), out_path("S9_sessionInfo.txt"))
  
  print(out, n = 50)
  if (have_human) { print(concordance); cat("one-sided Wilcoxon for human set shift:", shift_p, "\n\n") }
  cat("Output folder:", OUT_DIR, "\n")
}

# =============================================================================
# Values from the provisional re-implementation, for comparison only.
# Differences are expected; identify the cause before choosing a final value.
#   mouse Glutathione-related      Fisher 2.07e-04  Bonf 1.24e-03  cameraPR 1.8e-06
#   mouse Microsomal P450          Fisher 5.33e-07  Bonf 3.20e-06  cameraPR 0.096
#   mouse NRF2                     Fisher 0.145     Bonf 0.872     cameraPR 0.806
#   mouse Lipid peroxidation       Fisher 0.165     Bonf 0.988     cameraPR 0.943
#   mouse Iron homeostasis         Fisher 0.590     Bonf 1.00      cameraPR 0.126
#   mouse Fatty acid beta-oxid.    Fisher 1.00      Bonf 1.00      cameraPR 0.390
#   mouse glutathione sensitivity  full 2.07e-04 / Gstm collapsed 3.56e-04 /
#                                  Gstm removed 2.18e-03 / all transferases removed 0.374
#   human WT glutathione (26)      Fisher 0.752     cameraPR 0.337
#   human WT common 18             Fisher 1.00  (0 responsive)
#   concordance                    14 of 18, p vs 0.5 = 0.0154, background 0.629,
#                                  p vs background = 0.143, Wilcoxon one sided = 0.0108
#   median abs log2FC              mouse 0.695, human 0.327
#   (Note) The cameraPR values above were obtained with the signed score (STAT_INPUT = "signed"); the mouse signed values are reproduced by this code.
# =============================================================================