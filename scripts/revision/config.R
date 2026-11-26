# =============================================================================
# config.R  -  edit this file once; scripts 02-07, 10, 12 and 13 read it with source().
# (Script 01 only needs DEG_hepT3_raw.csv, redox_gene_sets_human.csv and,
#  optionally, mouse_T3_results.csv; its paths are set at the top of the script.)
# Run every script from this folder (scripts/revision): all paths are relative to it.
# =============================================================================

# ---- folder with the input files (counts, metadata, gene sets) ---------------
data_dir <- "data"
out_root <- "results"

# ---- gene sets (exported from Table S11; symbols are trimmed in the scripts) --
human_sets_file <- file.path(data_dir, "redox_gene_sets_human.csv")  # set, human_symbol
mouse_sets_file <- file.path(data_dir, "redox_gene_sets_mouse.csv")  # set, mouse_symbol

# ---- 02_biopsy_stage_analysis.R ---------------------------------------------------
# Raw count matrix: first column = gene symbol (HGNC), remaining columns = samples
# (column names must match the 'sample' column of the metadata file).
cohorts <- list(
  GSE135251 = list(counts = file.path(data_dir, "GSE135251_counts.csv"),
                   meta   = file.path(data_dir, "GSE135251_meta.csv")),
  GSE130970 = list(counts = file.path(data_dir, "GSE130970_counts.csv"),
                   meta   = file.path(data_dir, "GSE130970_meta.csv"))
)
# Metadata column names. Set a name to NA when the information does not exist
# (for example age and sex are not deposited for GSE135251).
meta_cols <- list(
  sample       = "sample",
  stage        = "fibrosis_stage",         # integer 0-4
  group        = "group",                  # "NAFLD" for patients (controls are dropped); NA = keep all
  group_keep   = "NAFLD",
  nas          = "nas",
  steatosis    = "steatosis",
  inflammation = "lobular_inflammation",
  ballooning   = "ballooning",
  age          = "age",
  sex          = "sex"
)
candidate_genes <- c("ABHD1", "F13B", "ADH4", "LGALS4")
# fibrosis marker genes used as positive controls
control_genes   <- c("COL1A1", "COL1A2", "COL3A1", "TIMP1", "ACTA2", "LUM", "PDGFRB", "THY1")
hepatocyte_markers <- c("ALB", "APOA1", "APOC3", "TTR", "SERPINA1", "AHSG", "FGB", "APOH")
min_count <- 10      # gene kept when count >= min_count in at least half of the samples
min_n_stratum <- 12  # minimum number of samples for the stage-stratified models

# ---- 03_mouse_T3_robustness.R -------------------------------------------------------
mouse_counts_file <- file.path(data_dir, "GSE85793_counts.csv")   # gene symbol (mouse) x samples
mouse_meta_file   <- file.path(data_dir, "GSE85793_meta.csv")     # columns: sample, group
mouse_group_col   <- "group"
mouse_treat_label <- "T3"
mouse_ctrl_label  <- "vehicle"
mouse_filter_min_count <- 10
mouse_padj_cut <- 0.05
mouse_lfc_cut  <- 1

# =============================================================================
# Scripts 04-07, 10 and 12
# =============================================================================
# ---- 04_human_final_tables.R ---------------------------------------------------
hep_deg_file    <- file.path(data_dir, "DEG_hepT3_raw.csv")      # DESeq2 output of GSE270233
mouse_deg_file  <- file.path(data_dir, "mouse_T3_results.csv")   # mouse T3 result with human orthologs
# optional: raw counts of GSE270233 for the maturity table (gene symbol x sample) + metadata
hep_counts_file <- file.path(data_dir, "GSE270233_counts.csv")
hep_meta_file   <- file.path(data_dir, "GSE270233_meta.csv")     # columns: sample, genotype, treatment
maturity_genes  <- c("AFP", "DLK1", "ALB", "ASGR1", "CYP3A4", "CYP2E1", "CYP1A2", "UGT1A1", "GSTM1")
# 26-gene panel of the submitted analysis. Script 12 recomputes the panel-size dependent statistics for the 22-gene
# panel of the IVW-centred Mendelian randomization definition (four genes removed).
panel_genes <- c("CUX2","BDH2","ADH4","LILRB3","SBK1","ABHD1","INCENP","CDC42EP5","CD19","GPCPD1",   # 10 positional
                 "BTC","CFHR1","CLIC5","COL11A2","COL15A1","DEF6","DNAJB4","F13B","LGALS4","PKLR",   # 16 protein-genetic
                 "POR","SPARCL1","TIMD4","GSTM1","HSPA1A",                                            # (ADH4 is shared)
                 "SLC39A8")                                                                            # colocalization lead gene

# ---- 06_biopsy_F4_exclusion.R (uses cohorts / meta_cols / hepatocyte_markers above) --------

# ---- 07_coloc_nonUKB_sensitivity.R -------------------------------------------------------
s2_file <- file.path(data_dir, "Antioxidant_ThyroidLiver_Supplementary_Tables_S1-S13.xlsx")       # Supplementary Tables workbook with sheet "Table S2" (window summary + pair table)

# ---- 10_TRb_ChIP_annotation.R ------------------------------------------------------------------
# Ayers et al. 2014 (PLoS ONE) HepG2 TR-beta BioChIP peaks, Supplementary Table S5 saved as CSV: chr,start,end
trb_peak_file <- file.path(data_dir, "Ayers2014_TRb_peaks.csv")
trb_peak_build <- "hg19"                                         # build of the peak coordinates
# TSS table of the same build: gene,chr,tss,strand (+/-); created by 10a_make_tss_hg19.R
tss_file <- file.path(data_dir, "tss_hg19.csv")

# ---- 13_recount_MR_support_IVW.R ---------------------------------------------------------------
# Revised Supplementary Tables workbook (sheets "Table S7", "Table S14", "Table S15", "Table S16"); the environment variable
# SUPP_XLSX overrides this path
supp_revised_file <- file.path(data_dir, "Antioxidant_ThyroidLiver_Supplementary_Tables_revised.xlsx")
