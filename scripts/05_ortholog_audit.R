################################################################################
#  05_ortholog_audit.R
#  (Original: 261021_ortholog_audit_biomaRt.R)
#  [1] Reads the Ensembl BioMart human -> mouse orthologue table (a) from a file downloaded by hand from the web or (b) downloads it with R,
#  [2] audits the mouse orthologues of the human genes in the COLOC-passing loci and evaluates their hepatic T3 response,
#  [3] saves the values specified for the download and everything produced during the run to files.
#
#  AUDIT STATES (order of the decision)
#     1. identifier conversion failed : BioMart has no row for this symbol (old names, clone names, ...)
#     2. no mouse ortholog            : the gene is in BioMart but has no mouse counterpart
#     3. ambiguous mapping            : not one-to-one (one2many, many2many, two or more mouse genes)
#                                       -> excluded from the primary analysis, but the result of evaluating all counterparts is kept in the permissive columns
#     4. absent from the mouse count matrix : the mouse gene name is not in the raw count matrix
#     5. excluded for low expression  : in the count matrix but missing from the DESeq2 result (pre-filter / independent filtering / Cook's outlier)
#     6. tested normally              : a padj exists. Responsive (padj < 0.05, |log2FC| >= 1) / below the response cut / not responsive
#
#  INPUT FILES (DATA_DIR)
#     MAPPING_XLSX  : 261019_gene_locus_full_mapping.xlsx  (sheet 3_Gene_in_locus; gene-window mapping table)
#     SYMBOLS_FILE  : genes_in_pass_loci.txt  (human symbols sent to BioMart; generated from MAPPING_XLSX if absent)
#     DEG_ALL       : 250912__DEG_all_genes.csv  (all DESeq2 results of 03_deseq2_mouse_liver_T3.R, 12,135 genes; also called 2_DEG_all_genes.csv)
#     COUNT_MATRIX  : 0_count_matrix.csv   (raw count matrix before the pre-filter; columns C1..C4, T1..T4 are used)
#     HUMAN_RESP_FILE (optional) : 5_LTl1_human_unique_genes.tsv  (earlier human list of responsive genes, for QC)
#     CROSSCHECK_MOUSE_FILE (optional) : 261019_ensembl_GRCm39.txt  (mouse-direction table; if present the two directions are compared)
#
#  OUTPUTS (in OUT_DIR; file names start with OUT_PREFIX)
#     _result.xlsx           : 0_download_record / 0b_settings / 0c_input_files / 1_summary / 1b_detail / 2_audit_table /
#                              3_candidates_responsive / 4_ambiguous_detail / 5_QC / 6_BioMart_table / 7_direction_comparison (optional)
#     _table.tsv             : the same content as 2_audit_table
#     _biomart_raw.tsv       : raw table received from BioMart (column names in the form recognised by the script)
#     _biomart_download_record.tsv : the values specified for the download (release, host, dataset, filter, attributes, date, ...) and the number of rows
#     _input_symbols.txt     : list of symbols sent to BioMart
#     _symbols_not_returned.txt : symbols for which BioMart returned no row
#     _input_files_manifest.tsv : path, size, modification time and md5 of the input files
#     _settings.tsv          : all settings of this run
#     _run_log.txt           : all screen output, warnings and messages
#     _sessionInfo.txt
#
#  Supplementary Table S6b is built from 2_audit_table (protein-coding genes). Three manual edits were made in the table: the aliases SPX/C12orf39 and
#  MYMK/TMEM8C are listed under both symbols and counted once, and four genes (ZNF358, ZNF512, ZNF513, ZNF687) that were tested under the mouse
#  symbols Zfp358, Zfp512, Zfp513 and Zfp687 are annotated; the audit rule itself matches human and mouse symbols directly.
################################################################################
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })

# ==============================================================================
# [SECTION 1]  Settings
# ==============================================================================
DATA_DIR     <- Sys.getenv("DATA_DIR", "data")
MAPPING_XLSX <- file.path(DATA_DIR, "261019_gene_locus_full_mapping.xlsx")
SHEET_GENES  <- "3_Gene_in_locus"
SCOPE        <- "pass"            # "pass": genes of the COLOC-passing loci / "all": genes of all 458 loci
SYMBOLS_FILE <- file.path(DATA_DIR, "genes_in_pass_loci.txt")   # NA, or a missing file: generated from MAPPING_XLSX according to SCOPE

# ---- BioMart download settings (kept in the record)
DOWNLOAD        <- FALSE          # TRUE: download with R / FALSE: use the file downloaded by hand from the web BioMart (BIOMART_MANUAL_FILE) or an already saved _biomart_raw.tsv
BIOMART_MANUAL_FILE  <- file.path(DATA_DIR, "mart_export.txt")   # path of the file downloaded from the web BioMart (used when DOWNLOAD <- FALSE)
MANUAL_DATASET_NOTE  <- "Ensembl Genes 116 / Human genes (GRCh38.p14)"   # Database / Dataset chosen on the web (for the record)
MANUAL_DOWNLOAD_DATE <- NA                  # download date, e.g. "2026-10-21". If NA the file modification time is recorded
DOWNLOAD_METHOD <- "rest"         # "rest": query BioMart directly at the release archive address (recommended; does not use the archive listing of biomaRt)
                                  # "biomaRt": use the biomaRt package (only possible while the former ensembl.org mirrors are open)
ENSEMBL_RELEASE <- 116            # fixed to the same release as the DEG remapping
ENSEMBL_HOST    <- "https://jun2026.archive.ensembl.org"   # archive address of release 116 (June 2026). The former www/asia/useast mirrors were shut down in July 2026
BM_RETRIES      <- 3              # number of retries after a connection failure
BM_BIOMART      <- "genes"
BM_DATASET      <- "hsapiens_gene_ensembl"      # human dataset
BM_FILTER       <- "hgnc_symbol"                # type of the input symbols
BM_BATCH        <- 200                          # number of symbols per request
BM_ATTRIBUTES   <- c("ensembl_gene_id", "external_gene_name",
                     "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name",
                     "mmusculus_homolog_orthology_type", "mmusculus_homolog_orthology_confidence")
BM_HEADERS      <- c("Gene stable ID", "Gene name", "Mouse gene stable ID", "Mouse gene name",
                     "Mouse homology type", "Mouse orthology confidence [0 low, 1 high]")  # in the same order as the attributes above

# ---- audit inputs
DEG_ALL       <- file.path(DATA_DIR, "250912__DEG_all_genes.csv")   # = 2_DEG_all_genes.csv of the analysis folder (same 12,135 genes); output of 03_deseq2_mouse_liver_T3.R
COUNT_MATRIX  <- file.path(DATA_DIR, "0_count_matrix.csv")
DESEQ_SAMPLES <- c("C1", "C2", "C3", "C4", "T1", "T2", "T3", "T4")   # sample columns that entered DESeq2
MIN_COUNT     <- 10
MIN_SAMPLES   <- ceiling(length(DESEQ_SAMPLES) / 2)
PADJ_CUT <- 0.05
LFC_CUT  <- 1
HUMAN_RESP_FILE       <- file.path(DATA_DIR, "5_LTl1_human_unique_genes.tsv")   # NA if absent
CROSSCHECK_MOUSE_FILE <- file.path(DATA_DIR, "261019_ensembl_GRCm39.txt")      # NA if absent

OUT_DIR    <- Sys.getenv("OUT_DIR", "results/05_ortholog_audit")
OUT_PREFIX <- "261021_ortholog_audit"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
P <- function(suffix) file.path(OUT_DIR, paste0(OUT_PREFIX, suffix))
RUN_START <- Sys.time()

# ==============================================================================
# [SECTION 2]  Log: all screen output, warnings and messages are saved to _run_log.txt
# ==============================================================================
con_log <- file(P("_run_log.txt"), open = "wt", encoding = "UTF-8")
sink(con_log, split = TRUE)
status <- tryCatch({
withCallingHandlers({

cat(sprintf("Start: %s\n", format(RUN_START, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Working directory: %s\n", getwd()))

# ==============================================================================
# [SECTION 3]  Helpers
# ==============================================================================
up       <- function(x) toupper(trimws(as.character(x)))
blank_na <- function(x) { x <- trimws(as.character(x)); x[x == "" | x == "NA"] <- NA_character_; x }
hdr_style <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", halign = "center", border = "Bottom")
add_sheet <- function(wb, name, df) {
  addWorksheet(wb, name); df <- as.data.frame(df)
  writeData(wb, name, df, headerStyle = hdr_style); freezePane(wb, name, firstRow = TRUE)
  if (ncol(df)) setColWidths(wb, name, cols = seq_len(ncol(df)),
                             widths = vapply(seq_len(ncol(df)), function(j) min(max(nchar(c(names(df)[j], as.character(head(df[[j]], 200)))), na.rm = TRUE) + 2, 45), numeric(1)))
}
pick_col <- function(nm, cands) { i <- match(tolower(cands), tolower(nm)); i <- i[!is.na(i)]; if (length(i)) nm[i[1]] else NA_character_ }

# ==============================================================================
# [SECTION 4]  Read the target genes
# ==============================================================================
cat("\n[1] Reading the target genes\n")
g3 <- as.data.table(read.xlsx(MAPPING_XLSX, sheet = SHEET_GENES))
if (SCOPE == "pass") g3 <- g3[locus_pass == TRUE]
genes <- g3[, .(biotype = biotype[1], loci = paste(sort(unique(Locus_ID)), collapse = ";"),
                tier = paste(sort(unique(locus_tier_manuscript)), collapse = ";"),
                min_dist_kb = min(dist_kb), GWAS_liver = max(GWAS_liver), GWAS_thyroid = max(GWAS_thyroid),
                GWAS_shared = max(GWAS_shared)), by = Gene]
genes[, gkey := up(Gene)]
cat(sprintf("   target genes %d (scope: %s) | coding %d, noncoding %d\n", nrow(genes), SCOPE,
            sum(genes$biotype == "coding"), sum(genes$biotype != "coding")))

if (!is.na(SYMBOLS_FILE) && file.exists(SYMBOLS_FILE)) {
  symbols <- unique(trimws(readLines(SYMBOLS_FILE, warn = FALSE))); symbols <- symbols[nzchar(symbols)]
  cat(sprintf("   symbol file %s: %d symbols | identical to the target genes: %s\n", SYMBOLS_FILE, length(symbols), setequal(up(symbols), genes$gkey)))
  if (!setequal(up(symbols), genes$gkey)) warning("The symbols in SYMBOLS_FILE differ from the target genes in MAPPING_XLSX.")
} else { symbols <- genes$Gene; cat("   symbol list generated from MAPPING_XLSX\n") }
writeLines(symbols, P("_input_symbols.txt"))

# ==============================================================================
# [SECTION 5]  BioMart download (R, biomaRt, fixed release) + record
# ==============================================================================
# ---- (a) REST method: send the BioMart XML query directly to the archive address (only httr is needed)
download_rest <- function(symbols) {
  if (!requireNamespace("httr", quietly = TRUE)) stop("The httr package is required: install.packages('httr')")
  base <- sub("/+$", "", ENSEMBL_HOST); url <- paste0(base, "/biomart/martservice")
  rec <- list(method = "BioMart XML query (httr POST)", url = url, host = base); t0 <- Sys.time()
  http <- function(fun, ...) {                                   # with retries
    last <- NULL
    for (a in seq_len(BM_RETRIES)) {
      r <- tryCatch(fun(...), error = function(e) e)
      if (!inherits(r, "error") && httr::status_code(r) == 200) return(httr::content(r, "text", encoding = "UTF-8"))
      last <- if (inherits(r, "error")) conditionMessage(r) else paste("HTTP", httr::status_code(r))
      cat(sprintf("   connection failed (%d/%d): %s\n", a, BM_RETRIES, last)); if (a < BM_RETRIES) Sys.sleep(5 * a)
    }
    stop("Could not connect to the BioMart server (", url, "): ", last,
         "\n  Check in a browser that ", url, "?type=registry opens. On a company or university network a proxy setting (Sys.setenv(https_proxy=...)) may be needed.")
  }
  reg <- http(httr::GET, url, query = list(type = "registry"), httr::timeout(60))
  m <- regmatches(reg, regexpr("ENSEMBL GENES [0-9]+", reg))
  rec$registry <- if (length(m)) m else paste("release number not found in the registry:", substr(gsub("\\s+", " ", reg), 1, 120))
  if (length(m) && as.integer(sub(".*GENES ", "", m)) != ENSEMBL_RELEASE)
    stop(sprintf("The release of the server (%s) differs from the requested release (%s). Check ENSEMBL_HOST.", m, ENSEMBL_RELEASE))
  if (!length(m)) warning("Could not verify the release number of the server. Check that ENSEMBL_HOST is the archive address of release ", ENSEMBL_RELEASE, ".")
  ds <- tryCatch(http(httr::GET, url, query = list(type = "datasets", mart = "ENSEMBL_MART_ENSEMBL"), httr::timeout(60)), error = function(e) NA_character_)
  rec$dataset_info <- if (!is.na(ds)) { l <- grep(BM_DATASET, strsplit(ds, "\n")[[1]], value = TRUE); if (length(l)) gsub("\t+", " | ", l[1]) else NA_character_ } else NA_character_
  attr_xml <- paste0('<Attribute name="', BM_ATTRIBUTES, '"/>', collapse = "")
  chunks <- split(symbols, ceiling(seq_along(symbols) / BM_BATCH)); out <- list()
  for (i in seq_along(chunks)) {
    cat(sprintf("   BioMart request %d/%d (%d symbols)\n", i, length(chunks), length(chunks[[i]])))
    xml <- sprintf('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="%s" interface="default"><Filter name="%s" value="%s"/>%s</Dataset></Query>',
                   BM_DATASET, BM_FILTER, paste(chunks[[i]], collapse = ","), attr_xml)
    txt <- http(httr::POST, url, body = list(query = xml), encode = "form", httr::timeout(300))
    if (grepl("^\\s*Query ERROR|^\\s*Problems with the Query", txt)) stop("BioMart query error: ", substr(txt, 1, 300))
    out[[i]] <- if (nzchar(trimws(txt))) utils::read.delim(text = txt, header = FALSE, sep = "\t", quote = "", colClasses = "character", na.strings = "", fill = TRUE)
                else as.data.frame(matrix(character(), 0, length(BM_HEADERS)))
    if (ncol(out[[i]]) != length(BM_HEADERS)) stop("The number of columns of the downloaded table (", ncol(out[[i]]), ") differs from the number of attributes (", length(BM_HEADERS), ").")
    names(out[[i]]) <- BM_HEADERS
  }
  bm <- unique(rbindlist(out, fill = TRUE))
  rec$seconds <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  list(bm = as.data.frame(bm), record = rec)
}

# ---- (b) biomaRt method (only while the former ensembl.org mirrors are open)
download_biomaRt <- function(symbols) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) stop("The biomaRt package is required: BiocManager::install('biomaRt')")
  rec <- list(method = "biomaRt package", url = ENSEMBL_HOST, host = ENSEMBL_HOST); t0 <- Sys.time()
  mart <- biomaRt::useEnsembl(biomart = BM_BIOMART, dataset = BM_DATASET, version = ENSEMBL_RELEASE)
  rec$host <- tryCatch(mart@host, error = function(e) NA_character_)
  chunks <- split(symbols, ceiling(seq_along(symbols) / BM_BATCH)); out <- list()
  for (i in seq_along(chunks)) {
    cat(sprintf("   BioMart request %d/%d (%d symbols)\n", i, length(chunks), length(chunks[[i]])))
    out[[i]] <- biomaRt::getBM(attributes = BM_ATTRIBUTES, filters = BM_FILTER, values = chunks[[i]], mart = mart)
  }
  bm <- unique(rbindlist(out, fill = TRUE)); setnames(bm, BM_ATTRIBUTES, BM_HEADERS)
  rec$biomaRt_version <- as.character(packageVersion("biomaRt")); rec$seconds <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  list(bm = as.data.frame(bm), record = rec)
}

if (!exists("download_biomart")) {            # if a function with the same name has been defined beforehand, that one is used
  download_biomart <- function(symbols) if (DOWNLOAD_METHOD == "rest") download_rest(symbols) else download_biomaRt(symbols)
}

BIOMART_RAW_FILE <- P("_biomart_raw.tsv")
dl_record <- list()
cat("\n[2] BioMart table\n")
if (DOWNLOAD) {
  cat(sprintf("   downloading with R (%s): dataset=%s, release=%s, filter=%s, %d symbols\n", DOWNLOAD_METHOD, BM_DATASET, ENSEMBL_RELEASE, BM_FILTER, length(symbols)))
  dl <- download_biomart(symbols); bm_raw <- as.data.table(dl$bm); dl_record <- dl$record
  fwrite(bm_raw, BIOMART_RAW_FILE, sep = "\t", quote = FALSE, na = "")
  cat(sprintf("   rows received %d | saved: %s\n", nrow(bm_raw), BIOMART_RAW_FILE))
} else {
  src <- if (!is.na(BIOMART_MANUAL_FILE)) BIOMART_MANUAL_FILE else BIOMART_RAW_FILE
  if (!file.exists(src)) stop("BioMart table file not found: ", src, "\n  Set BIOMART_MANUAL_FILE to the path of the file downloaded from the web BioMart.")
  l1 <- readLines(src, n = 1, warn = FALSE); sp <- if (grepl("\t", l1, fixed = TRUE)) "\t" else ","
  bm_raw <- data.table::fread(src, sep = sp, quote = "\"", fill = TRUE, data.table = TRUE, colClasses = "character", na.strings = "")
  cat(sprintf("   using the manually downloaded file: %s (%d rows)\n", src, nrow(bm_raw)))
  if (!is.na(BIOMART_MANUAL_FILE) && normalizePath(src, mustWork = FALSE) != normalizePath(BIOMART_RAW_FILE, mustWork = FALSE))
    fwrite(bm_raw, BIOMART_RAW_FILE, sep = "\t", quote = FALSE, na = "")                 # save a copy under the standard name (for re-runs)
  dl_record <- list(method = "downloaded manually from the web BioMart", url = NA, host = NA,
                    registry = paste0("entered manually (not verified on the server): ", MANUAL_DATASET_NOTE), dataset_info = paste0("entered manually: ", MANUAL_DATASET_NOTE),
                    manual_file = src, manual_md5 = unname(tools::md5sum(src)),
                    manual_date = if (!is.na(MANUAL_DOWNLOAD_DATE)) MANUAL_DOWNLOAD_DATE else format(file.mtime(src), "%Y-%m-%d (file modification time)"))
}
hcol <- pick_col(names(bm_raw), c("Gene name", "HGNC symbol", "Human gene name"))
mcol <- pick_col(names(bm_raw), c("Mouse gene name", "Mouse Gene name"))
tcol <- pick_col(names(bm_raw), c("Mouse homology type"))
if (anyNA(c(hcol, mcol, tcol))) { cat("   Columns of the BioMart file:\n"); print(names(bm_raw))
  stop("Required columns not found. Check that Gene stable ID, Gene name, Mouse gene stable ID, Mouse gene name, Mouse homology type and Mouse orthology confidence were selected under Attributes (Homologues) on the web BioMart.") }
bm_hsym <- unique(up(bm_raw[[hcol]]))
not_returned <- symbols[!(up(symbols) %in% bm_hsym)]
writeLines(not_returned, P("_symbols_not_returned.txt"))
g <- function(x) if (length(x)) as.character(x)[1] else NA_character_
dl_tab <- data.table(
  item = c("run mode", "download method", "Ensembl release (requested)", "server address used", "BioMart query address", "release confirmed in the server registry", "dataset", "dataset information (server response)",
          "manual download file", "manual file md5", "download date", "biomart", "filter", "attributes", "saved column names", "batch size", "number of retries", "number of input symbols", "input symbol file", "rows received (after de-duplication)",
          "symbols present in BioMart", "symbols with no row received", "download time (s)", "biomaRt version (biomaRt method only)", "R version", "record time"),
  value = c(if (DOWNLOAD) "downloaded with R" else "manually downloaded file used", g(dl_record$method), ENSEMBL_RELEASE, g(dl_record$host), g(dl_record$url), g(dl_record$registry), BM_DATASET, g(dl_record$dataset_info),
        g(dl_record$manual_file), g(dl_record$manual_md5), g(dl_record$manual_date), BM_BIOMART, BM_FILTER, paste(BM_ATTRIBUTES, collapse = ", "), paste(BM_HEADERS, collapse = " | "), BM_BATCH, BM_RETRIES, length(symbols),
        if (!is.na(SYMBOLS_FILE)) SYMBOLS_FILE else "MAPPING_XLSX (automatic)", nrow(unique(bm_raw)), length(bm_hsym), length(not_returned),
        g(dl_record$seconds), g(dl_record$biomaRt_version), R.version.string, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
dl_tab[, value := as.character(value)]
fwrite(dl_tab, P("_biomart_download_record.tsv"), sep = "\t")
print(dl_tab, nrows = 40)

# ==============================================================================
# [SECTION 6]  Cleaning the BioMart table (common to both directions)
# ==============================================================================
read_biomart <- function(bm, direction, cols = NULL) {
  nm <- names(bm)
  if (is.null(cols)) {
    cols <- if (direction == "mouse") {
      c(human = pick_col(nm, c("Human gene name", "Human Gene name")), mouse = pick_col(nm, c("Gene name", "Mouse gene name")),
        type = pick_col(nm, c("Human homology type")), mid = pick_col(nm, c("Gene stable ID", "Mouse gene stable ID")),
        conf = pick_col(nm, c("Human orthology confidence [0 low, 1 high]", "Human orthology confidence")))
    } else {
      c(human = pick_col(nm, c("Gene name", "HGNC symbol", "Human gene name")), mouse = pick_col(nm, c("Mouse gene name", "Mouse Gene name")),
        type = pick_col(nm, c("Mouse homology type")), mid = pick_col(nm, c("Mouse gene stable ID")),
        conf = pick_col(nm, c("Mouse orthology confidence [0 low, 1 high]", "Mouse orthology confidence")))
    }
  }
  if (anyNA(cols[c("human", "mouse", "type")])) { print(nm); stop("BioMart column names not recognised.") }
  t <- data.table(hsym = blank_na(bm[[cols["human"]]]), msym = blank_na(bm[[cols["mouse"]]]), mtype = blank_na(bm[[cols["type"]]]),
                  mid = if (!is.na(cols["mid"])) blank_na(bm[[cols["mid"]]]) else NA_character_,
                  conf = if (!is.na(cols["conf"])) blank_na(bm[[cols["conf"]]]) else NA_character_)
  t <- unique(t[!is.na(hsym)]); t[, hkey := up(hsym)]; t
}
bmt <- read_biomart(bm_raw, "human")
fwrite(bmt, P("_biomart_clean.tsv"), sep = "\t")
cat(sprintf("   cleaned BioMart table: %d rows | %d human symbols\n", nrow(bmt), uniqueN(bmt$hkey)))

# ==============================================================================
# [SECTION 7]  DEG, count matrix
# ==============================================================================
cat("\n[3] Reading the DEG results and the count matrix\n")
deg <- data.table::fread(DEG_ALL, data.table = TRUE); deg[, mkey := up(gene)]
raw_mkey <- up(data.table::fread(COUNT_MATRIX, select = 1, data.table = TRUE)[[1]])
cm_head <- names(data.table::fread(COUNT_MATRIX, nrows = 0))
miss_s <- setdiff(DESEQ_SAMPLES, cm_head); if (length(miss_s)) stop("Sample columns missing in COUNT_MATRIX: ", paste(miss_s, collapse = ", "))
cm <- data.table::fread(COUNT_MATRIX, select = c(cm_head[1], DESEQ_SAMPLES), data.table = TRUE); setnames(cm, 1, "gene")
pre_keep <- rowSums(as.matrix(cm[, -1, with = FALSE]) >= MIN_COUNT) >= MIN_SAMPLES
pre_pass_key <- up(cm$gene[pre_keep]); pre_ok <- setequal(pre_pass_key, deg$mkey)
cat(sprintf("   count matrix genes %d | passing the pre-filter %d | DEG results %d | (pre-filter pass = DEG results: %s)\n", length(raw_mkey), length(pre_pass_key), nrow(deg), pre_ok))
writeLines(sort(pre_pass_key), P("_prefilter_pass_genes.txt"))

trace_mouse <- function(mkey) {
  if (!(mkey %in% raw_mkey)) return(list(state = "absent from the mouse count matrix", detail = "mouse name absent from the count matrix (name mismatch possible)", resp = NA, lfc = NA_real_, padj = NA_real_))
  i <- match(mkey, deg$mkey)
  if (is.na(i)) return(list(state = "excluded for low expression", detail = "excluded by the pre-filter (count >= 10 in at least half of the samples)", resp = NA, lfc = NA_real_, padj = NA_real_))
  d <- deg[i]
  if (is.na(d$padj)) return(list(state = "excluded for low expression", detail = if (is.na(d$pvalue)) "padj NA (Cook's outlier)" else "padj NA (independent filtering)", resp = NA, lfc = d$log2FoldChange, padj = NA_real_))
  resp <- d$padj < PADJ_CUT && abs(d$log2FoldChange) >= LFC_CUT
  det <- if (resp) "responsive" else if (d$padj < PADJ_CUT) "below the response cut (only padj significant, |log2FC| < cut)" else "not responsive"
  list(state = "tested normally", detail = det, resp = resp, lfc = d$log2FoldChange, padj = d$padj)
}

# ==============================================================================
# [SECTION 8]  Audit per gene (function used for both directions)
# ==============================================================================
audit_genes <- function(bmt, direction) {
  res <- vector("list", nrow(genes)); amb <- list()
  for (i in seq_len(nrow(genes))) {
    k <- genes$gkey[i]; rows <- bmt[hkey == k]
    out <- list(Gene = genes$Gene[i], state = NA_character_, detail = NA_character_, mouse_genes = NA_character_, homology_types = NA_character_,
                n_mouse = 0L, strict_1to1 = FALSE, responsive = NA, log2FC = NA_real_, padj = NA_real_,
                permissive_n_tested = NA_integer_, permissive_any_responsive = NA, ensembl_conf = NA_character_)
    if (!nrow(rows)) {
      if (direction == "mouse") { out$state <- "undetermined (absent from table)"; out$detail <- "this human gene is not in the mouse-to-human table: no mouse ortholog / mouse gene not in the list (count matrix) / mouse name mismatch" }
      else { out$state <- "identifier conversion failed"; out$detail <- "no BioMart row for this symbol" }
    } else if (all(is.na(rows$msym))) {
      out$state <- "no mouse ortholog"; out$detail <- "the gene is in BioMart but has no mouse counterpart"
    } else {
      r2 <- unique(rows[!is.na(msym), .(msym, mtype, conf)])
      out$mouse_genes <- paste(sort(unique(r2$msym)), collapse = ";"); out$homology_types <- paste(sort(unique(r2$mtype)), collapse = ";")
      out$n_mouse <- uniqueN(r2$msym); out$ensembl_conf <- paste(sort(unique(r2$conf)), collapse = ";")
      out$strict_1to1 <- (out$n_mouse == 1L) && all(r2$mtype == "ortholog_one2one")
      tr <- lapply(unique(r2$msym), function(m) c(list(msym = m), trace_mouse(up(m))))
      if (out$strict_1to1) {
        t1 <- tr[[1]]; out$state <- t1$state; out$detail <- t1$detail; out$responsive <- t1$resp; out$log2FC <- t1$lfc; out$padj <- t1$padj
      } else {
        out$state <- "ambiguous mapping"; out$detail <- sprintf("not one-to-one (%s; %d mouse genes)", out$homology_types, out$n_mouse)
        tested <- Filter(function(t) t$state == "tested normally", tr)
        out$permissive_n_tested <- length(tested)
        out$permissive_any_responsive <- any(vapply(tested, function(t) isTRUE(t$resp), logical(1)))
        best <- if (length(tested)) tested[[which.min(vapply(tested, function(t) t$padj, numeric(1)))]] else tr[[1]]
        out$log2FC <- best$lfc; out$padj <- best$padj
        amb[[length(amb) + 1]] <- rbindlist(lapply(tr, function(t) data.table(Gene = genes$Gene[i], mouse = t$msym, homology = paste(unique(r2[msym == t$msym]$mtype), collapse = ";"),
                                                                              state = t$state, detail = t$detail, log2FC = t$lfc, padj = t$padj)), fill = TRUE)
      }
    }
    res[[i]] <- as.data.table(out)
  }
  aud <- rbindlist(res, fill = TRUE)
  aud <- merge(genes[, .(Gene, biotype, loci, tier, min_dist_kb, GWAS_liver, GWAS_thyroid, GWAS_shared)], aud, by = "Gene", sort = FALSE)
  aud[, state_detail := ifelse(state %in% c("tested normally", "excluded for low expression"), paste0(state, " - ", detail), state)]
  aud[, `responsive_primary_one2one` := fifelse(is.na(responsive), NA_character_, fifelse(responsive, "responsive", "not responsive/below cut"))]
  aud[, `responsive_permissive` := NA_character_]
  aud[state == "ambiguous mapping" & permissive_any_responsive %in% TRUE, `responsive_permissive` := "possibly responsive"]
  aud[state == "ambiguous mapping" & !(permissive_any_responsive %in% TRUE) & !is.na(permissive_n_tested), `responsive_permissive` := "not responsive"]
  list(aud = aud, amb = if (length(amb)) rbindlist(amb, fill = TRUE) else NULL)
}

cat("\n[4] Running the audit (human -> mouse table)\n")
A <- audit_genes(bmt, "human"); aud <- A$aud; amb <- A$amb

# ==============================================================================
# [SECTION 9]  Summary, candidates, QC
# ==============================================================================
summ <- dcast(aud[, .N, by = .(state, biotype)], state ~ biotype, value.var = "N", fill = 0L)
summ[, total := rowSums(.SD), .SDcols = setdiff(names(summ), "state")]
ord <- c("identifier conversion failed", "undetermined (absent from table)", "no mouse ortholog", "ambiguous mapping", "absent from the mouse count matrix", "excluded for low expression", "tested normally")
summ <- summ[order(match(state, ord))]
detail_sum <- aud[state %in% c("tested normally", "excluded for low expression"), .N, by = .(state = state, detail = detail, biotype)][order(state, detail)]
cand <- aud[`responsive_primary_one2one` == "responsive" | `responsive_permissive` == "possibly responsive"]
cand <- cand[, .(Gene, biotype, loci, tier, group = fifelse(`responsive_primary_one2one` %in% "responsive", "1_primary (one-to-one orthologue)", "2_permissive (ambiguous mapping)"),
                 mouse_genes, homology_types, log2FC, padj, GWAS_liver, GWAS_thyroid, GWAS_shared)]
setorder(cand, group, padj)

qc <- data.table(item = character(), value = character()); qadd <- function(a, b) qc <<- rbind(qc, data.table(item = a, value = as.character(b)))
qadd("target genes (scope)", sprintf("%d (%s)", nrow(aud), SCOPE))
qadd("Ensembl release (requested) / direction", sprintf("%s / human -> mouse", ENSEMBL_RELEASE))
qadd("BioMart cleaned rows / human symbols", sprintf("%d / %d", nrow(bmt), uniqueN(bmt$hkey)))
qadd("symbols for which BioMart returned no row", length(not_returned))
qadd("pre-filter reproduction: genes passing / DEG results / identical", sprintf("%d / %d / %s", length(pre_pass_key), nrow(deg), pre_ok))
qadd("genes with a strict one-to-one orthologue", sum(aud$strict_1to1))
qadd("primary responsive genes | genes possibly responsive in permissive only", sprintf("%d | %d", sum(aud$`responsive_primary_one2one` %in% "responsive"), sum(aud$`responsive_permissive` %in% "possibly responsive")))
if (!is.na(HUMAN_RESP_FILE) && file.exists(HUMAN_RESP_FILE)) {
  hr <- data.table::fread(HUMAN_RESP_FILE, data.table = TRUE); hcol <- pick_col(names(hr), c("Human gene name", "Human.gene.name"))
  if (!is.na(hcol)) {
    hs <- unique(up(hr[[hcol]])); mine <- up(aud[`responsive_primary_one2one` %in% "responsive" | `responsive_permissive` %in% "possibly responsive", Gene])
    qadd(sprintf("comparison with the earlier human list of responsive genes (%d genes): responsive candidates of this audit that are not in the list", length(hs)), paste(setdiff(mine, hs), collapse = ", "))
    in_list <- aud[up(Gene) %in% hs, .(Gene, state, detail)]
    qadd("genes in the earlier list that are not responsive in this audit (reason to be checked)",
         paste(sprintf("%s[%s]", in_list[!(Gene %in% cand$Gene)]$Gene, in_list[!(Gene %in% cand$Gene)]$state), collapse = ", "))
  }
}

# ==============================================================================
# [SECTION 10]  (optional) comparison with the mouse-direction table
# ==============================================================================
cmp_tab <- NULL
if (!is.na(CROSSCHECK_MOUSE_FILE) && file.exists(CROSSCHECK_MOUSE_FILE)) {
  cat("\n[5] Comparison with the mouse-direction table:", CROSSCHECK_MOUSE_FILE, "\n")
  l1 <- readLines(CROSSCHECK_MOUSE_FILE, n = 1, warn = FALSE); sp <- if (grepl("\t", l1, fixed = TRUE)) "\t" else ","
  bmm <- read_biomart(data.table::fread(CROSSCHECK_MOUSE_FILE, sep = sp, quote = "\"", fill = TRUE, data.table = TRUE, colClasses = "character"), "mouse")
  Am <- audit_genes(bmm, "mouse")$aud
  coarse <- function(s) ifelse(s %in% c("identifier conversion failed", "no mouse ortholog", "undetermined (absent from table)"), "absent from table or no ortholog", s)
  m <- merge(aud[, .(Gene, biotype, state_human = state, mouse_human = mouse_genes, strict_human = strict_1to1, resp_human = `responsive_primary_one2one`)],
             Am[, .(Gene, state_mouse = state, mouse_mouse = mouse_genes, strict_mouse = strict_1to1, resp_mouse = `responsive_primary_one2one`)], by = "Gene", sort = FALSE)
  m[, same_state := coarse(state_human) == coarse(state_mouse)]
  m[, same_result := same_state & (strict_human == strict_mouse) & (resp_human %in% resp_mouse | (is.na(resp_human) & is.na(resp_mouse)))]
  cmp_tab <- m[same_result == FALSE]
  cat(sprintf("   genes with different states in the two directions: %d / %d\n", nrow(cmp_tab), nrow(m)))
  qadd("direction comparison (human table vs mouse table): genes with different results", sprintf("%d / %d", nrow(cmp_tab), nrow(m)))
  xt <- as.data.table(table(human_table = m$state_human, mouse_table = m$state_mouse)); xt <- xt[N > 0]
  print(xt)
  cmp_tab <- list(diff = cmp_tab, cross = xt)
}

# ==============================================================================
# [SECTION 11]  Save: results, settings, input file record, session information
# ==============================================================================
cat("\n[6] Saving\n")
all_settings <- data.table(item = c("MAPPING_XLSX", "SHEET_GENES", "SCOPE", "SYMBOLS_FILE", "DOWNLOAD", "BIOMART_MANUAL_FILE", "MANUAL_DATASET_NOTE", "MANUAL_DOWNLOAD_DATE", "DOWNLOAD_METHOD", "ENSEMBL_RELEASE", "ENSEMBL_HOST", "BM_RETRIES", "BM_BIOMART", "BM_DATASET", "BM_FILTER",
                                    "BM_BATCH", "BM_ATTRIBUTES", "BM_HEADERS", "DEG_ALL", "COUNT_MATRIX", "DESEQ_SAMPLES", "MIN_COUNT", "MIN_SAMPLES", "PADJ_CUT", "LFC_CUT",
                                    "HUMAN_RESP_FILE", "CROSSCHECK_MOUSE_FILE", "OUT_DIR", "OUT_PREFIX", "run start", "working directory"),
                           value = c(MAPPING_XLSX, SHEET_GENES, SCOPE, SYMBOLS_FILE, DOWNLOAD, BIOMART_MANUAL_FILE, MANUAL_DATASET_NOTE, MANUAL_DOWNLOAD_DATE, DOWNLOAD_METHOD, ENSEMBL_RELEASE, if (is.null(ENSEMBL_HOST)) "NULL (default)" else ENSEMBL_HOST, BM_RETRIES, BM_BIOMART, BM_DATASET, BM_FILTER,
                                  BM_BATCH, paste(BM_ATTRIBUTES, collapse = ", "), paste(BM_HEADERS, collapse = " | "), DEG_ALL, COUNT_MATRIX, paste(DESEQ_SAMPLES, collapse = ", "),
                                  MIN_COUNT, MIN_SAMPLES, PADJ_CUT, LFC_CUT, HUMAN_RESP_FILE, CROSSCHECK_MOUSE_FILE, OUT_DIR, OUT_PREFIX, format(RUN_START, "%Y-%m-%d %H:%M:%S"), getwd()))
all_settings[, value := as.character(value)]
fwrite(all_settings, P("_settings.tsv"), sep = "\t")
inp <- unique(c(MAPPING_XLSX, SYMBOLS_FILE, DEG_ALL, COUNT_MATRIX, HUMAN_RESP_FILE, CROSSCHECK_MOUSE_FILE, if (!DOWNLOAD) c(BIOMART_MANUAL_FILE, BIOMART_RAW_FILE))); inp <- inp[!is.na(inp)]
man <- rbindlist(lapply(inp, function(f) if (file.exists(f)) data.table(file = f, size_bytes = file.size(f), modified = format(file.mtime(f), "%Y-%m-%d %H:%M:%S"), md5 = unname(tools::md5sum(f)))
                        else data.table(file = f, size_bytes = NA, modified = NA, md5 = "file not found")))
fwrite(man, P("_input_files_manifest.tsv"), sep = "\t")

wb <- createWorkbook()
add_sheet(wb, "0_download_record", dl_tab); add_sheet(wb, "0b_settings", all_settings); add_sheet(wb, "0c_input_files", man)
add_sheet(wb, "1_summary", summ); add_sheet(wb, "1b_detail", detail_sum); add_sheet(wb, "2_audit_table", aud); add_sheet(wb, "3_candidates_responsive", cand)
if (!is.null(amb)) add_sheet(wb, "4_ambiguous_detail", amb)
add_sheet(wb, "5_QC", qc); add_sheet(wb, "6_BioMart_table", bm_raw)
if (!is.null(cmp_tab)) { add_sheet(wb, "7_direction_comparison", if (nrow(cmp_tab$diff)) cmp_tab$diff else data.table(note = "the results of the two directions are identical")); add_sheet(wb, "7b_state_crosstab", cmp_tab$cross) }
saveWorkbook(wb, P("_result.xlsx"), overwrite = TRUE)
fwrite(aud, P("_table.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), P("_sessionInfo.txt"))

cat("\n===== Summary =====\n"); print(summ); cat("\nResponsive candidates:\n"); print(cand[, .(Gene, group, mouse_genes, log2FC, padj)])
cat(sprintf("\nSaved in: %s\nEnd: %s (elapsed %.1f s)\n", normalizePath(OUT_DIR), format(Sys.time(), "%Y-%m-%d %H:%M:%S"), as.numeric(difftime(Sys.time(), RUN_START, units = "secs"))))

}, warning = function(w) { cat("[warning]", conditionMessage(w), "\n"); invokeRestart("muffleWarning") },
   message = function(m) { cat("[message]", conditionMessage(m)); invokeRestart("muffleMessage") })
"OK"
}, error = function(e) { cat("\n[error]", conditionMessage(e), "\n"); paste("ERROR:", conditionMessage(e)) },
   finally = { sink(); close(con_log) })
cat("\nRun status:", status, "\nLog:", P("_run_log.txt"), "\n")
