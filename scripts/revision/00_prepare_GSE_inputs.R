#!/usr/bin/env Rscript
# =============================================================================
# 00_prepare_GSE_inputs.R
#
# Builds, from the ORIGINAL GEO files, the input files that config.R expects in data/.
# Run this once before 02, 03, 04 (maturity table), 05 and 06.
#
#   config.R entry                          file written to data/              original files it is built from
#   --------------------------------------  ---------------------------------  ----------------------------------------------------
#   mouse_counts_file / mouse_meta_file     GSE85793_counts.csv / _meta.csv    GSE85793_RAW.tar (or the per-sample htseq files)
#                                                                              + GSE85793_series_matrix.txt
#   cohorts$GSE130970$counts / $meta        GSE130970_counts.csv / _meta.csv   GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv
#                                                                              + GSE130970_series_matrix.txt
#   cohorts$GSE135251$counts / $meta        GSE135251_counts.csv / _meta.csv   GSE135251_RAW.tar (or the 216 per-sample count files)
#                                                                              + GSE135251_series_matrix.txt
#   hep_counts_file / hep_meta_file         GSE270233_counts.csv / _meta.csv   GSE270233_counts.csv + GSE270233_series_matrix.txt
#
# Which script reads what:   02, 06 -> GSE135251, GSE130970      03, 05 -> GSE85793      04 (table T3) -> GSE270233
#
# Processing rules (same as the submitted analysis, scripts/09_human_liver_fibrosis.R)
#   * GSE85793   : only the 8 samples used in the manuscript (4 Control = "vehicle", 4 T3 = "T3"); gene symbols as deposited.
#                  The Gluc, Gluc+T3 and Cp groups are dropped.
#                  Note: the deposited file of animal T3 (GSM2330539) is comma-separated while the other 19 are tab-separated;
#                  the delimiter is therefore detected per file.
#   * GSE130970  : tximport counts with Entrez IDs -> HGNC symbol (org.Hs.eg.db, first match). When several rows give the
#                  same symbol, the row with the largest total count is kept. Counts are rounded to integers.
#                  All 78 samples are written with group = "NAFLD" (the histologically normal samples cannot be told apart
#                  from one column, so none is removed, as in the submitted analysis).
#   * GSE135251  : per-sample htseq count files (Ensembl IDs) -> HGNC symbol, same duplicate rule. 216 samples are written;
#                  group = disease ("NAFLD" 206, "Control" 10). config.R keeps group_keep = "NAFLD".
#                  Age and sex are not deposited for this cohort, so those columns do not exist in its meta file.
#   * GSE270233  : the deposited count table is copied (first column renamed to "gene"); S1-S12 = GSM8337685-GSM8337696 in
#                  the order of the series matrix (WT/KO x vehicle/T3, 3 each). The order is checked against the series matrix.
#
# Working directory: scripts/revision (the folder that contains config.R and data/).
# Packages: data.table; org.Hs.eg.db and AnnotationDbi (Bioconductor) for the two human biopsy cohorts.
# Note: the symbol of an Entrez/Ensembl ID depends on the installed org.Hs.eg.db version. Gene lists therefore can differ by a
#       few renamed genes from a run with a different annotation version.
# Output: data/*_counts.csv, data/*_meta.csv, data/00_prepare_GSE_inputs_log.txt
# =============================================================================
suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

## ---- [settings] ------------------------------------------------------------------------------------------------------
data_dir <- "data"                       # same folder as data_dir in config.R

# Folders that contain the original GEO files. They are searched from top to bottom, so putting the files in
# data/raw_geo is the simplest option. Add other folders (full paths on your computer) to this vector if needed.
raw_dirs <- c(
  file.path(data_dir, "raw_geo")
)

# optional: the 8-sample count matrix written by scripts/02_build_count_matrix_GSE85793.R, copied into data/;
# used only for a comparison with the matrix built here. Set to NA_character_ to skip.
reference_mouse_matrix <- file.path(data_dir, "1_cleaned_count_matrix.csv")

# genes used for the coverage report at the end (same lists as config.R)
candidate_genes    <- c("ABHD1", "F13B", "ADH4", "LGALS4")
hepatocyte_markers <- c("ALB", "APOA1", "APOC3", "TTR", "SERPINA1", "AHSG", "FGB", "APOH")
control_genes      <- c("COL1A1", "COL1A2", "COL3A1", "TIMD4", "TIMP1", "ACTA2", "LUM", "PDGFRB", "THY1")
maturity_genes     <- c("AFP", "DLK1", "ALB", "ASGR1", "CYP3A4", "CYP2E1", "CYP1A2", "UGT1A1", "GSTM1")
mouse_sets_file    <- file.path(data_dir, "redox_gene_sets_mouse.csv")                  # optional, coverage check only

.cfg <- Sys.getenv("GSE_PREP_CFG", ""); if (nzchar(.cfg)) source(.cfg)   # optional override (used for testing)

## ---- helpers ---------------------------------------------------------------------------------------------------------
log_lines <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); log_lines <<- c(log_lines, s) }
as_num <- function(x) suppressWarnings(as.numeric(x))
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

find_file <- function(pattern) {
  for (d in raw_dirs) {
    if (!dir.exists(d)) next
    f <- list.files(d, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    f <- f[!dir.exists(f)]
    if (length(f)) return(f[1])
  }
  stop("File not found (pattern: ", pattern, ").\n  Folders searched:\n    ", paste(raw_dirs, collapse = "\n    "),
       "\n  Put the original file in one of them or edit raw_dirs at the top of this script.", call. = FALSE)
}

# series matrix -> one row per sample (gsm, title and every "key: value" characteristic, key made lower case with "_")
read_series_matrix <- function(f) {
  L <- readLines(f, warn = FALSE, encoding = "UTF-8")
  L <- L[startsWith(L, "!Sample_") & !startsWith(L, "!Sample_supplementary_file")]
  cut1 <- function(l) { p <- strsplit(l, "\t", fixed = TRUE)[[1]]; list(key = sub("^!", "", p[1]), val = gsub('^"|"$', "", p[-1])) }
  P <- lapply(L, cut1)
  get1 <- function(k) P[[which(vapply(P, function(x) x$key == k, TRUE))[1]]]$val
  df <- data.frame(gsm = get1("Sample_geo_accession"), title = get1("Sample_title"), stringsAsFactors = FALSE)
  for (x in P[vapply(P, function(x) x$key == "Sample_characteristics_ch1", TRUE)]) {
    k <- trimws(sub(":.*$", "", x$val[1])); k <- gsub("[^a-z0-9]+", "_", tolower(k))
    df[[make.unique(c(names(df), k))[ncol(df) + 1]]] <- trimws(sub("^[^:]*:", "", x$val))
  }
  df
}
pick_col <- function(df, key) if (key %in% names(df)) df[[key]] else rep(NA, nrow(df))

gsm_of <- function(x) {
  b <- basename(x); m <- regexpr("GSM[0-9]+", b)
  ifelse(m > 0, substr(b, m, m + attr(m, "match.length") - 1), NA_character_)
}

# per-sample files (a tar, an unpacked folder, or loose files) -> named vector of file paths, one per requested GSM id
collect_gsm_files <- function(tar_pattern, dir_pattern, gsms) {
  srcs <- character()
  for (d in raw_dirs) if (dir.exists(d)) {
    srcs <- c(srcs, list.files(d, pattern = tar_pattern, full.names = TRUE, ignore.case = TRUE))
    sub_d <- list.files(d, pattern = dir_pattern, full.names = TRUE, include.dirs = TRUE, ignore.case = TRUE)
    srcs <- c(srcs, sub_d[dir.exists(sub_d)], d)
  }
  for (s in srcs) {
    if (dir.exists(s)) {
      fs <- list.files(s, full.names = TRUE); fs <- fs[!dir.exists(fs)]
    } else {
      td <- file.path(tempdir(), paste0("untar_", gsub("[^A-Za-z0-9]", "_", basename(s))))
      if (!dir.exists(td)) { dir.create(td, recursive = TRUE); untar(s, exdir = td) }
      fs <- list.files(td, full.names = TRUE, recursive = TRUE)
    }
    g <- gsm_of(fs)
    if (all(gsms %in% g)) { cat("  per-sample files from:", s, "\n"); return(setNames(fs[match(gsms, g)], gsms)) }
  }
  stop("Could not find the per-sample count files for ", length(gsms), " samples (tar pattern ", tar_pattern, ").\n",
       "  Folders searched:\n    ", paste(raw_dirs, collapse = "\n    "), call. = FALSE)
}

# a 2-column count file (gene id, count), gz or plain, with or without header -> named numeric vector.
# The delimiter is detected per file: one of the deposited GSE85793 files (T3) is comma-separated, the other 19 are tab-separated.
read_two_col <- function(f, header) {
  first <- readLines(gzfile(f), n = 1, warn = FALSE)
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (grepl(",", first, fixed = TRUE)) "," else
    stop("Cannot detect the delimiter of ", basename(f), call. = FALSE)
  x <- scan(gzfile(f), what = list(gene = "", n = 0), sep = sep, skip = as.integer(header), quote = "",
            comment.char = "", flush = TRUE, quiet = TRUE)
  ok <- !is.na(x$n) & !grepl("^(__|N_)", x$gene)           # htseq / STAR summary rows
  g <- x$gene[ok]; v <- x$n[ok]
  if (anyDuplicated(g)) { s <- tapply(v, g, sum); g <- names(s); v <- as.numeric(s) }
  structure(setNames(v, g), sep = sep)
}
report_delimiters <- function(vecs) {
  s <- vapply(vecs, function(v) attr(v, "sep"), "")
  if (any(s != "\t")) say("  note: %d of %d files are not tab-separated (read correctly): %s", sum(s != "\t"), length(s), paste(names(s)[s != "\t"], collapse = ", "))
}

merge_sample_vectors <- function(vecs) {
  genes <- Reduce(intersect, lapply(vecs, names))
  lost <- 1 - length(genes) / length(unique(unlist(lapply(vecs, names))))
  if (lost > 0.01) warning(sprintf("%.1f%% of the gene ids are not shared by all samples and were dropped.", 100 * lost))
  do.call(cbind, lapply(vecs, function(v) v[genes]))
}

ids_to_symbol <- function(ids, keytype) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("Bioconductor packages org.Hs.eg.db and AnnotationDbi are required: BiocManager::install(c('org.Hs.eg.db','AnnotationDbi'))", call. = FALSE)
  ids <- as.character(ids)
  m <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ids), column = "SYMBOL",
                                              keytype = keytype, multiVals = "first"))
  unname(m[ids])
}
# several rows with the same symbol -> keep the row with the largest total count; drop rows without symbol
collapse_to_symbol <- function(mat, sym) {
  keep <- !is.na(sym) & sym != ""
  mat <- mat[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(rowSums(mat), decreasing = TRUE); mat <- mat[o, , drop = FALSE]; sym <- sym[o]
  d <- !duplicated(sym); mat <- mat[d, , drop = FALSE]; rownames(mat) <- sym[d]
  mat
}
to_integer <- function(mat) { m <- round(mat); storage.mode(m) <- "integer"; m }     # integers avoid 1e+05 style output

write_counts <- function(mat, name) {
  f <- file.path(data_dir, paste0(name, "_counts.csv"))
  fwrite(data.table(gene = rownames(mat), as.data.table(mat)), f)
  say("  wrote %s  (%d genes x %d samples)", f, nrow(mat), ncol(mat))
}
write_meta <- function(meta, name) {
  f <- file.path(data_dir, paste0(name, "_meta.csv")); fwrite(meta, f)
  say("  wrote %s  (%d samples; columns: %s)", f, nrow(meta), paste(names(meta), collapse = ", "))
}
coverage <- function(label, mat, genes) {
  have <- genes[genes %in% rownames(mat)]
  say("  %-18s %2d/%2d present%s", label, length(have), length(genes),
      if (length(have) < length(genes)) paste0("  | missing: ", paste(setdiff(genes, have), collapse = ", ")) else "")
}
expect <- function(what, got, want) say("  check %-34s %s (found %s, expected %s)", what, if (identical(as.integer(got), as.integer(want))) "OK  " else "DIFF", got, want)

## =============================================================================================================
## 1. GSE85793  mouse liver, vehicle vs T3 (4 vs 4)       -> used by 03 and 05
## =============================================================================================================
say("== GSE85793 (mouse liver) ==")
sm <- read_series_matrix(find_file("^GSE85793_series_matrix\\.txt$"))
sm <- sm[sm$treatment %in% c("Control", "T3"), ]
sm$sample <- sub("^Liver_", "", sm$title)
sm$group  <- ifelse(sm$treatment == "Control", "vehicle", "T3")
files <- collect_gsm_files("^GSE85793_RAW\\.tar$", "^GSE85793_RAW", sm$gsm)
vecs  <- lapply(files, read_two_col, header = TRUE)
report_delimiters(vecs)
mat   <- merge_sample_vectors(vecs); colnames(mat) <- sm$sample
mouse <- to_integer(mat)
# header of each per-sample file must carry the same animal name as the series matrix title
hdr <- vapply(files, function(f) strsplit(readLines(gzfile(f), n = 1), "[\t,]")[[1]][2], "")
if (!all(hdr == sm$sample)) warning("Column header of a per-sample file differs from the series matrix title: ",
                                    paste(paste(hdr, sm$sample, sep = " vs "), collapse = "; "))
write_counts(mouse, "GSE85793")
write_meta(data.table(sample = sm$sample, gsm = sm$gsm, title = sm$title, treatment = sm$treatment, group = sm$group), "GSE85793")
expect("samples (4 vehicle + 4 T3)", ncol(mouse), 8)
say("  group sizes: %s", paste(names(table(sm$group)), table(sm$group), sep = " = ", collapse = ", "))
if (!is.na(reference_mouse_matrix) && file.exists(reference_mouse_matrix)) {
  ref <- fread(reference_mouse_matrix, data.table = FALSE); rg <- ref[[1]]; rm <- as.matrix(ref[, -1]); rownames(rm) <- rg
  same <- identical(rownames(mouse), rg) && identical(colnames(mouse), colnames(rm)) && all(mouse == rm[rownames(mouse), colnames(mouse)])
  say("  comparison with the matrix of the submitted DEG (%s): %s", basename(reference_mouse_matrix),
      if (same) "IDENTICAL (genes, order, samples, values)" else "DIFFERENT, check the original file")
}
say("  genes with count >= 10 in at least half of the samples: %d", sum(rowSums(mouse >= 10) >= ncol(mouse) / 2))
if (file.exists(mouse_sets_file)) {
  ms <- fread(mouse_sets_file); ms[, mouse_symbol := trimws(mouse_symbol)]
  say("  redox sets: %s", paste(sprintf("%s %d/%d", names(split(ms$mouse_symbol, ms$set)),
        vapply(split(ms$mouse_symbol, ms$set), function(g) sum(unique(g) %in% rownames(mouse)), 1L),
        vapply(split(ms$mouse_symbol, ms$set), function(g) length(unique(g)), 1L)), collapse = "; "))
}

## =============================================================================================================
## 2. GSE130970  human liver biopsies (78)                 -> used by 02 and 06
## =============================================================================================================
say("\n== GSE130970 (human liver biopsies) ==")
sm <- read_series_matrix(find_file("^GSE130970_series_matrix\\.txt$"))
cnt <- fread(find_file("^GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID\\.csv$"), data.table = FALSE)
ids <- as.character(cnt[[1]]); mat <- as.matrix(cnt[, -1, drop = FALSE]); storage.mode(mat) <- "double"
idx <- match(colnames(mat), sm$title); if (anyNA(idx)) idx <- match(colnames(mat), sm$gsm)
if (anyNA(idx)) stop("Column names of the count file do not match the series matrix title or GSM id.", call. = FALSE)
sm <- sm[idx, ]; colnames(mat) <- sm$gsm
sym <- ids_to_symbol(ids, "ENTREZID")
say("  Entrez id -> symbol: %d of %d rows mapped", sum(!is.na(sym)), length(ids))
mat <- to_integer(collapse_to_symbol(mat, sym))
write_counts(mat, "GSE130970")
meta130 <- data.table(sample = sm$gsm, title = sm$title, group = "NAFLD",
                      fibrosis_stage = as.integer(as_num(pick_col(sm, "fibrosis_stage"))),
                      nas = as_num(pick_col(sm, "nafld_activity_score")),
                      steatosis = as_num(pick_col(sm, "steatosis_grade")),
                      lobular_inflammation = as_num(pick_col(sm, "lobular_inflammation_grade")),
                      ballooning = as_num(pick_col(sm, "cytological_ballooning_grade")),
                      age = as_num(pick_col(sm, "age_at_biopsy")), sex = pick_col(sm, "sex"))
write_meta(meta130, "GSE130970")
expect("samples", nrow(meta130), 78)
say("  fibrosis stage F0-F4: %s", paste(sprintf("F%s=%d", names(table(meta130$fibrosis_stage)), as.integer(table(meta130$fibrosis_stage))), collapse = ", "))
say("  missing values: stage %d, NAS %d, age %d, sex %d | NAS = steatosis + inflammation + ballooning in %d of %d samples",
    sum(is.na(meta130$fibrosis_stage)), sum(is.na(meta130$nas)), sum(is.na(meta130$age)), sum(is.na(meta130$sex)),
    sum(meta130$nas == meta130$steatosis + meta130$lobular_inflammation + meta130$ballooning, na.rm = TRUE), nrow(meta130))
coverage("candidates", mat, candidate_genes); coverage("hepatocyte markers", mat, hepatocyte_markers); coverage("fibrosis controls", mat, control_genes)
mat130 <- mat

## =============================================================================================================
## 3. GSE135251  human liver biopsies (206 NAFLD + 10 control)   -> used by 02 and 06
## =============================================================================================================
say("\n== GSE135251 (human liver biopsies) ==")
sm <- read_series_matrix(find_file("^GSE135251_series_matrix\\.txt$"))
files <- collect_gsm_files("^GSE135251_RAW\\.tar$", "^GSE135251_RAW", sm$gsm)
vecs  <- lapply(files, read_two_col, header = FALSE)
report_delimiters(vecs)
mat   <- merge_sample_vectors(vecs); colnames(mat) <- sm$gsm
sym   <- ids_to_symbol(sub("\\..*$", "", rownames(mat)), "ENSEMBL")        # drop an Ensembl version suffix if present
say("  Ensembl id -> symbol: %d of %d rows mapped", sum(!is.na(sym)), nrow(mat))
mat <- to_integer(collapse_to_symbol(mat, sym))
write_counts(mat, "GSE135251")
meta135 <- data.table(sample = sm$gsm, title = sm$title, group = pick_col(sm, "disease"),
                      group_in_paper = pick_col(sm, "group_in_paper"),
                      fibrosis_stage = as.integer(as_num(pick_col(sm, "fibrosis_stage"))),
                      nas = as_num(pick_col(sm, "nas_score")))
write_meta(meta135, "GSE135251")
expect("samples (206 NAFLD + 10 control)", nrow(meta135), 216)
expect("NAFLD samples", sum(meta135$group == "NAFLD"), 206)
nafld <- meta135[group == "NAFLD"]
say("  NAFLD fibrosis stage F0-F4: %s", paste(sprintf("F%s=%d", names(table(nafld$fibrosis_stage)), as.integer(table(nafld$fibrosis_stage))), collapse = ", "))
say("  missing values among NAFLD: stage %d, NAS %d | age and sex are not deposited for this cohort", sum(is.na(nafld$fibrosis_stage)), sum(is.na(nafld$nas)))
coverage("candidates", mat, candidate_genes); coverage("hepatocyte markers", mat, hepatocyte_markers); coverage("fibrosis controls", mat, control_genes)
mat135 <- mat

## =============================================================================================================
## 4. GSE270233  iPSC-derived hepatocytes, WT / THRB-KO x vehicle / T3 (3 each)   -> used by 04 (table T3)
## =============================================================================================================
say("\n== GSE270233 (human iPSC hepatocytes) ==")
sm  <- read_series_matrix(find_file("^GSE270233_series_matrix\\.txt$"))
cnt <- fread(find_file("^GSE270233_counts\\.csv$"), data.table = FALSE)
genes <- trimws(as.character(cnt[[1]])); mat <- as.matrix(cnt[, -1, drop = FALSE]); storage.mode(mat) <- "double"
if (ncol(mat) != nrow(sm)) stop("The count file has ", ncol(mat), " samples but the series matrix lists ", nrow(sm), ".", call. = FALSE)
nd <- sum(duplicated(genes)); mat <- to_integer(collapse_to_symbol(mat, genes))
if (nd) say("  %d duplicated gene symbols collapsed (row with the largest total kept)", nd)
meta270 <- data.table(sample = colnames(mat), gsm = sm$gsm, title = sm$title,
                      genotype  = ifelse(grepl("knockout|KO", pick_col(sm, "genotype"), ignore.case = TRUE), "KO", "WT"),
                      treatment = ifelse(grepl("T3", pick_col(sm, "treatment"), ignore.case = TRUE), "T3", "Veh"))
# S1-S12 are assigned to the series matrix samples by position; the design must be WT x6 then KO x6, Veh x3 / T3 x3 within each
ok <- identical(meta270$genotype, rep(c("WT", "KO"), each = 6)) && identical(meta270$treatment, rep(rep(c("Veh", "T3"), each = 3), 2))
say("  sample order check (S1-S12 vs series matrix): %s", if (ok) "OK (matches the design used for the submitted DESeq2)" else "DIFFERENT, check the order before using this file")
write_counts(mat, "GSE270233")
write_meta(meta270, "GSE270233")
coverage("maturity markers", mat, maturity_genes)
mat270 <- mat

## ---- summary ---------------------------------------------------------------------------------------------------------
say("\nDone. config.R can now read: mouse_counts_file, mouse_meta_file, cohorts (GSE135251, GSE130970), hep_counts_file, hep_meta_file.")
say("The files in data/ for 01 (DEG_hepT3_raw.csv, redox gene sets) and 04 (mouse_T3_results.csv) were not touched.")
writeLines(log_lines, file.path(data_dir, "00_prepare_GSE_inputs_log.txt"))
