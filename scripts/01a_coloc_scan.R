################################################################################
#  01a_coloc_scan.R
#  Colocalization scan of the liver MRI GWAS and the thyroid GWAS.
#  Output = the 458 windows and the 29,802 completed dataset-pair tests of Supplementary Table S2b.
#
#  METHODS (as in the manuscript)
#    Pass 1  Each of the 12 GWAS files is read with chromosome, position and P only. Variants with P < 1e-6 are collected.
#    Windows Distance-based greedy clumping over the union of all 12 datasets: variants are taken in order of P; each
#            variant not yet covered becomes a lead variant and defines a window of lead position +/- 500 kb (variants
#            within 500 kb are covered). This gives 458 windows (GRCh37), named L00001 ... in order of chromosome and position.
#    Pass 2  Each dataset is read again, filtered (valid beta, SE, P, allele frequency; MAF >= 0.01; one row per position)
#            and restricted to the windows.
#    Test    coloc.abf (R package coloc) for each window and each of the 66 dataset pairs (12 choose 2), with
#            p1 = p2 = 1e-4 and p12 = 1e-5. A pair is tested only if the two datasets share at least 50 variants (by position)
#            in the window and the smallest P of at least one dataset is <= 0.05. Variants are matched by position.
#            The PP.H4 >= 0.8 flag written here is the original screening flag; the tiers of the manuscript are
#            assigned by 01b_coloc_window_tiers.R from the posterior probabilities.
#
#  INPUTS (DATA_DIR, searched by file name, also in sub-folders; public data, see data/README.md)
#    GCST90446646.tsv ... GCST90446649.tsv   cT1 GWAS, 4 adjustment models (GWAS Catalog)
#    GCST90446650.tsv ... GCST90446653.tsv   PDFF GWAS, 4 adjustment models (GWAS Catalog)
#    GCST010653.txt                          TSH GWAS, dataset A (GWAS Catalog)
#    GCST90296333.tsv                        TSH GWAS, dataset B (GWAS Catalog)
#    FT3_summary_data.txt, FT4_summary_data.txt   FT3 and FT4 GWAS (ThyroidOmics Consortium)
#    All files are on genome build GRCh37. A file is also found if only its folder or file name contains the GWAS Catalog accession;
#    .gz files are read directly (R package R.utils).
#
#  OUTPUTS (OUT_DIR, default results/01a_coloc)
#    coloc_windows.tsv   window, chr, start, end, lead variant (chr:pos), lead position, lead P
#    coloc_pairs.tsv     one row per tested window and dataset pair: PP.H0 - PP.H4, number of variants, group labels
#    01a_coloc_scan_sessionInfo.txt
#
#  RUNTIME AND MEMORY  Reading the 12 files takes about 1 h; the 458 x 66 tests take about 1 h more. Roughly 16 GB of memory.
#
#  REPRODUCTION CHECK  The 29,802 pair rows of the original scan equal the rows of Supplementary Table S2b (all posterior
#  probabilities identical). 01b_coloc_window_tiers.R repeats this comparison when both are available (see docs/PROVENANCE.md).
#
#  USAGE  DATA_DIR=/path/to/data OUT_DIR=results/01a_coloc Rscript 01a_coloc_scan.R
################################################################################
suppressPackageStartupMessages({ library(data.table); library(coloc) })

## ---- settings ----------------------------------------------------------------
DATA_DIR <- Sys.getenv("DATA_DIR", "data")
OUT_DIR  <- Sys.getenv("OUT_DIR", file.path("results", "01a_coloc"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

WINDOW_KB       <- 500     # window = lead variant +/- WINDOW_KB
LEAD_P          <- 1e-6    # pass 1: a variant with P below this can become a lead variant
COLOC_P         <- 0.05    # a pair is tested only if min P of at least one dataset is <= this value in the window
MAF_MIN         <- 0.01
MIN_SHARED_SNPS <- 50      # minimum number of variants shared by the two datasets in a window
PPH4_THRESHOLD  <- 0.8     # screening flag only (tiers are assigned in 01b)
PRIOR_P1  <- 1e-4
PRIOR_P2  <- 1e-4
PRIOR_P12 <- 1e-5
set.seed(42)

## ---- datasets (file names, sample sizes and column names) ----------------------
liver_cols <- list(chr = "chromosome", pos = "base_pair_location", pval = "p_value",
                   beta = "beta", se = "standard_error", eaf = "effect_allele_frequency")
th_cols    <- list(chr = "chr", pos = "position", pval = "P.value", beta = "Effect", se = "StdErr", eaf = "Freq1")
GWAS_FILES <- list(
  PDFF_1 = list(file = "GCST90446650.tsv", alias = "PDFF_650", group = "PDFF", type = "quant", N = 36394, cols = liver_cols),
  PDFF_2 = list(file = "GCST90446651.tsv", alias = "PDFF_651", group = "PDFF", type = "quant", N = 35146, cols = liver_cols),
  PDFF_3 = list(file = "GCST90446652.tsv", alias = "PDFF_652", group = "PDFF", type = "quant", N = 34367, cols = liver_cols),
  PDFF_4 = list(file = "GCST90446653.tsv", alias = "PDFF_653", group = "PDFF", type = "quant", N = 33540, cols = liver_cols),
  cT1_1  = list(file = "GCST90446646.tsv", alias = "cT1_646",  group = "cT1",  type = "quant", N = 30481, cols = liver_cols),
  cT1_2  = list(file = "GCST90446647.tsv", alias = "cT1_647",  group = "cT1",  type = "quant", N = 29312, cols = liver_cols),
  cT1_3  = list(file = "GCST90446648.tsv", alias = "cT1_648",  group = "cT1",  type = "quant", N = 28638, cols = liver_cols),
  cT1_4  = list(file = "GCST90446649.tsv", alias = "cT1_649",  group = "cT1",  type = "quant", N = 27782, cols = liver_cols),
  FT3    = list(file = "FT3_summary_data.txt", alias = "FT3", group = "FT3", type = "quant", N = 57107,  cols = th_cols),
  FT4    = list(file = "FT4_summary_data.txt", alias = "FT4", group = "FT4", type = "quant", N = 117166, cols = th_cols),
  ## TSH dataset A has no variant ID column; variants are matched by position in all datasets
  TSH_A  =list(file = "GCST010653.txt", alias = "TSH_A", group = "TSH_A", type = "quant", N = 55342,
                cols = list(chr = "Chr", pos = "Pos", pval = "P-value", beta = "Effect", se = "SE", eaf = "Freq")),
  TSH_B  = list(file = "GCST90296333.tsv", alias = "TSH_B", group = "TSH_B", type = "quant", N = 124358,
                cols = list(chr = "chromosome", pos = "base_pair_location", pval = "p_value",
                            beta = "beta", se = "standard_error", eaf = "effect_allele_frequency"))
)

## ---- locate the input files by name ----------------------------------------------
## A file is found by its exact name; the GWAS Catalog files are also found when only the folder or the file name contains the
## accession (for example GCST90446650/GCST90446650_buildGRCh37.tsv.gz). Compressed files (.gz) are read with R.utils.
ROOTS <- unique(normalizePath(c(DATA_DIR, "data")[dir.exists(c(DATA_DIR, "data"))]))
ALL_FILES <- unique(normalizePath(unlist(lapply(ROOTS, function(r) list.files(r, recursive = TRUE, full.names = TRUE)))))
ALL_FILES <- ALL_FILES[!grepl("(^|/)results/", ALL_FILES) & grepl("\\.(tsv|txt|csv)(\\.gz)?$", ALL_FILES)]
find_file <- function(fname) {
  hits <- ALL_FILES[basename(ALL_FILES) == fname]
  acc  <- regmatches(fname, regexpr("^GCST[0-9]+", fname))
  if (!length(hits) && length(acc)) hits <- ALL_FILES[grepl(acc, ALL_FILES, fixed = TRUE)]
  if (!length(hits)) stop(sprintf("Input file not found: %s\n  Folders searched: %s\n  See data/README.md for the source of this file.",
                                  fname, paste(ROOTS, collapse = " | ")), call. = FALSE)
  if (length(hits) > 1) message("  more than one file for ", fname, "; using ", hits[1])
  hits[1]
}
for (k in names(GWAS_FILES)) GWAS_FILES[[k]]$path <- find_file(GWAS_FILES[[k]]$file)

## ---- helpers --------------------------------------------------------------------
col_map <- function(info) { cm <- unlist(info$cols); cm[!is.na(cm)] }
std_rename <- function(dt, cmap) {
  src <- cmap[cmap %in% names(dt)]
  if (length(src) > 0) setnames(dt, as.character(src), names(src))
  invisible(dt)
}

## Pass 1: significant variants (chromosome, position, P)
scan_sig_snps <- function(key, info, lead_p = LEAD_P) {
  message("  Scanning [", key, "]: ", info$file)
  cm <- col_map(info)
  need_src <- cm[names(cm) %in% c("chr", "pos", "pval")]
  dt <- fread(info$path, select = as.character(need_src), showProgress = FALSE)
  std_rename(dt, need_src)
  dt[, chr  := as.integer(gsub("chr|Chr|CHR", "", chr))]
  dt[, pos  := as.integer(pos)]
  dt[, pval := as.numeric(pval)]
  out <- dt[!is.na(chr) & !is.na(pos) & !is.na(pval) & pval > 0 & pval < lead_p, .(chr, pos, pval)]
  message("    ", format(nrow(out), big.mark = ","), " variants with P < ", lead_p)
  out
}

## Windows: distance-based greedy clumping (the most significant uncovered variant becomes the lead variant)
define_loci <- function(sig_all, window_kb = WINDOW_KB) {
  sig_all   <- unique(sig_all, by = c("chr", "pos"))
  window_bp <- as.integer(window_kb * 1000)
  message("  Unique significant positions: ", format(nrow(sig_all), big.mark = ","))
  loci <- rbindlist(lapply(split(sig_all, sig_all$chr), function(d) {
    d    <- d[order(pval)]
    used <- logical(nrow(d))
    out  <- list(); lid <- 1L
    for (i in seq_len(nrow(d))) {
      if (used[i]) next
      lp <- d$pos[i]
      used[abs(d$pos - lp) <= window_bp] <- TRUE
      out[[lid]] <- data.table(chr = d$chr[i], start = max(1L, lp - window_bp), end = lp + window_bp,
                               lead_snp = paste0(d$chr[i], ":", lp), lead_pos = lp, lead_pval = d$pval[i])
      lid <- lid + 1L
    }
    rbindlist(out)
  }))
  loci <- loci[order(chr, lead_pos)]
  loci[, locus_id := paste0("L", sprintf("%05d", .I))]
  setcolorder(loci, c("locus_id", "chr", "start", "end", "lead_snp", "lead_pos", "lead_pval"))
  message("  Windows defined: ", nrow(loci))
  loci
}

## Pass 2: variants inside the windows, all columns needed by coloc
load_gwas_in_loci <- function(key, info, loci) {
  message("  Loading [", key, "] restricted to the windows ...")
  cm <- col_map(info)
  header <- names(fread(info$path, nrows = 0L))
  src_avail <- cm[cm %in% header]
  dt <- fread(info$path, select = as.character(src_avail), showProgress = FALSE)
  std_rename(dt, src_avail)
  dt[, chr  := as.integer(gsub("chr|Chr|CHR", "", chr))]
  dt[, pos  := as.integer(pos)]
  dt[, pval := as.numeric(pval)]
  dt[, beta := as.numeric(beta)]
  dt[, se   := as.numeric(se)]
  dt[, eaf  := as.numeric(eaf)]
  dt <- dt[!is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(se) & se > 0 &
           !is.na(pval) & pval > 0 & pval <= 1 & !is.na(eaf)]
  dt[, maf     := ifelse(eaf > 0.5, 1 - eaf, eaf)]
  dt[, varbeta := se^2]
  dt <- dt[maf >= MAF_MIN]
  dt <- dt[order(chr, pos, pval)][!duplicated(paste(chr, pos))]
  dt[, pos_end := pos]
  setkey(dt, chr, pos, pos_end)
  loci_key <- loci[, .(chr, start, end)]
  setkey(loci_key, chr, start, end)
  dt_in <- foverlaps(dt, loci_key, by.x = c("chr", "pos", "pos_end"), by.y = c("chr", "start", "end"),
                     nomatch = 0L, type = "any")
  rm(dt); gc()
  dt_in <- unique(dt_in, by = c("chr", "pos"))
  dt_in <- dt_in[, .(chr, pos, pval, beta, varbeta, maf)]
  message("    ", format(nrow(dt_in), big.mark = ","), " variants in the windows")
  dt_in
}

## coloc.abf for one window and one dataset pair
make_ds <- function(dt, info) list(snp = as.character(dt$pos), beta = dt$beta, varbeta = dt$varbeta,
                                   N = info$N, MAF = dt$maf, type = info$type)
run_coloc_single <- function(dt1, dt2, info1, info2, locus) {
  s1 <- dt1[chr == locus$chr & pos >= locus$start & pos <= locus$end]
  s2 <- dt2[chr == locus$chr & pos >= locus$start & pos <= locus$end]
  common <- intersect(s1$pos, s2$pos)
  if (length(common) < MIN_SHARED_SNPS) return(NULL)
  s1 <- s1[pos %in% common][order(pos)]
  s2 <- s2[pos %in% common][order(pos)]
  if (min(s1$pval) > COLOC_P && min(s2$pval) > COLOC_P) return(NULL)
  res <- tryCatch(
    coloc.abf(dataset1 = make_ds(s1, info1), dataset2 = make_ds(s2, info2),
              p1 = PRIOR_P1, p2 = PRIOR_P2, p12 = PRIOR_P12),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      coloc.abf(make_ds(s1, info1), make_ds(s2, info2), p1 = PRIOR_P1, p2 = PRIOR_P2, p12 = PRIOR_P12)))
  if (is.null(res)) return(NULL)
  sm <- res$summary
  data.table(locus_id = locus$locus_id, chr = locus$chr, start_bp = locus$start, end_bp = locus$end,
             lead_snp = locus$lead_snp, lead_pval = locus$lead_pval,
             trait1 = info1$alias, trait2 = info2$alias, group1 = info1$group, group2 = info2$group,
             nSNPs = as.integer(sm["nsnps"]),
             PP.H0 = round(sm["PP.H0.abf"], 5), PP.H1 = round(sm["PP.H1.abf"], 5), PP.H2 = round(sm["PP.H2.abf"], 5),
             PP.H3 = round(sm["PP.H3.abf"], 5), PP.H4 = round(sm["PP.H4.abf"], 5),
             coloc = sm["PP.H4.abf"] >= PPH4_THRESHOLD)
}

## ---- main -------------------------------------------------------------------------
message("[1/4] Pass 1: collecting significant variants")
all_sig <- rbindlist(lapply(names(GWAS_FILES), function(k) scan_sig_snps(k, GWAS_FILES[[k]])))
gc()

message("\n[2/4] Defining windows")
loci <- define_loci(all_sig)
rm(all_sig); gc()
fwrite(loci[, .(window = locus_id, chr, start_GRCh37 = start, end_GRCh37 = end, lead_variant = lead_snp,
                lead_pos, lead_P = lead_pval)],
       file.path(OUT_DIR, "coloc_windows.tsv"), sep = "\t")

message("\n[3/4] Pass 2: loading the variants inside the windows")
gwas_data <- lapply(names(GWAS_FILES), function(k) { dt <- load_gwas_in_loci(k, GWAS_FILES[[k]], loci); gc(); dt })
names(gwas_data) <- names(GWAS_FILES)

message("\n[4/4] coloc.abf for all windows and dataset pairs")
all_pairs  <- combn(names(GWAS_FILES), 2, simplify = FALSE)
total_runs <- nrow(loci) * length(all_pairs)
message("  Windows: ", nrow(loci), "  |  dataset pairs: ", length(all_pairs), "  |  possible tests: ", format(total_runs, big.mark = ","))
all_results <- list(); run_idx <- 0L
for (pair in all_pairs) {
  k1 <- pair[1]; k2 <- pair[2]
  info1 <- GWAS_FILES[[k1]]; info2 <- GWAS_FILES[[k2]]
  for (i in seq_len(nrow(loci))) {
    run_idx <- run_idx + 1L
    if (run_idx %% 5000 == 0 || run_idx == 1L)
      cat(sprintf("  [%d/%d] %s x %s | window %s\n", run_idx, total_runs, info1$alias, info2$alias, loci$locus_id[i]))
    res <- run_coloc_single(gwas_data[[k1]], gwas_data[[k2]], info1, info2, loci[i])
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
}
raw_coloc <- rbindlist(Filter(Negate(is.null), all_results), fill = TRUE)
message("  Completed tests: ", format(nrow(raw_coloc), big.mark = ","))
fwrite(raw_coloc, file.path(OUT_DIR, "coloc_pairs.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "01a_coloc_scan_sessionInfo.txt"))
message("\nDone. Output: ", normalizePath(OUT_DIR))
