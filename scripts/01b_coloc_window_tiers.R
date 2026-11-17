################################################################################
#  01b_coloc_window_tiers.R
#  Window summary of the colocalization scan: highest cross-domain PP.H4, tier and prior sensitivity of every window.
#  Output = Supplementary Table S2a (458 windows; 11 high, 16 moderate, 431 below 0.5).
#
#  METHODS (as in the manuscript)
#    - Only cross-domain pairs are used (a liver MRI dataset against a thyroid dataset: 4 thyroid x 8 liver = 32 pairs).
#      The 34 within-domain pairs (6 thyroid-thyroid, 28 liver-liver) do not enter the window judgment.
#    - TSH has two datasets (A, B). For each liver dataset the lower PP.H4 of the two TSH datasets is used; both must have been
#      tested. FT3 and FT4 have one dataset each.
#    - For each trait-level pair (PDFF or cT1 against TSH, FT3 or FT4) the maximum over the 4 adjustment models of the liver
#      trait is taken; the window value is the maximum over the 6 trait-level pairs.
#    - tier: high if the window value is >= 0.8, moderate if 0.5 - 0.8, otherwise "no". supporting_pairs lists the trait-level
#      pairs with PP.H4 >= 0.5. A window without any cross-domain pair has an empty value and tier "no".
#    - Prior sensitivity (supported windows only): with p1 and p2 fixed, only H4 scales with p12, so PP.H4 is recalculated from the
#      baseline posterior (p12 = 1e-5) as PP.H4 * r / (1 - PP.H4 + PP.H4 * r), r = new p12 / 1e-5, for p12 = 5e-6, 1e-5, 5e-5, 1e-4;
#      the window rule above is then applied again.
#
#  INPUTS  (one of two sources; SOURCE=auto uses the scan output if it is found, otherwise the Supplementary Tables)
#    scan   coloc_windows.tsv and coloc_pairs.tsv from 01a_coloc_scan.R (searched by name in DATA_DIR, OUT_DIR and results/01a_coloc)
#    supp   Antioxidant_ThyroidLiver_Supplementary_Tables_S1-S13.xlsx (file name must contain "Supplementary_Tables_S1-S13"),
#           tables S2a (windows) and S2b (all dataset pairs); no other data are needed
#
#  OUTPUTS (OUT_DIR, default results/01b_tiers)
#    window_tiers.tsv                 window, chr, start_GRCh37, end_GRCh37, lead_variant, lead_P, tier, max_cross_domain_PP.H4,
#                                     supporting_pairs, PP.H4_at_p12_5e-06 ... PP.H4_at_p12_0.0001   (= Table S2a)
#    cross_domain_dataset_pairs.tsv   the 32 cross-domain pairs with the TSH rule applied (window, liver dataset, thyroid trait, PP.H4)
#    01b_coloc_window_tiers_sessionInfo.txt
#  Check printed to the console: the result is compared with Table S2a when the workbook is available, and the scan output is
#  compared with Table S2b when both are available.
#
#  REPRODUCTION CHECK  Run on Table S2b, the script gives the maximum PP.H4, tier and supporting pairs of all 458 windows of
#  Table S2a and the prior-sensitivity values of the 27 supported windows (see docs/PROVENANCE.md).
#
#  USAGE  DATA_DIR=/path/to/data Rscript 01b_coloc_window_tiers.R          (SOURCE=supp or SOURCE=scan to force one source)
################################################################################
invisible(try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(readxl) })

## ---- settings ----------------------------------------------------------------
DATA_DIR <- Sys.getenv("DATA_DIR", "data")
OUT_DIR  <- Sys.getenv("OUT_DIR", file.path("results", "01b_tiers"))
SOURCE   <- Sys.getenv("SOURCE", "auto")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
HIGH <- 0.8; MODERATE <- 0.5
BASE_P12 <- 1e-5
P12S <- c(5e-6, 1e-5, 5e-5, 1e-4)

## ---- locate input files by name ---------------------------------------------------
ROOTS <- unique(c(DATA_DIR, OUT_DIR, file.path("results", "01a_coloc"), "data", getwd())); ROOTS <- ROOTS[dir.exists(ROOTS)]
find_quiet <- function(pattern) {
  hits <- unlist(lapply(ROOTS, function(r) list.files(r, pattern = pattern, recursive = TRUE, full.names = TRUE)))
  hits <- hits[!grepl("^~\\$", basename(hits))]
  if (length(hits)) hits[which.max(file.mtime(hits))] else NA_character_
}
SUPP <- { e <- Sys.getenv("SUPP_XLSX", ""); if (nzchar(e)) e else find_quiet("Supplementary_Tables_S1-S13.*\\.xlsx$") }
F_WIN  <- find_quiet("^coloc_windows\\.tsv$")
F_PAIR <- find_quiet("^coloc_pairs\\.tsv$")
if (SOURCE == "auto") SOURCE <- if (!is.na(F_WIN) && !is.na(F_PAIR)) "scan" else "supp"
if (SOURCE == "supp" && is.na(SUPP))
  stop("No input found. Put the Supplementary Tables workbook (file name containing 'Supplementary_Tables_S1-S13') in DATA_DIR, ",
       "or run 01a_coloc_scan.R first.\n  Folders searched: ", paste(ROOTS, collapse = " | "), call. = FALSE)
message("Source: ", SOURCE)

## ---- Supplementary Tables reader (two layouts: one sheet per table, or one sheet per table number with lettered tables side by side)
if (!is.na(SUPP)) {
  SUPP_SHEETS <- excel_sheets(SUPP)
  locate_supp <- function(key) {
    sh <- SUPP_SHEETS[SUPP_SHEETS %in% c(key, paste("Table", key))]
    if (length(sh) == 1) return(list(sheet = sh, c0 = 1L, c1 = NA_integer_))
    sh <- SUPP_SHEETS[SUPP_SHEETS %in% paste("Table", sub("^(S[0-9]+).*$", "\\1", key))]
    if (length(sh) != 1) return(NULL)
    r1 <- unlist(read_excel(SUPP, sheet = sh, col_names = FALSE, n_max = 1, col_types = "text", .name_repair = "minimal"), use.names = FALSE)
    c0 <- which(startsWith(r1, paste0("Table ", key, ".")))[1]
    if (is.na(c0)) return(NULL)
    nxt <- which(grepl("^Table S[0-9]+[a-z]?\\.", r1)); nxt <- nxt[nxt > c0]
    list(sheet = sh, c0 = c0, c1 = if (length(nxt)) nxt[1] - 1L else NA_integer_)
  }
  read_supp <- function(key, first_col) {
    loc <- locate_supp(key)
    if (is.null(loc)) stop(sprintf("Table %s not found in %s", key, basename(SUPP)), call. = FALSE)
    rng <- function(r0, r1) readxl::cell_limits(c(r0, loc$c0), c(r1, loc$c1))
    top <- read_excel(SUPP, sheet = loc$sheet, range = rng(1, 25), col_names = FALSE, col_types = "text", .name_repair = "minimal")
    h <- which(trimws(as.character(top[[1]])) == first_col)[1]
    if (is.na(h)) stop(sprintf("Header row (first cell '%s') not found for Table %s", first_col, key), call. = FALSE)
    d <- as.data.table(read_excel(SUPP, sheet = loc$sheet, range = rng(h, NA), .name_repair = "minimal"))
    blank <- (is.na(names(d)) | names(d) == "" | grepl("^\\.\\.\\.[0-9]+$", names(d))) & vapply(d, function(x) all(is.na(x)), logical(1))
    d <- d[, which(!blank), with = FALSE]
    d[rowSums(!is.na(d)) > 0]
  }
}

## ---- windows and dataset pairs in a common format ---------------------------------
is_liver <- function(x) grepl("^(PDFF|cT1)", x)
is_thy   <- function(x) x %in% c("TSH_A", "TSH_B", "FT3", "FT4")
if (SOURCE == "scan") {
  W <- fread(F_WIN)[, .(window, chr = as.character(chr), start_GRCh37, end_GRCh37, lead_variant, lead_P)]
  P <- fread(F_PAIR)[, .(window = locus_id, dataset_1 = trait1, dataset_2 = trait2, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)]
} else {
  W <- read_supp("S2a", "window")[, .(window, chr = as.character(chr), start_GRCh37 = as.numeric(start_GRCh37), end_GRCh37 = as.numeric(end_GRCh37),
                                      lead_variant, lead_P = as.character(lead_P))]
  P <- read_supp("S2b", "window")[, .(window, dataset_1, dataset_2, PP.H0 = as.numeric(PP.H0), PP.H1 = as.numeric(PP.H1),
                                      PP.H2 = as.numeric(PP.H2), PP.H3 = as.numeric(PP.H3), PP.H4 = as.numeric(PP.H4))]
}
message(sprintf("  %d windows, %d tested dataset pairs", nrow(W), nrow(P)))
P[, pair_type := fifelse((is_liver(dataset_1) & is_thy(dataset_2)) | (is_liver(dataset_2) & is_thy(dataset_1)), "cross-domain", "within-domain")]
message(sprintf("  cross-domain pairs: %d, within-domain pairs: %d", sum(P$pair_type == "cross-domain"), sum(P$pair_type == "within-domain")))

## ---- window rule ---------------------------------------------------------------------
## Returns the window value (best), its trait-level pairs (tl) and the cross-domain table after the TSH rule (liver dataset level)
window_rule <- function(pairs) {
  cd <- pairs[pair_type == "cross-domain"]
  cd[, `:=`(liver = fifelse(is_liver(dataset_1), dataset_1, dataset_2), thy = fifelse(is_thy(dataset_1), dataset_1, dataset_2))]
  cd[, `:=`(h4 = PP.H4, lt = sub("_.*", "", liver), tt = fifelse(thy %in% c("TSH_A", "TSH_B"), "TSH", thy))]
  keep <- c("window", "liver", "lt", "tt", "h4")
  tsh <- cd[tt == "TSH"]; tsh[, n := .N, by = .(window, liver)]
  tsh <- tsh[n == 2][, .SD[which.min(h4)], by = .(window, liver)]            # lower of the two TSH datasets, both required
  both <- rbind(cd[tt != "TSH", ..keep], tsh[, ..keep])
  tl   <- both[, .SD[which.max(h4)], by = .(window, lt, tt)]                  # maximum over the 4 adjustment models
  best <- tl[, .SD[which.max(h4)], by = window]                                # maximum over the 6 trait-level pairs
  list(best = best, tl = tl, both = both)
}
base <- window_rule(P)
fwrite(base$both[order(window, liver, tt), .(window, liver_dataset = liver, thyroid_trait = tt, PP.H4 = h4)],
       file.path(OUT_DIR, "cross_domain_dataset_pairs.tsv"), sep = "\t")

tier_of <- function(x) fifelse(is.na(x), "no", fifelse(x >= HIGH, "high", fifelse(x >= MODERATE, "moderate", "no")))
sup_pairs <- base$tl[h4 >= MODERATE, .(supporting_pairs = paste(
  { v <- paste0(lt, "_x_", tt); v[order(tolower(v))] }, collapse = "; ")), by = window]
S <- merge(W, base$best[, .(window, max_cross_domain_PP.H4 = h4)], by = "window", all.x = TRUE, sort = FALSE)
S <- merge(S, sup_pairs, by = "window", all.x = TRUE, sort = FALSE)
S[, tier := tier_of(max_cross_domain_PP.H4)]

## prior sensitivity (supported windows)
sens <- rbindlist(lapply(P12S, function(p12) {
  d <- copy(P); r <- p12 / BASE_P12
  d[, PP.H4 := PP.H4 * r / (1 - PP.H4 + PP.H4 * r)]
  b <- window_rule(d)$best
  b[, .(window, p12 = p12, h4)]
}))
sens_w <- dcast(sens, window ~ p12, value.var = "h4")          # columns in increasing order of p12
stopifnot(!is.unsorted(P12S), ncol(sens_w) == length(P12S) + 1)
setnames(sens_w, 2:ncol(sens_w), paste0("PP.H4_at_p12_", c("5e-06", "1e-05", "5e-05", "0.0001")))   # names as in Table S2a
S <- merge(S, sens_w, by = "window", all.x = TRUE, sort = FALSE)
for (cn in grep("^PP.H4_at_p12_", names(S), value = TRUE)) { S[tier == "no", (cn) := NA_real_]; S[, (cn) := round(get(cn), 5)] }
S <- S[order(window)]
fwrite(S, file.path(OUT_DIR, "window_tiers.tsv"), sep = "\t", na = "")
message(sprintf("\nTiers: high %d, moderate %d, no %d (windows: %d)", sum(S$tier == "high"), sum(S$tier == "moderate"), sum(S$tier == "no"), nrow(S)))

## ---- checks against the Supplementary Tables ------------------------------------------
if (!is.na(SUPP)) {
  a <- read_supp("S2a", "window")
  m <- merge(S, a, by = "window", suffixes = c("", ".supp"))
  n_tier <- sum(m$tier != m$tier.supp)
  mx <- as.numeric(m$`max_cross_domain_PP.H4.supp`)
  ## Table S2a shows 4 decimals, the posterior probabilities of S2b have 5: a difference below 1e-4 is rounding of the last digit
  n_max <- sum(!((is.na(m$max_cross_domain_PP.H4) & is.na(mx)) | (!is.na(mx) & !is.na(m$max_cross_domain_PP.H4) & abs(m$max_cross_domain_PP.H4 - mx) < 1e-4)))
  n_sup <- sum(fifelse(is.na(m$supporting_pairs), "", m$supporting_pairs) != fifelse(is.na(m$supporting_pairs.supp), "", m$supporting_pairs.supp))
  message(sprintf("Check against Table S2a (%d windows): tier differs in %d, max cross-domain PP.H4 differs (beyond rounding) in %d, supporting pairs differ in %d",
                  nrow(m), n_tier, n_max, n_sup))
  for (cn in grep("^PP.H4_at_p12_", names(S), value = TRUE)) {
    x <- m[[cn]]; y <- as.numeric(m[[paste0(cn, ".supp")]])
    ok <- (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & abs(round(x, 4) - round(y, 4)) < 1.5e-4)
    message(sprintf("  %s: %d of %d windows differ", cn, sum(!ok), nrow(m)))
  }
  if (SOURCE == "scan") {
    b <- read_supp("S2b", "window")
    b[, key := paste(window, dataset_1, dataset_2)]; P[, key := paste(window, dataset_1, dataset_2)]
    mm <- merge(P, b, by = "key", suffixes = c("", ".supp"))
    dmax <- max(abs(mm$PP.H4 - as.numeric(mm$PP.H4.supp)), abs(mm$PP.H0 - as.numeric(mm$PP.H0.supp)))
    message(sprintf("Check of the scan against Table S2b: %d of %d scan rows found in S2b (S2b has %d); largest difference in PP.H0 / PP.H4: %g",
                    nrow(mm), nrow(P), nrow(b), dmax))
  }
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "01b_coloc_window_tiers_sessionInfo.txt"))
message("Done. Output: ", normalizePath(OUT_DIR))
