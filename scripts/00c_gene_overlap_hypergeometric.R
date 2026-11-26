################################################################################
#  00c_gene_overlap_hypergeometric.R
#  Exploratory test of the overlap between the thyroid-derived and the liver-derived GWAS gene sets (shared genes).
#
#  STATUS  This test was run during the analysis (meeting level) and is NOT part of the submitted Supplementary Tables:
#          Table S3a gives descriptive counts only ("no significance test").  The test reported in the manuscript for the
#          relation of the two gene sets is the competitive MAGMA gene-set analysis (00b_magma_gene_set_enrichment.R, Table S1).
#          Hypergeometric and permutation tests assume independent gene draws and ignore LD and gene density, so they are
#          anti-conservative; the numbers are given for orientation only.
#
#  TESTS
#    (1) hypergeometric test  P(X >= k) = phyper(k-1, K, N-K, n, lower.tail = FALSE)  for background sizes N = 19,000 / 20,000 / 25,000
#    (2) permutation test: sets of size K and n are drawn at random from N genes, 100,000 times; empirical p = (#overlap >= k + 1) / (iterations + 1)
#        Fixing one set and drawing only the other gives the same overlap distribution as drawing both, and is much faster.
#        If no draw reaches the observed overlap the empirical p equals 1 / 100,001 = 1e-5 and must be read as "p < 1e-5".
#
#  INPUTS  (searched by name in DATA_DIR; MODE=genes|supp forces one)
#    genes  gwas_gene_flags.tsv (data/gene_sets/, or the output of 00_gwas_gene_sets_from_FUMA.R):  K = thyroid-derived genes,
#           n = liver-derived genes, k = shared genes (1,709 / 1,120 / 199 symbols in the flags file; ABHD16B and C20ORF135 are one gene, i.e. 1,119 liver genes as in Table S3a)
#    supp   Antioxidant_ThyroidLiver_Supplementary_Tables_S1-S13.xlsx, Table S3a: K = 1,709, n = 1,119, k = 199 (as printed in the table)
#  OUTPUTS (OUT_DIR, default results/00c_overlap)
#    gene_overlap_hypergeometric.tsv      N, K, n, k_obs, expected, fold_enrichment, p_value (+ permutation p, same for every N row)
#    permutation_overlap_null_distribution.png
#    00c_gene_overlap_hypergeometric_sessionInfo.txt
#  EXPECTED (genes mode, flags of the manuscript): N = 19,000 -> expected 100.74, fold 1.975, p 1.03e-21;  N = 20,000 -> 95.70, 2.079, 2.46e-24;
#                                                  N = 25,000 -> 76.56, 2.599, 7.26e-37;  permutation p < 1e-5
#
#  USAGE   DATA_DIR=/path/to/data Rscript 00c_gene_overlap_hypergeometric.R
################################################################################
invisible(try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(readxl) })

DATA_DIR <- Sys.getenv("DATA_DIR", "data")
OUT_DIR  <- Sys.getenv("OUT_DIR", file.path("results", "00c_overlap"))
MODE     <- Sys.getenv("MODE", "auto")
N_PERM   <- as.integer(Sys.getenv("N_PERM", "100000"))
N_BACKGROUND <- c(19000, 20000, 25000)
N_PERM_BG    <- 20000                                   # background size of the permutation test
set.seed(20260916)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ROOTS <- unique(c(DATA_DIR, "data", getwd())); ROOTS <- ROOTS[dir.exists(ROOTS)]
find_quiet <- function(pattern) {
  hits <- unlist(lapply(ROOTS, function(r) list.files(r, pattern = pattern, recursive = TRUE, full.names = TRUE)))
  hits <- hits[!grepl("^~\\$", basename(hits)) & !grepl("(^|/)results/", hits)]
  if (length(hits)) hits[1] else NA_character_
}
F_FLAGS <- find_quiet("^gwas_gene_flags\\.tsv$")
SUPP    <- { e <- Sys.getenv("SUPP_XLSX", ""); if (nzchar(e)) e else find_quiet("Supplementary_Tables_S1-S13.*\\.xlsx$") }
if (MODE == "auto") MODE <- if (!is.na(F_FLAGS)) "genes" else "supp"

if (MODE == "genes") {
  if (is.na(F_FLAGS)) stop("gwas_gene_flags.tsv not found in ", paste(ROOTS, collapse = " | "), call. = FALSE)
  fl <- fread(F_FLAGS)
  thy <- unique(fl[GWAS_TSH == 1 | GWAS_FT3 == 1 | GWAS_FT4 == 1, Gene])
  liv <- unique(fl[GWAS_PDFF == 1 | GWAS_cT1 == 1, Gene])
  K <- length(thy); n <- length(liv); k_obs <- length(intersect(thy, liv))
  message("Source: ", F_FLAGS)
} else {
  if (is.na(SUPP)) stop("Neither gwas_gene_flags.tsv nor the Supplementary Tables workbook found in ", paste(ROOTS, collapse = " | "), call. = FALSE)
  sh <- excel_sheets(SUPP); sh <- sh[sh %in% c("S3", "Table S3")]
  s3 <- as.data.table(read_excel(SUPP, sheet = sh, range = readxl::cell_limits(c(4, 1), c(12, 3)), .name_repair = "minimal"))
  K <- s3[item == "Thyroid-derived GWAS gene set", genes]; n <- s3[item == "Liver-derived GWAS gene set", genes]; k_obs <- s3[item == "Shared (in both sets)", genes]
  message("Source: Table S3a of ", basename(SUPP))
}
cat(sprintf("Thyroid-derived genes (K) = %d\nLiver-derived genes (n)   = %d\nShared genes (k_obs)      = %d\n\n", K, n, k_obs))

## ---- (1) hypergeometric test ------------------------------------------------------
hyper <- rbindlist(lapply(N_BACKGROUND, function(N) {
  expected <- K * n / N
  data.table(N = N, K = K, n = n, k_obs = k_obs, expected = round(expected, 2), fold_enrichment = round(k_obs / expected, 3),
             p_value = phyper(k_obs - 1, K, N - K, n, lower.tail = FALSE))
}))
cat("(1) Hypergeometric test, sensitivity to the background gene number N\n"); print(hyper, row.names = FALSE)

## ---- (2) permutation test -------------------------------------------------------------
overlaps <- integer(N_PERM)
for (i in seq_len(N_PERM)) overlaps[i] <- sum(sample.int(N_PERM_BG, n) <= K)       # genes 1..K = the fixed thyroid set
emp_p <- (sum(overlaps >= k_obs) + 1) / (N_PERM + 1)
cat(sprintf("\n(2) Permutation test (N = %d, %d draws): mean overlap %.2f, maximum %d, observed %d, empirical p = %s%s\n",
            N_PERM_BG, N_PERM, mean(overlaps), max(overlaps), k_obs, format(emp_p, digits = 3),
            if (all(overlaps < k_obs)) sprintf("  (no draw reached the observed overlap: report as p < %s)", format(1 / (N_PERM + 1), digits = 1)) else ""))
hyper[, permutation_empirical_p := emp_p]
fwrite(hyper, file.path(OUT_DIR, "gene_overlap_hypergeometric.tsv"), sep = "\t")

png(file.path(OUT_DIR, "permutation_overlap_null_distribution.png"), width = 1600, height = 1000, res = 150)
hist(overlaps, breaks = 40, col = "#C8C8C8", border = "white", xlim = c(min(overlaps), max(k_obs, max(overlaps)) + 25),
     main = sprintf("Permutation null of the overlap (N = %s, %s draws)", format(N_PERM_BG, big.mark = ","), format(N_PERM, big.mark = ",")),
     xlab = "Number of shared genes in a random draw", ylab = "Frequency")
abline(v = k_obs, col = "#E67E22", lwd = 3, lty = 2)
text(k_obs, par("usr")[4] * 0.9, labels = sprintf("observed = %d", k_obs), col = "#E67E22", pos = 2, font = 2)
invisible(dev.off())
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "00c_gene_overlap_hypergeometric_sessionInfo.txt"))
message("Done. Output: ", normalizePath(OUT_DIR))
