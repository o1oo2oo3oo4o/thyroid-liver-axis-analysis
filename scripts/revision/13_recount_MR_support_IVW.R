#!/usr/bin/env Rscript
# =============================================================================
# 13_recount_MR_support_IVW.R
#
# Purpose
#   Recount, from the Supplementary Tables, every number of the manuscript that depends on the definition of
#   "genetic support" by Mendelian randomization (MR):
#
#       A protein-outcome pair has genetic support when its inverse-variance weighted (IVW) estimate has FDR P < 0.05.
#       The weighted-median (WM) and MR-Egger results are sensitivity analyses; they are not used for the definition.
#
#   Nothing is re-estimated (no MR, no colocalization): the script reads the tables that hold the estimates and counts.
#
#     [1] Table S7   IVW pairs with FDR P < 0.05: pairs, proteins, outcome groups, direction, Steiger P, MR-Egger intercept P
#     [2] Table S7   the five genes of the original submission whose support came only from the weighted median
#                    (ADH4, BTC, COL11A2, HSPA1A, LGALS4)
#     [3] Table S14  colocalization status and MR-RAPS status of the same pairs
#     [4] Table S16  instrument strength: F statistics of the unique instruments, instruments per pair
#     [5] Table S15  the summary table is compared with the values recomputed from S16 and S7
#
# Run from scripts/revision (the folder that contains config.R, data/ and results/):
#   Rscript 13_recount_MR_support_IVW.R
#   SUPP_XLSX="path/to/Supplementary_Tables.xlsx" Rscript 13_recount_MR_support_IVW.R
#
# Input   the revised Supplementary Tables workbook (config.R: supp_revised_file; or the environment variable SUPP_XLSX) with the
#         sheets "Table S7" (header row "pQTL source / Protein / MR method / Beta (95% CI) / FDR P / Pleiotropy P / Steiger P",
#         one block of four columns per outcome, the outcome name above the first column of each block), "Table S14",
#         "Table S15" and "Table S16"
# Output  results/13_MR_support_IVW/
#           1_ivw_significant_pairs.csv       one row per IVW pair with FDR P < 0.05 (669 rows)
#           2_protein_summary.csv             one row per protein: pairs, sources, outcome groups, direction
#           3_count_summary.csv               every recounted number next to the value printed in the manuscript
#           4_WM_only_genes.csv               ADH4, BTC, COL11A2, HSPA1A, LGALS4: IVW and WM significant pairs
#           5_F_statistics_recomputed.csv     F-statistic summary (overall, by source) recomputed from Table S16
#           sessionInfo.txt
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(readxl) })

cfg <- new.env(); sys.source("config.R", envir = cfg)
supp_file <- Sys.getenv("SUPP_XLSX", unset = cfg$supp_revised_file)
if (!file.exists(supp_file)) stop("Supplementary Tables workbook not found: ", supp_file,
                                  "\n  set supp_revised_file in config.R or the environment variable SUPP_XLSX", call. = FALSE)
out_dir <- file.path(cfg$out_root, "13_MR_support_IVW"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("workbook:", supp_file, "\n")

## ---- helpers -------------------------------------------------------------------------------------------
read_raw <- function(sheet) {                         # a sheet as text, no header
  as.data.frame(read_excel(supp_file, sheet = sheet, col_names = FALSE, col_types = "text", .name_repair = "minimal"),
                stringsAsFactors = FALSE)
}
fill_fwd <- function(x) { x <- trimws(x); x[!nzchar(x)] <- NA; for (i in seq_along(x)[-1]) if (is.na(x[i])) x[i] <- x[i - 1]; x }
pnum <- function(x) suppressWarnings(as.numeric(gsub("[*§<]", "", trimws(x))))   # "4.47E-02*" -> 0.0447, "<1.00E-16" -> 1e-16
beta_part <- function(x, which) {                    # "-0.05 (-0.09, -0.01)": which = 1 beta, 2 lower, 3 upper
  m <- regmatches(x, regexec("^\\s*(-?[0-9.]+(?:[eE][+-]?[0-9]+)?)\\s*\\(\\s*(-?[0-9.]+(?:[eE][+-]?[0-9]+)?)\\s*,\\s*(-?[0-9.]+(?:[eE][+-]?[0-9]+)?)\\s*\\)", x, perl = TRUE))
  vapply(m, function(v) if (length(v) == 4) as.numeric(v[which + 1]) else NA_real_, numeric(1))
}
round1 <- function(x) floor(x * 10 + 0.5 + 1e-8) / 10                      # one decimal, halves rounded up (as printed in Table S15)
quant <- function(x) c(min = min(x), p5 = unname(quantile(x, .05)), q1 = unname(quantile(x, .25)), median = median(x),
                       q3 = unname(quantile(x, .75)), p95 = unname(quantile(x, .95)), max = max(x))

STEATOTIC <- c("Liver PDFFnon", "Liver PDFFbmi", "Liver PDFFwfm", "Liver PDFFvat", "MASLD", "FLI")
FIBROTIC  <- c("Liver cT1non", "Liver cT1bmi", "Liver cT1wfm", "Liver cT1vat", "FIB4", "APRI")
FIVE      <- c("ADH4", "BTC", "COL11A2", "HSPA1A", "LGALS4")

## ---- [1] Table S7 --------------------------------------------------------------------------------------
M  <- read_raw("Table S7")
hr <- which(trimws(M[[1]]) == "pQTL source" & trimws(M[[2]]) == "Protein")[1]
if (is.na(hr)) stop("Table S7: header row ('pQTL source', 'Protein', 'MR method', ...) not found", call. = FALSE)
lab_row <- unlist(M[hr - 1, ], use.names = FALSE)
bstart  <- which(!is.na(lab_row) & nzchar(trimws(lab_row)) & seq_along(lab_row) >= 4)       # first column of each outcome block
stopifnot(all(trimws(unlist(M[hr, bstart + 1], use.names = FALSE)) == "FDR P"))
outcomes <- trimws(lab_row[bstart])
cat("outcomes in Table S7 (", length(outcomes), "): ", paste(outcomes, collapse = ", "), "\n", sep = "")
stopifnot(all(c(STEATOTIC, FIBROTIC) %in% outcomes))

B <- M[(hr + 1):nrow(M), , drop = FALSE]
B[[1]] <- fill_fwd(B[[1]]); B[[2]] <- fill_fwd(B[[2]]); B[[3]] <- trimws(B[[3]])
B <- B[B[[1]] %in% c("deCODE", "ARIC") & B[[3]] %in% c("IVW", "WM"), , drop = FALSE]
cat("Table S7 rows (source x protein x method):", nrow(B), "\n")

long <- rbindlist(lapply(seq_along(bstart), function(k) {
  j <- bstart[k]
  data.table(source = B[[1]], protein = B[[2]], method = B[[3]], outcome = outcomes[k],
             beta_ci = trimws(B[[j]]), fdr_txt = trimws(B[[j + 1]]), fdr = pnum(B[[j + 1]]),
             egger_p = pnum(B[[j + 2]]), steiger_p = pnum(B[[j + 3]]))
}))
long[, star := grepl("\\*\\s*$", fdr_txt)]
long[, sig  := !is.na(fdr) & fdr < 0.05]
cat(sprintf("FDR < 0.05 and asterisk agree in %d of %d cells\n", long[, sum(sig == star)], nrow(long)))
long[, beta := beta_part(beta_ci, 1)]; long[, lo := beta_part(beta_ci, 2)]; long[, hi := beta_part(beta_ci, 3)]
long[, dir_sign := ifelse(!is.na(beta) & beta != 0, sign(beta), sign((lo + hi) / 2))]    # beta printed as 0.00: midpoint of the CI
long[, group := fifelse(outcome %in% STEATOTIC, "steatotic", fifelse(outcome %in% FIBROTIC, "fibrotic", NA_character_))]

ivw <- long[method == "IVW"]
tested   <- ivw[, uniqueN(protein)]
sigp     <- ivw[sig == TRUE]
n_pairs  <- nrow(sigp)
prot_all <- unique(sigp$protein)
per_src  <- sigp[, .(pairs = .N, proteins = uniqueN(protein)), by = source]
cat(sprintf("\n[1] proteins with an IVW estimate: %d; IVW pairs with FDR < 0.05: %d (deCODE %d, ARIC %d); proteins with >= 1 such pair: %d\n",
            tested, n_pairs, per_src[source == "deCODE", pairs], per_src[source == "ARIC", pairs], length(prot_all)))

S_set <- unique(sigp[group == "steatotic", protein]); F_set <- unique(sigp[group == "fibrotic", protein])
same_dir <- function(grp, set) sum(vapply(set, function(p) uniqueN(sigp[protein == p & group == grp, dir_sign]) == 1, logical(1)))
cat(sprintf("    steatotic outcomes %d proteins, fibrotic outcomes %d proteins, both %d; one direction only: %d and %d\n",
            length(S_set), length(F_set), length(intersect(S_set, F_set)), same_dir("steatotic", S_set), same_dir("fibrotic", F_set)))

# MR-Egger rule of the earlier definition (not used any more; kept to state how many proteins would be lost)
egger_ok  <- sigp[!is.na(egger_p) & egger_p >= 0.05]
prot_egg  <- unique(egger_ok$protein)
lost_egg  <- setdiff(prot_all, prot_egg)
cat(sprintf("    with the additional MR-Egger intercept rule (P >= 0.05): %d pairs, %d proteins (lost: %s)\n",
            nrow(egger_ok), length(prot_egg), paste(sort(lost_egg), collapse = ", ")))

# Steiger test (P < 0.05: the protein-to-outcome direction is supported)
st_yes <- sigp[!is.na(steiger_p) & steiger_p < 0.05, .N]
st_no  <- sigp[!is.na(steiger_p) & steiger_p >= 0.05]
st_na  <- sigp[is.na(steiger_p), .N]
cat(sprintf("    Steiger among the %d pairs: P < 0.05 in %d (%.1f%%), P >= 0.05 in %d [%s], not computed in %d\n", n_pairs, st_yes,
            100 * st_yes / n_pairs, nrow(st_no), paste(sprintf("%s %s x %s, P = %.3f", st_no$source, st_no$protein, st_no$outcome, st_no$steiger_p), collapse = "; "), st_na))

prot_tab <- sigp[, .(n_pairs = .N, sources = paste(sort(unique(source)), collapse = "+"),
                     outcomes_steatotic = sum(group == "steatotic"), outcomes_fibrotic = sum(group == "fibrotic"),
                     direction = fcase(uniqueN(dir_sign) == 1 & dir_sign[1] > 0, "up (all pairs)",
                                       uniqueN(dir_sign) == 1 & dir_sign[1] < 0, "down (all pairs)", default = "mixed")),
                 by = protein][order(-n_pairs, protein)]
fwrite(sigp[, .(source, protein, outcome, group, beta_ci, fdr = fdr_txt, pleiotropy_P = egger_p, steiger_P = steiger_p)][order(source, protein, outcome)],
       file.path(out_dir, "1_ivw_significant_pairs.csv"))
fwrite(prot_tab, file.path(out_dir, "2_protein_summary.csv"))

## ---- [2] the five genes whose only support was the weighted median -------------------------------------
five <- long[protein %in% FIVE, .(IVW_sig_pairs = sum(method == "IVW" & sig), WM_sig_pairs = sum(method == "WM" & sig)), by = protein][order(protein)]
fwrite(five, file.path(out_dir, "4_WM_only_genes.csv"))
cat("\n[2] significant pairs of the five genes of the original submission (support by weighted median only):\n"); print(five)

## ---- [3] Table S14: colocalization and MR-RAPS ---------------------------------------------------------
M14 <- read_raw("Table S14")
h14 <- which(trimws(M14[[1]]) == "Outcome trait" & trimws(M14[[2]]) == "pQTL source")[1]
nm14 <- trimws(unlist(M14[h14, ], use.names = FALSE))
S14 <- M14[(h14 + 1):nrow(M14), , drop = FALSE]; S14 <- S14[trimws(S14[[2]]) %in% c("deCODE", "ARIC"), , drop = FALSE]
col14 <- function(name) which(nm14 == name)[1]
d14 <- data.table(outcome = trimws(S14[[1]]), source = trimws(S14[[2]]), protein = trimws(S14[[3]]),
                  coloc = trimws(S14[[col14("Coloc status")]]), n_snps = suppressWarnings(as.numeric(S14[[col14("N snps")]])),
                  raps = trimws(S14[[col14("RAPS status")]]))
cat(sprintf("\n[3] Table S14: %d pairs; same pair set as the IVW pairs of Table S7: %s\n", nrow(d14),
            identical(sort(paste(d14$source, d14$protein, d14$outcome)), sort(paste(sigp$source, sigp$protein, sigp$outcome)))))
cat("    colocalization status:", paste(names(table(d14$coloc)), table(d14$coloc), collapse = "; "), "\n")
cat("    MR-RAPS status:      ", paste(names(table(d14$raps)), table(d14$raps), collapse = "; "), "\n")
raps_cons <- d14[, sum(raps == "Consistent")]; raps_est <- d14[, sum(raps != "NA (nIV < 3)")]
cat(sprintf("    MR-RAPS Consistent: %d of %d pairs (%.1f%%), %d of %d estimable pairs (%.1f%%); pairs with < 3 variants: %d\n",
            raps_cons, nrow(d14), 100 * raps_cons / nrow(d14), raps_cons, raps_est, 100 * raps_cons / raps_est, d14[, sum(n_snps < 3, na.rm = TRUE)]))

## ---- [4] Table S16: instrument strength ----------------------------------------------------------------
M16 <- read_raw("Table S16")
h16 <- which(trimws(M16[[1]]) == "pQTL source" & trimws(M16[[2]]) == "Protein" & trimws(M16[[3]]) == "Outcome")[1]
nm16 <- trimws(unlist(M16[h16, ], use.names = FALSE))
S16 <- M16[(h16 + 1):nrow(M16), , drop = FALSE]; S16 <- S16[trimws(S16[[1]]) %in% c("deCODE", "ARIC"), , drop = FALSE]
d16 <- data.table(source = trimws(S16[[1]]), protein = trimws(S16[[2]]), outcome = trimws(S16[[3]]), variant = trimws(S16[[4]]),
                  F = as.numeric(S16[[which(nm16 == "F-statistic")[1]]]), n_final = as.numeric(S16[[which(nm16 == "N instruments (final)")[1]]]))
# One instrument (source x protein x variant) is counted once. A few instruments have different F values in the rows of different
# outcomes; their F is the median of those rows (the summary does not depend on this choice beyond the last printed decimal).
n_varF <- d16[, .(k = uniqueN(F)), by = .(source, protein, variant)][k > 1, .N]
uniq <- d16[, .(F = median(F)), by = .(source, protein, variant)]
per_pair <- unique(d16[, .(source, protein, outcome, n_final)], by = c("source", "protein", "outcome"))
cat(sprintf("\n[4] Table S16: %d rows, %d unique instruments (source x protein x variant; %d of them with different F values in different outcome rows), %d pairs\n",
            nrow(d16), nrow(uniq), n_varF, nrow(per_pair)))
fq <- function(x) { q <- round1(quant(x)); sprintf("min %.1f, Q1 %.1f, median %.1f, Q3 %.1f, max %.1f", q["min"], q["q1"], q["median"], q["q3"], q["max"]) }
cat("    all:   ", fq(uniq$F), "; N with F < 10:", sum(uniq$F < 10), "\n")
for (s in c("deCODE", "ARIC")) cat(sprintf("    %-6s %s (%d instruments)\n", s, fq(uniq[source == s, F]), uniq[source == s, .N]))
cat(sprintf("    pairs with fewer than 3 instruments in S16: %d\n", per_pair[n_final < 3, .N]))

fs <- rbindlist(list(
  data.table(Group = "All",    `N proteins` = uniqueN(sigp[, .(source, protein)]), `N instruments` = nrow(uniq), as.data.table(as.list(round1(quant(uniq$F))))),
  data.table(Group = "deCODE", `N proteins` = uniqueN(sigp[source == "deCODE", protein]), `N instruments` = uniq[source == "deCODE", .N], as.data.table(as.list(round1(quant(uniq[source == "deCODE", F]))))),
  data.table(Group = "ARIC",   `N proteins` = uniqueN(sigp[source == "ARIC", protein]),   `N instruments` = uniq[source == "ARIC", .N],   as.data.table(as.list(round1(quant(uniq[source == "ARIC", F])))))))
fs[, `N with F < 10` := c(sum(uniq$F < 10), uniq[source == "deCODE", sum(F < 10)], uniq[source == "ARIC", sum(F < 10)])]
fwrite(fs, file.path(out_dir, "5_F_statistics_recomputed.csv"))

## ---- [5] Table S15 against the recomputed values -------------------------------------------------------
M15 <- read_raw("Table S15")
h15 <- which(trimws(M15[[1]]) == "Set")[1]
nm15 <- trimws(unlist(M15[h15, ], use.names = FALSE))
S15 <- M15[(h15 + 1):nrow(M15), , drop = FALSE]
s15 <- function(group, col) as.numeric(S15[[which(nm15 == col)[1]]][trimws(S15[[2]]) == group][1])
s15_diff <- 0L
s15_cols <- list(c("N proteins", "N proteins"), c("N instruments", "N instruments"), c("F min", "min"), c("F Q1", "q1"),
                 c("F median", "median"), c("F Q3", "q3"), c("F max", "max"))
for (g in c("All", "deCODE", "ARIC")) for (pr in s15_cols) {
  a <- s15(g, pr[1]); b <- fs[Group == g][[pr[2]]]
  if (is.na(a) || abs(a - b) > 0.051) {                       # the table prints one decimal (counts are integers)
    s15_diff <- s15_diff + 1L
    cat(sprintf("    S15 differs: %s / %s: table %s, recomputed %s\n", g, pr[1], a, signif(b, 6)))
  }
}
cat(sprintf("[5] Table S15 (rows All, deCODE, ARIC; 7 columns each) differs from the recomputed values in %d of 21 cells (tolerance 0.05, the table prints one decimal)\n", s15_diff))

## ---- count summary: recomputed value next to the value printed in the manuscript ----------------------
chk <- function(what, got, ms) data.table(quantity = what, recomputed = as.character(got), manuscript = as.character(ms),
                                          agrees = if (is.na(ms)) NA else isTRUE(all.equal(as.numeric(got), as.numeric(ms), tolerance = 1e-6)))
summ <- rbindlist(list(
  chk("proteins tested (IVW rows)", tested, 262),
  chk("IVW pairs with FDR < 0.05", n_pairs, 669),
  chk("  deCODE pairs", per_src[source == "deCODE", pairs], 411),
  chk("  ARIC pairs", per_src[source == "ARIC", pairs], 258),
  chk("proteins with genetic support (IVW)", length(prot_all), 224),
  chk("  deCODE proteins", per_src[source == "deCODE", proteins], 172),
  chk("  ARIC proteins", per_src[source == "ARIC", proteins], 100),
  chk("proteins with support for steatotic outcomes", length(S_set), 139),
  chk("proteins with support for fibrotic outcomes", length(F_set), 128),
  chk("proteins with support for both", length(intersect(S_set, F_set)), 43),
  chk("steatotic: one direction only", same_dir("steatotic", S_set), 136),
  chk("fibrotic: one direction only", same_dir("fibrotic", F_set), 122),
  chk("proteins left when the MR-Egger rule is added", length(prot_egg), 222),
  chk("Steiger P < 0.05", st_yes, 654),
  chk("Steiger P >= 0.05", nrow(st_no), 1),
  chk("Steiger not computed", st_na, 14),
  chk("MR-RAPS Consistent", raps_cons, 643),
  chk("MR-RAPS estimable pairs", raps_est, 644),
  chk("pairs with < 3 MR-RAPS variants (S14)", d14[, sum(n_snps < 3, na.rm = TRUE)], 25),
  chk("pairs with < 3 instruments (S16)", per_pair[n_final < 3, .N], 24),
  chk("colocalization strong (PP.H4 >= 0.8)", d14[, sum(coloc == "H4 strong")], 51),
  chk("colocalization suggestive", d14[, sum(coloc == "H4 suggestive")], 20),
  chk("unique instruments", nrow(uniq), 9356),
  chk("F median, all instruments", round1(median(uniq$F)), 72.3),
  chk("F Q1, all instruments", round1(unname(quantile(uniq$F, .25))), 44.2),
  chk("F Q3, all instruments", round1(unname(quantile(uniq$F, .75))), 155.6),
  chk("F minimum", round1(min(uniq$F)), 29.7),
  chk("instruments with F < 10", sum(uniq$F < 10), 0),
  chk("F median, deCODE", round1(median(uniq[source == "deCODE", F])), 74.5),
  chk("F median, ARIC", round1(median(uniq[source == "ARIC", F])), 61.8),
  chk("IVW pairs of the five WM-only genes", five[, sum(IVW_sig_pairs)], 0),
  chk("Table S15 cells that differ from the recomputed values", s15_diff, 0)))
fwrite(summ, file.path(out_dir, "3_count_summary.csv"))
options(width = 200)
cat("\n== recomputed value versus the value printed in the manuscript ==\n"); print(summ, row.names = FALSE)
if (any(!summ$agrees, na.rm = TRUE)) warning("some recomputed values differ from the manuscript; see 3_count_summary.csv") else cat("\nall values agree with the manuscript\n")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\noutput:", normalizePath(out_dir), "\n")
