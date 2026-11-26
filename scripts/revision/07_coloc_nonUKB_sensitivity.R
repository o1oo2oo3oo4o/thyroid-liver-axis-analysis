#!/usr/bin/env Rscript
# =============================================================================
# 07_coloc_nonUKB_sensitivity.R        (revision analysis; Table S2)
#
# Reviewer 2 major comment 1: TSH dataset B contains UK Biobank participants, who also
# contribute to the liver MRI GWAS. How robust is the colocalization result when TSH
# dataset B is excluded?
#
# Rule of the manuscript (Methods 2.4), reproduced here from the pair-level table of Table S2:
#   * a liver phenotype (PDFF or cT1) supports a window if ANY of its 4 adjustment models does
#   * TSH must reach the threshold in BOTH TSH datasets, so the LOWER value decides
#   * a window is supported when the highest cross-domain PP.H4 is >= 0.5 (high: >= 0.8)
# Sensitivity rule: TSH dataset B removed, only TSH dataset A, FT4 and FT3 are used
# (TSH A is nested within TSH B, so this is a non-UKB thyroid-input sensitivity analysis,
#  not an independent replication).
#
# Input : Table S2 workbook (sheet "Table S2"; window summary on the left, pair table on the right)
# Output: results/coloc_nonUKB/*.csv
# Expected with the submitted Table S2: original rule 27 windows (11 high, 16 moderate);
#   without TSH-B 28 windows (13 high, 15 moderate); L00359 newly supported, L00414 moderate -> high.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(readxl) })
source(Sys.getenv("STAGE_CFG", "config.R"))
out_dir <- file.path(out_root, "07_coloc_nonUKB"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- as.data.table(read_excel(s2_file, sheet = "Table S2", col_names = FALSE, .name_repair = "minimal"))
hdr_row <- which(sapply(seq_len(min(15, nrow(raw))), function(i) "dataset_1" %in% as.character(unlist(raw[i]))))[1]
hdr <- as.character(unlist(raw[hdr_row]))
c_ds1 <- which(hdr == "dataset_1")
pair <- raw[(hdr_row + 1):nrow(raw), (c_ds1 - 1):(c_ds1 + 7), with = FALSE]
setnames(pair, c("window", "ds1", "ds2", "H0", "H1", "H2", "H3", "H4", "pair_type"))
pair <- pair[grepl("^L0", window)]; pair[, H4 := as.numeric(H4)]
# window summary (reported tiers) on the left of the same sheet
summ_hdr <- which(sapply(seq_len(min(15, nrow(raw))), function(i) identical(as.character(raw[[1]][i]), "window")))[1]
win <- raw[(summ_hdr + 1):nrow(raw), 1:9, with = FALSE]
setnames(win, c("window", "chr", "start", "end", "lead_variant", "lead_P", "tier", "maxPP", "supporting_pairs"))
win <- win[grepl("^L0", window)]; win[, maxPP := as.numeric(maxPP)]

cd <- pair[pair_type == "cross-domain"]
thy <- c("FT3", "FT4", "TSH_A", "TSH_B")
cd[, `:=`(thy = fifelse(ds1 %in% thy, ds1, ds2), liv = fifelse(ds1 %in% thy, ds2, ds1))]
wide <- dcast(cd, window + liv ~ thy, value.var = "H4", fun.aggregate = max)
for (cc in thy) if (!cc %in% names(wide)) wide[, (cc) := NA_real_]
best_pp <- function(wide, tsh_rule) {
  w <- copy(wide)
  w[, TSH := if (tsh_rule == "both") pmin(TSH_A, TSH_B) else TSH_A]       # NA when a dataset is missing
  w[, best := pmax(FT3, FT4, TSH, na.rm = TRUE)]
  w[is.infinite(best), best := NA_real_]
  w[, .(maxPP = suppressWarnings(max(best, na.rm = TRUE))), by = window][is.infinite(maxPP), maxPP := NA_real_]
}
orig <- best_pp(wide, "both"); sens <- best_pp(wide, "A")
res <- merge(win[, .(window, chr, start, end, lead_variant, reported_tier = tier, reported_maxPP = maxPP)],
             orig[, .(window, maxPP_original_rule = maxPP)], by = "window", all.x = TRUE)
res <- merge(res, sens[, .(window, maxPP_without_TSH_B = maxPP)], by = "window", all.x = TRUE)
tier_of <- function(x) fifelse(is.na(x), "no", fifelse(x >= 0.8, "high", fifelse(x >= 0.5, "moderate", "no")))
res[, `:=`(tier_original = tier_of(maxPP_original_rule), tier_without_TSH_B = tier_of(maxPP_without_TSH_B))]
res[, change := fifelse(tier_original == tier_without_TSH_B, "", paste(tier_original, "->", tier_without_TSH_B))]
fwrite(res, file.path(out_dir, "windows_original_vs_without_TSH_B.csv"))

cnt <- function(x) c(supported = sum(x != "no"), high = sum(x == "high"), moderate = sum(x == "moderate"))
cat("[reproduction] reported tiers in Table S2 :", cnt(res$reported_tier), "\n")
cat("[reproduction] original rule, recomputed  :", cnt(res$tier_original),
    "| max |difference| of maxPP vs Table S2 =", signif(max(abs(res$maxPP_original_rule - res$reported_maxPP), na.rm = TRUE), 2), "\n")
cat("[sensitivity ] without TSH dataset B      :", cnt(res$tier_without_TSH_B), "\n\n")
cat("windows whose tier changes (add the genes of newly supported windows to the response):\n")
print(res[change != "", .(window, chr, start, end, lead_variant, maxPP_original_rule, maxPP_without_TSH_B, change)])
cat("\nlost windows (should be none):", res[tier_original != "no" & tier_without_TSH_B == "no", .N], "\n")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
