## =====================================================================
## Supplementary Figure S1 (drawn from the Supplementary Tables). Competitive gene-set (MAGMA) results for the 12 tests
##   Thyroid-derived gene set (1,473 genes) x 8 liver datasets; liver-derived gene sets (958-965 genes) x 4 thyroid datasets
##
## Input (Supplementary Tables workbook): Table S1 (MAGMA competitive gene-set output)
## Thresholds: Bonferroni P = 0.05/12 = 0.0042, nominal P = 0.05
## Output: figures_out_supp_S1/ (FigureS1 as .pdf and .png, FigureS1_legend.txt, sessionInfo.txt)
## Usage : Rscript FigS1_MAGMA_gene_set_supp.R
## =====================================================================
## ---------------------------------------------------------------------------------------------
## [0] Setup
##   * Standalone script: it does not source any other file.
##   * Main input: the Supplementary Tables workbook (file name containing "Supplementary_Tables_S1-S13", .xlsx).
##     It is searched, recursively, in the folder given by the DATA_DIR environment variable (if set), in ./data
##     and in the current working directory; or set SUPP_XLSX to the full path of the workbook.
##     Two layouts are accepted: one sheet per table ("Table S2a" or "S2a"), or one sheet per table number with the
##     lettered tables side by side ("Table S2" holding S2a and S2b), as in the submitted workbook.
##   * Public-data inputs are used only where stated in the header; they are searched by name in the same folders.
##   * Input files are only read, never modified. Figures are written to OUT_DIR in the current working directory.
##   * R packages: ggplot2, data.table, readxl, patchwork, ggrepel, scales, ragg
## ---------------------------------------------------------------------------------------------
invisible(try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE))   # needed for non-ASCII file names
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(readxl); library(patchwork); library(ggrepel); library(scales) })
ROOTS <- unique(c(Sys.getenv("DATA_DIR", unset = ""), "data", getwd()))
ROOTS <- ROOTS[nzchar(ROOTS) & dir.exists(ROOTS)]
if (!length(ROOTS)) stop("No input folder found. Create ./data or set the DATA_DIR environment variable.")
find_file <- function(stem, ext) {                       # a copy suffix after the stem (_1 or _<6 hex characters>) is treated as the same file
  pat  <- paste0("^", stem, "(_[0-9a-f]{6}|_[0-9]+)?\\.", ext, "$")
  hits <- unlist(lapply(ROOTS, function(r) list.files(r, pattern = pat, recursive = TRUE, full.names = TRUE)))
  if (!length(hits)) stop(sprintf("File not found: %s.%s\n  Folders searched: %s", gsub("\\(20\\)\\?", "", stem), ext, paste(ROOTS, collapse = " | ")), call. = FALSE)
  hit <- hits[which.max(file.mtime(hits))]
  message(sprintf("  Using file: %s  (modified %s)", hit, format(file.mtime(hit), "%Y-%m-%d %H:%M")))
  hit
}
find_quiet <- function(stem, ext) tryCatch(suppressMessages(find_file(stem, ext)), error = function(e) NA_character_)   # NA when the file is absent
SUPP <- local({
  e <- Sys.getenv("SUPP_XLSX", unset = "")
  if (nzchar(e)) return(e)
  hits <- unlist(lapply(ROOTS, function(r) list.files(r, pattern = "Supplementary_Tables_S1-S13.*\\.xlsx$", recursive = TRUE, full.names = TRUE)))
  hits <- hits[!grepl("^~\\$", basename(hits))]
  if (!length(hits)) stop("Supplementary Tables workbook not found (file name must contain 'Supplementary_Tables_S1-S13'). Folders searched: ", paste(ROOTS, collapse = " | "), call. = FALSE)
  hits[which.max(file.mtime(hits))]
})
message("  Supplementary Tables: ", SUPP)
SUPP_SHEETS <- excel_sheets(SUPP)
## Locates one table of the workbook. Two layouts are accepted: one sheet per table ("Table S2a" or "S2a"), or one sheet per table
## number ("Table S2") with the lettered tables side by side, each starting with a title cell "Table S2a. ..." in row 1.
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
has_supp <- function(key) !is.null(locate_supp(key))
## Reads one table as a data.table. The header row is the first row of the table whose first cell equals `first_col`
## (the rows above it hold the table title and notes). Empty spacer columns and empty rows are dropped.
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
OUT_DIR <- "figures_out_supp_S1"; dir.create(OUT_DIR, showWarnings = FALSE)
P <- function(key) find_file(FILES[[key]][1], FILES[[key]][2])

## [1] All inputs come from the Supplementary Tables workbook (sheets named in the header)
FILES <- list()

## [2] Common settings (colours, fonts, save function)
## ---- Colours and theme (palette distinguishable with colour-vision deficiency) ---------------------------
COL <- list(up = "#C0392B", down = "#2E6DA4", ns = "grey75", high = "#B2182B", moderate = "#EF8A62",
            below = "grey70", thy = "#E08A00", liv = "#3B7DA6", cand = "#C0392B", shared = "#1F4E79",
            wt = "#2A9D8F", ko = "#E9C46A", int = "#6A4C93", green = "#2E8B57")
theme_fig <- function(base = 7) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          plot.tag = element_text(face = "bold", size = base + 3),
          legend.key.size = unit(3, "mm"), legend.title = element_text(size = base), legend.text = element_text(size = base - 0.5),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = base),
          plot.title = element_text(size = base + 0.5, face = "bold"))
}
theme_set(theme_fig())
MM <- 1 / 25.4
save_fig <- function(p, name, w_mm = 180, h_mm = 120) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w_mm * MM, height = h_mm * MM, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w_mm * MM, height = h_mm * MM, dpi = 300, device = ragg::agg_png)
  message("Saved: ", file.path(OUT_DIR, paste0(name, ".pdf / .png")))
  invisible(p)
}
sig_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))
nf <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
## Writes the figure legend to OUT_DIR/<name>_legend.txt (the `notes` argument is accepted but not written)
write_legend <- function(name, legend, notes = character()) {
  f <- file.path(OUT_DIR, paste0(name, "_legend.txt"))
  con <- file(f, open = "w", encoding = "UTF-8"); writeLines(legend, con); close(con)
  message("Saved: ", f); invisible(f)
}


## ---------------------------------------------------------------------------------------------
## [5] Drawing the figures
## ---------------------------------------------------------------------------------------------
s1 <- read_supp("S1", "gene_set")   # Supplementary Table S1
TSH_NAME <- c("TSH dataset B (GCST90296333)" = "TSH_study1", "TSH dataset A (GCST010653)" = "TSH_study2")   # names used in the MAGMA input files (note of Table S1)
s1[, ds := ifelse(tested_in %in% names(TSH_NAME), unname(TSH_NAME[tested_in]), tested_in)]
a <- s1[grepl("^Thyroid-derived", gene_set), .(liver_dataset = ds, NGENES = n_genes_in_set, BETA = beta, SE, P)]
b <- s1[grepl("^Liver-derived", gene_set), .(thyroid_dataset = ds, NGENES = n_genes_in_set, BETA = beta, SE, P)]
A <- a[, .(panel = "Thyroid-derived set (1,473 genes) tested in liver datasets", dataset = liver_dataset, NGENES, BETA, SE, P)]
B <- b[, .(panel = "Liver-derived set (958\u2013965 genes) tested in thyroid datasets", dataset = thyroid_dataset, NGENES, BETA, SE, P)]
d <- rbind(A, B); BONF <- 0.05 / nrow(d); stopifnot(nrow(d) == 12)
d[, cls := factor(fifelse(P < BONF, "Bonferroni significant", fifelse(P < 0.05, "nominal (P < 0.05)", "not significant")), levels = c("Bonferroni significant", "nominal (P < 0.05)", "not significant"))]
message(sprintf("Bonferroni threshold %.4f | significant %d: %s", BONF, sum(d$P < BONF), paste(d[P < BONF]$dataset, collapse = ", ")))
## Dataset names: the names of the MAGMA input (TSH_study1, TSH_study2) are changed to the names of manuscript Table 1. The name in square brackets is the original name.
## To be verified: TSH_study1 = GCST90296333 (TSH dataset B), TSH_study2 = GCST010653 (TSH dataset A)
##   -> Check against the GWAS input definition of the script that ran MAGMA (260916_2_magma_gene_set_enrichment.R). If it differs, only LBL below needs to be changed.
d[, panel := factor(panel, levels = unique(panel))]; PL <- levels(d$panel)
LBL <- c(TSH_study1 = "TSH dataset B (GCST90296333) [study1]", TSH_study2 = "TSH dataset A (GCST010653) [study2]", FT3 = "FT3", FT4 = "FT4")
nice <- function(x) { y <- unname(LBL[x]); ifelse(is.na(y), sub("_unadjusted", ", unadjusted", sub("_(BMI|WFM|VAT)adjusted", ", \\1-adjusted", x)), y) }
d[, label := nice(as.character(dataset))]
xr <- c(min(d$BETA - 1.96 * d$SE), max(d$BETA + 1.96 * d$SE) + 0.09)
mk <- function(dd, ttl, last = FALSE) {
  dd <- copy(dd); dd[, label := factor(label, levels = rev(unique(label)))]
  ggplot(dd, aes(BETA, label, colour = cls)) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_pointrange(aes(xmin = BETA - 1.96 * SE, xmax = BETA + 1.96 * SE), size = 0.4, linewidth = 0.45) +
    geom_text(aes(x = BETA + 1.96 * SE + 0.004, label = paste0("P = ", fmt_p(P))), hjust = 0, size = 2, show.legend = FALSE) +
    scale_colour_manual(values = c("Bonferroni significant" = COL$high, "nominal (P < 0.05)" = COL$moderate, "not significant" = "grey55"), name = NULL, drop = FALSE) +
    coord_cartesian(xlim = xr) +
    labs(title = ttl, x = if (last) "Competitive gene-set effect (\u03b2 \u00b1 95% CI)" else NULL, y = NULL,
         caption = if (last) sprintf("Bonferroni threshold P = %.4f (12 tests)", BONF) else NULL) +
    theme(plot.title = element_text(size = 7, face = "bold"), plot.caption = element_text(size = 5.5, hjust = 0))
}
p <- (mk(d[panel == PL[1]], "Thyroid-derived gene set (1,473 genes) tested in the 8 liver datasets") /
        mk(d[panel == PL[2]], "Liver-derived gene set (958\u2013965 genes) tested in the 4 thyroid datasets", last = TRUE)) +
  plot_layout(heights = c(8, 4), guides = "collect") & theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))
save_fig(p, "FigS1_MAGMA_gene_set", 180, 100)

## ---- Legend txt (the legend of this figure; the numbers are the values computed above) -----------------------
sig_names <- nice(as.character(d[P < BONF]$dataset))
leg <- sprintf("Supplementary Figure S1. Competitive gene-set analysis with MAGMA. The thyroid-derived gene set (1,473 genes) was tested in the 8 liver datasets and the liver-derived gene set (958 to 965 genes) in the 4 thyroid datasets. Points show the competitive gene-set effect (\u03b2) with its 95%% confidence interval; red marks Bonferroni significance (P < %.4f, 12 tests; %d tests: %s) and orange nominal significance (P < 0.05).", BONF, sum(d$P < BONF), paste(sig_names, collapse = "; "))
write_legend("FigureS1", leg)


## ---- Record R and package versions used to draw the figures --------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
