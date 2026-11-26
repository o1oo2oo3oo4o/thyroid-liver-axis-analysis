## =====================================================================
## Supplementary Figure S3 (drawn from the Supplementary Tables). Regional layouts of all 27 supported windows
##   Windows (coloured by PP.H4), lead variants (diamonds), coding genes
##   (shared-domain genes in blue, T3-responsive candidates in red)
##
## Inputs: Supplementary Tables workbook: S2a, S2b (colocalization), S4 (gene-to-window assignment)
##   Optional (regional association plots only; if absent a window/gene layout is drawn instead): public GWAS summary statistics listed in [1-b]
##   -- GWAS Catalog GCST90446646-GCST90446653, GCST010653, GCST90296333, and the ThyroidOmics Consortium
##   files FT3_summary_data.txt and FT4_summary_data.txt. They are not part of the Supplementary Tables and are not redistributed.
## Output: figures_out_supp_S3/ (FigureS3 as .pdf and .png, FigureS3_legend.txt, sessionInfo.txt)
## Usage : Rscript FigS3_window_layouts_supp.R
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
##   * R packages: ggplot2, data.table, readxl, patchwork, ggrepel, scales, ragg, cowplot
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
OUT_DIR <- "figures_out_supp_S3"; dir.create(OUT_DIR, showWarnings = FALSE)
P <- function(key) find_file(FILES[[key]][1], FILES[[key]][2])

## [1] All inputs come from the Supplementary Tables workbook (sheets named in the header)
FILES <- list()

## [1-b] GWAS summary statistics used for the regional association plots (Supplementary Figure S3) (optional: if the files are not found, a window/gene layout is drawn instead)
##   File names and column names are those of the public GWAS files used by the colocalization scan (step 1 in README.md).
##   Liver MRI GWAS : GCST90446646-49 (cT1), GCST90446650-53 (PDFF)  .tsv
##   Thyroid        : FT3_summary_data.txt, FT4_summary_data.txt, GCST010653.txt (TSH A), GCST90296333.tsv (TSH B)
##   each entry = c(file name, extension, chromosome column, position column, P-value column)
GW <- list(
  cT1_646  = c("GCST90446646", "tsv", "chromosome", "base_pair_location", "p_value"),
  cT1_647  = c("GCST90446647", "tsv", "chromosome", "base_pair_location", "p_value"),
  cT1_648  = c("GCST90446648", "tsv", "chromosome", "base_pair_location", "p_value"),
  cT1_649  = c("GCST90446649", "tsv", "chromosome", "base_pair_location", "p_value"),
  PDFF_650 = c("GCST90446650", "tsv", "chromosome", "base_pair_location", "p_value"),
  PDFF_651 = c("GCST90446651", "tsv", "chromosome", "base_pair_location", "p_value"),
  PDFF_652 = c("GCST90446652", "tsv", "chromosome", "base_pair_location", "p_value"),
  PDFF_653 = c("GCST90446653", "tsv", "chromosome", "base_pair_location", "p_value"),
  FT3      = c("FT3_summary_data", "txt", "chr", "position", "P.value"),
  FT4      = c("FT4_summary_data", "txt", "chr", "position", "P.value"),
  TSH_A    = c("GCST010653", "txt", "Chr", "Pos", "P-value"),
  TSH_B    = c("GCST90296333", "tsv", "chromosome", "base_pair_location", "p_value")
)
GW_LABEL <- c(cT1_646 = "cT1 (GCST90446646)", cT1_647 = "cT1 (GCST90446647)", cT1_648 = "cT1 (GCST90446648)", cT1_649 = "cT1 (GCST90446649)",
              PDFF_650 = "PDFF (GCST90446650)", PDFF_651 = "PDFF (GCST90446651)", PDFF_652 = "PDFF (GCST90446652)", PDFF_653 = "PDFF (GCST90446653)",
              FT3 = "FT3", FT4 = "FT4", TSH_A = "TSH dataset A (GCST010653)", TSH_B = "TSH dataset B (GCST90296333)")
find_quiet <- function(stem, ext) tryCatch(suppressMessages(find_file(stem, ext)), error = function(e) NA_character_)

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


## [3] Candidate gene lists (as in manuscript Tables 3 and 4)
## ---- Candidate gene lists (as in manuscript Tables 3 and 4) ------------------------------
POS10   <- c("CUX2","BDH2","ADH4","LILRB3","SBK1","ABHD1","INCENP","CDC42EP5","CD19","GPCPD1")
PROT16  <- c("ADH4","BTC","CFHR1","CLIC5","COL11A2","COL15A1","DEF6","DNAJB4","F13B","LGALS4","PKLR","POR","SPARCL1","TIMD4","GSTM1","HSPA1A")
PERMISSIVE <- c("CFHR1","GSTM1","HSPA1A")                      # candidates included through a non-one-to-one mapping
MOUSE_OF   <- c(GSTM1 = "Gstm6", HSPA1A = "Hspa1b", F13B = "F13b")   # all others: only the first letter in upper case
mouse_symbol <- function(g) ifelse(g %in% names(MOUSE_OF), MOUSE_OF[g], paste0(substr(g,1,1), tolower(substr(g,2,nchar(g)))))
PANEL26 <- sort(unique(c(POS10, PROT16, "SLC39A8")))


## [4] COLOC classification and gene layout functions
## ---- COLOC ---------------------------------------------------------------
HG19_LEN <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,
              133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560)
chr_num <- function(x) { x <- sub("^chr","",as.character(x)); x[x=="X"] <- "23"; as.integer(x) }
load_coloc <- function() {
  a <- read_supp("S2a", "window")   # window summary: window, chr, start_GRCh37, end_GRCh37, lead_variant, lead_P, tier, ...
  b <- read_supp("S2b", "window")   # all 29,802 dataset pairs: window, dataset_1, dataset_2, PP.H0-PP.H4, pair_type
  S1 <- a[, .(Locus_ID = window, Chr = chr, Start_bp = start_GRCh37, End_bp = end_GRCh37, Lead_SNP = lead_variant, Lead_pval = lead_P, Final_tier = tier)]
  S2 <- b[, .(locus_id = window, trait1 = dataset_1, trait2 = dataset_2, PP.H3, PP.H4)]
  S1[, chr := chr_num(Chr)]
  S1[, lead_pos := as.numeric(sub("^.*:", "", sub("_.*$", "", Lead_SNP)))]
  list(S1 = S1, S2 = S2)
}
## Cross-domain classification rule: liver traits take the maximum over the 4 adjustment models, TSH takes the lower value of the two datasets (both must be available),
## and the window value is the maximum over the 6 trait-level pairs (PDFF/cT1 x TSH/FT3/FT4)
window_best <- function(S2, by = "locus_id") {
  d <- as.data.table(S2)
  isL <- function(x) grepl("^(PDFF|cT1)", x); isT <- function(x) x %in% c("TSH_A","TSH_B","FT3","FT4")
  o1 <- d[isL(trait1) & isT(trait2)][, `:=`(liver = trait1, thy = trait2)]
  o2 <- d[isL(trait2) & isT(trait1)][, `:=`(liver = trait2, thy = trait1)]
  cd <- rbind(o1, o2, fill = TRUE)
  cd[, `:=`(h3 = PP.H3, h4 = PP.H4, lt = sub("_.*", "", liver), tt = ifelse(thy %in% c("TSH_A","TSH_B"), "TSH", thy))]
  keep <- c(by, "liver", "lt", "tt", "h3", "h4")
  tsh <- cd[tt == "TSH"]; tsh[, n := .N, by = c(by, "liver")]
  tsh <- tsh[n == 2][, .SD[which.min(h4)], by = c(by, "liver")]
  oth <- cd[tt != "TSH"]
  both <- rbind(oth[, ..keep], tsh[, ..keep])
  tl   <- both[, .SD[which.max(h4)], by = c(by, "lt", "tt")]
  best <- tl[, .SD[which.max(h4)], by = by]
  best[, pair := paste(lt, "\u00d7", tt)]
  list(trait_level = tl, best = best)
}

## ---- Genes and windows ---------------------------------------------------------
load_genes <- function() {
  s4 <- read_supp("S4", "window")   # gene-to-window assignment of the 27 supported windows
  G <- s4[, .(Locus_ID = window, Gene = gene, gene_start = gene_start_GRCh37, gene_end = gene_end_GRCh37, biotype, GWAS_shared = in_shared_199, locus_pass = TRUE)]
  G[, gene_cat := fcase(Gene %in% POS10, "T3-responsive positional candidate",
                        GWAS_shared == 1, "shared-domain gene", default = "other")]
  G[]
}
## Assign rows so that labels do not overlap
pack_rows <- function(start, end, pad = 0) {
  o <- order(start); row <- integer(length(start)); ends <- numeric(0)
  for (i in o) { r <- which(ends + pad < start[i])[1]; if (is.na(r)) { ends <- c(ends, end[i]); r <- length(ends) } else ends[r] <- end[i]; row[i] <- r }
  row
}
## Layout of one region (several windows): windows (PP.H4 colour), lead variants, genes
region_layout <- function(loci, S1, W, G, title = NULL, label_genes = TRUE, gene_text = 1.9, win_text = 1.9, always_label = character()) {
  w <- merge(S1[Locus_ID %in% loci, .(locus_id = Locus_ID, Start_bp, End_bp, lead_pos, Lead_pval)], W[, .(locus_id, h4, pair)], by = "locus_id", all.x = TRUE)
  w[, `:=`(xmin = Start_bp / 1e6, xmax = End_bp / 1e6, lead = lead_pos / 1e6)]
  w[, row := pack_rows(xmin, xmax, 0.02)]
  g <- unique(G[Locus_ID %in% loci & biotype == "coding", .(Gene, gene_start, gene_end, gene_cat)])
  g[, gmin := gene_start / 1e6]; g[, gmax := gene_end / 1e6]
  g[, grow := pack_rows(gmin, gmax, 0.05)]
  nwr <- max(w$row); base <- -0.4
  p <- ggplot() +
    geom_rect(data = w, aes(xmin = xmin, xmax = xmax, ymin = nwr - row + 0.1, ymax = nwr - row + 0.8, fill = h4), colour = "grey30", linewidth = 0.25) +
    geom_text(data = w, aes(x = (xmin + xmax)/2, y = nwr - row + 0.45, label = sprintf("%s\nPP.H4 %.2f", locus_id, h4)), size = win_text, lineheight = 0.9) +
    geom_point(data = w, aes(x = lead, y = nwr + 0.35), shape = 18, size = 2.4, colour = "black") +
    geom_vline(data = w, aes(xintercept = lead), linetype = "dotted", linewidth = 0.25, colour = "grey40") +
    geom_segment(data = g, aes(x = gmin, xend = gmax, y = base - grow * 0.5, yend = base - grow * 0.5, colour = gene_cat), linewidth = 1.4, lineend = "butt") +
    scale_fill_gradient(low = "#FEE8C8", high = "#B30000", limits = c(0.5, 1), oob = squish, name = "PP.H4") +
    scale_colour_manual(values = c("T3-responsive positional candidate" = COL$cand, "shared-domain gene" = COL$shared, "other" = "grey70"), limits = c("T3-responsive positional candidate", "shared-domain gene", "other"), drop = FALSE, name = "Gene") +
    scale_y_continuous(NULL, breaks = NULL, expand = expansion(mult = c(0.12, 0.05))) + labs(x = "Position (Mb, GRCh37)", title = title) +
    theme(axis.line.y = element_blank())
  if (label_genes) {
    gl <- g[gene_cat != "other" | Gene %in% always_label]
    p <- p + geom_text_repel(data = gl, aes(x = (gmin + gmax)/2, y = base - grow * 0.5, label = Gene, colour = gene_cat),
                             size = gene_text, fontface = "italic", segment.size = 0.15, min.segment.length = 0, box.padding = 0.3, point.padding = 0.1, max.overlaps = Inf, show.legend = FALSE, direction = "both", nudge_y = -0.35, force = 3)
  }
  p
}

## ---- Reading GWAS intervals ------------------------------------------------------
## An interval that has been read is saved in cache_gwas_*.rds in OUT_DIR and loaded directly in later runs (the original GWAS files are only read)
read_gwas <- function(key, wins) {          # wins: data.table(region, chr, start, end)
  spec <- GW[[key]]; cache <- file.path(OUT_DIR, paste0("cache_gwas_", key, ".rds"))
  sig <- paste(wins$region, wins$chr, wins$start, wins$end, collapse = "|")
  if (file.exists(cache)) { z <- readRDS(cache); if (identical(z$sig, sig)) { message("  using cache: ", key); return(z$dat) } }
  f <- find_file(spec[1], spec[2])
  d <- fread(f, select = spec[3:5], showProgress = FALSE); setnames(d, spec[3:5], c("chr", "pos", "p"))
  d[, chr := chr_num(chr)]; d[, p := pmax(as.numeric(p), 1e-300)]; d <- d[!is.na(p) & !is.na(pos)]
  out <- rbindlist(lapply(seq_len(nrow(wins)), function(i) d[chr == wins$chr[i] & pos >= wins$start[i] & pos <= wins$end[i]][, region := wins$region[i]]))
  saveRDS(list(sig = sig, dat = out), cache); out
}

## ---------------------------------------------------------------------------------------------
## [5] Drawing the figures
## ---------------------------------------------------------------------------------------------
co <- load_coloc(); S1 <- co$S1; S2 <- co$S2; G <- load_genes()
wb <- window_best(S2)$best
sup <- S1[Final_tier %in% c("high", "moderate")][order(chr, lead_pos)]
stopifnot(nrow(sup) == 27)
## For each window the pair (thyroid trait x liver trait) with the highest PP.H4 is chosen; -log10 P is drawn at the top for the thyroid trait and at the bottom (upside down) for the liver trait
## Diamond = lead variant, lines at the bottom = shared-domain genes (blue) and positional candidates responsive to T3 (red). For TSH pairs dataset B (GCST90296333) is drawn
sup <- merge(sup, wb[, .(Locus_ID = locus_id, h4, pair, tt, liver)], by = "Locus_ID", all.x = TRUE, sort = FALSE)[order(chr, lead_pos)]
sup[, thy := fifelse(tt == "TSH", "TSH_B", tt)]
need <- unique(c(sup$thy, sup$liver))
have <- sapply(need, function(k) !is.na(find_quiet(GW[[k]][1], GW[[k]][2])))
if (all(have)) {
  message("  traits used in the regional association plots: ", paste(need, collapse = ", "))
  wins <- sup[, .(region = Locus_ID, chr, start = Start_bp, end = End_bp, thy, liv = liver)]
  gdat <- setNames(lapply(need, function(k) read_gwas(k, wins[thy == k | liv == k, .(region, chr, start, end)])), need)
  CAP <- 40; GRP <- c("Thyroid trait (up)", "Liver trait (down)", "T3-responsive positional candidate", "shared-domain gene")
  mini <- function(i) {
    r <- sup[i]; l <- r$Locus_ID
    d <- rbind(gdat[[r$thy]][region == l][, .(x = pos / 1e6, y =  pmin(-log10(p), CAP), grp = GRP[1])],
               gdat[[r$liver]][region == l][, .(x = pos / 1e6, y = -pmin(-log10(p), CAP), grp = GRP[2])])
    g <- unique(G[Locus_ID == l & biotype == "coding" & gene_cat != "other", .(Gene, gene_start, gene_end, gene_cat)])
    g[, `:=`(grp = fifelse(gene_cat == "T3-responsive positional candidate", GRP[3], GRP[4]), gmin = gene_start / 1e6, gmax = gene_end / 1e6)]
    p <- ggplot() +
      geom_hline(yintercept = c(-1, 1) * -log10(1e-6), linetype = "dashed", linewidth = 0.2, colour = "grey50") + geom_hline(yintercept = 0, linewidth = 0.25) +
      geom_point(data = d, aes(x, y, colour = grp), size = 0.2, alpha = 0.7) +
      geom_point(aes(x = r$lead_pos / 1e6, y = 0), shape = 18, size = 1.8)
    if (nrow(g)) p <- p + geom_segment(data = g, aes(x = gmin, xend = gmax, y = -CAP * 1.18, yend = -CAP * 1.18, colour = grp), linewidth = 1.2) +
      geom_text_repel(data = g[grp == GRP[3]], aes(x = (gmin + gmax) / 2, y = -CAP * 1.18, label = Gene, colour = grp), size = 1.5, fontface = "italic", nudge_y = -CAP * 0.3,
                      segment.size = 0.15, min.segment.length = 0, show.legend = FALSE)
    p + scale_colour_manual(values = setNames(c(COL$thy, COL$liv, COL$cand, COL$shared), GRP), limits = GRP, drop = FALSE, name = NULL) +
      scale_y_continuous(labels = function(v) abs(v), limits = c(-CAP * 1.6, CAP * 1.05), breaks = c(-40, -20, 0, 20, 40)) +
      labs(x = NULL, y = NULL, title = sprintf("%s  chr%d:%.1f Mb  PP.H4 %.2f (%s)", l, r$chr, r$lead_pos / 1e6, r$h4, r$pair)) +
      theme(plot.title = element_text(size = 5.5), axis.text = element_text(size = 5), legend.position = "none")
  }
  plots <- lapply(seq_len(nrow(sup)), mini)
  leg <- cowplot::get_legend(plots[[1]] + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2, linewidth = 1.2))) + theme(legend.position = "bottom"))
  p <- wrap_plots(plots, ncol = 3) / wrap_elements(leg) + plot_layout(heights = c(40, 1)) +
    plot_annotation(caption = "Position on x axis (Mb, GRCh37), window = lead variant \u00b1 500 kb.  y axis: \u2212log10 P capped at 40; dashed lines P = 1e-6.", theme = theme(plot.caption = element_text(size = 6, hjust = 0)))
} else {
  message("  GWAS summary statistics not found; a window/gene layout is drawn instead. Files not found: ", paste(need[!have], collapse = ", "))
  plots <- lapply(sup$Locus_ID, function(l) {
    r <- sup[Locus_ID == l]
    region_layout(l, S1, wb, G, title = sprintf("%s  chr%d:%.1f Mb  (%s)", l, r$chr, r$lead_pos / 1e6, ifelse(r$Final_tier == "nominal", "moderate", r$Final_tier)), gene_text = 1.5, win_text = 1.7) +
      theme(legend.position = "none", plot.title = element_text(size = 6), axis.text = element_text(size = 5), axis.title.x = element_blank())
  })
  leg <- cowplot::get_legend(region_layout(sup$Locus_ID[1], S1, wb, G) + theme(legend.position = "bottom", legend.box = "horizontal"))
  p <- wrap_plots(plots, ncol = 3) / wrap_elements(leg) + plot_layout(heights = c(40, 1))
}
save_fig(p, "FigS3_window_layouts", 180, 265)

## ---- Legend txt (the legend of this figure) -----------------------------------------------
leg <- if (all(have)) sprintf("Supplementary Figure S3. Regional association plots of the %d supported windows (lead variant \u00b1 500 kb). For each window, the \u2212log10 P of the thyroid trait is plotted upward and that of the liver trait downward, for the pair with the highest PP.H4 (given with its PP.H4 in each panel title; TSH pairs show dataset B). Diamonds mark lead variants, dashed lines mark P = 1 \u00d7 10\u207b\u2076 and the lower bars mark coding genes (shared-domain genes in blue; T3-responsive positional candidates in red and labelled). \u2212log10 P is capped at 40.", nrow(sup)) else
  sprintf("Supplementary Figure S3. Genomic layout of the %d supported windows (lead variant \u00b1 500 kb), showing the colocalization window coloured by PP.H4, the lead variant (diamond) and coding genes (shared-domain genes in blue; T3-responsive positional candidates in red).", nrow(sup))
write_legend("FigureS3", leg)


## ---- Record R and package versions used to draw the figures --------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
