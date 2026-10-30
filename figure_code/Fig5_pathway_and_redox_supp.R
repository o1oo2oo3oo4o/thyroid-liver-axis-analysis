## =====================================================================
## Figure 5 (drawn from the Supplementary Tables). Functional composition of the hepatic T3 response in mouse liver
##   (A) Leading GO / KEGG / Reactome over-representation terms for the 448 human ortholog entries
##       (background = tested genes)
##   (B) Observed vs expected responsive genes in the 6 redox gene sets (Fisher, raw and Bonferroni)
##   (C) Gene-level fold change of the Gstm family, Gclc, Por, Gpx4, Nqo1 and Hmox1
##   (D) Sensitivity of the glutathione set to inter-gene correlation and to exclusion of transferases
##       (raw vs adjusted)
##
## Inputs (Supplementary Tables workbook): S10b (over-representation results), S11a (redox gene-set results),
##   S11b (cameraPR correlation sensitivity), S11d (gene-level values of the redox sets)
## Output: figures_out_supp_5/ (Fig5 panels and Figure5 as .pdf and .png, Figure5_legend.txt, sessionInfo.txt)
## Usage : Rscript Fig5_pathway_and_redox_supp.R
## =====================================================================
## ---------------------------------------------------------------------------------------------
## [0] Setup
##   * Standalone script: it does not source any other file.
##   * Main input: the Supplementary Tables workbook (file name containing "Supplementary_Tables_S1-S13", .xlsx).
##     It is searched, recursively, in the folder given by the DATA_DIR environment variable (if set), in ./data
##     and in the current working directory; or set SUPP_XLSX to the full path of the workbook.
##     Sheet names with or without the prefix "Table " (e.g. "Table S2a" or "S2a") are accepted.
##   * Public-data inputs are used only where stated in the header; they are searched by name in the same folders.
##   * Input files are only read, never modified. Figures are written to OUT_DIR in the current working directory.
##   * R packages: ggplot2, data.table, readxl, patchwork, ggrepel, scales, ragg, stringr
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
has_supp <- function(key) any(SUPP_SHEETS %in% c(key, paste("Table", key)))
## Reads one sheet of the workbook as a data.table. The header row is the row whose first cell equals `first_col`
## (the rows above it hold the table title and notes).
read_supp <- function(key, first_col) {
  sh <- SUPP_SHEETS[SUPP_SHEETS %in% c(key, paste("Table", key))]
  if (length(sh) != 1) stop(sprintf("Sheet %s not found in %s", key, basename(SUPP)), call. = FALSE)
  top <- read_excel(SUPP, sheet = sh, col_names = FALSE, n_max = 25, col_types = "text", .name_repair = "minimal")
  h <- which(trimws(as.character(top[[1]])) == first_col)[1]
  if (is.na(h)) stop(sprintf("Header row (first cell '%s') not found in sheet %s", first_col, sh), call. = FALSE)
  as.data.table(read_excel(SUPP, sheet = sh, skip = h - 1, .name_repair = "minimal"))
}
OUT_DIR <- "figures_out_supp_5"; dir.create(OUT_DIR, showWarnings = FALSE)
P <- function(key) find_file(FILES[[key]][1], FILES[[key]][2])

## [1] All inputs come from the Supplementary Tables workbook (sheets named in the header)
FILES <- list()

## [2] 공통 설정(색, 글꼴, 저장 함수)
## ---- 색과 테마 (색각 이상자도 구분되는 팔레트) ---------------------------
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
  message("저장: ", file.path(OUT_DIR, paste0(name, ".pdf / .png")))
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
## [5] 그림 만들기
## ---------------------------------------------------------------------------------------------
## ---- (A) ORA ----------------------------------------------------------------------------------
## Supplementary Table S10b: full over-representation output (ontology, ID, Description, Count, FoldEnrichment, p.adjust, geneID)
ONT <- c("GO BP" = "GO biological process", "GO MF" = "GO molecular function", "GO CC" = "GO cellular component", "KEGG" = "KEGG", "Reactome" = "Reactome")
s10b <- read_supp("S10b", "ontology")
stopifnot(all(names(ONT) %in% s10b$ontology))
ora <- s10b[ontology %in% names(ONT), .(ontology = unname(ONT[ontology]), Description, Count, FoldEnrichment, p.adjust)]
top <- ora[p.adjust < 0.05][order(p.adjust), head(.SD, 5), by = ontology]
top[, ontology := factor(ontology, levels = c("GO biological process", "GO molecular function", "GO cellular component", "KEGG", "Reactome"))]
top[, lab := stringr::str_trunc(Description, 52)]
top[, gsh := grepl("glutathione|xenobiotic|oxidoreductase", Description, ignore.case = TRUE)]
top[, lab := paste(ontology, lab, sep = "|")]
top[, lab := factor(lab, levels = rev(top[order(ontology, p.adjust)]$lab))]
pA <- ggplot(top, aes(-log10(p.adjust), lab, fill = gsh)) +
  geom_col(width = 0.7) + geom_text(aes(label = Count), hjust = -0.25, size = 2) +
  facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^|]*\\|", "", x)) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = COL$int), labels = c("other", "xenobiotic / glutathione / oxidoreductase"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = expression(-log[10]~adjusted~P), y = NULL, tag = "A", caption = "Numbers: genes in the term") +
  theme(strip.text.y = element_text(angle = 0, size = 5.5, hjust = 0), axis.text.y = element_text(size = 5.5), legend.position = "top", legend.margin = margin(0, 0, 0, 0), plot.caption = element_text(size = 5.5, hjust = 0))

## ---- (B) 관찰 대 기대 --------------------------------------------------------------------------
S <- read_supp("S11a", "dataset"); m <- S[dataset == "mouse"]   # Supplementary Table S11a
m[, set_short := set]
m[grepl("beta-oxidation", set), set_short := "Fatty acid \u03b2-oxidation"]
m[grepl("P450", set), set_short := "Microsomal P450 oxidoreduction"]
m[grepl("Lipid", set), set_short := "Lipid peroxidation / ferroptosis"]
m[, set_short := factor(set_short, levels = rev(m[order(fisher_p)]$set_short))]
pB <- ggplot(m, aes(y = set_short)) +
  geom_col(aes(x = n_responsive), fill = COL$int, width = 0.6) +
  geom_point(aes(x = expected, shape = "expected by chance"), size = 2.2, colour = "black") +
  geom_text(aes(x = pmax(n_responsive, expected) + 0.3, label = sprintf("P=%s\nBonf. %s", fmt_p(fisher_p), fmt_p(fisher_bonf))), hjust = 0, size = 1.7, lineheight = 0.85) +
  scale_shape_manual(values = c("expected by chance" = 18), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 1.0))) +
  labs(x = "Responsive genes in set (n)", y = NULL, tag = "B") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))

## ---- (C) 유전자 수준 ---------------------------------------------------------------------------
D <- unique(read_supp("S11d", "dataset")[dataset == "mouse"][, .(gene, log2FoldChange = log2FC, padj)], by = "gene")   # Supplementary Table S11d
GC <- c(paste0("Gstm", 1:7), "Gclc", "Por", "Gpx4", "Nqo1", "Hmox1")
C <- D[gene %in% GC][, gene := factor(gene, levels = GC)]
message("(C) 찾은 유전자: ", paste(C$gene, collapse = ", "))
pC <- ggplot(C, aes(gene, log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col(width = 0.7) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = sig_stars(padj), y = log2FoldChange + sign(log2FoldChange) * 0.12), size = 3, vjust = 0.75) +
  scale_fill_manual(values = c("TRUE" = COL$up, "FALSE" = COL$down), guide = "none") +
  labs(x = NULL, y = expression(log[2]~fold~change), tag = "C") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, face = "italic"))

## ---- (D) 민감도 --------------------------------------------------------------------------------
NSET <- 6
rho <- read_supp("S11b", "dataset")[dataset == "mouse" & analysis == "cameraPR inter-gene correlation"][, .(dataset, rho = as.numeric(sub("rho = ", "", setting)), p = camera_p)][order(rho)]   # Supplementary Table S11b
dA <- rho[, .(panel = "cameraPR vs \u03c1", cond = factor(sprintf("\u03c1 = %s", format(rho)), levels = sprintf("\u03c1 = %s", format(sort(rho)))), raw = p)]
sens <- S[dataset == "mouse_glutathione_sensitivity"]
nm <- c(full = "full set (28)", "Gstm collapsed" = "Gstm collapsed", "Gstm removed" = "Gstm removed", "all transferases removed" = "all transferases removed")
dB <- rbind(sens[, .(panel = "Fisher exact", cond = nm[set], raw = fisher_p)], sens[, .(panel = "cameraPR (\u03c1 = 0.01)", cond = nm[set], raw = camera_p)])
dB[, cond := factor(cond, levels = nm)]
dd <- rbind(dA, dB, fill = TRUE); dd[, panel := factor(panel, levels = c("cameraPR vs \u03c1", "Fisher exact", "cameraPR (\u03c1 = 0.01)"))]
dd[, adj := pmin(1, raw * NSET)]
dl <- melt(dd, id.vars = c("panel", "cond"), measure.vars = c("raw", "adj"), variable.name = "type", value.name = "p")
dl[, type := factor(c(raw = "raw P", adj = "adjusted (\u00d76 sets)")[as.character(type)], levels = c("raw P", "adjusted (\u00d76 sets)"))]
pD <- ggplot(dl, aes(cond, -log10(p), shape = type, colour = type)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_point(size = 2, position = position_dodge(0.4)) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  scale_shape_manual(values = c(16, 1), name = NULL) + scale_colour_manual(values = c("black", "grey30"), name = NULL) +
  labs(x = NULL, y = expression(-log[10]~P), tag = "D", caption = "Dashed line: P = 0.05. Adjusted = Bonferroni across the 6 curated sets") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 5.5), legend.position = "top", legend.margin = margin(0, 0, 0, 0), plot.caption = element_text(size = 5.5, hjust = 0))

save_fig(pA, "Fig5A_ORA", 110, 110); save_fig(pB, "Fig5B_redox_sets", 85, 60); save_fig(pC, "Fig5C_gene_level", 85, 60); save_fig(pD, "Fig5D_sensitivity", 150, 65)
save_fig(wrap_plots(A = pA, B = pB, C = pC, D = pD, design = "AB\nAC\nDD", widths = c(1, 1), heights = c(1, 1, 1.15)), "Figure5", 180, 215)

## ---- 범례 txt (이 그림에 들어갈 범례. 숫자는 위에서 계산한 값) -----------------------
nresp <- D[!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1, .N]; ntest <- D[!is.na(padj), .N]
legA <- "(A) Leading over-representation terms (top 5 by adjusted P in each of Gene Ontology biological process, molecular function and cellular component, KEGG and Reactome) for the 448 human orthologue entries of the T3-responsive set, against the tested genes as background. Numbers give the genes in each term and purple bars mark xenobiotic, glutathione or oxidoreductase terms."
legB <- sprintf("(B) Observed (bars) and expected (diamonds) numbers of responsive genes in each of the %d curated redox sets, using %s responsive genes among %s testable genes. One-sided Fisher exact P values and Bonferroni-adjusted values across the %d sets are given.", nrow(m), nf(nresp), nf(ntest), NSET)
legC <- "(C) Gene-level log2 fold changes for the Gstm family, Gclc, Por, Gpx4, Nqo1 and Hmox1 (* adjusted P < 0.05, ** < 0.01, *** < 0.001)."
legD <- sprintf("(D) Sensitivity of the glutathione set to the assumed inter-gene correlation in cameraPR (\u03c1 = %s) and to excluding transferases (Fisher exact test and cameraPR at \u03c1 = 0.01). Filled symbols are raw P values and open symbols Bonferroni-adjusted values (\u00d7 %d sets); the dashed line marks P = 0.05.", paste(format(sort(rho$rho)), collapse = ", "), NSET)
write_legend("Figure5", paste("Figure 5. Functional composition of the hepatic T3 response in mouse liver.", legA, legB, legC, legD))


## ---- Record R and package versions used to draw the figures --------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
