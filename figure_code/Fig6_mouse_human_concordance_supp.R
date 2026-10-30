## =====================================================================
## Figure 6 (drawn from the Supplementary Tables). Direction of the glutathione-related response in mouse and human
##   (A) Mouse vs human log2 fold change for the 18 genes testable in both species (concordant / discordant)
##   (B) Human log2 fold change of set members vs all other genes
##   (C) Human estimates with 95% confidence intervals for the 18 genes
##
## Inputs (Supplementary Tables workbook): S9c (mouse-human concordance), S9a (wild-type human hepatocytes, T3 vs vehicle,
##   GEO GSE270233) and S11d (Wald statistic; the standard error is log2 fold change / Wald statistic).
##   Optional: DEG_WT_T3_effect.csv in the data folder (column lfcSE); if present, its standard errors are used instead.
## Output: figures_out_supp_6/ (Fig6 panels and Figure6 as .pdf and .png, Figure6_legend.txt, sessionInfo.txt)
## Usage : Rscript Fig6_mouse_human_concordance_supp.R
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
OUT_DIR <- "figures_out_supp_6"; dir.create(OUT_DIR, showWarnings = FALSE)
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
cg <- read_supp("S9c", "KEY")   # Supplementary Table S9c
cg[, concord := as.logical(concordant)]
if (anyNA(cg$concord)) cg[, concord := sign(mouse_lfc) == sign(human_lfc)]
W <- read_supp("S9a", "gene")[, .(gene, log2FoldChange = WT_T3_log2FC, padj = WT_T3_adj_P)][!is.na(padj)]   # wild-type human hepatocytes, T3 vs vehicle (S9a)
se <- read_supp("S11d", "dataset")[dataset == "human_WT" & !is.na(wald_stat) & wald_stat != 0, .(KEY = gene, human_se = abs(log2FC / wald_stat))]
f_wt <- find_quiet("DEG_WT_T3_effect", "csv")
if (!is.na(f_wt)) { se <- fread(f_wt)[, .(KEY = gene, human_se = lfcSE)]; message("  lfcSE from ", f_wt) }   # exact standard errors when the DESeq2 file is available
cg <- merge(cg, unique(se, by = "KEY"), by = "KEY", all.x = TRUE)   # otherwise: standard error = log2 fold change / Wald statistic (S11d)
n_conc <- sum(cg$concord)
## 일치 개수의 단측 이항 검정(0.5 기준)
bt <- binom.test(n_conc, nrow(cg), 0.5, alternative = "greater")
message(sprintf("유전자 %d개, 방향 일치 %d개, 단측 이항 P = %.4f (원고: 18, 14, 0.015)", nrow(cg), n_conc, bt$p.value))
cg[, cls := factor(ifelse(concord, "concordant", "discordant"), levels = c("concordant", "discordant"))]

## ---- (A) -------------------------------------------------------------------------------------
pA <- ggplot(cg, aes(mouse_lfc, human_lfc, colour = cls)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_point(aes(size = -log10(human_padj))) +
  geom_text_repel(aes(label = KEY), size = 2.2, fontface = "italic", segment.size = 0.2, box.padding = 0.3, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = c(concordant = COL$down, discordant = COL$up), name = NULL) +
  scale_size_continuous(range = c(1, 3.5), name = expression(-log[10]~human~adj.~P)) +
  labs(x = expression(Mouse~log[2]~fold~change), y = expression(Human~log[2]~fold~change), tag = "A") +
  theme(legend.position = "top", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))

## ---- (B) -------------------------------------------------------------------------------------
set_lab <- "Glutathione-set members (n = 18)"; oth_lab <- sprintf("All other genes (n = %s)", comma(sum(!W$gene %in% cg$KEY)))
W[, grp := factor(ifelse(gene %in% cg$KEY, set_lab, oth_lab), levels = c(oth_lab, set_lab))]
mw <- wilcox.test(W[gene %in% cg$KEY]$log2FoldChange, W[!gene %in% cg$KEY]$log2FoldChange, alternative = "less")
med <- W[, .(med = median(log2FoldChange)), by = grp]
message(sprintf("(B) 중앙값: 세트 %.3f, 나머지 %.3f, 단측 Mann-Whitney P = %.4f (원고: -0.27, -0.02, 0.011)", med[grp == set_lab]$med, med[grp == oth_lab]$med, mw$p.value))
pB <- ggplot(W, aes(grp, log2FoldChange)) +
  geom_violin(data = W[grp == oth_lab], fill = "grey85", colour = NA, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  geom_jitter(data = W[grp == set_lab], width = 0.12, size = 1.2, colour = COL$down) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  annotate("text", x = 1.5, y = max(W$log2FoldChange) * 0.92, label = sprintf("one-sided Mann\u2013Whitney P = %s", fmt_p(mw$p.value)), size = 2.3) +
  scale_x_discrete(labels = function(x) sub(" \\(", "\n(", x)) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = NULL, y = expression(Human~log[2]~fold~change), tag = "B") + theme(axis.text.x = element_text(size = 6))

## ---- (C) -------------------------------------------------------------------------------------
cg[, KEY := factor(KEY, levels = cg[order(human_lfc)]$KEY)]
pC <- ggplot(cg, aes(human_lfc, KEY, colour = cls)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_pointrange(aes(xmin = human_lfc - 1.96 * human_se, xmax = human_lfc + 1.96 * human_se), size = 0.3, linewidth = 0.4) +
  geom_text(aes(x = human_lfc + 1.96 * human_se + 0.03, label = sig_stars(human_padj)), size = 3, hjust = 0, vjust = 0.75, show.legend = FALSE) +
  scale_colour_manual(values = c(concordant = COL$down, discordant = COL$up), name = NULL) +
  labs(x = expression(Human~log[2]~fold~change~"(95% CI)"), y = NULL, tag = "C", caption = "* adj. P < 0.05, ** < 0.01") +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "top", legend.margin = margin(0, 0, 0, 0), plot.caption = element_text(size = 5.5, hjust = 0))

save_fig(pA, "Fig6A_mouse_vs_human", 85, 80); save_fig(pB, "Fig6B_distribution", 70, 80); save_fig(pC, "Fig6C_human_estimates", 75, 95)
save_fig((pA | pB) / pC + plot_layout(heights = c(1, 1.1)), "Figure6", 180, 190)

## ---- 범례 txt (이 그림에 들어갈 범례. 숫자는 위에서 계산한 값) -----------------------
legA <- sprintf("(A) Mouse against human log2 fold change for the %d genes testable in both species. %d genes changed in the same direction in human (one-sided binomial P = %s) and are shown in blue; discordant genes are red and point size shows the human \u2212log10 adjusted P.", nrow(cg), n_conc, fmt_p(bt$p.value))
mn <- function(x) sub("^-", "\u2212", sprintf("%.2f", x))
legB <- sprintf("(B) Distribution of human log2 fold change for the %d glutathione-set members (median %s) and for all other tested genes (n = %s, median %s); one-sided Mann\u2013Whitney P = %s.", nrow(cg), mn(med[grp == set_lab]$med), nf(sum(!W$gene %in% cg$KEY)), mn(med[grp == oth_lab]$med), fmt_p(mw$p.value))
legC <- "(C) Human log2 fold change with 95% confidence interval for the same genes (* adjusted P < 0.05, ** < 0.01)."
write_legend("Figure6", paste("Figure 6. Direction of the glutathione-related response in mouse and human.", legA, legB, legC))


## ---- Record R and package versions used to draw the figures --------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
