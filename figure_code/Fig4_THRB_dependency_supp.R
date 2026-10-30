## =====================================================================
## Figure 4 (drawn from the Supplementary Tables). Effect of receptor beta loss on the hepatic T3 response in human hepatocytes
##   (A) Wild-type T3 effect vs the genotype x treatment interaction term (all tested genes; altered responses highlighted)
##   (B) Canonical target genes and immediate-early control genes, both genotypes            (needs the standard errors)
##   (C) Genetically prioritized candidates: wild-type response, knockout response, interaction (needs the standard errors)
##   (D) Fraction of altered responses: all wild-type responsive genes vs candidates
##
## Inputs: Supplementary Tables workbook, sheet S9a (DESeq2 results of GEO GSE270233: log2 fold changes and adjusted P of the
##   wild-type effect, the knockout effect and the interaction). The standard errors (lfcSE) needed for panels B and C are read from
##   sheet S9d of the workbook if it exists, otherwise from the three DESeq2 result files DEG_WT_T3_effect.csv,
##   DEG_KO_T3_effect.csv and DEG_Interaction.csv (columns gene, lfcSE) placed in the data folder.
##   If neither is found, panels B and C are replaced by a note and the other panels are drawn.
## Response criterion: adjusted P < 0.05 and |log2 fold change| >= 1 (same as for mouse)
## Output: figures_out_supp_4/ (Fig4 panels and Figure4 as .pdf and .png, Figure4_legend.txt, sessionInfo.txt)
## Usage : Rscript Fig4_THRB_dependency_supp.R
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
OUT_DIR <- "figures_out_supp_4"; dir.create(OUT_DIR, showWarnings = FALSE)
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


## [3] 후보 유전자 목록 (원고 표 3, 표 4 기준)
## ---- 후보 유전자 목록 (원고 표 3, 표 4 기준) ------------------------------
POS10   <- c("CUX2","BDH2","ADH4","LILRB3","SBK1","ABHD1","INCENP","CDC42EP5","CD19","GPCPD1")
PROT16  <- c("ADH4","BTC","CFHR1","CLIC5","COL11A2","COL15A1","DEF6","DNAJB4","F13B","LGALS4","PKLR","POR","SPARCL1","TIMD4","GSTM1","HSPA1A")
PERMISSIVE <- c("CFHR1","GSTM1","HSPA1A")                      # 일대일이 아닌 매핑으로 포함된 후보
MOUSE_OF   <- c(GSTM1 = "Gstm6", HSPA1A = "Hspa1b", F13B = "F13b")   # 그 외는 첫 글자만 대문자
mouse_symbol <- function(g) ifelse(g %in% names(MOUSE_OF), MOUSE_OF[g], paste0(substr(g,1,1), tolower(substr(g,2,nchar(g)))))
PANEL26 <- sort(unique(c(POS10, PROT16, "SLC39A8")))


## ---------------------------------------------------------------------------------------------
## [5] 그림 만들기
## ---------------------------------------------------------------------------------------------
s9a <- read_supp("S9a", "gene")
dt <- s9a[, .(gene, wt_lfc = WT_T3_log2FC, wt_p = WT_T3_adj_P, ko_lfc = KO_T3_log2FC, ko_p = KO_T3_adj_P, int_lfc = Interaction_log2FC, int_p = Interaction_adj_P)]
HAVE_SE <- has_supp("S9d")
if (HAVE_SE) {
  s9d <- read_supp("S9d", "gene")
  dt <- merge(dt, s9d[, .(gene, wt_se = WT_T3_lfcSE, ko_se = KO_T3_lfcSE, int_se = Interaction_lfcSE)], by = "gene", all.x = TRUE)
} else {
  csv_se <- c(wt = "DEG_WT_T3_effect", ko = "DEG_KO_T3_effect", int = "DEG_Interaction")
  f_se <- sapply(csv_se, function(k) find_quiet(k, "csv"))
  if (!anyNA(f_se)) {          # fallback: the three DESeq2 result files of GEO GSE270233 (columns gene, lfcSE)
    for (k in names(f_se)) { d <- fread(f_se[[k]]); message("  lfcSE from ", f_se[[k]]); d <- d[, .(gene, se = lfcSE)]; setnames(d, "se", paste0(k, "_se")); dt <- merge(dt, d, by = "gene", all.x = TRUE) }
    HAVE_SE <- TRUE
  } else {
    message("  Neither Table S9d nor the three DESeq2 files (DEG_WT_T3_effect.csv, DEG_KO_T3_effect.csv, DEG_Interaction.csv) were found: panels B and C are replaced by a note")
    dt[, `:=`(wt_se = NA_real_, ko_se = NA_real_, int_se = NA_real_)]
  }
}
note_panel <- function(tag) ggplot() + annotate("text", x = 0, y = 0, label = "Requires Supplementary Table S9d\n(standard errors of the log2 fold changes)", size = 2.4) + theme_void() + labs(tag = tag)
isresp <- function(p, l) !is.na(p) & p < 0.05 & abs(l) >= 1
dt[, wt_resp := isresp(wt_p, wt_lfc)][, ko_resp := isresp(ko_p, ko_lfc)][, int_sig := isresp(int_p, int_lfc)][, altered := wt_resp & int_sig]
tested <- dt[!is.na(wt_p) & !is.na(int_p)]
message(sprintf("WT 검정 %d, WT 반응 %d, KO 반응 %d, 변화한 반응 %d (원고: 15,297 / 799 / 35 / 194)", sum(!is.na(dt$wt_p)), sum(dt$wt_resp), sum(dt$ko_resp), sum(dt$altered)))
CANON <- c("SLC5A5", "KLF9", "CPM", "MMP11", "EGR1", "FOS"); CAND3 <- c("ADH4", "F13B", "GPCPD1")

## ---- (A) -------------------------------------------------------------------------------------
tested[, cls := factor(fifelse(altered, "altered in KO (n = 194)", fifelse(wt_resp, "WT-responsive, not altered", "other")),
                       levels = c("other", "WT-responsive, not altered", "altered in KO (n = 194)"))]
setorder(tested, cls)
lab <- tested[gene %in% c(CANON, CAND3)]
pA <- ggplot(tested, aes(wt_lfc, int_lfc)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_point(aes(colour = cls, size = cls, alpha = cls)) +
  geom_text_repel(data = lab, aes(label = gene), size = 2.2, fontface = "italic", segment.size = 0.2, box.padding = 0.4, max.overlaps = Inf) +
  scale_colour_manual(values = c("grey80", "#7FB8AE", COL$int), name = NULL) +
  scale_size_manual(values = c(0.3, 0.9, 1.1), guide = "none") + scale_alpha_manual(values = c(0.5, 0.9, 0.95), guide = "none") +
  labs(x = expression(WT~T3~effect~(log[2]~fold~change)), y = expression(Interaction~(log[2]~fold~change)), tag = "A",
       caption = "Dashed line: complete loss of the wild-type response") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), plot.caption = element_text(size = 5.5, hjust = 0))

pA <- pA + guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2)))

## ---- (B) -------------------------------------------------------------------------------------
if (HAVE_SE) {
B <- melt(dt[gene %in% CANON, .(gene, wt_lfc, wt_se, ko_lfc, ko_se)], id.vars = "gene", measure.vars = patterns(lfc = "_lfc$", se = "_se$"), variable.name = "geno")
B[, geno := factor(c("Wild type", "THRB KO")[geno], levels = c("Wild type", "THRB KO"))]
B[, gene := factor(gene, levels = CANON)]
pB <- ggplot(B, aes(gene, lfc, colour = geno)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_pointrange(aes(ymin = lfc - 1.96 * se, ymax = lfc + 1.96 * se), position = position_dodge(0.5), size = 0.35, linewidth = 0.4) +
  scale_colour_manual(values = c("Wild type" = COL$wt, "THRB KO" = COL$ko), name = NULL) +
  scale_x_discrete(labels = c(CANON[1:5], "FOS\n(IEG control)")) +
  labs(x = NULL, y = expression(T3~effect~(log[2]~fold~change)), tag = "B") +
  theme(axis.text.x = element_text(face = "italic", size = 6), legend.position = "top", legend.margin = margin(0, 0, 0, 0))
} else pB <- note_panel("B")

## ---- (C) -------------------------------------------------------------------------------------
cand <- dt[gene %in% PANEL26 & !is.na(wt_p) & !is.na(int_p)]   # KO 가 검정되지 않은 유전자는 KO 칸만 비움
message(sprintf("(C) 패널 26개 중 검정 가능 %d개 (원고: 23개; 검정 불가: %s)", nrow(cand), paste(setdiff(PANEL26, cand$gene), collapse = ", ")))
if (HAVE_SE) {
C <- melt(cand[, .(gene, wt_lfc, wt_se, ko_lfc, ko_se, int_lfc, int_se, hit = gene %in% CAND3)], id.vars = c("gene", "hit"), measure.vars = patterns(lfc = "_lfc$", se = "_se$"), variable.name = "term")
C[, term := factor(c("WT response", "KO response", "Interaction")[term], levels = c("WT response", "KO response", "Interaction"))]
ord <- cand[order(int_lfc)]$gene; C[, gene := factor(gene, levels = ord)]
pC <- ggplot(C, aes(lfc, gene, colour = term, alpha = hit)) +
  geom_vline(xintercept = 0, linewidth = 0.25) + geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 0.25, colour = "grey60") +
  geom_pointrange(aes(xmin = lfc - 1.96 * se, xmax = lfc + 1.96 * se), size = 0.25, linewidth = 0.3, na.rm = TRUE) +
  facet_grid(. ~ term) +
  scale_colour_manual(values = c("WT response" = COL$wt, "KO response" = COL$ko, "Interaction" = COL$int), guide = "none") +
  scale_alpha_manual(values = c("FALSE" = 0.45, "TRUE" = 1), guide = "none") +
  labs(x = expression(log[2]~fold~change~"(95% CI)"), y = NULL, tag = "C") +
  theme(axis.text.y = element_text(face = "italic", size = 6, colour = "black"))
} else pC <- note_panel("C")

## ---- (D) -------------------------------------------------------------------------------------
nW <- sum(dt$wt_resp); nA <- sum(dt$altered)
cw <- dt[gene %in% PANEL26 & wt_resp]; kc <- nrow(cw); ka <- sum(cw$altered)
fp <- fisher.test(matrix(c(ka, kc - ka, nA - ka, nW - kc - (nA - ka)), 2), alternative = "greater")$p.value
message(sprintf("(D) 전체 %d/%d (%.1f%%), 후보 %d/%d, 조건부 Fisher P = %.4f (원고 0.014)", nA, nW, 100 * nA / nW, ka, kc, fp))
Dd <- data.table(group = factor(c("All WT-responsive genes", "Prioritized candidates"), levels = c("All WT-responsive genes", "Prioritized candidates")),
                 altered = c(nA, ka), total = c(nW, kc))[, prop := altered / total]
## 비율의 95% 신뢰구간(Clopper-Pearson). 후보는 n = 3 이라 3/3 이어도 하한이 낮다는 점이 그림에 보이도록 한다
Dd[, c("lo", "hi") := { z <- mapply(function(x, n) binom.test(x, n)$conf.int, altered, total); list(z[1, ], z[2, ]) }]
pD <- ggplot(Dd, aes(group, prop, colour = group)) +
  geom_hline(yintercept = nA / nW, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.7, linewidth = 0.7) +
  geom_text(aes(y = hi + 0.03, label = sprintf("%d / %d (%.1f%%)\n95%% CI %.0f\u2013%.0f%%", altered, total, 100 * prop, 100 * lo, 100 * hi)), vjust = 0, size = 2.3, lineheight = 0.95, show.legend = FALSE) +
  annotate("segment", x = 1, xend = 2, y = 1.2, yend = 1.2, linewidth = 0.3) + annotate("text", x = 1.5, y = 1.23, label = sprintf("one-sided Fisher P = %s", fmt_p(fp)), size = 2.3) +
  scale_colour_manual(values = c("grey45", COL$int), guide = "none") +
  scale_y_continuous(labels = percent, limits = c(0, 1.32), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Proportion with altered response", tag = "D", caption = "Dashed line: proportion among all WT-responsive genes.\nError bars: Clopper-Pearson 95% CI.") +
  theme(axis.text.x = element_text(size = 6), plot.caption = element_text(size = 5.5, hjust = 0))

save_fig(pA, "Fig4A_WT_vs_interaction", 95, 85); save_fig(pB, "Fig4B_canonical_targets", 85, 70)
save_fig(pC, "Fig4C_candidates", 120, 100); save_fig(pD, "Fig4D_altered_proportion", 75, 75)
save_fig((pA | pB) / (pC | pD) + plot_layout(heights = c(1, 1.15)), "Figure4", 180, 190)

## ---- 범례 txt (이 그림에 들어갈 범례. 숫자는 위에서 계산한 값) -----------------------
legA <- sprintf("(A) Wild-type T3 treatment effect against the genotype \u00d7 treatment interaction for the %s genes tested in both analyses. Genes with a significantly altered response in knockout cells (%d of the %d wild-type responsive genes) are highlighted; the dashed line marks complete loss of the wild-type response. Responsive and altered responses were defined as adjusted P < 0.05 and |log2 fold change| \u2265 1.", nf(nrow(tested)), nA, nW)
legB <- "(B) T3 effect (log2 fold change with 95% confidence interval) at 5 canonical targets and the immediate-early control FOS in wild-type and THRB knockout cells."
legC <- sprintf("(C) Wild-type response, knockout response and interaction (log2 fold change with 95%% confidence interval) for the %d genetically prioritized candidates that were expressed in hepatocytes; the %d responsive in wild-type cells are emphasized.", nrow(cand), kc)
legD <- sprintf("(D) Proportion of wild-type responsive genes with a significantly altered response: %d of %d (%.1f%%) among all genes and %d of %d (Clopper-Pearson 95%% confidence interval %.0f to %.0f%%) among the genetically prioritized candidates; one-sided Fisher P = %s. The dashed line marks the proportion among all wild-type responsive genes.", nA, nW, 100 * nA / nW, ka, kc, 100 * Dd$lo[2], 100 * Dd$hi[2], fmt_p(fp))
write_legend("Figure4", paste("Figure 4. Effect of receptor beta loss on the hepatic T3 response in human hepatocytes.", legA, legB, legC, legD))


## ---- Record R and package versions used to draw the figures --------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
