## =====================================================================
## Figure 3 (drawn from the Supplementary Tables). The protein-level candidate set
##   (A) Mendelian randomization (MR) support: outcome x estimator x instrument source x adjustment model
##       (13 primary + 3 permissive candidates)
##   (B) Hepatic T3 response (log2 fold change, mouse) of the same 16 genes
##   (C) Plasma associations of the measurable candidates with the 6 liver phenotypes
##
## Inputs (all in the Supplementary Tables workbook):
##   S7   MR matrix. An estimate supports a direction when FDR P < 0.05 and the MR-Egger pleiotropy P is not < 0.05
##        (this reproduces the supporting-estimate counts n of Table 4); the direction is the sign of beta
##   S6c  DESeq2 results (GEO GSE85793) of the mouse genes matched to the 16 candidates
##   S8   plasma proteomic associations (beta or OR, FDR P)
## Output: figures_out_supp_3/ (Fig3 panels and Figure3 as .pdf and .png, Figure3_legend.txt, sessionInfo.txt)
## Usage : Rscript Fig3_protein_candidates_supp.R
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
OUT_DIR <- "figures_out_supp_3"; dir.create(OUT_DIR, showWarnings = FALSE)
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
## ---- (A) sum_MR 시트 파싱: 머리글 4줄(출처, 추정법, 결과)을 앞으로 채워 열 의미를 만든다 ----------
## Supplementary Table S7: one row per protein x method (IVW, WM); per outcome 3 columns (Beta (95% CI), FDR P, Pleiotropy P);
## the instrument source (deCODE or ARIC) is written above the outcome label.
M <- as.data.frame(read_excel(SUPP, sheet = SUPP_SHEETS[SUPP_SHEETS %in% c("S7", "Table S7")], col_names = FALSE, col_types = "text", .name_repair = "minimal"))
hr <- which(trimws(M[[1]]) == "Protein")[1]
fill_fwd <- function(v) { for (i in seq_along(v)) if (is.na(v[i]) && i > 1) v[i] <- v[i - 1]; v }
src_row <- fill_fwd(as.character(unlist(M[hr - 2, ]))); out_row <- fill_fwd(as.character(unlist(M[hr - 1, ])))
bcols <- which(trimws(as.character(unlist(M[hr, ]))) == "Beta (95% CI)")
body <- M[(hr + 1):nrow(M), ]; body[[1]] <- fill_fwd(as.character(body[[1]]))
num1  <- function(x) suppressWarnings(as.numeric(sub("^\\s*(-?[0-9.]+).*$", "\\1", x)))
pnum  <- function(x) suppressWarnings(as.numeric(gsub("[*\u00a7<]", "", x)))
OUT_LAB <- c("Liver PDFFnon" = "PDFF unadj.", "Liver PDFFbmi" = "PDFF BMI", "Liver PDFFwfm" = "PDFF WFM", "Liver PDFFvat" = "PDFF VAT",
             "Liver cT1non" = "cT1 unadj.", "Liver cT1bmi" = "cT1 BMI", "Liver cT1wfm" = "cT1 WFM", "Liver cT1vat" = "cT1 VAT",
             "FLI" = "FLI", "FIB4" = "FIB-4", "APRI" = "APRI", "MASLD" = "MASLD")
stopifnot(all(out_row[bcols] %in% names(OUT_LAB)))
long <- rbindlist(lapply(PROT16, function(g) {
  r <- body[which(body[[1]] == g), , drop = FALSE]
  if (!nrow(r)) return(NULL)
  rbindlist(lapply(seq_len(nrow(r)), function(i) rbindlist(lapply(bcols, function(j) {
    b <- num1(r[i, j]); p <- pnum(r[i, j + 1]); ple <- pnum(r[i, j + 2])
    ok <- !is.na(b) && !is.na(p) && p < 0.05 && !(!is.na(ple) && ple < 0.05)          # FDR P < 0.05 and no pleiotropy signal
    data.table(gene = g, source = sub("pQTL source: ", "", src_row[j]), est = as.character(r[i, 2]), outcome = unname(OUT_LAB[out_row[j]]),
               dir = if (ok) ifelse(b > 0, "up", "down") else "no support")
  }))))
}))
long[, panel := factor(paste(source, est), levels = c("deCODE IVW", "deCODE WM", "ARIC IVW", "ARIC WM"))]
long[, outcome := factor(outcome, levels = unname(OUT_LAB))]
long[, dir := factor(dir, levels = c("down", "up", "no support"))]
long[, gene := factor(gene, levels = rev(PROT16))]
ylab <- function(x) ifelse(x %in% PERMISSIVE, paste0(x, "*"), x)
n_est <- long[dir != "no support", .N, by = gene]; message("(A) 지지 추정치 수 (표 4의 n과 같아야 함):"); print(n_est[match(PROT16, gene)])
## 한 칸 = 유전자 x 결과. 점 크기와 숫자 = 그 결과를 지지하는 추정치 수(최대 4 = deCODE/ARIC x IVW/WM), 색 = 방향
## 한 칸에 위/아래가 모두 있으면 "mixed". 유전자 이름 옆 n = 지지 추정치 총수(표 4의 n과 같다)
agg <- long[dir != "no support", .(n = .N, nup = sum(dir == "up"), ndown = sum(dir == "down")), by = .(gene, outcome)]
agg[, cls := factor(fifelse(nup > 0 & ndown > 0, "mixed", fifelse(nup > 0, "up", "down")), levels = c("down", "up", "mixed"))]
gene_n <- long[dir != "no support", .N, by = gene]
ylab_n <- setNames(sprintf("%s  (n = %d)", ylab(as.character(gene_n$gene)), gene_n$N), as.character(gene_n$gene))
pA <- ggplot(agg, aes(outcome, gene)) +
  geom_vline(xintercept = c(4.5, 8.5), colour = "grey75", linewidth = 0.3) +
  geom_point(aes(size = n, colour = cls)) +
  geom_text(aes(label = n), colour = "white", size = 2.1, fontface = "bold") +
  scale_size_continuous(range = c(3.4, 7.2), breaks = 1:4, limits = c(1, 4), name = "Supporting estimates") +
  scale_colour_manual(values = c(down = COL$down, up = COL$up, mixed = COL$int), name = "Genetically predicted\nprotein effect", drop = FALSE) +
  scale_x_discrete(drop = FALSE) + scale_y_discrete(labels = ylab_n, drop = FALSE) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 4)), size = guide_legend(order = 2)) +
  labs(x = NULL, y = NULL, tag = "A", caption = "* permissive (non one-to-one) orthologue mapping.  Estimates: deCODE and ARIC x IVW and weighted median (up to 4 per outcome).") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 6), axis.text.y = element_text(face = "italic", size = 6.5), panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
        plot.caption = element_text(size = 5.5, hjust = 0), legend.position = "right", axis.line = element_blank(), axis.ticks = element_blank())

## ---- (B) 마우스 반응 ----------------------------------------------------------------------
D <- read_supp("S6c", "human_gene")[, .(gene = mouse_gene, log2FoldChange, padj)]; D <- D[!is.na(padj)]   # DESeq2 of the mouse genes matched to the 16 candidates
B <- data.table(gene = PROT16, mouse = mouse_symbol(PROT16))
B <- merge(B, D[, .(mouse = gene, lfc = log2FoldChange, padj)], by = "mouse", all.x = TRUE)
stopifnot(!anyNA(B$lfc))
B[, dir := ifelse(lfc > 0, "up", "down")]; B[, gene := factor(gene, levels = B[order(lfc)]$gene)]
pB <- ggplot(B, aes(lfc, gene, fill = dir)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = paste0(fmt_p(padj)), hjust = ifelse(lfc > 0, -0.1, 1.1)), size = 1.9) +
  scale_fill_manual(values = c(down = COL$down, up = COL$up), guide = "none") +
  scale_y_discrete(labels = ylab) + scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(x = expression(log[2]~fold~change~(T3~vs~vehicle)), y = NULL, tag = "B") +
  theme(axis.text.y = element_text(face = "italic", size = 6.5))
message(sprintf("(B) 감소 %d, 증가 %d", sum(B$lfc < 0), sum(B$lfc > 0)))

## ---- (C) 혈장 단백체 연관 -----------------------------------------------------------------
## Supplementary Table S8: wide table, one block of 2 columns (beta or OR with 95% CI, FDR P) per phenotype
s8 <- as.data.frame(read_excel(SUPP, sheet = SUPP_SHEETS[SUPP_SHEETS %in% c("S8", "Table S8")], col_names = FALSE, col_types = "text", .name_repair = "minimal"))
h8 <- which(trimws(s8[[1]]) == "Protein")[1]; ph8 <- as.character(unlist(s8[h8 - 1, ]))
TRAITCODE <- c("MRI-PDFF" = "PDFF", "MRI-cT1" = "cT1", "MASLD" = "v56_wDM2", "FLI" = "fli", "FIB4" = "fib4_0", "APRI" = "apri_0")
stopifnot(all(ph8[!is.na(ph8)] %in% names(TRAITCODE)))
rows8 <- (h8 + 1):nrow(s8)
O <- rbindlist(lapply(which(!is.na(ph8)), function(j) data.table(protein = as.character(s8[rows8, 1]), trait = unname(TRAITCODE[ph8[j]]),
        analysis = ifelse(ph8[j] == "MASLD", "glm", "lm"), Effect = num1(s8[rows8, j]), adjP_bh = pnum(s8[rows8, j + 1]))))
TR <- c(PDFF = "PDFF", cT1 = "cT1", fli = "FLI", fib4_0 = "FIB-4", apri_0 = "APRI", v56_wDM2 = "MASLD")
measured <- PROT16[PROT16 %in% O$protein]
C <- O[protein %in% measured & trait %in% names(TR)]
C[, sig := adjP_bh < 0.05]
C[, dir := ifelse(analysis == "glm", ifelse(Effect > 1, "up", "down"), ifelse(Effect > 0, "up", "down"))]   # 로지스틱은 Effect=OR(>1이 증가)
C[, cls := factor(ifelse(sig, dir, "n.s."), levels = c("down", "up", "n.s."))]
C[, trait := factor(TR[trait], levels = unname(TR))]; C[, protein := factor(protein, levels = rev(measured))]
pC <- ggplot(C, aes(trait, protein)) +
  geom_point(aes(colour = cls, size = pmin(-log10(adjP_bh), 8))) +
  scale_colour_manual(values = c(down = COL$down, up = COL$up, "n.s." = "grey75"), name = "Association\n(BH P < 0.05)") +
  scale_size_continuous(range = c(0.8, 4), breaks = c(2, 4, 8), labels = c("2", "4", "\u22658"), name = expression(-log[10]~adj.~P)) +
  scale_y_discrete(labels = ylab) + labs(x = NULL, y = NULL, tag = "C") +
  theme(axis.text.y = element_text(face = "italic", size = 6.5), axis.text.x = element_text(angle = 45, hjust = 1))
message(sprintf("(C) 측정된 단백질 %d개 (원고: 10개)", length(measured)))

save_fig(pA, "Fig3A_MR_support", 150, 85)
save_fig(pB, "Fig3B_mouse_response", 85, 85)
save_fig(pC, "Fig3C_plasma_associations", 95, 85)
save_fig(pA / (pB | pC) + plot_layout(heights = c(1.1, 1)), "Figure3", 180, 185)

## ---- 범례 txt (이 그림에 들어갈 범례. 숫자는 위에서 계산한 값) -----------------------
tot_est <- sum(gene_n$N); n_dec <- sum(B$lfc < 0); n_inc <- sum(B$lfc > 0)
legA <- sprintf("(A) Mendelian randomization support for the %d protein-level candidates, %d primary and %d with permissive orthologue mapping (asterisks). Each dot is one outcome (proton density fat fraction, PDFF, and corrected T1, cT1, each unadjusted and with the 3 adiposity adjustments given in Table 1, plus the fatty liver index, FIB-4, APRI and MASLD). Dot size and number show the supporting estimates (up to 4: deCODE and ARIC instruments \u00d7 inverse-variance weighted and weighted median), colour shows the direction of the genetically predicted protein effect, and the number beside each gene gives its total supporting estimates.", length(PROT16), length(PROT16) - length(PERMISSIVE), length(PERMISSIVE))
legB <- sprintf("(B) Hepatic T3 response in mouse liver for the same %d genes (%d decreased, %d increased); labels give adjusted P values.", length(PROT16), n_dec, n_inc)
legC <- sprintf("(C) Associations of the %d measurable candidates with 6 liver phenotypes (PDFF, cT1, fatty liver index, FIB-4, APRI and MASLD) in plasma proteomic data. Colour shows the direction of associations with Benjamini-Hochberg adjusted P < 0.05 (grey, not significant) and size shows \u2212log10 adjusted P (capped at 8). For the binary MASLD outcome, direction is taken from the odds ratio.", length(measured))
write_legend("Figure3", paste("Figure 3. The protein-level candidate set.", legA, legB, legC))


## ---- Record R and package versions used to draw the figures --------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
