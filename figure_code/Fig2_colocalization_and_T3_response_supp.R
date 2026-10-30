## =====================================================================
## Figure 2 (drawn from the Supplementary Tables). Shared thyroid-liver genetic signal and hormone responsiveness
##   (A) Cross-domain colocalization PP.H4 of the 458 locus windows (Manhattan-style plot)
##   (B) Regional association plots of three regions (9q34, 4q24, 12q24); windows, lead variants and genes
##   (C) Volcano plot of the hepatic T3 response of the 11,819 testable mouse genes
##
## Inputs:
##   Supplementary Tables workbook: S2a (window summary), S2b (colocalization of all dataset pairs), S4 (gene-to-window assignment)
##   GSE85793_mouse_liver_T3_DESeq2_all_genes.csv   DESeq2 results of the public dataset GEO GSE85793 (mouse liver, T3 vs vehicle,
##       4 mice per group), all tested genes; columns gene, baseMean, log2FoldChange, lfcSE, pvalue, padj (Methods, section 2.4).
##       Required for panel C only (if the file is absent, panel C is replaced by a note and A and B are drawn). The file name
##       is a label chosen for this repository; it is generated from the GEO counts.
##   Optional (panel B only; if absent a window/gene layout is drawn instead): public GWAS summary statistics listed in [1-b]
##   -- GWAS Catalog GCST90446646-GCST90446653, GCST010653, GCST90296333, and the ThyroidOmics Consortium
##   files FT3_summary_data.txt and FT4_summary_data.txt. They are not part of the Supplementary Tables and are not redistributed.
## Output: figures_out_supp_2/ (Fig2A/Fig2B/Fig2C and Figure2 as .pdf and .png, Figure2_legend.txt, sessionInfo.txt)
## Usage : Rscript Fig2_colocalization_and_T3_response_supp.R
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
OUT_DIR <- "figures_out_supp_2"; dir.create(OUT_DIR, showWarnings = FALSE)
P <- function(key) find_file(FILES[[key]][1], FILES[[key]][2])

## [1] Public-data input of panel C (all other inputs come from the Supplementary Tables workbook)
FILES <- list(
  deg_mouse = c("GSE85793_mouse_liver_T3_DESeq2_all_genes", "csv")
)

## [1-b] 지역 연관 플롯(그림 2B)에 쓰는 GWAS 요약통계 (선택 사항: 파일을 못 찾으면 창/유전자 배치도로 자동 대체)
##   파일 이름과 열 이름은 COLOC 분석 스크립트(260518_1_COLOC_thyroid_liver.R)의 GWAS_FILES와 같다.
##   간 MRI GWAS : GCST90446646~49(cT1), GCST90446650~53(PDFF)  .tsv
##   갑상선      : FT3_summary_data.txt, FT4_summary_data.txt, GCST010653.txt(TSH A), GCST90296333.tsv(TSH B)
##   각 항목 = c(파일 이름, 확장자, 염색체 열, 위치 열, P값 열)
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


## [4] COLOC 판정과 유전자 배치 함수
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
## 교차 도메인 판정 규칙: 간 형질은 보정 모델 4개 중 최댓값, TSH는 두 데이터셋 중 낮은 값(둘 다 있어야 함),
## 창의 값은 6개 trait 수준 쌍(PDFF/cT1 x TSH/FT3/FT4) 중 최댓값
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

## ---- 유전자와 창 ---------------------------------------------------------
load_genes <- function() {
  s4 <- read_supp("S4", "window")   # gene-to-window assignment of the 27 supported windows
  G <- s4[, .(Locus_ID = window, Gene = gene, gene_start = gene_start_GRCh37, gene_end = gene_end_GRCh37, biotype, GWAS_shared = in_shared_199, locus_pass = TRUE)]
  G[, gene_cat := fcase(Gene %in% POS10, "T3-responsive positional candidate",
                        GWAS_shared == 1, "shared-domain gene", default = "other")]
  G[]
}
## 겹치지 않게 행(row) 배정
pack_rows <- function(start, end, pad = 0) {
  o <- order(start); row <- integer(length(start)); ends <- numeric(0)
  for (i in o) { r <- which(ends + pad < start[i])[1]; if (is.na(r)) { ends <- c(ends, end[i]); r <- length(ends) } else ends[r] <- end[i]; row[i] <- r }
  row
}
## 한 구역(여러 창) 배치도: 창(PP.H4 색), lead 변이, 유전자
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

## ---- GWAS 구간 읽기와 지역 연관 플롯 ------------------------------------
## 한 번 읽은 구간은 OUT_DIR 의 cache_gwas_*.rds 에 저장해 다음 실행부터 바로 불러온다(원본 GWAS 파일은 읽기만 한다)
read_gwas <- function(key, wins) {          # wins: data.table(region, chr, start, end)
  spec <- GW[[key]]; cache <- file.path(OUT_DIR, paste0("cache_gwas_", key, ".rds"))
  sig <- paste(wins$region, wins$chr, wins$start, wins$end, collapse = "|")
  if (file.exists(cache)) { z <- readRDS(cache); if (identical(z$sig, sig)) { message("  캐시 사용: ", key); return(z$dat) } }
  f <- find_file(spec[1], spec[2])
  d <- fread(f, select = spec[3:5], showProgress = FALSE); setnames(d, spec[3:5], c("chr", "pos", "p"))
  d[, chr := chr_num(chr)]; d[, p := pmax(as.numeric(p), 1e-300)]; d <- d[!is.na(p) & !is.na(pos)]
  out <- rbindlist(lapply(seq_len(nrow(wins)), function(i) d[chr == wins$chr[i] & pos >= wins$start[i] & pos <= wins$end[i]][, region := wins$region[i]]))
  saveRDS(list(sig = sig, dat = out), cache); out
}
gene_track <- function(loci, G, xlim, always_label = character(), gene_text = 1.9) {
  g <- unique(G[Locus_ID %in% loci & biotype == "coding", .(Gene, gene_start, gene_end, gene_cat)])
  g[, `:=`(gmin = gene_start / 1e6, gmax = gene_end / 1e6)]; g[, grow := pack_rows(gmin, gmax, 0.05)]
  p <- ggplot(g) +
    geom_segment(aes(x = gmin, xend = gmax, y = -grow, yend = -grow, colour = gene_cat), linewidth = 1.3, lineend = "butt") +
    scale_colour_manual(values = c("T3-responsive positional candidate" = COL$cand, "shared-domain gene" = COL$shared, "other" = "grey70"),
                        limits = c("T3-responsive positional candidate", "shared-domain gene", "other"), drop = FALSE, name = "Gene") +
    coord_cartesian(xlim = xlim) + scale_y_continuous(NULL, breaks = NULL, expand = expansion(mult = c(0.2, 0.1))) +
    labs(x = "Position (Mb, GRCh37)") + theme(axis.line.y = element_blank())
  gl <- g[gene_cat != "other" | Gene %in% always_label]
  if (nrow(gl)) p <- p + geom_text_repel(data = gl, aes(x = (gmin + gmax) / 2, y = -grow, label = Gene, colour = gene_cat), size = gene_text, fontface = "italic",
                                         segment.size = 0.15, min.segment.length = 0, box.padding = 0.25, max.overlaps = Inf, direction = "both", nudge_y = -0.35, force = 3, show.legend = FALSE)
  p
}
## 위: 갑상선 형질의 -log10 P, 아래: 간 형질의 -log10 P(뒤집어 그림). 창은 PP.H4 색의 배경, 마름모는 lead 변이
region_assoc <- function(r_name, loci, pair, thy, liv, S1, wb, G, cap = 50, always_label = character()) {
  w <- merge(S1[Locus_ID %in% loci, .(locus_id = Locus_ID, Start_bp, End_bp, lead_pos)], wb[, .(locus_id, h4)], by = "locus_id", all.x = TRUE)
  xlim <- c(min(w$Start_bp), max(w$End_bp)) / 1e6
  w <- w[order(Start_bp)]; w[, ypos := cap * c(1.20, 1.08)[(seq_len(.N) %% 2) + 1]]     # 겹치는 창의 라벨은 높이를 번갈아 배치
  d <- rbind(thy[, .(x = pos / 1e6, y =  pmin(-log10(p), cap), trait = "thyroid")], liv[, .(x = pos / 1e6, y = -pmin(-log10(p), cap), trait = "liver")])
  p <- ggplot() +
    geom_rect(data = w, aes(xmin = Start_bp / 1e6, xmax = End_bp / 1e6, ymin = -Inf, ymax = Inf, fill = h4), alpha = 0.16) +
    geom_hline(yintercept = c(-1, 1) * -log10(1e-6), linetype = "dashed", linewidth = 0.25, colour = "grey40") + geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_point(data = d, aes(x, y, colour = trait), size = 0.3, alpha = 0.7) +
    geom_point(data = w, aes(lead_pos / 1e6, 0), shape = 18, size = 2.6) +
    geom_text(data = w, aes((Start_bp + End_bp) / 2e6, ypos, label = sprintf("%s (PP.H4 %.2f)", locus_id, h4)), size = 1.7) +
    scale_fill_gradient(low = "#FEE8C8", high = "#B30000", limits = c(0.5, 1), oob = squish, name = "PP.H4") +
    scale_colour_manual(values = c(thyroid = COL$thy, liver = COL$liv), labels = c(thyroid = "Thyroid trait (up)", liver = "Liver trait (down)"), name = NULL) +
    scale_y_continuous(labels = function(v) abs(v), limits = c(-cap * 1.05, cap * 1.30)) + coord_cartesian(xlim = xlim) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2)) +
    labs(x = NULL, y = expression(-log[10]~italic(P)), title = r_name, subtitle = sprintf("Up: %s
Down: %s", pair$thy_lab, pair$liv_lab)) + theme(plot.subtitle = element_text(size = 5.8, lineheight = 1)) + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
  (p / gene_track(loci, G, xlim, always_label)) + plot_layout(heights = c(3.2, 1))
}

## ---------------------------------------------------------------------------------------------
## [5] 그림 만들기
## ---------------------------------------------------------------------------------------------
co <- load_coloc(); S1 <- co$S1; S2 <- co$S2
G  <- load_genes()
wb <- window_best(S2)$best
W  <- merge(S1[, .(locus_id = Locus_ID, chr, lead_pos, Final_tier, Lead_pval)], wb[, .(locus_id, h4, pair, lt, tt, liver)], by = "locus_id", all.x = TRUE)
stopifnot(nrow(W) == 458)

## ---- (A) 창별 교차 도메인 PP.H4 (Manhattan 형태) ---------------------------
## 점 색 = 판정 등급, 점 모양 = 신호를 만든 쌍의 간 형질(PDFF 원, cT1 삼각형), x = 창 안 lead 변이의 위치
## 교차 도메인 결과가 없는 창(결과 없음 1개 + 간 데이터셋 쌍만 계산된 6개)은 y = 0 에 x 로 표시한다
off <- c(0, cumsum(as.numeric(HG19_LEN)))[1:23]
W[, x := (off[chr] + lead_pos) / 1e6]
TIERS <- c("high (PP.H4 \u2265 0.8)", "moderate (0.5\u20130.8)", "below 0.5", "no cross-domain result")
W[, tier := factor(fcase(is.na(h4), TIERS[4], h4 >= 0.8, TIERS[1], h4 >= 0.5, TIERS[2], default = TIERS[3]), levels = TIERS)]
W[, shp := factor(fifelse(is.na(h4), "no cross-domain result", as.character(lt)), levels = c("PDFF", "cT1", "no cross-domain result"))]
W[, y := fifelse(is.na(h4), 0, h4)]
chr_mid <- data.table(chr = 1:23, a = off / 1e6, b = c(off[-1], sum(HG19_LEN)) / 1e6)[, mid := (a + b) / 2]
LAB <- c(L00114 = "SLC39A8", L00237 = "ABO", L00304 = "ATXN2", L00056 = "GCKR", L00303 = "CUX2")
pA <- ggplot(W, aes(x, y)) +
  geom_rect(data = chr_mid[chr %% 2 == 0], aes(xmin = a, xmax = b, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey95") +
  geom_hline(yintercept = c(0.5, 0.8), linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_point(aes(colour = tier, shape = shp, size = tier), alpha = 0.9, stroke = 0.6) +
  geom_text_repel(data = W[locus_id %in% names(LAB)], aes(label = sprintf("%s (%s)", LAB[locus_id], pair)), size = 2.2, fontface = "italic", min.segment.length = 0, segment.size = 0.2, box.padding = 0.4, max.overlaps = Inf) +
  scale_colour_manual(values = c(COL$high, COL$moderate, COL$below, "grey35"), name = "Tier") +
  scale_shape_manual(values = c(PDFF = 16, cT1 = 17, "no cross-domain result" = 4), name = "Liver trait of the best pair") +
  scale_size_manual(values = c(2.1, 1.8, 1.1, 1.3), guide = "none") +
  scale_x_continuous(breaks = chr_mid$mid, labels = c(1:22, "X"), expand = expansion(mult = 0.005)) +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.25)) +
  guides(colour = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 2, nrow = 1, override.aes = list(colour = "grey30"))) +
  labs(x = "Chromosome", y = "Cross-domain PP.H4 (per window)", tag = "A") +
  theme(legend.position = "top", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0), panel.grid = element_blank(), axis.text.x = element_text(size = 5.5))
message(sprintf("(A) 창 %d개 중 PP.H4 계산 %d개 (결과 없음 %d개), high %d, moderate %d", nrow(W), sum(!is.na(W$h4)), sum(is.na(W$h4)), sum(W$h4 >= 0.8, na.rm = TRUE), sum(W$h4 >= 0.5 & W$h4 < 0.8, na.rm = TRUE)))

## ---- (B) 세 구역의 지역 연관 플롯 ---------------------------------------
## 구역마다 PP.H4 가 가장 높은 창의 쌍(갑상선 형질 x 간 형질)을 골라, 갑상선은 위쪽, 간은 아래쪽(뒤집어서)에 -log10 P 를 그린다
## TSH 쌍은 데이터셋 B(GCST90296333)를 그린다. 점선 = P 1e-6(창 정의 기준), 마름모 = lead 변이
REG <- list("9q34" = "L00237", "4q24" = c("L00113", "L00114", "L00115"), "12q24" = c("L00303", "L00304", "L00305"))
reg_ext <- rbindlist(lapply(names(REG), function(r) { s <- S1[Locus_ID %in% REG[[r]]]; data.table(region = r, chr = s$chr[1], start = min(s$Start_bp), end = max(s$End_bp)) }))
pairs <- lapply(REG, function(l) { b <- wb[locus_id %in% l][which.max(h4)]; list(thy = if (b$tt == "TSH") "TSH_B" else b$tt, liv = b$liver, win = b$locus_id, pair = b$pair) })
need <- unique(unlist(lapply(pairs, function(p) c(p$thy, p$liv))))
have <- sapply(need, function(k) !is.na(find_quiet(GW[[k]][1], GW[[k]][2])))
if (all(have)) {
  message("  지역 연관 플롯에 쓰는 형질: ", paste(need, collapse = ", "))
  gdat <- setNames(lapply(need, function(k) { use <- names(Filter(function(p) k %in% c(p$thy, p$liv), pairs)); read_gwas(k, reg_ext[region %in% use]) }), need)
  pB <- lapply(names(REG), function(r) { p <- pairs[[r]]; p$thy_lab <- GW_LABEL[[p$thy]]; p$liv_lab <- GW_LABEL[[p$liv]]
  region_assoc(r, REG[[r]], p, gdat[[p$thy]][region == r], gdat[[p$liv]][region == r], S1, wb, G, always_label = c("ATXN2", "ABO")) })
} else {
  message("  GWAS 요약통계를 찾지 못해 창/유전자 배치도로 대신합니다. 못 찾은 파일: ", paste(need[!have], collapse = ", "))
  pB <- lapply(names(REG), function(r) region_layout(REG[[r]], S1, wb, G, title = r, always_label = c("ATXN2", "ABO")) + theme(legend.position = if (r == "12q24") "right" else "none"))
}
if (all(have)) if (all(have)) pB[[1]][[1]] <- pB[[1]][[1]] + labs(tag = "B") else pB[[1]] <- pB[[1]] + labs(tag = "B")   # regional plots are patchworks, layouts are single plots else pB[[1]] <- pB[[1]] + labs(tag = "B")   # regional plots are patchworks, layouts are single plots
pBc <- wrap_plots(pB, nrow = 1, widths = c(1, 1.1, 1.3), guides = "collect") & theme(legend.position = "bottom")

## ---- (C) 마우스 간 T3 반응 volcano plot -------------------------------
## 반응 기준(adjusted P < 0.05 와 |log2FC| >= 1)이 두 선으로 모두 보인다
HAVE_DEG <- !is.na(find_quiet(FILES$deg_mouse[1], FILES$deg_mouse[2]))
if (HAVE_DEG) {
D <- fread(P("deg_mouse")); D <- D[!is.na(padj)]
D[, resp := padj < 0.05 & abs(log2FoldChange) >= 1]
D[, nlp := -log10(pmax(padj, 1e-300))]
## 위치 유전자(지지 창 안의 코딩 유전자)는 사람 symbol 을 대문자로 바꿔 마우스 symbol 과 직접 일치시킨다(원고와 같은 방식)
pos_genes <- unique(toupper(G[locus_pass == TRUE & biotype == "coding", Gene]))
D[, human := toupper(gene)]
D[, evaluated := human %in% pos_genes]
D[, resp_ev := evaluated & resp]
D[evaluated == FALSE, human := NA_character_]
pC <- ggplot(D, aes(log2FoldChange, nlp)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_point(data = D[evaluated == FALSE], colour = "grey82", size = 0.35) +
  geom_point(data = D[evaluated == TRUE & !(resp_ev %in% TRUE)], aes(colour = "evaluated (symbol-matched), not responsive"), size = 0.9) +
  geom_point(data = D[resp_ev %in% TRUE], aes(colour = "evaluated and responsive"), size = 1.7) +
  geom_text_repel(data = D[resp_ev %in% TRUE], aes(label = human), size = 2.3, fontface = "italic", segment.size = 0.2, box.padding = 0.4, max.overlaps = Inf) +
  scale_colour_manual(values = c("evaluated (symbol-matched), not responsive" = COL$liv, "evaluated and responsive" = COL$cand), name = NULL) +
  labs(x = expression(log[2]~fold~change~(T3~vs~vehicle)), y = expression(-log[10]~adjusted~italic(P)), tag = "C",
       caption = "Dashed lines: adjusted P = 0.05 and |log2 fold change| = 1 (response criteria)") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), plot.caption = element_text(size = 5.5, hjust = 0))
message(sprintf("(C) 검정 가능 %d개, 평가(이름 일치) %d개, 반응 %d개  (원고: 11,819 / 194 / 10)", nrow(D), sum(D$evaluated), sum(D$resp_ev)))
} else {
  message("  Mouse DESeq2 file (GSE85793_mouse_liver_T3_DESeq2_all_genes.csv) not found: panel C is replaced by a note")
  D <- data.table(gene = character(), evaluated = logical(), resp_ev = logical())
  pC <- ggplot() + annotate("text", x = 0, y = 0, label = "Requires the DESeq2 results of GEO GSE85793\n(GSE85793_mouse_liver_T3_DESeq2_all_genes.csv)", size = 2.4) + theme_void() + labs(tag = "C")
}

## ---- 저장 -----------------------------------------------------------------
save_fig(pA,  "Fig2A_colocalization_PPH4", 180, 66)
save_fig(pBc, "Fig2B_regional_association", 180, 100)
save_fig(pC,  "Fig2C_hepatic_response_volcano", 110, 90)
save_fig((pA / pBc / pC) + plot_layout(heights = c(1, 1.6, 1.35)), "Figure2", 180, 255)

## ---- 범례 txt (이 그림에 들어갈 범례. 숫자는 위에서 계산한 값) -----------------------
n_high <- sum(W$h4 >= 0.8, na.rm = TRUE); n_mod <- sum(W$h4 >= 0.5 & W$h4 < 0.8, na.rm = TRUE); n_nores <- sum(is.na(W$h4))
legA <- sprintf("(A) Cross-domain colocalization posterior probability (PP.H4) for each of the %d windows, taken as the highest value across the 6 trait-level pairs (3 thyroid traits \u00d7 2 liver traits; for TSH, the lower of the 2 datasets). Point colour shows the tier (%d high, PP.H4 \u2265 0.8; %d moderate, 0.5 to 0.8) and point shape shows the liver trait of the best pair; crosses mark the %d windows without a cross-domain result.", nrow(W), n_high, n_mod, n_nores)
pr <- if (all(have)) paste(sapply(names(REG), function(r) sprintf("%s, window %s, %s and %s", r, pairs[[r]]$win, GW_LABEL[[pairs[[r]]$thy]], GW_LABEL[[pairs[[r]]$liv]])), collapse = "; ") else ""
legB <- if (all(have)) sprintf("(B) Regional association plots of the 3 most pleiotropic regions, 9q34, 4q24 and 12q24. For the pair with the highest PP.H4 in each region (%s), the \u2212log10 P of the thyroid trait is plotted upward and that of the liver trait downward. Shading marks colocalization windows coloured by PP.H4, diamonds mark lead variants, dashed lines mark P = 1 \u00d7 10\u207b\u2076 and the lower track shows coding genes (shared-domain genes in blue, T3-responsive positional candidates in red). \u2212log10 P is capped at 50.", pr) else
  "(B) Genomic layout of the 3 most pleiotropic regions, 9q34, 4q24 and 12q24, showing colocalization windows coloured by PP.H4, lead variants (diamonds) and coding genes (shared-domain genes in blue, T3-responsive positional candidates in red)."
legC <- if (!HAVE_DEG) "(C) Volcano plot of the hepatic response to T3 in mouse liver (not drawn: the DESeq2 results file was not found)." else sprintf("(C) Volcano plot of the hepatic response to T3 across %s testable mouse genes. Dashed lines mark the response criteria (adjusted P < 0.05 and |log2 fold change| \u2265 1). Genes evaluated within supported windows by matching human and mouse gene symbols (%d) are blue and the %d meeting the response criteria are labelled.", nf(nrow(D)), sum(D$evaluated), sum(D$resp_ev))
write_legend("Figure2", paste("Figure 2. Shared thyroid-liver genetic signal and hormone responsiveness within supported windows.", legA, legB, legC, "Competitive gene-set results, prior sensitivity and full window layouts are given in Supplementary Figures S1 to S3."))


## ---- Record R and package versions used to draw the figures --------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
