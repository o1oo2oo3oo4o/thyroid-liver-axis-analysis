## =====================================================================
## Supplementary Figure S2 (drawn from the Supplementary Tables). Sensitivity to the shared-signal prior p12 (5e-6 to 1e-4)
##   Does the highest cross-domain PP.H4 of the 27 supported windows hold when p12 is changed?
##
## The coloc.abf posterior is proportional to prior x likelihood, so with p1 and p2 fixed and only p12 changed,
## only H4 scales with p12:
##     PP.H4(new) = PP.H4 * r / (1 - PP.H4 + PP.H4 * r),   r = new p12 / baseline p12
## Baseline p12 = 1e-5 (p1 = p2 = 1e-4). The values for each p12 are recomputed exactly from the baseline posteriors
## of the dataset pairs (sheet S2b), and the decision rule of the main text is applied (liver trait: maximum over the
## 4 adjustment models; TSH: lower of the 2 datasets).
##
## Inputs (Supplementary Tables workbook): S2a (window summary and tier), S2b (colocalization of all dataset pairs)
## Output: figures_out_supp_S2/ (FigureS2 as .pdf and .png, FigureS2_legend.txt, sessionInfo.txt)
## Usage : Rscript FigS2_prior_sensitivity_supp.R
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
OUT_DIR <- "figures_out_supp_S2"; dir.create(OUT_DIR, showWarnings = FALSE)
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

## ---------------------------------------------------------------------------------------------
## [5] 그림 만들기
## ---------------------------------------------------------------------------------------------
co <- load_coloc(); S1 <- co$S1; S2 <- co$S2
BASE_P12 <- 1e-5
P12S <- c(5e-6, 1e-5, 5e-5, 1e-4)
res <- rbindlist(lapply(P12S, function(p12) {
  d <- copy(S2); r <- p12 / BASE_P12
  d[, PP.H4 := PP.H4 * r / (1 - PP.H4 + PP.H4 * r)]
  b <- window_best(d)$best; b[, .(locus_id, h4, p12 = p12)]
}))
sup <- S1[Final_tier %in% c("high", "moderate"), .(locus_id = Locus_ID)]
stopifnot(nrow(sup) == 27)
res <- merge(res, sup, by = "locus_id")
base <- res[abs(p12 - BASE_P12) < 1e-12, .(locus_id, base_h4 = h4)]
res <- merge(res, base, by = "locus_id")
res[, grp := factor(ifelse(base_h4 >= 0.8, "high at p12 = 1e-5", "moderate at p12 = 1e-5"), levels = c("high at p12 = 1e-5", "moderate at p12 = 1e-5"))]
## 요약 출력 (본문의 "사전확률을 바꿔도 안정적" 문장과 대조용)
for (p in P12S) message(sprintf("  p12 = %-6g : PP.H4 >= 0.5 인 창 %2d / 27,  >= 0.8 인 창 %2d / 27", p, sum(res[p12 == p]$h4 >= 0.5), sum(res[p12 == p]$h4 >= 0.8)))
stay <- res[, .(ok = all(h4 >= 0.5)), by = locus_id]
message(sprintf("  모든 p12 에서 PP.H4 >= 0.5 를 유지하는 창: %d / 27", sum(stay$ok)))
message("  p12 = 5e-6 에서 0.5 아래로 내려가는 창: ", paste(res[p12 == 5e-6 & h4 < 0.5]$locus_id, collapse = ", "))
p <- ggplot(res, aes(p12, h4, group = locus_id, colour = grp)) +
  geom_hline(yintercept = c(0.5, 0.8), linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_line(linewidth = 0.4, alpha = 0.8) + geom_point(size = 1) +
  scale_x_log10(breaks = P12S, labels = c("5e-6", "1e-5", "5e-5", "1e-4")) +
  scale_colour_manual(values = c(COL$high, COL$moderate), name = NULL) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Prior probability of a shared association (p12)", y = "Cross-domain PP.H4 (per window)") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))
save_fig(p, "FigS2_prior_sensitivity", 110, 85)

## ---- 범례 txt (이 그림에 들어갈 범례. 숫자는 위에서 계산한 값) -----------------------
cnt <- sapply(P12S, function(p) c(ge05 = sum(res[p12 == p]$h4 >= 0.5), ge08 = sum(res[p12 == p]$h4 >= 0.8)))
colnames(cnt) <- format(P12S); n_win <- uniqueN(res$locus_id); n_stay <- sum(stay$ok)
P12LAB <- c("5 \u00d7 10\u207b\u2076", "1 \u00d7 10\u207b\u2075", "5 \u00d7 10\u207b\u2075", "1 \u00d7 10\u207b\u2074")
txt_p <- paste(sprintf("p12 = %s: %d windows with PP.H4 \u2265 0.5 and %d with PP.H4 \u2265 0.8", P12LAB, cnt["ge05", ], cnt["ge08", ]), collapse = "; ")
leg <- sprintf("Supplementary Figure S2. Sensitivity of the cross-domain colocalization signal to the prior probability of a shared causal variant. For each of the %d supported windows, the highest cross-domain PP.H4 is shown at p12 = 5 \u00d7 10\u207b\u2076, 1 \u00d7 10\u207b\u2075 (baseline), 5 \u00d7 10\u207b\u2075 and 1 \u00d7 10\u207b\u2074, with p1 = p2 = 1 \u00d7 10\u207b\u2074 held fixed; values were recalculated from the baseline posterior probabilities. Lines are coloured by the tier at baseline and dashed lines mark PP.H4 = 0.5 and 0.8. %s. %d of %d windows kept PP.H4 \u2265 0.5 at every prior.", n_win, txt_p, n_stay, n_win)
write_legend("FigureS2", leg)


## ---- Record R and package versions used to draw the figures --------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
