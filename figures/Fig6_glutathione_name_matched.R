################################################################################
#  Fig6_glutathione_name_matched.R
#
#  Figure 6 (direction agreement of the glutathione set between mouse and human), drawn for the name-matched genes (20).
#  Standalone script: no other R file is needed.
#
#  Inputs
#    GENES_FILE    : 11b_concordance_genes_name_matched.csv   <- output of revision/11_glutathione_name_matched_concordance.R
#    SUMMARY_FILE  : 11d_concordance_summary.csv               <- same script (numbers printed for the legend)
#    HUMAN_WT_FILE : DEG_WT_T3_effect.csv                      (DESeq2 wild-type T3 effect, GSE270233; for the "all other genes" distribution in panel B)
#
#  Panels
#    A  mouse log2FC against human log2FC (point size = human -log10 adjusted P; blue = same direction, red = opposite direction)
#    B  human log2FC distribution: all other tested genes (violin + box) against the glutathione-set members (points + box)
#    C  human log2FC with 95% confidence interval (log2FC +/- 1.96 x lfcSE), ordered by log2FC; * adjusted P < 0.05, ** adjusted P < 0.01
#       (the interval is the unadjusted DESeq2 interval; the asterisks refer to the adjusted P)
#
#  Output (OUT_DIR)
#    Figure6_name_matched_glutathione.png (300 dpi) and Figure6_name_matched_glutathione.pdf; the legend and the numbers for the
#    text are printed to the console.
#
#  Usage:  Rscript figures/Fig6_glutathione_name_matched.R   (run from the repository root; every path can be set with the
#          environment variables RES_DIR, GENES_FILE, SUMMARY_FILE, HUMAN_WT_FILE, OUT_DIR)
#  R packages: ggplot2, dplyr, readr, ggrepel, patchwork; optional: ragg
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
  library(patchwork)
})

# ==============================================================================
# Settings
# ==============================================================================
RES_DIR       <- Sys.getenv("RES_DIR", "results/11_glutathione_name_matched")
GENES_FILE    <- Sys.getenv("GENES_FILE",   file.path(RES_DIR, "11b_concordance_genes_name_matched.csv"))
SUMMARY_FILE  <- Sys.getenv("SUMMARY_FILE", file.path(RES_DIR, "11d_concordance_summary.csv"))
HUMAN_WT_FILE <- Sys.getenv("HUMAN_WT_FILE", "data/DEG_WT_T3_effect.csv")
OUT_DIR       <- Sys.getenv("OUT_DIR", "figures_out_6")
FONT          <- Sys.getenv("FIG_FONT", "Arial")     # use another font name if Arial is not installed

COL_CONC <- "#4269A8"      # same direction
COL_DISC <- "#B23434"      # opposite direction
Y_LIM_B  <- c(-3.3, 3.3)   # y range of panel B (the violin is cut outside this range)
FIG_W <- 10.5; FIG_H <- 12 # inch

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(1)   # ggrepel places the gene labels of panel A with a random start; the seed makes the layout identical in every run

# Device used to measure text sizes. The same font engine as for saving (ragg) avoids font warnings and gives the exact layout.
USE_RAGG <- requireNamespace("ragg", quietly = TRUE)
if (USE_RAGG) ragg::agg_png(tempfile(fileext = ".png"), width = FIG_W, height = FIG_H, units = "in", res = 300)

# ==============================================================================
# Inputs
# ==============================================================================
for (f in c(GENES_FILE, SUMMARY_FILE, HUMAN_WT_FILE)) if (!file.exists(f)) stop("Input file not found: ", f)
g <- read_csv(GENES_FILE, show_col_types = FALSE)
sm <- read_csv(SUMMARY_FILE, show_col_types = FALSE) %>% filter(set == "name_matched")
h <- read_csv(HUMAN_WT_FILE, show_col_types = FALSE) %>% filter(!is.na(padj)) %>% distinct(gene, .keep_all = TRUE)
stopifnot(nrow(g) == sm$n_genes, nrow(sm) == 1)

g <- g %>%
  mutate(direction = factor(ifelse(concordant, "Concordant", "Discordant"), levels = c("Concordant", "Discordant")),
         nlogp = -log10(human_padj),
         lo = human_lfc - 1.96 * human_lfcSE, hi = human_lfc + 1.96 * human_lfcSE,
         star = case_when(human_padj < 0.01 ~ "**", human_padj < 0.05 ~ "*", TRUE ~ ""))
n_set <- nrow(g)
others <- h %>% filter(!(gene %in% g$KEY))
n_other <- nrow(others)
cols <- c(Concordant = COL_CONC, Discordant = COL_DISC)

# ==============================================================================
# Common theme
# ==============================================================================
th <- theme_classic(base_size = 12, base_family = FONT) +
  theme(axis.text = element_text(colour = "black", size = 11),
        axis.title = element_text(colour = "black", size = 12.5),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.tag = element_text(size = 22, face = "bold", family = FONT),
        plot.tag.position = c(0.005, 0.995),
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", fill = "white", linewidth = 0.6),
        legend.box.margin = margin(4, 8, 4, 8))

ylab_h <- expression("Human log"[2]*" FoldChange")

# Extracts the legend only (works with ggplot2 3.4 and with 3.5 or later)
extract_legend <- function(p) {
  gt <- ggplotGrob(p)
  i  <- grep("^guide-box", gt$layout$name)
  i  <- i[!vapply(gt$grobs[i], function(x) inherits(x, "zeroGrob"), logical(1))][1]
  if (is.na(i)) stop("legend not found")
  gt$grobs[[i]]
}

# ==============================================================================
# Panel A: mouse versus human log2FC
# ==============================================================================
pA_base <- ggplot(g, aes(mouse_lfc, human_lfc)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_point(aes(colour = direction, size = nlogp)) +
  geom_text_repel(aes(label = KEY, colour = direction), family = FONT, fontface = "bold.italic", size = 3.6,
                  min.segment.length = 0.15, segment.size = 0.3, box.padding = 0.55, point.padding = 0.45,
                  force = 2, max.overlaps = Inf, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_size_continuous(range = c(1.5, 6.5), breaks = c(0.5, 1, 1.5, 2),
                        name = expression("Human -log"[10]*" adjusted "*italic(P))) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.10))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(x = expression("Mouse log"[2]*" FoldChange"), y = ylab_h, tag = "A") +
  guides(colour = guide_legend(order = 1, nrow = 2, override.aes = list(size = 5)),
         size = guide_legend(order = 2, nrow = 1, title.position = "top", override.aes = list(colour = "black"))) +
  th
legA <- extract_legend(pA_base + theme(legend.position = "bottom", legend.box = "horizontal",
                                            legend.key = element_blank(), legend.key.width = unit(14, "pt")))
pA <- pA_base + theme(legend.position = "none")

# ==============================================================================
# Panel B: distribution of the human log2FC
# ==============================================================================
lab_other <- sprintf("All other genes\n(n = %s)", format(n_other, big.mark = ","))
lab_set   <- sprintf("Glutathione-set members\n(n = %d)", n_set)
dB_other <- tibble(x = lab_other, y = others$log2FoldChange)
dB_set   <- tibble(x = lab_set,   y = g$human_lfc)
pB <- ggplot() +
  geom_violin(data = dB_other, aes(x, y), fill = "grey80", colour = NA, width = 0.95, scale = "width") +
  geom_boxplot(data = dB_other, aes(x, y), width = 0.20, outlier.shape = NA, fill = "white", colour = "black",
               linewidth = 0.5) +
  geom_boxplot(data = dB_set, aes(x, y), width = 0.28, outlier.shape = NA, fill = "white", colour = "black",
               linewidth = 0.5) +
  geom_point(data = dB_set, aes(x, y), position = position_jitter(width = 0.12, height = 0, seed = 1),
             colour = COL_CONC, size = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.4) +
  scale_x_discrete(limits = c(lab_other, lab_set)) +
  coord_cartesian(ylim = Y_LIM_B) +
  labs(x = NULL, y = ylab_h, tag = "B") +
  th + theme(axis.text.x = element_text(size = 11, lineheight = 1.0))

# ==============================================================================
# Panel C: human log2FC with 95% CI
# ==============================================================================
g_c <- g %>% arrange(human_lfc) %>% mutate(KEY = factor(KEY, levels = KEY))
x_pad <- 0.025 * diff(range(c(g_c$lo, g_c$hi)))
pC_base <- ggplot(g_c, aes(y = KEY, colour = direction)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.7) +
  geom_point(aes(x = human_lfc), size = 2.8) +
  geom_text(aes(x = hi + x_pad, label = star), family = FONT, size = 5, vjust = 0.75, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = expression("Human log"[2]*" FoldChange (95% CI)"), y = NULL, tag = "C",
       caption = "* adjusted P < 0.05, ** adjusted P < 0.01 (DESeq2, Benjamini-Hochberg). Intervals are unadjusted 95% confidence intervals (log2FC \u00b1 1.96 \u00d7 lfcSE).") +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.8, size = 3.5))) +
  th + theme(axis.text.y = element_text(face = "italic", size = 12.5),
             plot.caption = element_text(size = 10, hjust = 0, colour = "black", margin = margin(t = 6)))
legC <- extract_legend(pC_base + theme(legend.position = "bottom", legend.key = element_blank(),
                                            legend.key.width = unit(26, "pt"), legend.text = element_text(size = 12)))
pC <- pC_base + theme(legend.position = "none")

# ==============================================================================
# Assembly and saving
# ==============================================================================
top   <- (pA | pB) + plot_layout(widths = c(1, 0.95))
leg_r <- (wrap_elements(full = legA) | plot_spacer()) + plot_layout(widths = c(1, 0.95))
fig <- top / leg_r / pC / wrap_elements(full = legC) +
  plot_layout(heights = c(1, 0.17, 1.2, 0.12))

png_dev <- if (USE_RAGG) ragg::agg_png else "png"
f_png <- file.path(OUT_DIR, "Figure6_name_matched_glutathione.png")
f_pdf <- file.path(OUT_DIR, "Figure6_name_matched_glutathione.pdf")
ggsave(f_png, fig, device = png_dev, width = FIG_W, height = FIG_H, units = "in", dpi = 300, bg = "white")
ggsave(f_pdf, fig, device = cairo_pdf, width = FIG_W, height = FIG_H, units = "in", bg = "white")
if (USE_RAGG) invisible(dev.off())
cat("Saved:", f_png, "\n     ", f_pdf, "\n")

# ==============================================================================
# Numbers for the legend and the text
# ==============================================================================
cat("\n---------- numbers for legend and text (11d_concordance_summary.csv, name_matched) ----------\n")
cat(sprintf("set genes %d, same direction %d (one-sided binomial test against 0.5: P = %.3f; against the background rate of %.1f%%: P = %.2f)\n",
            sm$n_genes, sm$n_concordant, sm$p_vs_half, 100 * sm$background_rate, sm$p_vs_background))
cat(sprintf("human log2FC median: set %.2f, all other genes %.2f (n = %s), one-sided Mann-Whitney P = %.3f\n",
            sm$median_human_lfc, sm$median_human_lfc_others, format(sm$n_other_genes, big.mark = ","), sm$wilcox_one_sided_p))
cat(sprintf("median |log2FC|: mouse %.2f, human %.2f\n", sm$median_abs_mouse, sm$median_abs_human))
cat("genes with human adjusted P < 0.05:", paste(as.character(g_c$KEY[g_c$human_padj < 0.05]), collapse = ", "), "\n")
cat("genes not in the same direction:", paste(as.character(g_c$KEY[!g_c$concordant]), collapse = ", "), "\n")
