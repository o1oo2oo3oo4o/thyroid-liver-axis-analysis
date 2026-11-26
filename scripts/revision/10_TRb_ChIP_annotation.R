#!/usr/bin/env Rscript
# =============================================================================
# 10_TRb_ChIP_annotation.R        (revision analysis; Reviewer 2 major comment 3)
#
# "Could the authors use published ChIP-seq data for THRB to confirm direct binding at the regulatory
#  regions of these candidate genes?"
#
# Approach: annotate the published human HepG2 TR-beta BioChIP-seq peaks (Ayers et al. 2014, PLoS ONE,
# doi:10.1371/journal.pone.0081186, Supplementary Table S5) over the 26 prioritized genes and the curated
# redox sets with one uniform rule (no re-alignment of raw reads):
#   promoter window : TSS -2 kb to +1 kb          proximal windows: +-10 kb and +-50 kb of the TSS
# The result is DESCRIPTIVE context in an engineered cell model: a peak nearby is not proof of causal
# regulation, and the absence of a nearby peak does not exclude distal or indirect regulation.
#
# Inputs (config.R)
#   trb_peak_file : peaks as CSV with chr,start,end (the build of the file is given in trb_peak_build)
#   tss_file      : gene,chr,tss,strand for the SAME build (strand "+" or "-"); ideally for ALL genes,
#                   which also gives a genome-wide background fraction for context.
#   How to build the TSS table (hg19 example, run once):
#     library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db)
#     g <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#     tss <- data.frame(gene = mapIds(org.Hs.eg.db, g$gene_id, "SYMBOL", "ENTREZID"), chr = as.character(seqnames(g)),
#                       tss = ifelse(strand(g) == "+", start(g), end(g)), strand = as.character(strand(g)))
#     write.csv(na.omit(tss), "data/tss_hg19.csv", row.names = FALSE)
#   If the peaks are in hg19 and the TSS table in hg38 (or the reverse), lift one of them first
#   (rtracklayer::liftOver); the two files must use the same build.
# Output: results/TRb_ChIP/*.csv
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
source(Sys.getenv("STAGE_CFG", "config.R"))
out_dir <- file.path(out_root, "10_TRb_ChIP"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

norm_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
peaks <- fread(trb_peak_file)
setnames(peaks, 1:3, c("chr", "start", "end")); peaks[, `:=`(chr = norm_chr(chr), start = as.integer(start), end = as.integer(end))]
tss <- fread(tss_file); tss[, `:=`(gene = trimws(gene), chr = norm_chr(chr), tss = as.integer(tss))]
tss <- tss[!duplicated(gene)]
sets <- fread(human_sets_file); sets[, human_symbol := trimws(human_symbol)]
cat(sprintf("peaks: %d | genes with a TSS: %d | build: %s\n", nrow(peaks), nrow(tss), trb_peak_build))

# strand-aware windows around the TSS
windows <- function(t, up, down) {                 # up = bp upstream, down = bp downstream of the TSS
  w <- copy(t)
  w[, `:=`(wstart = fifelse(strand == "+", tss - up, tss - down), wend = fifelse(strand == "+", tss + down, tss + up))]
  w[, .(gene, chr, wstart, wend)]
}
peak_count <- function(w) {
  setkey(peaks, chr, start, end)
  ov <- foverlaps(w[, .(gene, chr, start = wstart, end = wend)], peaks, by.x = c("chr", "start", "end"),
                  by.y = c("chr", "start", "end"), type = "any")
  # in the result the columns start/end belong to the peaks (y); they are NA when a window has no overlapping peak
  ov[, .(n_peaks = sum(!is.na(start))), by = gene]
}
win <- list(promoter_m2kb_p1kb = windows(tss, 2000, 1000), proximal_10kb = windows(tss, 10000, 10000),
            proximal_50kb = windows(tss, 50000, 50000))
res <- tss[, .(gene, chr, tss, strand)]
for (nm in names(win)) { pc <- peak_count(win[[nm]]); res <- merge(res, setnames(pc, "n_peaks", nm), by = "gene", all.x = TRUE) }
for (nm in names(win)) res[is.na(get(nm)), (nm) := 0L]
# nearest peak distance to the TSS
nearest <- function(g, ch, ts) { p <- peaks[chr == ch]; if (!nrow(p)) return(NA_integer_)
  d <- pmin(abs(p$start - ts), abs(p$end - ts)); d[ts >= p$start & ts <= p$end] <- 0L; min(d) }
res[, nearest_peak_distance_bp := mapply(nearest, gene, chr, tss)]
res[, any_promoter_peak := promoter_m2kb_p1kb > 0]

groups <- list(panel_26 = panel_genes,
               three_THRB_sensitive = c("GPCPD1", "ADH4", "F13B"))
for (s in unique(sets$set)) groups[[paste0("redox: ", s)]] <- sets[set == s, human_symbol]
out_tab <- rbindlist(lapply(names(groups), function(g) {
  r <- res[gene %in% groups[[g]]]
  data.table(group = g, genes_requested = length(groups[[g]]), genes_with_TSS = nrow(r),
             with_promoter_peak = sum(r$promoter_m2kb_p1kb > 0), with_peak_10kb = sum(r$proximal_10kb > 0),
             with_peak_50kb = sum(r$proximal_50kb > 0))
}))
bg <- res[, .(group = "all genes in the TSS table (background)", genes_requested = .N, genes_with_TSS = .N,
              with_promoter_peak = sum(promoter_m2kb_p1kb > 0), with_peak_10kb = sum(proximal_10kb > 0),
              with_peak_50kb = sum(proximal_50kb > 0))]
out_tab <- rbind(out_tab, bg)
out_tab[, `:=`(frac_promoter = round(with_promoter_peak / genes_with_TSS, 3), frac_10kb = round(with_peak_10kb / genes_with_TSS, 3),
               frac_50kb = round(with_peak_50kb / genes_with_TSS, 3))]
fwrite(res[gene %in% unique(c(panel_genes, sets$human_symbol))][order(gene)], file.path(out_dir, "1_gene_level_peak_annotation.csv"))
fwrite(out_tab, file.path(out_dir, "2_group_summary.csv"))
missing <- setdiff(unique(c(panel_genes, sets$human_symbol)), res$gene)
if (length(missing)) message("genes without a TSS in the table: ", paste(missing, collapse = ", "))
cat("\n[group summary] (descriptive; compare each group with the background row)\n"); print(out_tab)
cat("\nthe three THRB-sensitive candidates:\n")
print(res[gene %in% c("GPCPD1", "ADH4", "F13B"), .(gene, promoter_m2kb_p1kb, proximal_10kb, proximal_50kb, nearest_peak_distance_bp)])
cat("\nreport the numbers as they are; state that a nearby peak is context, not proof of regulation.\n")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
