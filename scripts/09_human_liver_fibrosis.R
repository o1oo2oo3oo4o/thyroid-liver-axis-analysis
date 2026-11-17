################################################################################
#  09_human_liver_fibrosis.R
#  (Original: 261021_human_liver_fibrosis_analysis.R)
#  Expression change with fibrosis stage in human liver tissue (GSE135251, GSE130970)
#
#  QUESTION
#    How do the glutathione set (33 human genes, down-regulated by T3 in mouse) and the prioritised candidate genes
#    change with fibrosis stage (F0-F4) in human NAFLD liver biopsies, and is the direction the same in the two cohorts?
#    (This is a trend test for gene sets and genes fixed in advance, not a two-group DEG comparison.)
#
#  INPUTS (DATA_DIR)
#    GSE130970_series_matrix.txt, GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv
#    GSE135251_series_matrix.txt, GSE135251_RAW.tar   (the tar file itself or an unpacked folder)
#    S9_gene_level_all_sets.csv  (output of 07_redox_cameraPR.R; the set definitions are taken from the rows with dataset == "human_WT")
#
#  ANALYSIS RULES (fixed before the results were seen)
#    1. Samples: GSE135251: 206 NAFLD samples (the 10 controls are added in a sensitivity analysis only);
#                GSE130970: all 78 samples (6 histologically normal samples are not identified by a single column and are not removed)
#    2. Outcome: fibrosis stage 0-4 as a number (per stage). Effect = log2 fold change per stage
#    3. Gene filter: genes with count >= 10 in at least half of the samples (the same rule as for the mouse DEG analysis)
#    4. Main model (DESeq2, Wald):
#         GSE130970 : ~ sex + age_z + fibrosis
#         GSE135251 : ~ fibrosis   (no age or sex in GEO; give CLINICAL_FILE_135251 if a clinical file is available)
#       cameraPR input statistic = Wald statistic of the fibrosis coefficient (not shrunken)
#    5. Sensitivity models: adjusted for NAS / adjusted for a hepatocyte marker score / F4 excluded / F3-4 vs F0-1 / (GSE135251) controls included
#    6. Set level: cameraPR (use.ranks = FALSE, inter-gene correlation 0.01; sensitivity 0, 0.05, 0.10),
#                  Bonferroni x6 over the six redox sets (glutathione is the primary set, the others are comparison sets)
#                  + a glutathione set without transferases (GST, MGST) as sensitivity analysis
#                  + relation between the set score (mean z-score of the member genes) and fibrosis (lm, Spearman)
#    7. Candidate genes: effect per cohort, BH within the candidate list, fixed-effect meta-analysis across cohorts
#    8. Between cohorts: number of genes with the same direction (binomial test against 50% and against the genome-wide agreement rate);
#       "replicated" = same direction in both cohorts and nominal P < 0.05 in each
#    9. Direction is tested two-sided. The expected direction was fixed before the results were seen and is reported together with them
#   10. Cell composition: results adjusted for the hepatocyte marker score are reported as well (marker based, not a deconvolution)
#
#  OUTPUTS (OUT_DIR; file names start with OUT_PREFIX)
#    _results.xlsx  : settings, samples, per-model summary, sets (cameraPR), set scores, candidates, meta-analysis, concordance
#    _gene_results_<cohort>_<model>.csv : all-gene results
#    _plots.pdf, _sessionInfo.txt
#  (Supplementary Tables S12 and S13 of the manuscript are built from these outputs.)
################################################################################
suppressPackageStartupMessages({
  library(DESeq2); library(limma); library(AnnotationDbi); library(org.Hs.eg.db)
  library(data.table); library(openxlsx); library(ggplot2)
})
options(stringsAsFactors = FALSE)
set.seed(1)

# ==============================================================================
# [Settings]
# ==============================================================================
DATA_DIR <- Sys.getenv("DATA_DIR", "data")
OUT_DIR  <- Sys.getenv("OUT_DIR", "results/09_human_liver_fibrosis")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
GSE130970_COUNTS <- file.path(DATA_DIR, "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv")
GSE130970_MATRIX <- file.path(DATA_DIR, "GSE130970_series_matrix.txt")
GSE135251_RAW    <- file.path(DATA_DIR, "GSE135251_RAW.tar")         # tar file or folder
GSE135251_MATRIX <- file.path(DATA_DIR, "GSE135251_series_matrix.txt")
REDOX_SET_FILE   <- Sys.getenv("REDOX_SET_FILE", "results/07_redox_cameraPR/signed/S9_gene_level_all_sets.csv")
REDOX_DATASET    <- "human_WT"                  # the set and gene columns of the rows of this dataset define the sets
PRIMARY_SET      <- "Glutathione-related"
CLINICAL_FILE_135251 <- NA                      # (optional) columns: gsm or title + age, sex, bmi, t2dm, ...
CLINICAL_COVARS_135251 <- c("age", "sex", "bmi", "t2dm")
OUT_PREFIX       <- file.path(OUT_DIR, "261021_human_liver_fibrosis")
CHECK_ONLY       <- identical(toupper(Sys.getenv("CHECK_ONLY", "FALSE")), "TRUE")   # TRUE: only read the files, match the samples and check the genes, then stop (no model is fitted)
MIN_COUNT <- 10; MIN_FRAC <- 0.5
RHO_MAIN  <- 0.01; RHO_SENS <- c(0, 0.05, 0.10)
ALPHA     <- 0.05
# 26-gene panel of the manuscript: 10 locus-based candidates (POS10) + 16 protein-based candidates (PROT16), ADH4 in both, + SLC39A8
# (the same definition as in Fig4_THRB_dependency_supp.R)
POS10   <- c("CUX2","BDH2","ADH4","LILRB3","SBK1","ABHD1","INCENP","CDC42EP5","CD19","GPCPD1")
PROT16  <- c("ADH4","BTC","CFHR1","CLIC5","COL11A2","COL15A1","DEF6","DNAJB4","F13B","LGALS4","PKLR","POR","SPARCL1","TIMD4","GSTM1","HSPA1A")
CANDIDATES <- sort(unique(c(POS10, PROT16, "SLC39A8")))
HEP_MARKERS <- c("ALB","APOA1","APOC3","TTR","SERPINA1","AHSG","FGB","APOH")

# If an input file does not have the configured name, look for a similar name (pattern) in the same folder. If none is found, list the folder and stop
find_input <- function(path, pattern) {
  if (!is.na(path) && file.exists(path)) return(path)
  d <- if (!is.na(path) && dirname(path) != "") dirname(path) else "."
  hits <- list.files(d, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(hits)) { cat(sprintf("   '%s' not found; using a file with a similar name: %s\n", path, hits[1])); return(hits[1]) }
  cat(sprintf("   '%s' not found. Files in the folder (%s):\n", path, normalizePath(d, mustWork = FALSE))); print(list.files(d))
  stop("Input file not found: ", path)
}
GSE130970_COUNTS <- find_input(GSE130970_COUNTS, "GSE130970.*counts.*entrez.*\\.csv(\\.gz)?$")
GSE130970_MATRIX <- find_input(GSE130970_MATRIX, "GSE130970.*series_matrix")
GSE135251_MATRIX <- find_input(GSE135251_MATRIX, "GSE135251.*series_matrix")
GSE135251_RAW    <- find_input(GSE135251_RAW,    "GSE135251_RAW")
REDOX_SET_FILE   <- find_input(REDOX_SET_FILE,   "S9_gene_level_all_sets.*\\.csv$")

# ==============================================================================
# [1] Read the inputs
# ==============================================================================
read_series_matrix <- function(f) {
  L <- readLines(f, warn = FALSE, encoding = "UTF-8")
  L <- L[startsWith(L, "!Sample_") & !startsWith(L, "!Sample_supplementary_file")]
  cut <- function(l) { p <- strsplit(l, "\t", fixed = TRUE)[[1]]; list(key = sub("^!", "", p[1]), val = gsub('^"|"$', "", p[-1])) }
  P <- lapply(L, cut)
  get1 <- function(k) P[[which(vapply(P, function(x) x$key == k, TRUE))[1]]]$val
  df <- data.frame(gsm = get1("Sample_geo_accession"), title = get1("Sample_title"), stringsAsFactors = FALSE)
  for (x in P[vapply(P, function(x) x$key == "Sample_characteristics_ch1", TRUE)]) {
    k <- trimws(sub(":.*$", "", x$val[1])); k <- gsub("[^a-z0-9]+", "_", tolower(k))
    df[[k]] <- trimws(sub("^[^:]*:", "", x$val))
  }
  df
}

match_samples <- function(cn, meta) {
  tries <- list(match(cn, meta$title), match(cn, meta$gsm), match(make.names(cn), make.names(meta$title)))
  for (idx in tries) if (mean(!is.na(idx)) >= 0.9) return(idx)
  stop("The column names of the count file do not match the title/GSM of the series matrix.\n  count column examples: ",
       paste(head(cn, 4), collapse = ", "), "\n  title examples: ", paste(head(meta$title, 4), collapse = ", "))
}

ids_to_symbol <- function(ids) {
  ids <- as.character(ids)
  type <- if (mean(grepl("^ENSG", ids)) > 0.8) "ENSEMBL" else if (mean(grepl("^[0-9]+$", ids)) > 0.8) "ENTREZID" else "SYMBOL"
  if (type == "ENSEMBL") ids <- sub("\\..*$", "", ids)
  if (type == "SYMBOL") return(list(sym = ids, type = type))
  m <- suppressMessages(mapIds(org.Hs.eg.db, keys = unique(ids), column = "SYMBOL", keytype = type, multiVals = "first"))
  list(sym = unname(m[ids]), type = type)
}

collapse_to_symbol <- function(mat, sym) {   # if several rows share a symbol, keep only the row with the largest total count
  keep <- !is.na(sym) & sym != ""
  mat <- mat[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(rowSums(mat), decreasing = TRUE); mat <- mat[o, , drop = FALSE]; sym <- sym[o]
  d <- !duplicated(sym); mat <- mat[d, , drop = FALSE]; rownames(mat) <- sym[d]
  mat
}

load_gse130970 <- function() {
  meta <- read_series_matrix(GSE130970_MATRIX)
  cnt <- fread(GSE130970_COUNTS, data.table = FALSE)
  ids <- as.character(cnt[[1]]); mat <- as.matrix(cnt[, -1, drop = FALSE]); storage.mode(mat) <- "double"
  idx <- match_samples(colnames(mat), meta)
  meta <- meta[idx, ]; colnames(mat) <- meta$gsm
  s <- ids_to_symbol(ids)
  cat(sprintf("GSE130970: gene ID type %s, %d samples, symbol conversion %d/%d rows\n", s$type, ncol(mat), sum(!is.na(s$sym)), length(ids)))
  list(mat = collapse_to_symbol(mat, s$sym), meta = meta)
}

load_gse135251 <- function() {
  meta <- read_series_matrix(GSE135251_MATRIX)
  dir <- GSE135251_RAW
  if (file.exists(dir) && !dir.exists(dir)) {
    td <- file.path(tempdir(), "GSE135251_RAW"); dir.create(td, showWarnings = FALSE); untar(dir, exdir = td); dir <- td
  }
  files <- list.files(dir, full.names = TRUE, recursive = TRUE)
  files <- files[grepl("GSM[0-9]+", basename(files)) & grepl("\\.(txt|tsv|csv)(\\.gz)?$", files, ignore.case = TRUE)]
  if (length(files) == 0) stop("GSE135251 count files not found: ", dir)
  gsm <- sub("^.*?(GSM[0-9]+).*$", "\\1", basename(files))
  rd <- function(f) {
    x <- tryCatch(fread(f, data.table = FALSE), error = function(e) read.delim(gzfile(f), header = FALSE))
    g <- as.character(x[[1]]); v <- suppressWarnings(as.numeric(x[[ncol(x)]]))
    ok <- !is.na(v) & !grepl("^(__|N_)", g); g <- g[ok]; v <- v[ok]
    if (anyDuplicated(g)) { v <- tapply(v, g, sum); g <- names(v); v <- as.numeric(v) }
    setNames(v, g)
  }
  vecs <- lapply(files, rd); names(vecs) <- gsm
  genes <- Reduce(intersect, lapply(vecs, names))
  lost <- 1 - length(genes) / length(unique(unlist(lapply(vecs, names))))
  if (lost > 0.01) warning(sprintf("Gene IDs differ between samples; %.1f%% of the genes were lost.", 100 * lost))
  mat <- do.call(cbind, lapply(vecs, function(v) v[genes]))
  mat <- mat[, intersect(meta$gsm, colnames(mat)), drop = FALSE]; meta <- meta[match(colnames(mat), meta$gsm), ]
  s <- ids_to_symbol(rownames(mat))
  cat(sprintf("GSE135251: gene ID type %s, %d samples, symbol conversion %d/%d rows\n", s$type, ncol(mat), sum(!is.na(s$sym)), nrow(mat)))
  list(mat = collapse_to_symbol(mat, s$sym), meta = meta)
}

# ==============================================================================
# [2] Sample information per cohort (rules 1, 2)
# ==============================================================================
as_num <- function(x) suppressWarnings(as.numeric(x))

make_coldata_130970 <- function(meta) {
  cd <- data.frame(row.names = meta$gsm, sample = meta$gsm, fibrosis = as_num(meta$fibrosis_stage),
                   nas = as_num(meta$nafld_activity_score), sex = factor(meta$sex), age = as_num(meta$age_at_biopsy))
  cd$age_z <- as.numeric(scale(cd$age)); cd
}

make_coldata_135251 <- function(meta) {
  cd <- data.frame(row.names = meta$gsm, sample = meta$gsm, fibrosis = as_num(meta$fibrosis_stage),
                   nas = as_num(meta$nas_score), disease = meta$disease, group = meta$group_in_paper)
  covars <- character()
  if (!is.na(CLINICAL_FILE_135251) && file.exists(CLINICAL_FILE_135251)) {
    cl <- fread(CLINICAL_FILE_135251, data.table = FALSE); names(cl) <- tolower(names(cl))
    key <- if ("gsm" %in% names(cl)) "gsm" else "title"; m <- match(if (key == "gsm") meta$gsm else meta$title, cl[[key]])
    for (v in intersect(CLINICAL_COVARS_135251, names(cl))) {
      x <- cl[[v]][m]
      if (is.numeric(x)) { cd[[v]] <- as.numeric(scale(x)); covars <- c(covars, v) } else { cd[[v]] <- factor(x); covars <- c(covars, v) }
    }
    cat("GSE135251 clinical covariates used:", paste(covars, collapse = ", "), "\n")
  } else cat("GSE135251: no clinical covariate file; the main model uses the fibrosis stage only.\n")
  list(cd = cd, covars = covars)
}

# ==============================================================================
# [3] DESeq2 (rules 3, 4, 5)
# ==============================================================================
fit_deseq <- function(mat, cd, terms, coef) {
  cd <- droplevels(cd[colnames(mat), , drop = FALSE])
  for (v in terms) if (is.factor(cd[[v]]) && nlevels(cd[[v]]) < 2) stop("Covariate ", v, " has only one level.")
  m <- round(mat); storage.mode(m) <- "integer"
  dds <- suppressMessages(DESeqDataSetFromMatrix(m, cd, as.formula(paste("~", paste(terms, collapse = " + ")))))
  keep <- rowSums(counts(dds) >= MIN_COUNT) >= ceiling(ncol(dds) * MIN_FRAC)
  dds <- suppressMessages(DESeq(dds[keep, ], quiet = TRUE))
  r <- results(dds, name = coef)
  out <- data.frame(gene = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
                    stat = r$stat, pvalue = r$pvalue, padj = r$padj, row.names = NULL)
  list(res = out, dds = dds, n_samples = ncol(dds), n_genes = nrow(dds))
}

hep_score <- function(dds) {
  v <- tryCatch(vst(dds, blind = FALSE), error = function(e) varianceStabilizingTransformation(dds, blind = FALSE))
  mk <- intersect(HEP_MARKERS, rownames(v))
  if (length(mk) < 3) { warning("Fewer than 3 hepatocyte markers; the hepatocyte adjustment is skipped."); return(NULL) }
  z <- t(scale(t(assay(v)[mk, , drop = FALSE])))
  as.numeric(scale(colMeans(z))) -> s; names(s) <- colnames(v); s
}

run_cohort <- function(name, mat, cd, covars, primary_subset, extra_subsets = list()) {
  prim_ids <- rownames(cd)[primary_subset]
  M <- list()
  base <- c(covars, "fibrosis")
  M$primary <- fit_deseq(mat[, prim_ids, drop = FALSE], cd, base, "fibrosis")
  hs <- hep_score(M$primary$dds)
  if (!is.null(hs)) cd[names(hs), "hep_z"] <- hs
  cdp <- cd[prim_ids, , drop = FALSE]
  M$nas_adjusted <- fit_deseq(mat[, prim_ids, drop = FALSE], cd, c(covars, "nas", "fibrosis"), "fibrosis")
  if (!is.null(hs)) M$hepatocyte_adjusted <- fit_deseq(mat[, prim_ids, drop = FALSE], cd, c(covars, "hep_z", "fibrosis"), "fibrosis")
  ex <- prim_ids[cdp$fibrosis <= 3]; M$excl_F4 <- fit_deseq(mat[, ex, drop = FALSE], cd, base, "fibrosis")
  ad <- prim_ids[cdp$fibrosis != 2]; cd$stage_group <- factor(ifelse(cd$fibrosis >= 3, "advanced", "early"), levels = c("early", "advanced"))
  M$advF3F4_vs_F0F1 <- fit_deseq(mat[, ad, drop = FALSE], cd, c(covars, "stage_group"), "stage_group_advanced_vs_early")
  for (nm in names(extra_subsets)) {
    ids <- rownames(cd)[extra_subsets[[nm]]]; M[[nm]] <- fit_deseq(mat[, ids, drop = FALSE], cd, base, "fibrosis")
  }
  list(models = M, cd = cd, prim_ids = prim_ids, covars = covars, vst = {
    v <- tryCatch(vst(M$primary$dds, blind = FALSE), error = function(e) varianceStabilizingTransformation(M$primary$dds, blind = FALSE)); assay(v) })
}

# ==============================================================================
# [4] Set level (rule 6)
# ==============================================================================
load_sets <- function() {
  g <- read.csv(REDOX_SET_FILE, check.names = FALSE)
  need <- c("dataset", "set", "gene"); if (!all(need %in% names(g))) stop("Required columns are missing in the set file: ", paste(names(g), collapse = ", "))
  g <- g[g$dataset == REDOX_DATASET, ]
  sets <- lapply(split(g$gene, g$set), unique)
  sets[["Glutathione-related (transferases removed)"]] <- setdiff(sets[[PRIMARY_SET]], grep("^M?GST", sets[[PRIMARY_SET]], value = TRUE))
  sets
}

run_camera <- function(res, sets, rho) {
  st <- setNames(res$stat, res$gene); st <- st[is.finite(st)]
  idx <- lapply(sets, function(g) which(names(st) %in% g)); idx <- idx[lengths(idx) >= 2]
  cp <- cameraPR(st, idx, inter.gene.cor = rho, use.ranks = FALSE, sort = FALSE)
  data.frame(set = rownames(cp), n_tested = cp$NGenes, direction = cp$Direction, camera_p = cp$PValue,
             mean_stat = vapply(idx, function(i) mean(st[i]), 0)[rownames(cp)], row.names = NULL)
}

set_score_trend <- function(vstm, cd, ids, covars, sets) {
  rows <- lapply(names(sets), function(s) {
    g <- intersect(sets[[s]], rownames(vstm)); if (length(g) < 2) return(NULL)
    z <- t(scale(t(vstm[g, ids, drop = FALSE]))); sc <- colMeans(z); d <- cbind(cd[ids, , drop = FALSE], score = sc)
    f <- lm(as.formula(paste("score ~", paste(c(covars, "fibrosis"), collapse = " + "))), data = d)
    co <- summary(f)$coefficients["fibrosis", ]; ct <- suppressWarnings(cor.test(d$score, d$fibrosis, method = "spearman"))
    data.frame(set = s, n_genes = length(g), n_samples = length(ids), slope_per_stage_SD = co[1], slope_se = co[2], lm_p = co[4],
               spearman_rho = unname(ct$estimate), spearman_p = ct$p.value)
  })
  do.call(rbind, rows)
}

# ==============================================================================
# [5] Candidate genes, meta-analysis, concordance between cohorts (rules 7, 8)
# ==============================================================================
meta_two <- function(b1, s1, b2, s2) {
  w1 <- 1 / s1^2; w2 <- 1 / s2^2; b <- (w1 * b1 + w2 * b2) / (w1 + w2); se <- sqrt(1 / (w1 + w2)); z <- b / se
  Q <- (b1 - b2)^2 / (s1^2 + s2^2)
  data.frame(meta_log2FC = b, meta_se = se, meta_p = 2 * pnorm(-abs(z)), Q = Q, het_p = pchisq(Q, 1, lower.tail = FALSE))
}

candidate_tables <- function(R1, R2, model, n1, n2) {
  a <- R1$models[[model]]$res; b <- R2$models[[model]]$res
  a <- a[a$gene %in% CANDIDATES, ]; b <- b[b$gene %in% CANDIDATES, ]
  a$padj_candidates <- p.adjust(a$pvalue, "BH"); b$padj_candidates <- p.adjust(b$pvalue, "BH")
  long <- rbind(cbind(cohort = n1, model = model, a), cbind(cohort = n2, model = model, b))
  m <- merge(a[, c("gene", "log2FC", "lfcSE", "pvalue")], b[, c("gene", "log2FC", "lfcSE", "pvalue")], by = "gene", suffixes = c("_1", "_2"))
  if (nrow(m)) {
    mt <- meta_two(m$log2FC_1, m$lfcSE_1, m$log2FC_2, m$lfcSE_2)
    m <- cbind(model = model, m, mt); m$same_direction <- sign(m$log2FC_1) == sign(m$log2FC_2)
    m$replicated <- m$same_direction & m$pvalue_1 < ALPHA & m$pvalue_2 < ALPHA
    m$meta_padj_candidates <- p.adjust(m$meta_p, "BH")
    names(m) <- sub("_1$", paste0("_", n1), names(m)); names(m) <- sub("_2$", paste0("_", n2), names(m))
  }
  list(long = long, meta = m)
}

concordance <- function(R1, R2, sets, model = "primary") {
  a <- R1$models[[model]]$res; b <- R2$models[[model]]$res
  m <- merge(a[, c("gene", "stat", "log2FC")], b[, c("gene", "stat", "log2FC")], by = "gene", suffixes = c("_1", "_2"))
  m <- m[is.finite(m$stat_1) & is.finite(m$stat_2), ]
  bg <- mean(sign(m$stat_1) == sign(m$stat_2))
  rows <- lapply(names(sets), function(s) {
    d <- m[m$gene %in% sets[[s]], ]; if (nrow(d) < 3) return(NULL)
    k <- sum(sign(d$stat_1) == sign(d$stat_2)); n <- nrow(d)
    data.frame(set = s, n_common = n, n_same_direction = k, p_vs_half = binom.test(k, n, 0.5)$p.value,
               background_rate = bg, p_vs_background = binom.test(k, n, bg)$p.value,
               spearman_stat = unname(cor(d$stat_1, d$stat_2, method = "spearman")))
  })
  list(sets = do.call(rbind, rows), genomewide_spearman = cor(m$stat_1, m$stat_2, method = "spearman"), n_common_genes = nrow(m))
}

# ==============================================================================
# [6] Run
# ==============================================================================
cat("== Reading data ==\n")
D130 <- load_gse130970(); D135 <- load_gse135251()
cd130 <- make_coldata_130970(D130$meta); cl135 <- make_coldata_135251(D135$meta); cd135 <- cl135$cd
cat(sprintf("GSE130970: %d samples | fibrosis distribution %s\n", nrow(cd130), paste(names(table(cd130$fibrosis)), table(cd130$fibrosis), sep = ":", collapse = " ")))
cat(sprintf("GSE135251: %d samples (NAFLD %d, control %d)\n", nrow(cd135), sum(cd135$disease == "NAFLD"), sum(cd135$disease != "NAFLD")))
sets <- load_sets()
cat("Sets:", paste(sprintf("%s(%d)", names(sets), lengths(sets)), collapse = "; "), "\n")

# ---- Check report: were genes and samples read correctly (checked before the models are fitted)
cat("\n== Check report ==\n")
cover <- function(nm, mat) {
  cat(sprintf("[%s] %d genes x %d samples | total count %.3g | mean library size per sample %.3g\n", nm, nrow(mat), ncol(mat), sum(mat), mean(colSums(mat))))
  keep <- rowSums(round(mat) >= MIN_COUNT) >= ceiling(ncol(mat) * MIN_FRAC); cat(sprintf("   genes passing the pre-filter (count >= %d in at least half of the samples): %d\n", MIN_COUNT, sum(keep)))
  for (s in names(sets)[!grepl("transferases removed", names(sets))]) cat(sprintf("   set %-58s present in data %2d/%2d | passing filter %2d\n", s, sum(sets[[s]] %in% rownames(mat)), length(sets[[s]]), sum(sets[[s]] %in% rownames(mat)[keep])))
  cat(sprintf("   candidate genes %d/%d present | passing filter %d | absent: %s\n", sum(CANDIDATES %in% rownames(mat)), length(CANDIDATES), sum(CANDIDATES %in% rownames(mat)[keep]), paste(setdiff(CANDIDATES, rownames(mat)), collapse = ", ")))
  cat(sprintf("   hepatocyte markers %d/%d present\n", sum(HEP_MARKERS %in% rownames(mat)), length(HEP_MARKERS)))
}
cover("GSE130970", D130$mat); cover("GSE135251", D135$mat)
cat(sprintf("Sample check: GSE130970 count columns %d = series matrix samples %d | GSE135251 count columns %d = series matrix samples %d\n", ncol(D130$mat), nrow(D130$meta), ncol(D135$mat), nrow(D135$meta)))
if (CHECK_ONLY) { cat("\nCHECK_ONLY=TRUE: stopping here. If there is no problem, run again without CHECK_ONLY to run the full analysis.\n"); quit(save = "no", status = 0) }

cat("\n== Fitting models (this takes a while) ==\n")
R130 <- run_cohort("GSE130970", D130$mat, cd130, c("sex", "age_z"), rep(TRUE, nrow(cd130)))
R135 <- run_cohort("GSE135251", D135$mat, cd135, cl135$covars, cd135$disease == "NAFLD",
                   extra_subsets = list(with_controls = rep(TRUE, nrow(cd135))))
COH <- list(GSE130970 = R130, GSE135251 = R135)

# Summary per model
msum <- do.call(rbind, lapply(names(COH), function(c) do.call(rbind, lapply(names(COH[[c]]$models), function(m) {
  x <- COH[[c]]$models[[m]]; data.frame(cohort = c, model = m, n_samples = x$n_samples, n_genes_tested = x$n_genes,
                                         n_padj_lt_0.05 = sum(x$res$padj < ALPHA, na.rm = TRUE)) }))))
# Sets (cameraPR)
cam <- do.call(rbind, lapply(names(COH), function(c) do.call(rbind, lapply(names(COH[[c]]$models), function(m) {
  r <- run_camera(COH[[c]]$models[[m]]$res, sets, RHO_MAIN); r$cohort <- c; r$model <- m
  main <- r$set %in% setdiff(names(sets), "Glutathione-related (transferases removed)")
  r$camera_bonf <- NA_real_; r$camera_bonf[main] <- pmin(1, r$camera_p[main] * sum(main)); r }))))
cam <- cam[, c("cohort", "model", "set", "n_tested", "direction", "camera_p", "camera_bonf", "mean_stat")]
camrho <- do.call(rbind, lapply(names(COH), function(c) do.call(rbind, lapply(RHO_SENS, function(rho) {
  r <- run_camera(COH[[c]]$models$primary$res, sets, rho); r$cohort <- c; r$rho <- rho; r[, c("cohort", "rho", "set", "n_tested", "direction", "camera_p")] }))))
# Set scores
sst <- do.call(rbind, lapply(names(COH), function(c) { x <- set_score_trend(COH[[c]]$vst, COH[[c]]$cd, COH[[c]]$prim_ids, COH[[c]]$covars, sets); x$cohort <- c; x }))
# Candidates, meta-analysis
cand_long <- list(); cand_meta <- list()
for (m in intersect(names(R130$models), names(R135$models))) {
  ct <- candidate_tables(R130, R135, m, "GSE130970", "GSE135251"); cand_long[[m]] <- ct$long; cand_meta[[m]] <- ct$meta
}
cand_long <- do.call(rbind, cand_long); cand_meta <- do.call(rbind, cand_meta)
conc <- concordance(R130, R135, sets)

# ==============================================================================
# [7] Save
# ==============================================================================
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", halign = "center", border = "Bottom")
wb <- createWorkbook()
put <- function(name, df) { if (is.null(df) || !nrow(df)) return(invisible()); addWorksheet(wb, name); writeData(wb, name, df, headerStyle = hs); freezePane(wb, name, firstActiveRow = 2) }
settings <- data.frame(item = c("run date and time","minimum count / sample fraction","cameraPR correlation (main)","cameraPR correlation (sensitivity)","main model GSE130970","main model GSE135251","primary set","number of candidate genes","hepatocyte markers"),
  value = c(format(Sys.time(), "%Y-%m-%d %H:%M"), sprintf("%d / %.2f", MIN_COUNT, MIN_FRAC), RHO_MAIN, paste(RHO_SENS, collapse = ", "),
            paste("~", paste(c(R130$covars, "fibrosis"), collapse = " + ")), paste("~", paste(c(R135$covars, "fibrosis"), collapse = " + ")),
            PRIMARY_SET, length(CANDIDATES), paste(HEP_MARKERS, collapse = ", ")))
put("settings", settings); put("model_summary", msum)
put("sets_cameraPR", cam); put("sets_cameraPR_rho", camrho); put("sets_score_trend", sst)
put("candidates_by_cohort", cand_long); put("candidates_meta", cand_meta)
put("concordance_sets", conc$sets)
put("concordance_overall", data.frame(n_common_genes = conc$n_common_genes, genomewide_spearman_of_stat = conc$genomewide_spearman))
put("samples_GSE130970", cbind(cd130, hep_z = R130$cd$hep_z)); put("samples_GSE135251", cbind(cd135, hep_z = R135$cd$hep_z))
saveWorkbook(wb, paste0(OUT_PREFIX, "_results.xlsx"), overwrite = TRUE)
for (c in names(COH)) for (m in names(COH[[c]]$models)) fwrite(COH[[c]]$models[[m]]$res, sprintf("%s_gene_results_%s_%s.csv", OUT_PREFIX, c, m))

# Figures: set score vs fibrosis, candidate forest plot
pdf(paste0(OUT_PREFIX, "_plots.pdf"), width = 11, height = 7)
main_sets <- setdiff(names(sets), "Glutathione-related (transferases removed)")
pl <- do.call(rbind, lapply(names(COH), function(c) { v <- COH[[c]]$vst; ids <- COH[[c]]$prim_ids
  do.call(rbind, lapply(main_sets, function(s) { g <- intersect(sets[[s]], rownames(v)); if (length(g) < 2) return(NULL)
    data.frame(cohort = c, set = s, fibrosis = COH[[c]]$cd[ids, "fibrosis"], score = colMeans(t(scale(t(v[g, ids, drop = FALSE]))))) })) }))
print(ggplot(pl, aes(factor(fibrosis), score)) + geom_boxplot(outlier.shape = NA, fill = "grey90") + geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
        facet_grid(cohort ~ set, scales = "free_y", labeller = label_wrap_gen(18)) + labs(x = "Fibrosis stage", y = "Set score (mean z of vst)") + theme_bw(base_size = 9))
if (!is.null(cand_meta) && nrow(cand_meta[cand_meta$model == "primary", ])) {
  cm <- cand_meta[cand_meta$model == "primary", ]
  fd <- rbind(data.frame(gene = cm$gene, cohort = "GSE130970", b = cm$log2FC_GSE130970, se = cm$lfcSE_GSE130970),
              data.frame(gene = cm$gene, cohort = "GSE135251", b = cm$log2FC_GSE135251, se = cm$lfcSE_GSE135251),
              data.frame(gene = cm$gene, cohort = "meta", b = cm$meta_log2FC, se = cm$meta_se))
  print(ggplot(fd, aes(b, reorder(gene, b), colour = cohort)) + geom_vline(xintercept = 0, linetype = 2) +
          geom_pointrange(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), position = position_dodge(width = 0.6), size = 0.3) +
          labs(x = "log2 fold change per fibrosis stage (primary model)", y = NULL) + theme_bw(base_size = 9))
}
dev.off()
writeLines(capture.output(sessionInfo()), paste0(OUT_PREFIX, "_sessionInfo.txt"))

# Console summary
cat("\n================ Summary ================\n")
print(cam[cam$model == "primary" & cam$set == PRIMARY_SET, ], row.names = FALSE)
print(cam[cam$model == "hepatocyte_adjusted" & cam$set == PRIMARY_SET, ], row.names = FALSE)
cat("\n[Direction agreement between cohorts]\n"); print(conc$sets, row.names = FALSE)
cat("\n[Candidates replicated in both cohorts (main model)]\n"); print(cand_meta[cand_meta$model == "primary" & cand_meta$replicated, c("gene", "meta_log2FC", "meta_p")], row.names = FALSE)
cat("\nSaved:", paste0(OUT_PREFIX, "_results.xlsx"), "\nDone.\n")
