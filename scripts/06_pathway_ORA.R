################################################################################
#  06_pathway_ORA.R
#  Over-representation analysis (ORA) of the human counterparts of the T3-responsive mouse liver genes,
#  run against a background of the genes that were actually tested ("tested" background, main analysis)
#  or against the whole annotated genome ("genome" background, comparison).
#  (Original: 261021_pathway_ORA_tested_background.R, a generalised version of 261020_background_pathway.R, which produced the main analysis)
#
#  WHAT THE TWO BACKGROUNDS MEAN
#    genome : no universe is given -> every gene annotated in each ontology is the background (clusterProfiler default)
#    tested : universe = human counterparts of the mouse genes for which DESeq2 returned an adjusted P value
#             (11,819 mouse genes -> 10,085 human genes for the 448-gene analysis).
#
#  "Tested genes" are the genes that had a chance to become responsive (adjusted P < 0.05, |log2FC| >= 1).
#  Genes removed by the pre-filter (count >= 10 in fewer than half of the samples) or with padj NA (independent
#  filtering, Cook's outliers) could not be judged and are therefore excluded from the background. The human background is built
#  by 04_reortholog_synonym.R with the same mapping rules as the foreground, and the foreground must always be a subset of the background.
#
#  INPUTS (DATA_DIR)
#    GENE_FILE     : foreground, one human symbol per line (liver_T3_human_responsive_genes_448.txt)
#    UNIVERSE_FILE : background, one human symbol per line (human_universe_all.txt = universe_tested_10085.txt)
#    COMPARE_FILE  : (optional) workbook of the genome-background run; if present the change per term is added as sheets
#
#  PACKAGES: clusterProfiler, org.Hs.eg.db, ReactomePA, msigdbr, openxlsx, dplyr
#    KEGG needs an internet connection (enrichKEGG downloads from the KEGG server).
#
#  OUTPUT (OUT_XLSX): Used_genes, GO_BP, GO_CC, GO_MF, KEGG, Reactome, Hallmark
#      + Compare_summary, Key_terms, Compare_changed when COMPARE_FILE exists
#      + a *_sessionInfo.txt file next to the workbook
#
#  USAGE   Rscript 06_pathway_ORA.R               # main analysis: 448 genes, tested background
#          UNIVERSE_MODE=genome Rscript 06_pathway_ORA.R   # genome background (comparison)
################################################################################
suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
})

# ==============================================================================
# [Settings]
# ==============================================================================
DATA_DIR      <- Sys.getenv("DATA_DIR", "data")
OUT_DIR       <- Sys.getenv("OUT_DIR", "results/06_pathway_ORA")
GENE_FILE     <- file.path(DATA_DIR, "liver_T3_human_responsive_genes_448.txt")
UNIVERSE_FILE <- file.path(DATA_DIR, "human_universe_all.txt")
UNIVERSE_MODE <- Sys.getenv("UNIVERSE_MODE", "tested")   # "tested": background of tested genes / "genome": whole-genome background (comparison)
COMPARE_FILE  <- file.path(DATA_DIR, "261020_3_DEG-448_pathway_results.xlsx")   # NA if absent. Run first with UNIVERSE_MODE=genome to create this file, then with "tested"
OUT_XLSX      <- NA                # NA = automatic name according to the mode
P_CUT <- 1; Q_CUT <- 1             # keep all terms, as in the original workbooks (filter on p.adjust afterwards)
MIN_GS <- 10; MAX_GS <- 500        # clusterProfiler defaults
ALPHA <- 0.05

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(UNIVERSE_MODE %in% c("tested", "genome"))

suppressPackageStartupMessages({ library(clusterProfiler); library(ReactomePA); library(msigdbr) })
suppressPackageStartupMessages(library(org.Hs.eg.db))

# ==============================================================================
# [1] Symbol -> Entrez conversion (same rules for foreground and background)
# ==============================================================================
read_symbols <- function(f) { x <- trimws(readLines(f, warn = FALSE, encoding = "UTF-8")); unique(x[nzchar(x)]) }
map_symbols <- function(sym) {
  m1 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = sym, keytype = "SYMBOL", columns = "ENTREZID"))
  m1 <- m1[!is.na(m1$ENTREZID), , drop = FALSE]; m1 <- m1[!duplicated(m1$SYMBOL), , drop = FALSE]
  miss <- setdiff(sym, m1$SYMBOL)                      # symbols not found as current symbols are tried once more as aliases (ALIAS)
  if (length(miss)) {
    m2 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = miss, keytype = "ALIAS", columns = "ENTREZID"))
    m2 <- m2[!is.na(m2$ENTREZID), , drop = FALSE]; m2 <- m2[!duplicated(m2$ALIAS), , drop = FALSE]
    if (nrow(m2)) { names(m2)[names(m2) == "ALIAS"] <- "SYMBOL"; m1 <- rbind(m1[, c("SYMBOL", "ENTREZID")], m2[, c("SYMBOL", "ENTREZID")]) }
  }
  out <- data.frame(Input_gene = sym, stringsAsFactors = FALSE)
  out$ENTREZID <- m1$ENTREZID[match(sym, m1$SYMBOL)]
  out$Canonical_SYMBOL <- NA_character_
  ok <- !is.na(out$ENTREZID)
  if (any(ok)) out$Canonical_SYMBOL[ok] <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = out$ENTREZID[ok], column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
  out$Mapped <- ok
  out[, c("Input_gene", "Canonical_SYMBOL", "ENTREZID", "Mapped")]
}

# If the input list is not found under the configured name (e.g. a different gene count), use the most recent file of the same kind in that folder
pick_latest <- function(path, pattern) {
  if (file.exists(path)) return(path)
  hits <- list.files(dirname(path), pattern = pattern, full.names = TRUE)
  if (!length(hits)) stop("Input file not found: ", path)
  hit <- hits[which.max(file.mtime(hits))]; cat(sprintf("   '%s' not found; using the most recent file: %s\n", path, hit)); hit
}
GENE_FILE <- pick_latest(GENE_FILE, "^liver_T3_human_responsive_genes_[0-9]+\\.txt$")
if (UNIVERSE_MODE == "tested") UNIVERSE_FILE <- pick_latest(UNIVERSE_FILE, "^liver_T3_human_universe_tested_[0-9]+\\.txt$")
fg_sym <- read_symbols(GENE_FILE); fg_map <- map_symbols(fg_sym)
gene_entrez <- unique(fg_map$ENTREZID[fg_map$Mapped])
if (is.na(OUT_XLSX)) OUT_XLSX <- file.path(OUT_DIR, sprintf("261021_3_DEG-%d_pathway_results_%s_background.xlsx", length(fg_sym), UNIVERSE_MODE))
cat(sprintf("Foreground: %d input genes -> %d Entrez IDs (%d failed to map)\n", length(fg_sym), length(gene_entrez), sum(!fg_map$Mapped)))

universe_entrez <- NULL; uni_info <- list(n_input = NA, n_entrez = NA, n_added = 0L)
if (UNIVERSE_MODE == "tested") {
  uni_sym <- read_symbols(UNIVERSE_FILE); uni_map <- map_symbols(uni_sym)
  universe_entrez <- unique(uni_map$ENTREZID[uni_map$Mapped])
  add <- setdiff(gene_entrez, universe_entrez)           # the foreground must be contained in the background
  if (length(add)) { warning(length(add), " foreground genes were not in the background and were added to it: ", paste(fg_map$Input_gene[fg_map$ENTREZID %in% add], collapse = ", ")); universe_entrez <- union(universe_entrez, add) }
  uni_info <- list(n_input = length(uni_sym), n_entrez = length(universe_entrez), n_added = length(add))
  cat(sprintf("Background: %d input genes -> %d Entrez IDs (%d added from the foreground)\n", uni_info$n_input, uni_info$n_entrez, uni_info$n_added))
}

# ==============================================================================
# [2] Pathway analysis
# ==============================================================================
to_df <- function(x) { if (is.null(x)) return(data.frame()); d <- as.data.frame(x); rownames(d) <- NULL; d[order(d$p.adjust), , drop = FALSE] }
notes <- character()
run_enrich <- function() {
  res <- list()
  for (ont in c("BP", "CC", "MF")) {
    cat("  GO", ont, "...\n")
    res[[paste0("GO_", ont)]] <- to_df(enrichGO(gene = gene_entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
                                                  pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, minGSSize = MIN_GS, maxGSSize = MAX_GS, readable = TRUE))
  }
  cat("  KEGG ...\n")
  res$KEGG <- tryCatch({
    k <- enrichKEGG(gene = gene_entrez, universe = universe_entrez, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, minGSSize = MIN_GS, maxGSSize = MAX_GS)
    to_df(setReadable(k, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
  }, error = function(e) { notes <<- c(notes, paste("KEGG failed (check the internet connection):", conditionMessage(e))); warning("KEGG failed: ", conditionMessage(e)); data.frame() })
  cat("  Reactome ...\n")
  res$Reactome <- tryCatch(to_df(enrichPathway(gene = gene_entrez, universe = universe_entrez, organism = "human", pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT,
                                               minGSSize = MIN_GS, maxGSSize = MAX_GS, readable = TRUE)),
                           error = function(e) { notes <<- c(notes, paste("Reactome failed:", conditionMessage(e))); data.frame() })
  cat("  Hallmark ...\n")
  res$Hallmark <- tryCatch({
    h <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"), error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
    gcol <- intersect(c("ncbi_gene", "entrez_gene"), names(h))[1]
    t2g <- unique(data.frame(term = h$gs_name, gene = as.character(h[[gcol]])))
    e <- enricher(gene = gene_entrez, universe = universe_entrez, TERM2GENE = t2g, pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, minGSSize = MIN_GS, maxGSSize = MAX_GS)
    to_df(setReadable(e, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
  }, error = function(e) { notes <<- c(notes, paste("Hallmark failed:", conditionMessage(e))); data.frame() })
  res
}
sheets <- c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "Hallmark")
cat("\nRunning pathway analysis (background:", UNIVERSE_MODE, ")\n"); res <- run_enrich()

# ==============================================================================
# [3] Comparison with the genome-background result (optional)
# ==============================================================================
key_terms <- data.frame(
  ont = c("GO_BP", "GO_BP", "GO_BP", "GO_BP", "GO_BP", "GO_MF", "GO_MF", "KEGG", "KEGG", "KEGG", "Reactome", "Reactome", "Reactome"),
  term = c("xenobiotic metabolic process", "cellular response to xenobiotic stimulus", "detoxification", "sulfur compound metabolic process", "lipid transport",
           "glutathione transferase activity", "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen",
           "Metabolism of xenobiotics by cytochrome P450", "Drug metabolism - cytochrome P450", "Glutathione metabolism",
           "Cytochrome P450 - arranged by substrate type", "Biological oxidations", "Glutathione conjugation"), stringsAsFactors = FALSE)
cmp <- NULL
if (!is.na(COMPARE_FILE) && file.exists(COMPARE_FILE)) {
  cat("\nComparison with the earlier result:", COMPARE_FILE, "\n")
  pick <- c("ID", "Description", "p.adjust", "FoldEnrichment", "Count")
  cmp_list <- lapply(sheets, function(s) {
    o <- as.data.frame(openxlsx::read.xlsx(COMPARE_FILE, sheet = s)); n <- res[[s]]
    if (!nrow(o) || !nrow(n)) return(NULL)
    m <- merge(o[, pick], n[, pick], by = "ID", all = TRUE, suffixes = c("_old", "_new"))
    m$Description <- ifelse(is.na(m$Description_new), m$Description_old, m$Description_new)
    so <- !is.na(m$p.adjust_old) & m$p.adjust_old < ALPHA; sn <- !is.na(m$p.adjust_new) & m$p.adjust_new < ALPHA
    m$status <- ifelse(so & sn, "remains significant", ifelse(so & !sn, "loses significance", ifelse(!so & sn, "newly significant", "not significant")))
    m$ontology <- s
    m[, c("ontology", "ID", "Description", "p.adjust_old", "p.adjust_new", "FoldEnrichment_old", "FoldEnrichment_new", "Count_old", "Count_new", "status")]
  })
  cmp <- do.call(rbind, cmp_list)
  summ <- as.data.frame(table(cmp$ontology, cmp$status)); names(summ) <- c("ontology", "status", "n")
  summ <- reshape(summ, idvar = "ontology", timevar = "status", direction = "wide"); names(summ) <- sub("^n\\.", "", names(summ)); summ[is.na(summ)] <- 0
  cat("\n[Change in significance versus the genome-background run]\n"); print(summ, row.names = FALSE)
  kt <- do.call(rbind, lapply(seq_len(nrow(key_terms)), function(i) {
    d <- cmp[cmp$ontology == key_terms$ont[i] & (tolower(cmp$Description) == tolower(key_terms$term[i])), , drop = FALSE]
    if (!nrow(d)) d <- data.frame(ontology = key_terms$ont[i], ID = NA, Description = key_terms$term[i], p.adjust_old = NA, p.adjust_new = NA, FoldEnrichment_old = NA, FoldEnrichment_new = NA, Count_old = NA, Count_new = NA, status = "term not found")
    d[1, ]
  }))
  ktp <- kt[, c("ontology", "Description", "p.adjust_old", "p.adjust_new", "status")]; ktp$Description <- substr(ktp$Description, 1, 45)
  cat("\n[Key terms of the manuscript]\n"); print(ktp, row.names = FALSE)
  changed <- cmp[cmp$status %in% c("loses significance", "newly significant"), ]; changed <- changed[order(changed$status, changed$p.adjust_new), ]
}

# ==============================================================================
# [4] Save
# ==============================================================================
used <- data.frame(
  item = c("input foreground file", "number of input genes", "number of genes after de-duplication", "ENTREZID mapped", "ENTREZID not mapped",
           "background (universe)", "number of background input genes", "number of background ENTREZIDs", "foreground genes added to the background",
           "GO ontology", "KEGG organism", "Reactome organism", "Hallmark (MSigDB) collection", "p.adjust method / minimum, maximum term size", "run time", "note"),
  value = c(basename(GENE_FILE), length(fg_sym), length(fg_sym), sum(fg_map$Mapped), sum(!fg_map$Mapped),
          if (UNIVERSE_MODE == "tested") sprintf("human counterparts of the tested genes (%s)", basename(UNIVERSE_FILE)) else "whole genome (default)",
          uni_info$n_input, uni_info$n_entrez, uni_info$n_added,
          "BP, CC, MF", "hsa", "human", "H (Hallmark, 50 gene sets)", sprintf("BH / %d, %d", MIN_GS, MAX_GS), format(Sys.time(), "%Y-%m-%d %H:%M"),
          paste(c(notes, sprintf("clusterProfiler %s", as.character(packageVersion("clusterProfiler")))), collapse = " | ")),
  check.names = FALSE, stringsAsFactors = FALSE)
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", halign = "center", border = "Bottom")
put <- function(wb, name, df, start = 1) { if (is.null(df)) return(invisible()); if (!(name %in% names(wb))) addWorksheet(wb, name); writeData(wb, name, as.data.frame(df), startRow = start, headerStyle = hs); freezePane(wb, name, firstActiveRow = start + 1) }
wb <- createWorkbook()
addWorksheet(wb, "Used_genes"); writeData(wb, "Used_genes", used, headerStyle = hs); writeData(wb, "Used_genes", fg_map, startRow = nrow(used) + 3, headerStyle = hs)
setColWidths(wb, "Used_genes", cols = 1:4, widths = c(34, 60, 14, 10))
for (s in sheets) put(wb, s, res[[s]])
if (!is.null(cmp)) { put(wb, "Compare_summary", summ); put(wb, "Key_terms", kt); put(wb, "Compare_changed", changed) }
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
writeLines(capture.output(sessionInfo()), sub("\\.xlsx$", "_sessionInfo.txt", OUT_XLSX))
cat("\nSaved:", OUT_XLSX, "\nDone.\n")
