################################################################################
#  01c_window_gene_mapping.R
#  Assignment of genes to the colocalization windows (GRCh37). Output = Supplementary Table S4 (gene to window assignment of
#  the 27 supported windows) and the gene-to-window table that is read by 05_ortholog_audit.R.
#
#  METHODS (as in the manuscript)
#    - Gene coordinates: UCSC RefSeq genes (hg19 refGene). Transcripts of the same symbol on the same chromosome are merged
#      into one gene position when they overlap or lie within 100 kb of each other; a symbol with positions further apart (snRNA,
#      miRNA and duplicated genes) keeps one position per cluster. A gene is "coding" if one of its transcripts is an NM_ record,
#      otherwise "noncoding" (NR_ records only).
#    - Optional supplement from the FUMA gene files (genes.txt of the 12 GWAS datasets): symbols that are absent from refGene
#      (e.g. Ensembl clone-based names and former symbols) are added with the FUMA coordinates; the biotype is "coding" if the
#      FUMA type is protein_coding. A FUMA symbol whose position overlaps a refGene gene by >= 80 % (of both genes, same strand,
#      compatible biotype) is a former name of that gene and is not added. An overlap of 50 - 80 % is not merged, but if both
#      symbols fall into the same window they are marked in the column same_gene_renamed_as. These marks are candidates: the
#      manuscript lists two genes with a former and a current symbol (SPX = C12orf39, MYMK = TMEM8C); all pairs with a partial
#      overlap are written to fuma_possible_aliases.tsv for review.  With the real FUMA files four pairs are marked: the two of the
#      manuscript and, in addition, C19orf45 / TEX45 and RP11-45M22.4 / FLCN (Table S4 lists these as separate genes; see docs/PROVENANCE.md).
#    - A gene is assigned to a window if its interval overlaps the window (same chromosome). overlap_with_window_bp is the
#      length of the overlap, distance_to_lead_kb is the distance from the lead variant to the nearest gene end (0 if the lead
#      variant lies inside the gene). A gene that occurs twice in one window keeps the position nearest to the lead variant.
#    - n_windows_containing_gene counts the windows (of all 458, supported or not) that contain the gene. The representative window
#      of a gene is the window with the highest tier and, among equal tiers, the nearest lead variant.
#
#  INPUTS (DATA_DIR, searched by file name, also in sub-folders)
#    window_tiers.tsv                  windows and tiers from 01b_coloc_window_tiers.R; if absent, Table S2a of the
#                                      Supplementary Tables workbook (Antioxidant_ThyroidLiver_Supplementary_Tables_S1-S13.xlsx)
#    hg19_refgene.txt                  UCSC Table Browser, assembly hg19 (GRCh37), group "Genes and Gene Predictions", track
#                                      "NCBI RefSeq", table refGene, output format "all fields from selected table" or CSV with the
#                                      fields name, chrom, strand, txStart, txEnd, name2. The goldenPath file refGene.txt(.gz) of
#                                      hg19 (no header) is also accepted if it is renamed to hg19_refgene.txt.
#    gwas_gene_flags.tsv               (optional) Gene, GWAS_PDFF, GWAS_cT1, GWAS_TSH, GWAS_FT3, GWAS_FT4 (0/1); in data/gene_sets/.
#                                      Adds the GWAS columns and the shared-gene flag (in_shared_199).
#    FUMA/<dataset>/genes.txt          (optional, folder FUMA_DIR, default DATA_DIR/FUMA) FUMA SNP2GENE output of the GWAS datasets;
#                                      all genes.txt files below FUMA_DIR are used. Without them 509 of the 526 gene-window rows of
#                                      Table S4 are found; the 17 others have symbols that exist only in the FUMA files.
#
#  OUTPUTS (OUT_DIR, default results/01c_gene_mapping)
#    window_gene_mapping.tsv             all windows (supported and not), one row per window and gene
#    window_gene_mapping_supported.tsv   the 27 supported windows (columns of Table S4 except the mouse-retrieval columns)
#    shared_gene_coloc_status.tsv        the 199 shared genes with colocalization status (Table S3b; needs gwas_gene_flags.tsv)
#    261019_gene_locus_full_mapping.xlsx sheet 3_Gene_in_locus: all rows in the layout read by 05_ortholog_audit.R
#    fuma_possible_aliases.tsv           FUMA symbols with a partial overlap to a refGene gene (for review)
#    01c_window_gene_mapping_sessionInfo.txt
#
#  REPRODUCTION CHECK  With the Supplementary Tables workbook available the script compares its output for the 27 supported
#  windows with Table S4 (genes, coordinates, overlap, distance, biotype, number of windows, representative window).
#  See docs/PROVENANCE.md for the result.
#
#  USAGE  DATA_DIR=/path/to/data Rscript 01c_window_gene_mapping.R
################################################################################
invisible(try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(readxl); library(openxlsx) })

## ---- settings ----------------------------------------------------------------
DATA_DIR  <- Sys.getenv("DATA_DIR", "data")
OUT_DIR   <- Sys.getenv("OUT_DIR", file.path("results", "01c_gene_mapping"))
FUMA_DIR  <- Sys.getenv("FUMA_DIR", file.path(DATA_DIR, "FUMA"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
CLUSTER_GAP_BP    <- 100000   # transcripts of one symbol further apart than this are separate positions
ALIAS_MIN_OVERLAP <- 0.8      # FUMA symbol is a former name of a refGene gene
ALIAS_REVIEW_MIN  <- 0.5      # partial overlap: kept, both symbols marked
VALID_CHR <- c(as.character(1:22), "X", "Y")
SUPPORTED <- c("high", "moderate")

ROOTS <- unique(c(DATA_DIR, file.path("results", "01b_tiers"), "data", getwd())); ROOTS <- ROOTS[dir.exists(ROOTS)]
find_quiet <- function(pattern) {
  hits <- unlist(lapply(ROOTS, function(r) list.files(r, pattern = pattern, recursive = TRUE, full.names = TRUE)))
  hits <- hits[!grepl("^~\\$", basename(hits))]
  if (length(hits)) hits[which.max(file.mtime(hits))] else NA_character_
}
F_TIERS <- find_quiet("^window_tiers\\.tsv$")
F_REF   <- find_quiet("^hg19_refgene\\.txt(\\.gz)?$")
F_FLAGS <- find_quiet("^gwas_gene_flags\\.tsv$")
SUPP    <- { e <- Sys.getenv("SUPP_XLSX", ""); if (nzchar(e)) e else find_quiet("Supplementary_Tables_S1-S13.*\\.xlsx$") }
if (is.na(F_REF)) stop("hg19_refgene.txt not found in: ", paste(ROOTS, collapse = " | "),
                       "\n  Download it from the UCSC Table Browser (hg19, NCBI RefSeq, refGene); see data/README.md.", call. = FALSE)
if (is.na(F_TIERS) && is.na(SUPP)) stop("No window table found: run 01b_coloc_window_tiers.R or put the Supplementary Tables workbook in DATA_DIR.", call. = FALSE)

## ---- Supplementary Tables reader (only needed when 01b output is absent, and for the check) -----
if (!is.na(SUPP)) {
  SUPP_SHEETS <- excel_sheets(SUPP)
  locate_supp <- function(key) {
    sh <- SUPP_SHEETS[SUPP_SHEETS %in% c(key, paste("Table", key))]
    if (length(sh) == 1) return(list(sheet = sh, c0 = 1L, c1 = NA_integer_))
    sh <- SUPP_SHEETS[SUPP_SHEETS %in% paste("Table", sub("^(S[0-9]+).*$", "\\1", key))]
    if (length(sh) != 1) return(NULL)
    r1 <- unlist(read_excel(SUPP, sheet = sh, col_names = FALSE, n_max = 1, col_types = "text", .name_repair = "minimal"), use.names = FALSE)
    c0 <- which(startsWith(r1, paste0("Table ", key, ".")))[1]
    if (is.na(c0)) return(NULL)
    nxt <- which(grepl("^Table S[0-9]+[a-z]?\\.", r1)); nxt <- nxt[nxt > c0]
    list(sheet = sh, c0 = c0, c1 = if (length(nxt)) nxt[1] - 1L else NA_integer_)
  }
  read_supp <- function(key, first_col) {
    loc <- locate_supp(key)
    if (is.null(loc)) stop(sprintf("Table %s not found in %s", key, basename(SUPP)), call. = FALSE)
    rng <- function(r0, r1) readxl::cell_limits(c(r0, loc$c0), c(r1, loc$c1))
    top <- read_excel(SUPP, sheet = loc$sheet, range = rng(1, 25), col_names = FALSE, col_types = "text", .name_repair = "minimal")
    h <- which(trimws(as.character(top[[1]])) == first_col)[1]
    if (is.na(h)) stop(sprintf("Header row (first cell '%s') not found for Table %s", first_col, key), call. = FALSE)
    d <- as.data.table(read_excel(SUPP, sheet = loc$sheet, range = rng(h, NA), .name_repair = "minimal"))
    blank <- (is.na(names(d)) | names(d) == "" | grepl("^\\.\\.\\.[0-9]+$", names(d))) & vapply(d, function(x) all(is.na(x)), logical(1))
    d <- d[, which(!blank), with = FALSE]
    d[rowSums(!is.na(d)) > 0]
  }
}

## ---- helpers --------------------------------------------------------------------
norm_chr <- function(x) { x <- sub("^CHR", "", toupper(trimws(as.character(x)))); x[x == "23"] <- "X"; x[x == "24"] <- "Y"; x[x %in% c("M", "MT")] <- "MT"; x }
key_of <- function(x) toupper(trimws(as.character(x)))
nz <- function(x) !is.na(x) & x != ""

## UCSC refGene, one row per transcript (UCSC starts are 0-based, hence +1)
read_refgene <- function(f) {
  d <- fread(f, sep = "auto", quote = "", fill = TRUE, header = "auto", data.table = TRUE, colClasses = "character")
  setnames(d, gsub('"', "", names(d), fixed = TRUE))
  for (j in names(d)) set(d, j = j, value = gsub('"', "", d[[j]], fixed = TRUE))
  nm <- tolower(sub("^.*\\.", "", sub("^#", "", names(d))))
  if (all(c("name", "chrom", "txstart", "txend", "name2") %in% nm)) {                       # table with a header
    g <- function(n) d[[match(n, nm)]]
    out <- data.table(id = g("name"), symbol = g("name2"), chr = norm_chr(g("chrom")),
                      strand = if ("strand" %in% nm) g("strand") else NA_character_,
                      txStart = as.numeric(g("txstart")), txEnd = as.numeric(g("txend")))
  } else if (ncol(d) >= 13) {                                                              # goldenPath refGene.txt without header
    out <- data.table(id = d[[2]], symbol = d[[13]], chr = norm_chr(d[[3]]), strand = d[[4]],
                      txStart = as.numeric(d[[5]]), txEnd = as.numeric(d[[6]]))
  } else stop("refGene file not recognised: it needs the columns name, chrom, txStart, txEnd, name2 (see data/README.md).", call. = FALSE)
  out <- out[chr %in% VALID_CHR & nz(symbol) & !is.na(txStart) & !is.na(txEnd)]
  out[, .(id, symbol, chr, strand, gs = txStart + 1, ge = txEnd)]
}
## transcripts of one symbol and chromosome -> position clusters
cluster_intervals <- function(d, gap = CLUSTER_GAP_BP) {
  d <- copy(d); setorder(d, key, chr, gs, ge)
  d[, .prev := data.table::shift(cummax(ge)), by = .(key, chr)]
  d[, .newc := is.na(.prev) | gs > .prev + gap]
  d[, cl := cumsum(.newc), by = .(key, chr)]
  d[, c(".prev", ".newc") := NULL]
  d[]
}
read_fuma_genes <- function(f) {
  d <- fread(f, sep = "\t", quote = "", fill = TRUE, data.table = TRUE)
  if (!all(c("symbol", "chr", "start", "end") %in% names(d))) return(NULL)
  data.table(symbol = as.character(d$symbol), chr = norm_chr(d$chr), start = as.numeric(d$start), end = as.numeric(d$end),
             strand = if ("strand" %in% names(d)) as.character(d$strand) else NA_character_,
             type = if ("type" %in% names(d)) as.character(d$type) else NA_character_)[nz(symbol) & chr %in% VALID_CHR & !is.na(start) & !is.na(end)]
}

## ---- [1] windows -------------------------------------------------------------------
if (!is.na(F_TIERS)) {
  W <- fread(F_TIERS)[, .(window, chr = norm_chr(chr), ws = as.numeric(start_GRCh37), we = as.numeric(end_GRCh37), lead_variant, tier, maxPP4 = as.numeric(`max_cross_domain_PP.H4`),
                          supporting_pairs = as.character(supporting_pairs))]
} else {
  a <- read_supp("S2a", "window")
  W <- a[, .(window, chr = norm_chr(chr), ws = as.numeric(start_GRCh37), we = as.numeric(end_GRCh37), lead_variant, tier, maxPP4 = as.numeric(`max_cross_domain_PP.H4`),
             supporting_pairs = as.character(supporting_pairs))]
}
W[, lead_pos := as.numeric(sub("^.*:", "", sub("_.*$", "", lead_variant)))]
W[, rk := fcase(tier == "high", 3L, tier == "moderate", 2L, default = 1L)]
W[, supported := tier %in% SUPPORTED]
message(sprintf("[1] %d windows, %d supported (%s)", nrow(W), sum(W$supported), paste(names(table(W$tier)), table(W$tier), collapse = ", ")))

## ---- [2] gene positions: refGene (+ FUMA supplement) ----------------------------------
ref_tx <- read_refgene(F_REF); ref_tx[, key := key_of(symbol)]
ref_tx <- cluster_intervals(ref_tx)
ref <- ref_tx[, .(symbol = symbol[1], gs = min(gs), ge = max(ge), is_coding = any(grepl("^NM_", id)),
                  strand = strand[1], coord_src = "refGene"), by = .(key, chr, cl)]
message(sprintf("[2] refGene: %d symbols, %d positions", uniqueN(ref$key), nrow(ref)))

fuma_files <- if (dir.exists(FUMA_DIR)) list.files(FUMA_DIR, pattern = "^genes\\.txt$", recursive = TRUE, full.names = TRUE) else character(0)
alias_tab <- data.table(); possible_tab <- data.table()
fuma_add <- data.table(symbol = character(), chr = character(), cl = integer(), gs = numeric(), ge = numeric(),
                       is_coding = logical(), coord_src = character())
fuma_add[, key := character()]          # (data.table(key = ...) would be read as the key argument)
setcolorder(fuma_add, c("key", "symbol", "chr", "cl", "gs", "ge", "is_coding", "coord_src"))
if (length(fuma_files)) {
  fuma_raw <- rbindlist(lapply(fuma_files, read_fuma_genes), fill = TRUE)
  fuma_raw[, key := key_of(symbol)]
  fz <- cluster_intervals(fuma_raw[, .(key, symbol, chr, gs = start, ge = end, strand, type)])
  fuma_g <- fz[, .(symbol = symbol[1], gs = min(gs), ge = max(ge), strand = strand[1], type = type[which(!is.na(type))[1]]), by = .(key, chr, cl)]
  fuma_only <- fuma_g[!(key %in% ref$key)]
  message(sprintf("    FUMA genes.txt: %d files, %d symbols, %d of them absent from refGene", length(fuma_files), uniqueN(fuma_g$key), uniqueN(fuma_only$key)))
  if (nrow(fuma_only)) {
    r_i <- ref[, .(chr, rs = gs, re = ge, ref_symbol = symbol, ref_key = key, ref_strand = strand,
                   ref_biotype = fifelse(is_coding, "coding", "noncoding"))]
    setkey(r_i, chr, rs, re)
    f_i <- fuma_only[, .(fuma_key = key, fuma_cl = cl, fuma_symbol = symbol, chr, fs = gs, fe = ge, fuma_strand = strand, fuma_type = type)]
    ov <- suppressWarnings(foverlaps(f_i, r_i, by.x = c("chr", "fs", "fe"), by.y = c("chr", "rs", "re"), type = "any", nomatch = 0L))
    if (nrow(ov)) {
      ov[, ov_bp := pmin(fe, re) - pmax(fs, rs) + 1]
      ov[, `:=`(frac_fuma = ov_bp / (fe - fs + 1), frac_ref = ov_bp / (re - rs + 1))]
      ov[, score := pmin(frac_fuma, frac_ref)]
      known <- ov$fuma_strand %in% c("+", "-") & ov$ref_strand %in% c("+", "-")
      ov <- ov[!(known & ov$fuma_strand != ov$ref_strand)]                    # opposite strands: a different gene
      ov[, conflict := (!is.na(fuma_type) & fuma_type != "protein_coding" & ref_biotype == "coding") |
                       (!is.na(fuma_type) & fuma_type == "protein_coding" & ref_biotype == "noncoding")]
      setorder(ov, fuma_key, fuma_cl, conflict, -score)
      bestm <- ov[, .SD[1], by = .(fuma_key, fuma_cl)]
      alias_tab    <- bestm[score >= ALIAS_MIN_OVERLAP & !conflict]
      possible_tab <- bestm[(score >= ALIAS_REVIEW_MIN & score < ALIAS_MIN_OVERLAP) | (score >= ALIAS_MIN_OVERLAP & conflict)]
    }
  }
  fuma_add <- if (nrow(alias_tab)) fuma_only[!alias_tab, on = .(key = fuma_key, cl = fuma_cl)] else copy(fuma_only)
  fuma_add <- fuma_add[, .(key, symbol, chr, cl, gs, ge, is_coding = ifelse(is.na(type), NA, type == "protein_coding"), coord_src = "FUMA")]
  message(sprintf("    former names merged into refGene genes: %d | partial overlaps kept (listed for review): %d | FUMA genes added: %d",
                  nrow(alias_tab), nrow(possible_tab), nrow(fuma_add)))
  if (nrow(possible_tab)) fwrite(possible_tab[, .(fuma_symbol, ref_symbol, fuma_start = fs, fuma_end = fe, ref_start = rs, ref_end = re,
                                                   overlap_fraction_of_shorter = round(score, 3), fuma_strand, ref_strand, fuma_type, ref_biotype)],
                                  file.path(OUT_DIR, "fuma_possible_aliases.tsv"), sep = "\t")
} else {
  message("    No FUMA genes.txt found in ", FUMA_DIR, ": gene positions come from refGene only")
}
univ <- rbind(ref[, .(key, symbol, chr, cl, gs, ge, is_coding, coord_src)], fuma_add)
univ[, n_sites_on_chr := .N, by = .(key, chr)]
alias_map <- if (nrow(alias_tab)) setNames(alias_tab[, .SD[1], by = fuma_key]$ref_key, alias_tab[, .SD[1], by = fuma_key]$fuma_key) else character(0)

## ---- [3] genes in windows -----------------------------------------------------------
loci_key <- W[, .(chr, ls = ws, le = we, window)]; setkey(loci_key, chr, ls, le)
univ_k <- copy(univ); setkey(univ_k, chr, gs, ge)
hits <- foverlaps(univ_k, loci_key, by.x = c("chr", "gs", "ge"), by.y = c("chr", "ls", "le"), type = "any", nomatch = 0L)
hits <- merge(hits, W[, .(window, lead_variant, tier, rk, maxPP4, lead_pos, ws, we, supported)], by = "window")
hits[, overlap_bp := pmin(ge, we) - pmax(gs, ws)]
hits[, dist_bp := fifelse(lead_pos >= gs & lead_pos <= ge, 0, pmin(abs(lead_pos - gs), abs(lead_pos - ge)))]
hits[, distance_to_lead_kb := round(dist_bp / 1000, 1)]
setorder(hits, window, key, dist_bp)
hits <- hits[, .SD[1], by = .(window, key)]                                   # nearest position of a gene within a window
hits[, biotype := fifelse(is.na(is_coding), "unknown", fifelse(is_coding, "coding", "noncoding"))]
hits[, tier_manuscript := tier]
n_win <- hits[, .(n_windows_containing_gene = uniqueN(window)), by = key]          # all 458 windows
hits <- merge(hits, n_win, by = "key", all.x = TRUE, sort = FALSE)
rep_tab <- copy(hits); setorder(rep_tab, key, -rk, dist_bp, -maxPP4, window)
rep_tab <- rep_tab[, .SD[1], by = key][, .(key, rep_window = window)]
hits <- merge(hits, rep_tab, by = "key", all.x = TRUE, sort = FALSE)
hits[, in_representative_window := !is.na(rep_window) & window == rep_window]

## symbols marked as the same gene (partial overlap): both symbols must occur in the same window
hits[, same_gene_renamed_as := NA_character_]
pm <- data.table(fk = character(), rk_ = character(), fs_ = character(), rs_ = character())
if (nrow(possible_tab)) {
  ## partial overlap (50 - 80 % of the shorter gene), compatible biotype and strand; pairs with a biotype conflict are only listed in fuma_possible_aliases.tsv
  pm <- unique(possible_tab[score < ALIAS_MIN_OVERLAP & !conflict, .(fk = fuma_key, rk_ = key_of(ref_symbol), fs_ = fuma_symbol, rs_ = ref_symbol)])
  for (i in seq_len(nrow(pm))) {
    wa <- hits[key == pm$fk[i], window]; wb <- hits[key == pm$rk_[i], window]
    both <- intersect(wa, wb)
    if (length(both)) {
      hits[key == pm$fk[i] & window %in% both, same_gene_renamed_as := pm$rs_[i]]
      hits[key == pm$rk_[i] & window %in% both, same_gene_renamed_as := pm$fs_[i]]
    }
  }
}
hits[, distinct_gene := symbol]
for (i in seq_len(nrow(pm))) {
  keep <- hits[key == pm$rk_[i] & !is.na(same_gene_renamed_as), unique(symbol)]
  if (length(keep)) hits[key == pm$fk[i] & !is.na(same_gene_renamed_as), distinct_gene := keep[1]]
}

## ---- [4] GWAS flags (optional) ----------------------------------------------------------
fl_cols <- c("GWAS_PDFF", "GWAS_cT1", "GWAS_TSH", "GWAS_FT3", "GWAS_FT4")
if (!is.na(F_FLAGS)) {
  fl <- fread(F_FLAGS); fl[, key := key_of(Gene)]
  fl <- fl[, lapply(.SD, max), by = key, .SDcols = fl_cols]
  if (length(alias_map)) {                                                    # former FUMA names carry their flags to the refGene name
    add <- fl[key %in% names(alias_map)]; if (nrow(add)) { add[, key := unname(alias_map[key])]; fl <- rbind(fl, add)[, lapply(.SD, max), by = key, .SDcols = fl_cols] }
  }
  ii <- match(hits$key, fl$key)
  for (cn in fl_cols) set(hits, j = cn, value = ifelse(is.na(ii), 0L, fl[[cn]][ii]))
  hits[, `:=`(GWAS_liver = as.integer(GWAS_PDFF == 1 | GWAS_cT1 == 1), GWAS_thyroid = as.integer(GWAS_TSH == 1 | GWAS_FT3 == 1 | GWAS_FT4 == 1))]
  hits[, GWAS_shared := as.integer(GWAS_liver == 1 & GWAS_thyroid == 1)]
  hits[, in_shared_199 := GWAS_shared]
} else {
  message("[4] gwas_gene_flags.tsv not found: the GWAS columns are left empty")
  for (cn in c(fl_cols, "GWAS_liver", "GWAS_thyroid", "GWAS_shared", "in_shared_199")) hits[, (cn) := NA_integer_]
}
hits[, locus_pass := supported]

## ---- [5] output ---------------------------------------------------------------------------
setorder(hits, window, gs)
tab <- hits[, .(window, tier, lead_variant, chr, window_start_GRCh37 = ws, window_end_GRCh37 = we, gene = symbol,
                gene_start_GRCh37 = gs, gene_end_GRCh37 = ge, overlap_with_window_bp = overlap_bp, distance_to_lead_kb, biotype,
                in_shared_199, n_windows_containing_gene, in_representative_window, same_gene_renamed_as, distinct_gene,
                coord_src, supported)]
fwrite(tab[, !"supported"], file.path(OUT_DIR, "window_gene_mapping.tsv"), sep = "\t", na = "")
sup_tab <- tab[supported == TRUE, !"supported"]
fwrite(sup_tab, file.path(OUT_DIR, "window_gene_mapping_supported.tsv"), sep = "\t", na = "")
message(sprintf("[5] %d window-gene rows (%d genes) in all windows; %d rows (%d genes, %d coding) in the %d supported windows",
                nrow(tab), uniqueN(hits$key), nrow(sup_tab), uniqueN(hits[supported == TRUE, key]),
                uniqueN(hits[supported == TRUE & biotype == "coding", key]), sum(W$supported)))

## table in the layout read by 05_ortholog_audit.R
g3 <- hits[, .(Locus_ID = window, locus_tier = tier, locus_tier_manuscript = tier_manuscript, locus_pass = locus_pass, Lead_SNP = lead_variant,
               Chr = chr, Gene = symbol, gene_start = gs, gene_end = ge, dist_kb = distance_to_lead_kb, coord_src, n_sites_on_chr, biotype,
               GWAS_PDFF, GWAS_cT1, GWAS_TSH, GWAS_FT3, GWAS_FT4, GWAS_liver, GWAS_thyroid, GWAS_shared,
               is_representative = in_representative_window, n_loci_for_gene = n_windows_containing_gene, locus_maxPP4 = maxPP4)]
wb <- createWorkbook(); addWorksheet(wb, "3_Gene_in_locus"); writeData(wb, "3_Gene_in_locus", as.data.frame(g3))
saveWorkbook(wb, file.path(OUT_DIR, "261019_gene_locus_full_mapping.xlsx"), overwrite = TRUE)

## ---- check against Table S4 --------------------------------------------------------------------
if (!is.na(SUPP) && !is.null(locate_supp("S4"))) {
  s4 <- read_supp("S4", "window")
  s4[, key := paste(window, key_of(gene))]; sup_tab[, key := paste(window, key_of(gene))]
  mm <- merge(sup_tab, s4, by = "key", suffixes = c("", ".s4"))
  message(sprintf("\nCheck against Table S4: %d of %d S4 rows found (%d rows of the script output are not in S4)",
                  nrow(mm), nrow(s4), nrow(sup_tab) - nrow(mm)))
  cmp <- function(a, b) sum(as.numeric(a) != as.numeric(b), na.rm = TRUE)
  message(sprintf("  among the found rows: gene start differs in %d, gene end %d, overlap %d, distance %d, biotype %d, n windows %d",
                  cmp(mm$gene_start_GRCh37, mm$gene_start_GRCh37.s4), cmp(mm$gene_end_GRCh37, mm$gene_end_GRCh37.s4),
                  cmp(mm$overlap_with_window_bp, mm$overlap_with_window_bp.s4), cmp(mm$distance_to_lead_kb, mm$distance_to_lead_kb.s4),
                  sum(mm$biotype != mm$biotype.s4), cmp(mm$n_windows_containing_gene, mm$n_windows_containing_gene.s4)))
  message(sprintf("  representative window differs in %d, in_shared_199 differs in %d (NA = flags not provided), renamed-as differs in %d",
                  sum(as.logical(mm$in_representative_window) != as.logical(mm$in_representative_window.s4)),
                  if (all(is.na(mm$in_shared_199))) NA_integer_ else cmp(mm$in_shared_199, mm$in_shared_199.s4),
                  sum(fifelse(is.na(mm$same_gene_renamed_as), "", mm$same_gene_renamed_as) != fifelse(is.na(mm$same_gene_renamed_as.s4), "", as.character(mm$same_gene_renamed_as.s4)))))
  miss <- s4[!(key %in% sup_tab$key), paste(window, gene)]
  if (length(miss)) message("  S4 rows not found: ", paste(miss, collapse = ", "))
}
## ---- [6] colocalization status of the shared genes (Table S3b) ------------------------------------
## Shared genes = genes of both the thyroid-derived and the liver-derived GWAS gene set (199). Status of a gene: the tier of its representative window
## (high; "nominal" = the manuscript tier "moderate"; "no" = a tested window without PP.H4 >= 0.5); "outside_tested_loci" if the gene lies in no window;
## "coordinates_missing" if the gene has no position (not in refGene and no FUMA genes.txt supplied: former symbols such as HIST1H1A need the FUMA files).
if (!is.na(F_FLAGS)) {
  fl2 <- fread(F_FLAGS); fl2[, key := key_of(Gene)]
  shg <- fl2[(GWAS_PDFF == 1 | GWAS_cT1 == 1) & (GWAS_TSH == 1 | GWAS_FT3 == 1 | GWAS_FT4 == 1)]
  pos_key <- function(k) ifelse(k %in% univ$key, k, ifelse(k %in% names(alias_map), unname(alias_map[k]), k))
  shg[, pkey := pos_key(key)]
  shg[, has_pos := pkey %in% univ$key]
  shg[, rep_window := rep_tab$rep_window[match(pkey, rep_tab$key)]]
  shg[, tier_rep := W$tier[match(rep_window, W$window)]]
  shg[, COLOC_status := fcase(!has_pos, "coordinates_missing", is.na(rep_window), "outside_tested_loci", tier_rep == "high", "high",
                              tier_rep == "moderate", "nominal", default = "no")]
  shg[, supporting_pairs := fifelse(COLOC_status %in% c("high", "nominal"), W$supporting_pairs[match(rep_window, W$window)], NA_character_)]
  shg[, thyroid_traits := apply(.SD, 1, function(r) paste(c("TSH", "FT3", "FT4")[r == 1], collapse = "; ")), .SDcols = c("GWAS_TSH", "GWAS_FT3", "GWAS_FT4")]
  shg[, liver_traits := apply(.SD, 1, function(r) paste(c("PDFF", "cT1")[r == 1], collapse = "; ")), .SDcols = c("GWAS_PDFF", "GWAS_cT1")]
  s3b_out <- shg[, .(Gene, thyroid_traits, liver_traits, COLOC_status, nearest_tested_window = rep_window, supporting_pairs)]
  s3b_out <- s3b_out[order(factor(COLOC_status, c("high", "nominal", "no", "outside_tested_loci", "coordinates_missing")), Gene)]
  fwrite(s3b_out, file.path(OUT_DIR, "shared_gene_coloc_status.tsv"), sep = "\t", na = "")
  message(sprintf("\n[6] %d shared genes: %s", nrow(s3b_out), paste(sprintf("%s %d", names(table(s3b_out$COLOC_status)), table(s3b_out$COLOC_status)), collapse = ", ")))
  if (any(s3b_out$COLOC_status == "coordinates_missing"))
    message("    ", sum(s3b_out$COLOC_status == "coordinates_missing"), " shared genes have no position (former symbols, e.g. ",
            paste(head(s3b_out[COLOC_status == "coordinates_missing", Gene], 3), collapse = ", "), "): supply the FUMA genes.txt files (FUMA_DIR) to place them")
  if (!is.na(SUPP) && !is.null(locate_supp("S3b"))) {
    s3b <- read_supp("S3b", "Gene")
    mm <- merge(s3b_out, s3b, by = "Gene", suffixes = c("", ".S3b"))
    ok <- mm$COLOC_status != "coordinates_missing"
    chk <- function(a, b) sum(fifelse(is.na(a), "", as.character(a))[ok] != fifelse(is.na(b), "", as.character(b))[ok])
    message(sprintf("    Check against Table S3b: %d of %d genes found; among the %d genes with a position: status differs in %d, nearest window in %d, supporting pairs in %d, traits in %d",
                    nrow(mm), nrow(s3b), sum(ok), chk(mm$COLOC_status, mm$COLOC_status.S3b), chk(mm$nearest_tested_window, mm$nearest_tested_window.S3b),
                    chk(mm$supporting_pairs, mm$supporting_pairs.S3b), chk(paste(mm$thyroid_traits, mm$liver_traits), paste(mm$thyroid_traits.S3b, mm$liver_traits.S3b))))
  }
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "01c_window_gene_mapping_sessionInfo.txt"))
message("Done. Output: ", normalizePath(OUT_DIR))
