################################################################################
#  00b_magma_gene_set_enrichment.R
#  MAGMA competitive gene-set analysis between the thyroid-derived and the liver-derived GWAS gene sets (Supplementary Table S1).
#
#  DESIGN (as in the manuscript)
#    gene-level analysis of 12 GWAS datasets (SNP-wise mean model, 1000 Genomes Phase 3 EUR reference, build 37, gene window 10 kb up / 10 kb down)
#    competitive gene-set analysis in both directions, 12 tests, Bonferroni threshold P = 0.05 / 12 = 0.0042
#      thyroid-derived GWAS gene set tested in the 8 liver datasets (PDFF x 4 and cT1 x 4 adjustment models)
#      liver-derived GWAS gene set   tested in the 4 thyroid datasets (TSH A, TSH B, FT3, FT4)
#    SNP identifiers: the summary statistics are read by chromosome and position only; the rsID of the reference panel (.bim) at the same position is
#    used, so that all 12 files and the annotation use the same SNP IDs.  Rows without p-value (or N) are dropped before MAGMA is called.
#    Gene sets: gene symbols of data/gene_sets/gwas_gene_flags.tsv converted to Entrez IDs with org.Hs.eg.db (set names THYROID_GWAS_GENES and LIVER_GWAS_GENES).
#
#  INPUTS  (searched by name in DATA_DIR, and in OUT_DIR for finished MAGMA output)
#    MAGMA v1.10 program                 env MAGMA_BIN, or `magma` on the PATH, or a file named magma below DATA_DIR
#    NCBI37.3.gene.loc                   MAGMA gene locations, build 37
#    g1000_eur.bed / .bim / .fam         1000 Genomes Phase 3 EUR reference panel for MAGMA (build 37)
#    12 GWAS summary statistics (all GRCh37), found by the GWAS Catalog ID in the file or folder name or by the file name:
#      GCST90446650-653 (PDFF unadjusted, BMI-, WFM-, VAT-adjusted), GCST90446646-649 (cT1, same order), GCST90296333 (TSH B),
#      GCST010653 (TSH A), FT3_summary_data, FT4_summary_data          (.tsv / .txt, plain or .gz)
#    gwas_gene_flags.tsv                 data/gene_sets/ of this repository, or the output of 00_gwas_gene_sets_from_FUMA.R
#    R packages: data.table, readxl, org.Hs.eg.db (Bioconductor; without it NCBI37.3.gene.loc is used to convert symbols and set sizes may differ slightly)
#    These large files are not in the repository; data/README.md says where to download them.
#    Without MAGMA / the reference files the script only assembles the table from MAGMA output (*.gsa.out) that is already in OUT_DIR.
#
#  OUTPUTS (OUT_DIR, default results/00b_magma)
#    gene_annot.genes.annot, geneLevel_<dataset>.genes.out/.genes.raw/.log, gsa_<set>_in_<dataset>.gsa.out   MAGMA files (existing gene-level results are reused)
#    thyroid_geneset.txt, liver_geneset.txt                                                                   MAGMA set files
#    magma_competitive_gene_set_results.tsv                                                                   = Supplementary Table S1 (12 rows)
#    00b_magma_gene_set_enrichment_sessionInfo.txt
#  Check printed to the console: the table is compared with Table S1 when the Supplementary Tables workbook is found.
#
#  MAGMA runtime: one gene-level run takes 0.5 - 1.5 h and several GB of memory; runs are resumed (finished .genes.out files are skipped).
#
#  USAGE   DATA_DIR=/path/to/data  MAGMA_BIN=/path/to/magma  Rscript 00b_magma_gene_set_enrichment.R        (RUN_MAGMA=0 to only assemble the table)
################################################################################
invisible(try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(readxl) })

DATA_DIR <- Sys.getenv("DATA_DIR", "data")
OUT_DIR  <- Sys.getenv("OUT_DIR", file.path("results", "00b_magma"))
RUN_MAGMA <- Sys.getenv("RUN_MAGMA", "auto")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_DIR <- normalizePath(OUT_DIR)
N_TESTS <- 12; ALPHA_BONF <- 0.05 / N_TESTS

## ---- the 12 datasets --------------------------------------------------------------------------
## N_fixed: sample size of the liver GWAS models (Nature Medicine 2024, Table 1); the thyroid files carry N in their own column.
DS <- rbindlist(list(
  data.table(dataset = paste0("PDFF_", c("unadjusted", "BMIadjusted", "WFMadjusted", "VATadjusted")), domain = "liver", trait_group = "PDFF",
             file_id = paste0("GCST90446", 650:653), N_fixed = c(36394, 35146, 34367, 33540), chr_col = "chromosome", bp_col = "base_pair_location", p_col = "p_value", n_col = "n"),
  data.table(dataset = paste0("cT1_", c("unadjusted", "BMIadjusted", "WFMadjusted", "VATadjusted")), domain = "liver", trait_group = "cT1",
             file_id = paste0("GCST90446", 646:649), N_fixed = c(30481, 29312, 28638, 27782), chr_col = "chromosome", bp_col = "base_pair_location", p_col = "p_value", n_col = "n"),
  data.table(dataset = c("TSH_B", "TSH_A", "FT3", "FT4"), domain = "thyroid", trait_group = c("TSH dataset B (GCST90296333)", "TSH dataset A (GCST010653)", "FT3", "FT4"),
             file_id = c("GCST90296333", "GCST010653", "FT3_summary_data", "FT4_summary_data"), N_fixed = NA_real_,
             chr_col = c("chromosome", "Chr", "chr", "chr"), bp_col = c("base_pair_location", "Pos", "position", "position"),
             p_col = c("p_value", "P-value", "P.value", "P.value"), n_col = c("n", "N", "N", "N"))))
## label of the tested dataset in Table S1
DS[, tested_in := fifelse(domain == "liver", dataset, trait_group)]
LABEL_SET <- c(thyroid = "Thyroid-derived GWAS genes (tested in liver datasets)", liver = "Liver-derived GWAS genes (tested in thyroid datasets)")

## ---- locate files by name ----------------------------------------------------------------------------
ROOTS <- unique(normalizePath(c(DATA_DIR, OUT_DIR, "data")[dir.exists(c(DATA_DIR, OUT_DIR, "data"))]))
all_files <- unique(normalizePath(unlist(lapply(ROOTS, function(r) list.files(r, recursive = TRUE, full.names = TRUE)))))
all_files <- all_files[!grepl("^~\\$", basename(all_files))]
find_name <- function(pattern) { h <- all_files[grepl(pattern, basename(all_files))]; if (length(h)) h[1] else NA_character_ }
SUPP    <- { e <- Sys.getenv("SUPP_XLSX", ""); if (nzchar(e)) e else find_name("Supplementary_Tables_S1-S13.*\\.xlsx$") }
F_FLAGS <- { h <- all_files[basename(all_files) == "gwas_gene_flags.tsv" & !grepl("/results/", all_files)]; if (length(h)) h[1] else NA_character_ }
GENE_LOC <- { e <- Sys.getenv("GENE_LOC", ""); if (nzchar(e)) e else find_name("^NCBI37\\.3\\.gene\\.loc$") }
BFILE <- { e <- Sys.getenv("REF_BFILE", ""); if (nzchar(e)) e else { b <- find_name("^g1000_eur\\.bim$"); if (is.na(b)) NA_character_ else sub("\\.bim$", "", b) } }
MAGMA <- { e <- Sys.getenv("MAGMA_BIN", ""); if (nzchar(e)) path.expand(e) else { w <- unname(Sys.which("magma")); if (nzchar(w)) w else { h <- all_files[basename(all_files) %in% c("magma", "magma.exe")]; if (length(h)) h[1] else NA_character_ } } }
find_sumstats <- function(key) {
  h <- all_files[grepl(key, all_files, fixed = TRUE) & grepl("\\.(tsv|txt|csv)(\\.gz)?$", all_files) & !grepl("/results/", all_files)]
  h <- h[!grepl("\\.(out|raw|log|annot|loc)$", h)]
  if (length(h) > 1) { message("  more than one file for '", key, "', the first is used: ", h[1]); h <- h[1] }
  if (length(h)) h else NA_character_
}

## ---- gene sets as MAGMA set files ---------------------------------------------------------------------------
if (is.na(F_FLAGS)) stop("gwas_gene_flags.tsv not found in ", paste(ROOTS, collapse = " | "), call. = FALSE)
fl <- fread(F_FLAGS)
sym_thy <- unique(fl[GWAS_TSH == 1 | GWAS_FT3 == 1 | GWAS_FT4 == 1, Gene])
sym_liv <- unique(fl[GWAS_PDFF == 1 | GWAS_cT1 == 1, Gene])
symbol_to_entrez <- function(symbols) {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
    m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(symbols), keytype = "SYMBOL", columns = "ENTREZID"))
    return(unique(m$ENTREZID[!is.na(m$ENTREZID)]))
  }
  if (is.na(GENE_LOC)) stop("Neither org.Hs.eg.db nor NCBI37.3.gene.loc is available to convert gene symbols to Entrez IDs.", call. = FALSE)
  message("  org.Hs.eg.db not installed: symbols are converted with ", basename(GENE_LOC), " (set sizes may differ slightly from Table S1)")
  gl <- fread(GENE_LOC, header = FALSE, select = c(1, 6), col.names = c("entrez", "symbol"), colClasses = "character")
  unique(gl$entrez[gl$symbol %in% symbols])
}
ent_thy <- symbol_to_entrez(sym_thy); ent_liv <- symbol_to_entrez(sym_liv)
message(sprintf("Gene sets: thyroid-derived %d symbols -> %d Entrez IDs; liver-derived %d symbols -> %d Entrez IDs", length(sym_thy), length(ent_thy), length(sym_liv), length(ent_liv)))
SETFILE <- c(thyroid = file.path(OUT_DIR, "thyroid_geneset.txt"), liver = file.path(OUT_DIR, "liver_geneset.txt"))
SETNAME <- c(thyroid = "THYROID_GWAS_GENES", liver = "LIVER_GWAS_GENES")
writeLines(paste(c(SETNAME["thyroid"], ent_thy), collapse = "\t"), SETFILE["thyroid"])
writeLines(paste(c(SETNAME["liver"],   ent_liv), collapse = "\t"), SETFILE["liver"])

## ---- MAGMA runs ---------------------------------------------------------------------------------------------------
gsa_file <- function(set, ds) file.path(OUT_DIR, sprintf("gsa_%sset_in_%s.gsa.out", toupper(set), ds))
tests <- rbind(data.table(set = "thyroid", DS[domain == "liver", .(dataset)]), data.table(set = "liver", DS[domain == "thyroid", .(dataset)]))
runnable <- !is.na(MAGMA) && !is.na(GENE_LOC) && !is.na(BFILE) && all(file.exists(paste0(BFILE, c(".bed", ".bim", ".fam"))))
do_run <- if (RUN_MAGMA == "0") FALSE else if (RUN_MAGMA == "1") TRUE else runnable
if (do_run && !runnable)
  stop("MAGMA, NCBI37.3.gene.loc and the g1000_eur reference (.bed/.bim/.fam) are needed to run MAGMA. Found: magma = ", MAGMA, ", gene.loc = ", GENE_LOC, ", reference = ", BFILE, call. = FALSE)

if (do_run) {
  if (file.access(MAGMA, 1) != 0) stop("MAGMA is not executable: ", MAGMA, " (chmod +x)", call. = FALSE)
  message("MAGMA: ", MAGMA); message("gene.loc: ", GENE_LOC); message("reference: ", BFILE)
  run <- function(args) system2(MAGMA, args = args)

  ## SNP position file and chr:bp -> rsID lookup from the reference .bim
  bim <- fread(paste0(BFILE, ".bim"), header = FALSE, col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  snp_loc <- file.path(OUT_DIR, "snp_loc_from_bim.txt")
  if (!file.exists(snp_loc)) fwrite(bim[, .(SNP, CHR, BP)], snp_loc, sep = "\t", col.names = FALSE)     # MAGMA reads no header here
  bim[, CHR := as.character(CHR)]
  lookup <- bim[!duplicated(bim, by = c("CHR", "BP")), .(CHR, BP, SNP)]                  # one rsID per position (first one)
  setkey(lookup, CHR, BP); rm(bim)

  ## annotation, once for all datasets
  annot_prefix <- file.path(OUT_DIR, "gene_annot")
  if (!file.exists(paste0(annot_prefix, ".genes.annot")))
    run(c("--annotate", "window=10,10", "--snp-loc", shQuote(snp_loc), "--gene-loc", shQuote(GENE_LOC), "--out", shQuote(annot_prefix)))
  if (!file.exists(paste0(annot_prefix, ".genes.annot"))) stop("MAGMA annotation failed (no gene_annot.genes.annot).", call. = FALSE)

  ## rsID of the reference panel by chromosome (one chromosome at a time to limit memory)
  attach_rsid <- function(ss, chr_col, bp_col) {
    ss[, .chr := sub("^chr", "", as.character(get(chr_col)), ignore.case = TRUE)]
    ss[, .bp := as.integer(get(bp_col))]
    out <- lapply(unique(ss$.chr), function(ch) {
      s <- ss[.chr == ch]; l <- lookup[CHR == ch]
      if (!nrow(s) || !nrow(l)) return(NULL)
      merge(s, l, by.x = c(".chr", ".bp"), by.y = c("CHR", "BP"), all = FALSE)
    })
    out <- rbindlist(out, use.names = TRUE, fill = TRUE)
    out[, c(".chr", ".bp") := NULL]
    out
  }

  gene_level <- function(d) {
    prefix <- file.path(OUT_DIR, paste0("geneLevel_", d$dataset))
    if (file.exists(paste0(prefix, ".genes.out"))) { message(sprintf("[%s] gene-level result exists, skipped", d$dataset)); return(prefix) }
    f <- find_sumstats(d$file_id)
    if (is.na(f)) stop(sprintf("[%s] summary statistics not found (looked for a file or folder with '%s' in its name below %s)", d$dataset, d$file_id, paste(ROOTS, collapse = " | ")), call. = FALSE)
    cols <- unique(c(d$chr_col, d$bp_col, d$p_col, d$n_col)); if (!is.na(d$N_fixed)) cols <- setdiff(cols, d$n_col)
    ss <- fread(f, select = cols)
    n0 <- nrow(ss); ss <- attach_rsid(ss, d$chr_col, d$bp_col)
    message(sprintf("[%s] %s: %d of %d rows matched to a reference rsID (%.1f%%)", d$dataset, basename(f), nrow(ss), n0, 100 * nrow(ss) / n0))
    if (!nrow(ss)) stop(sprintf("[%s] no row matched by chromosome and position: check the column names and that the file is GRCh37", d$dataset), call. = FALSE)
    use_cols <- if (is.na(d$N_fixed)) c(d$p_col, d$n_col) else d$p_col
    ss <- ss[stats::complete.cases(ss[, ..use_cols])]                            # an empty field would shift the columns for MAGMA
    tmp <- file.path(OUT_DIR, paste0(d$dataset, "_pval.tmp"))
    fwrite(ss[, c("SNP", use_cols), with = FALSE], tmp, sep = "\t"); rm(ss); invisible(gc(verbose = FALSE))
    pval_args <- if (is.na(d$N_fixed)) c(sprintf("use=SNP,%s", d$p_col), sprintf("ncol=%s", d$n_col)) else c(sprintf("use=SNP,%s", d$p_col), sprintf("N=%d", as.integer(d$N_fixed)))
    st <- run(c("--bfile", shQuote(BFILE), "--pval", shQuote(tmp), shQuote(pval_args), "--gene-annot", shQuote(paste0(annot_prefix, ".genes.annot")), "--out", shQuote(prefix)))
    if (st != 0 || !file.exists(paste0(prefix, ".genes.out"))) stop(sprintf("[%s] MAGMA gene-level analysis failed (exit code %s); see %s.log", d$dataset, st, prefix), call. = FALSE)
    unlink(tmp)
    prefix
  }
  gl <- setNames(lapply(seq_len(nrow(DS)), function(i) gene_level(DS[i])), DS$dataset)

  for (i in seq_len(nrow(tests))) {
    set <- tests$set[i]; ds <- tests$dataset[i]
    out <- sub("\\.gsa\\.out$", "", gsa_file(set, ds))
    run(c("--gene-results", shQuote(paste0(gl[[ds]], ".genes.raw")), "--set-annot", shQuote(SETFILE[[set]]), "--out", shQuote(out)))
    if (!file.exists(gsa_file(set, ds))) stop(sprintf("MAGMA gene-set analysis failed for the %s set in %s", set, ds), call. = FALSE)
  }
} else {
  message("MAGMA is not run (RUN_MAGMA=", RUN_MAGMA, if (!runnable) "; MAGMA, gene.loc or the reference panel not found", "): the table is assembled from existing *.gsa.out files")
}

## ---- table S1 ------------------------------------------------------------------------------------------------------------
read_gsa <- function(f) fread(f, skip = "VARIABLE")
miss <- tests[!mapply(function(s, d) file.exists(gsa_file(s, d)), set, dataset)]
if (nrow(miss))
  stop("No MAGMA result for: ", paste(sprintf("%s set in %s", miss$set, miss$dataset), collapse = "; "),
       "\n  Provide MAGMA, NCBI37.3.gene.loc, g1000_eur and the 12 summary statistics (see data/README.md) or put the gsa_*.gsa.out files into ", OUT_DIR, call. = FALSE)
res <- rbindlist(lapply(seq_len(nrow(tests)), function(i) {
  r <- read_gsa(gsa_file(tests$set[i], tests$dataset[i]))[VARIABLE == SETNAME[[tests$set[i]]]]
  data.table(set = tests$set[i], dataset = tests$dataset[i], type = r$TYPE, n_genes_in_set = r$NGENES, beta = r$BETA, beta_std = r$BETA_STD, SE = r$SE, P = r$P)
}))
res[, gene_set := LABEL_SET[set]]
res[, tested_in := DS$tested_in[match(dataset, DS$dataset)]]
res[, trait_group := fifelse(set == "thyroid", DS$trait_group[match(dataset, DS$dataset)], tested_in)]
res[, significance := fifelse(P < ALPHA_BONF, sprintf("Bonferroni (P < %.4f)", ALPHA_BONF), fifelse(P < 0.05, "nominal (P < 0.05)", "not significant"))]
res[, `:=`(models_in_trait_group = .N, trait_group_min_P = min(P), trait_group_n_nominal = sum(P < 0.05)), by = .(set, trait_group)]
S1 <- res[, .(gene_set, tested_in, trait_group, type, n_genes_in_set, beta, beta_std, SE, P, significance, models_in_trait_group, trait_group_min_P, trait_group_n_nominal)]
fwrite(S1, file.path(OUT_DIR, "magma_competitive_gene_set_results.tsv"), sep = "\t")
cat("\nCompetitive gene-set results (Bonferroni threshold P =", format(ALPHA_BONF, digits = 2), ")\n")
print(S1[, .(gene_set = substr(gene_set, 1, 12), tested_in, n_genes_in_set, beta = signif(beta, 4), SE = signif(SE, 4), P = signif(P, 4), significance)], row.names = FALSE)

## ---- check against Table S1 -----------------------------------------------------------------------------------------------------
if (!is.na(SUPP)) {
  sh <- excel_sheets(SUPP); sh <- sh[sh %in% c("S1", "Table S1")]
  if (length(sh) == 1) {
    top <- read_excel(SUPP, sheet = sh, range = readxl::cell_limits(c(1, 1), c(25, 1)), col_names = FALSE, col_types = "text", .name_repair = "minimal")
    h <- which(trimws(as.character(top[[1]])) == "gene_set")[1]
    ref <- as.data.table(read_excel(SUPP, sheet = sh, range = readxl::cell_limits(c(h, 1), c(NA, 13)), .name_repair = "minimal"))
    ref <- ref[!is.na(gene_set)]
    mm <- merge(S1, ref, by = c("gene_set", "tested_in"), suffixes = c("", ".S1"))
    rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
    cat(sprintf("\nCheck against Table S1 (%s): %d of %d rows matched; n_genes_in_set differs in %d, max relative difference beta %.2g, beta_std %.2g, SE %.2g, P %.2g; significance differs in %d; trait-group columns differ in %d\n",
                basename(SUPP), nrow(mm), nrow(ref), sum(mm$n_genes_in_set != mm$n_genes_in_set.S1), rel(mm$beta, mm$beta.S1), rel(mm$beta_std, mm$beta_std.S1),
                rel(mm$SE, mm$SE.S1), rel(mm$P, mm$P.S1), sum(mm$significance != mm$significance.S1),
                sum(mm$models_in_trait_group != mm$models_in_trait_group.S1 | mm$trait_group_n_nominal != mm$trait_group_n_nominal.S1 | abs(mm$trait_group_min_P - mm$trait_group_min_P.S1) > 1e-3 * mm$trait_group_min_P.S1)))
  }
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "00b_magma_gene_set_enrichment_sessionInfo.txt"))
message("Done. Output: ", OUT_DIR)
