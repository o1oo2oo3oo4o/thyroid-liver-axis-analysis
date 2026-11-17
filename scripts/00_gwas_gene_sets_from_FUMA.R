################################################################################
#  00_gwas_gene_sets_from_FUMA.R
#  GWAS gene mapping, post-processing of FUMA output:  FUMA genes.txt / eqtl.txt / GenomicRiskLoci.txt  ->  cis-eQTL PASS genes of every GWAS dataset
#  ->  gene set of every trait  ->  thyroid-derived set, liver-derived set and shared set.
#
#  FUMA itself (SNP2GENE) is a web tool and is run outside this repository; its settings are recorded in docs/FUMA_settings.md.
#  This script only reads the downloaded FUMA result files.  It replaces the two interactive scripts used during the analysis
#  (FUMA_cis-eQTL_mapping and fuma_pass_union_intersection_integrator): all questions that were asked at the keyboard are now fixed rules.
#
#  RULES (as used for the manuscript)
#    1. genes.txt: keep genes with eqtlMapminQ < 0.05                                   (eQTL-mapped genes, FUMA eQTL mapping)
#    2. the gene must lie on the chromosome of its genomic risk locus                   (genes.txt GenomicLocus -> GenomicRiskLoci.txt chr)
#    3. every eQTL pair of the gene in eqtl.txt must be on the chromosome of the gene   (cis only); a gene without eQTL pair fails
#       -> cis PASS = 1 & 2 & 3.  [In the first, interactive version a gene that failed only because of a missing eQTL pair could leave an
#          empty row in the PASS file; here it is simply a fail, and the number of such genes is printed.]
#    4. trait gene set = union over the datasets of the trait (PDFF: 4 models, cT1: 4 models, FT3: 1, FT4: 1);
#       TSH = INTERSECTION of the two TSH datasets (A and B), i.e. a gene must pass in both
#    5. liver-derived set = PDFF or cT1;  thyroid-derived set = TSH or FT3 or FT4;  shared set = in both
#    Genes are counted by SYMBOL IN UPPER CASE (the key of all later tables).  FUMA lists C6orf165 and C6ORF165 as two different Ensembl genes;
#    after upper-casing they are one gene, as in the sets of the manuscript.  Ensembl IDs that share a symbol (e.g. PIK3R2, SOX7) are merged too.
#
#  INPUTS  (FUMA_DIR, default <DATA_DIR>/FUMA; the files are searched by name)
#    one folder per GWAS dataset, folder name = dataset name, each with  genes.txt, GenomicRiskLoci.txt, eqtl.txt
#      PDFF_650 PDFF_651 PDFF_652 PDFF_653   (GCST90446650-653)      cT1_646 cT1_647 cT1_648 cT1_649   (GCST90446646-649)
#      FT3   FT4   TSH_A (GCST010653)   TSH_B (GCST90296333)
#    A folder is also accepted when its path contains the GWAS Catalog ID.  If your folders are named differently, give a table
#    fuma_datasets.tsv (columns: dataset, folder) with the folder (or part of its path) of every dataset.
#    The FUMA result files themselves are not part of this repository (re-run FUMA with the settings of docs/FUMA_settings.md).
#
#  OUTPUTS (OUT_DIR, default results/00_gwas_gene_sets)
#    fuma_cis_genes_<dataset>.tsv   all genes of the dataset with the result of every rule and the final cis_pass
#    fuma_step_counts.tsv           genes after every rule, per dataset
#    trait_gene_sets.tsv            trait, gene, number of datasets in which the gene passed
#    gwas_gene_flags.tsv            Gene, GWAS_PDFF, GWAS_cT1, GWAS_TSH, GWAS_FT3, GWAS_FT4 (0/1), ensg, GWAS_liver, GWAS_thyroid, GWAS_shared
#                                   -> input of 00b, 00c and 01c (a copy of the flags used for the manuscript is in data/gene_sets/)
#    gene_set_counts.tsv            sizes of the thyroid, liver, shared and union sets
#    00_gwas_gene_sets_sessionInfo.txt
#  Check printed to the console: set sizes against Table S3a and, if data/gene_sets/gwas_gene_flags.tsv exists, gene by gene against that file.
#  EXPECTED with the FUMA results of the manuscript (checked on the real files): thyroid 1,709, liver 1,120, shared 199, union 2,630 gene symbols; per trait
#  PDFF 786, cT1 486, TSH 828, FT3 149, FT4 999 -- identical to data/gene_sets/gwas_gene_flags.tsv.  Table S3a (and the Results text) give 1,119 liver genes
#  and a union of 2,629: FUMA lists ABHD16B and C20ORF135 (an old name of ABHD16B) as two genes and the manuscript counts them as one, so the check
#  against Table S3a printed below shows a difference of 1 in the liver and union rows and 0 in the thyroid and shared rows (docs/PROVENANCE.md, "Gene counts").
#  RESUME=1  reuse fuma_cis_genes_<dataset>.tsv already in OUT_DIR (skips re-reading very large eqtl.txt files, e.g. the 670 MB file of TSH_B).
#
#  USAGE   DATA_DIR=/path/to/data  Rscript 00_gwas_gene_sets_from_FUMA.R
################################################################################
invisible(try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(readxl) })

## ---- settings ----------------------------------------------------------------
DATA_DIR <- Sys.getenv("DATA_DIR", "data")
FUMA_DIR <- Sys.getenv("FUMA_DIR", file.path(DATA_DIR, "FUMA"))
OUT_DIR  <- Sys.getenv("OUT_DIR", file.path("results", "00_gwas_gene_sets"))
REF_FLAGS <- Sys.getenv("REF_FLAGS", file.path(DATA_DIR, "gene_sets", "gwas_gene_flags.tsv"))
Q_THRESH <- 0.05
RESUME   <- Sys.getenv("RESUME", "0") == "1"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DATASETS <- data.table(
  dataset      = c(paste0("PDFF_", 650:653), paste0("cT1_", 646:649), "FT3", "FT4", "TSH_A", "TSH_B"),
  trait        = c(rep("PDFF", 4), rep("cT1", 4), "FT3", "FT4", "TSH", "TSH"),
  gwas_catalog = c(paste0("GCST90446", 650:653), paste0("GCST90446", 646:649), NA, NA, "GCST010653", "GCST90296333"))

## ---- helpers -------------------------------------------------------------------
norm_chr <- function(x) gsub("^chr", "", as.character(x), ignore.case = TRUE)
strip_version <- function(x) sub("\\..*$", "", as.character(x))
read_fuma <- function(f) fread(f, sep = "\t", quote = "", header = TRUE, data.table = TRUE, fill = TRUE, check.names = FALSE)
need <- function(d, cols, what) { m <- setdiff(cols, names(d)); if (length(m)) stop("Missing column(s) in ", what, ": ", paste(m, collapse = ", "), call. = FALSE) }

## ---- locate the FUMA folders ---------------------------------------------------------
if (!dir.exists(FUMA_DIR))
  stop("FUMA_DIR not found: ", FUMA_DIR, "\n  Put the FUMA result folders (genes.txt, eqtl.txt, GenomicRiskLoci.txt) of the 12 GWAS datasets there ",
       "(see header of this script and docs/FUMA_settings.md), or set FUMA_DIR.", call. = FALSE)
gfiles <- list.files(FUMA_DIR, pattern = "^genes\\.txt$", recursive = TRUE, full.names = TRUE)
fdirs  <- unique(dirname(gfiles))
fdirs  <- fdirs[file.exists(file.path(fdirs, "GenomicRiskLoci.txt")) & (RESUME | file.exists(file.path(fdirs, "eqtl.txt")))]
if (!length(fdirs)) stop("No folder with genes.txt + GenomicRiskLoci.txt + eqtl.txt under ", FUMA_DIR, call. = FALSE)
relp <- substring(fdirs, nchar(FUMA_DIR) + 2)

map_file <- list.files(c(FUMA_DIR, DATA_DIR), pattern = "^fuma_datasets\\.tsv$", full.names = TRUE)[1]
user_map <- if (!is.na(map_file)) fread(map_file) else NULL
if (!is.null(user_map)) need(user_map, c("dataset", "folder"), "fuma_datasets.tsv")

pick_dir <- function(ds, gcst) {
  if (!is.null(user_map) && ds %in% user_map$dataset) {
    pat <- user_map$folder[match(ds, user_map$dataset)]
    hit <- which(grepl(pat, relp, fixed = TRUE))
  } else {
    comps <- strsplit(relp, "[/\\\\]")
    hit <- which(vapply(comps, function(cc) any(tolower(cc) == tolower(ds)), logical(1)))
    if (!length(hit) && !is.na(gcst)) hit <- which(grepl(gcst, relp, fixed = TRUE))
  }
  if (length(hit) > 1) stop("More than one FUMA folder matches dataset ", ds, ": ", paste(relp[hit], collapse = " | "),
                            "\n  Use fuma_datasets.tsv to say which one to use.", call. = FALSE)
  if (length(hit) == 0) NA_character_ else fdirs[hit]
}
DATASETS[, fuma_folder := mapply(pick_dir, dataset, gwas_catalog)]
if (anyNA(DATASETS$fuma_folder))
  stop("No FUMA result folder found for: ", paste(DATASETS$dataset[is.na(DATASETS$fuma_folder)], collapse = ", "),
       "\n  All 12 datasets are needed (TSH is the intersection of A and B).\n  Folders found: ", paste(relp, collapse = " | "), call. = FALSE)

## ---- rules 1-3 per dataset -----------------------------------------------------------
process_dataset <- function(ds, dir) {
  genes <- read_fuma(file.path(dir, "genes.txt"))
  loci  <- read_fuma(file.path(dir, "GenomicRiskLoci.txt"))
  if (!file.exists(file.path(dir, "eqtl.txt"))) stop("eqtl.txt not found in ", dir, call. = FALSE)
  eqtl  <- read_fuma(file.path(dir, "eqtl.txt"))
  for (a in list(c("gene", "Gene"), c("Symbol", "symbol"), c("p", "P")))
    if (!(a[2] %in% names(eqtl)) && (a[1] %in% names(eqtl))) setnames(eqtl, a[1], a[2])
  need(genes, c("ensg", "symbol", "chr", "start", "end", "GenomicLocus", "eqtlMapminQ"), paste(ds, "genes.txt"))
  need(loci,  c("GenomicLocus", "chr"), paste(ds, "GenomicRiskLoci.txt"))
  need(eqtl,  c("Gene", "chr", "FDR", "db", "tissue"), paste(ds, "eqtl.txt"))
  genes[, ensg := strip_version(ensg)]; genes[, gene_chr := norm_chr(chr)]
  loci[, locus_chr := norm_chr(chr)]
  eqtl[, Gene := strip_version(Gene)]; eqtl[, eqtl_chr := norm_chr(chr)]

  genes[, rule1_q := !is.na(eqtlMapminQ) & eqtlMapminQ < Q_THRESH]
  lk <- unique(loci[, .(GenomicLocus, locus_chr)], by = "GenomicLocus")
  genes[, locus_chr := lk$locus_chr[match(GenomicLocus, lk$GenomicLocus)]]
  genes[, rule2_chr_match := !is.na(gene_chr) & !is.na(locus_chr) & gene_chr == locus_chr]

  gchr <- unique(genes[, .(Gene = ensg, gene_chr)], by = "Gene")
  eqtl[, gene_chr := gchr$gene_chr[match(Gene, gchr$Gene)]]
  es <- eqtl[, .(eqtl_pairs_n = .N,
                 eqtl_chr_unique = paste(sort(unique(eqtl_chr)), collapse = ";"),
                 eqtl_minFDR = suppressWarnings(min(FDR, na.rm = TRUE)),
                 eqtl_db_set = paste(sort(unique(db)), collapse = ";"),
                 eqtl_tissue_n = uniqueN(tissue[!is.na(tissue) & tissue != ""]),
                 rule3_eqtl_on_gene_chr = !is.na(gene_chr[1]) && all(eqtl_chr == gene_chr[1])), by = Gene]
  setnames(es, "Gene", "ensg")
  genes <- merge(genes, es, by = "ensg", all.x = TRUE, sort = FALSE)
  genes[is.na(eqtl_pairs_n), eqtl_pairs_n := 0L]
  genes[, rule3_eqtl_on_gene_chr := rule3_eqtl_on_gene_chr %in% TRUE]
  genes[, cis_pass := rule1_q & rule2_chr_match & rule3_eqtl_on_gene_chr]

  cnt <- count_steps(ds, genes)
  keep <- intersect(c("ensg", "symbol", "type", "pLI", "gene_chr", "start", "end", "GenomicLocus", "locus_chr", "eqtlMapminP", "eqtlMapminQ",
                      "rule1_q", "rule2_chr_match", "eqtl_pairs_n", "eqtl_chr_unique", "eqtl_minFDR", "eqtl_db_set", "eqtl_tissue_n",
                      "rule3_eqtl_on_gene_chr", "cis_pass"), names(genes))
  list(genes = genes[, ..keep], counts = cnt)
}

count_steps <- function(ds, genes)
  data.table(dataset = ds,
             genes_in_genes_txt   = uniqueN(genes$ensg),
             rule1_eqtlMapminQ    = uniqueN(genes[rule1_q == TRUE, ensg]),
             rule2_chr_match      = uniqueN(genes[rule1_q & rule2_chr_match, ensg]),
             rule3_eqtl_on_chr    = uniqueN(genes[rule1_q & rule2_chr_match & rule3_eqtl_on_gene_chr, ensg]),
             cis_pass             = uniqueN(genes[cis_pass == TRUE, ensg]),
             failed_only_no_eqtl_pair = uniqueN(genes[rule1_q & rule2_chr_match & eqtl_pairs_n == 0, ensg]))
get_dataset <- function(ds, dir) {
  f <- file.path(OUT_DIR, sprintf("fuma_cis_genes_%s.tsv", ds))
  if (RESUME && file.exists(f)) {
    message("  ", ds, ": RESUME, using ", f)
    g <- fread(f); return(list(genes = g, counts = count_steps(ds, g)))
  }
  message("  ", ds, ": reading ", dir)
  process_dataset(ds, dir)
}
message(sprintf("FUMA folders: %d datasets found under %s", nrow(DATASETS), FUMA_DIR))
res <- Map(get_dataset, DATASETS$dataset, DATASETS$fuma_folder)
counts <- rbindlist(lapply(res, `[[`, "counts"))
counts[, trait := DATASETS$trait[match(dataset, DATASETS$dataset)]]
setcolorder(counts, c("dataset", "trait"))
fwrite(counts, file.path(OUT_DIR, "fuma_step_counts.tsv"), sep = "\t")
for (ds in names(res)) fwrite(res[[ds]]$genes, file.path(OUT_DIR, sprintf("fuma_cis_genes_%s.tsv", ds)), sep = "\t", na = "")
print(counts, row.names = FALSE)
if (any(counts$failed_only_no_eqtl_pair > 0))
  message("Note: ", sum(counts$failed_only_no_eqtl_pair), " gene-dataset combinations failed only because no eQTL pair of the gene is in eqtl.txt (counted as fail).")

## ---- rules 4-5: trait sets, thyroid / liver / shared -------------------------------------
pass_long <- rbindlist(lapply(names(res), function(ds) {
  g <- res[[ds]]$genes[cis_pass == TRUE, .(ensg, symbol)]
  if (nrow(g)) g[, `:=`(dataset = ds, trait = DATASETS$trait[match(ds, DATASETS$dataset)])]
  g
}), fill = TRUE)
trait_sets <- pass_long[, .(symbol = symbol[1], n_datasets_passed = uniqueN(dataset), datasets = paste(sort(unique(dataset)), collapse = ",")), by = .(trait, ensg)]
n_per_trait <- DATASETS[, .(n_datasets = .N), by = trait]
trait_sets[, n_datasets_of_trait := n_per_trait$n_datasets[match(trait, n_per_trait$trait)]]
trait_sets[, in_trait_set := fifelse(trait == "TSH", n_datasets_passed == n_datasets_of_trait, n_datasets_passed >= 1L)]   # TSH = A and B
setorder(trait_sets, trait, -n_datasets_passed, symbol)
fwrite(trait_sets, file.path(OUT_DIR, "trait_gene_sets.tsv"), sep = "\t")

TRAITS <- c("PDFF", "cT1", "TSH", "FT3", "FT4")
sets <- lapply(TRAITS, function(tr) trait_sets[trait == tr & in_trait_set == TRUE])
names(sets) <- TRAITS
message(sprintf("Trait sets (genes by Ensembl ID): %s", paste(sprintf("%s %d", TRAITS, vapply(sets, nrow, 1L)), collapse = ", ")))

sym_all <- unique(rbindlist(lapply(sets, function(s) s[, .(ensg, symbol)])), by = "ensg")
flags <- data.table(ensg = sym_all$ensg, Gene = toupper(sym_all$symbol), symbol_fuma = sym_all$symbol)
for (tr in TRAITS) set(flags, j = paste0("GWAS_", tr), value = as.integer(flags$ensg %in% sets[[tr]]$ensg))
flags[, GWAS_liver := as.integer(GWAS_PDFF == 1 | GWAS_cT1 == 1)]
flags[, GWAS_thyroid := as.integer(GWAS_TSH == 1 | GWAS_FT3 == 1 | GWAS_FT4 == 1)]
flags[, GWAS_shared := as.integer(GWAS_liver == 1 & GWAS_thyroid == 1)]
fcols <- paste0("GWAS_", TRAITS)
## one row per gene symbol in upper case (the key of all later tables); Ensembl IDs that share a symbol (or differ only in case) are merged
by_sym <- flags[, c(list(symbol_fuma = paste(unique(symbol_fuma), collapse = ";"), ensg = paste(unique(ensg), collapse = ";")),
                    lapply(.SD, max)), by = Gene, .SDcols = c(fcols, "GWAS_liver", "GWAS_thyroid", "GWAS_shared")]
setorder(by_sym, Gene)
setcolorder(by_sym, c("Gene", fcols, "ensg", "symbol_fuma", "GWAS_liver", "GWAS_thyroid", "GWAS_shared"))
fwrite(by_sym, file.path(OUT_DIR, "gwas_gene_flags.tsv"), sep = "\t")

set_counts <- data.table(set = c("thyroid-derived (TSH or FT3 or FT4)", "liver-derived (PDFF or cT1)", "shared", "union"),
                         genes = c(sum(by_sym$GWAS_thyroid), sum(by_sym$GWAS_liver), sum(by_sym$GWAS_shared), nrow(by_sym)))
fwrite(set_counts, file.path(OUT_DIR, "gene_set_counts.tsv"), sep = "\t")
cat("\nGene sets (by gene symbol)\n"); print(set_counts, row.names = FALSE)
cat("Per trait: ", paste(sprintf("%s %d", TRAITS, colSums(by_sym[, ..fcols])), collapse = ", "), "\n", sep = "")

## ---- checks ------------------------------------------------------------------------------
## (a) Table S3a, when the Supplementary Tables workbook is found
SUPP <- { e <- Sys.getenv("SUPP_XLSX", ""); if (nzchar(e)) e else { h <- list.files(c(DATA_DIR, "data", getwd()), pattern = "Supplementary_Tables_S1-S13.*\\.xlsx$", recursive = TRUE, full.names = TRUE); h <- h[!grepl("^~\\$", basename(h))]; if (length(h)) h[1] else NA_character_ } }
if (!is.na(SUPP)) {
  sh <- excel_sheets(SUPP); sh <- sh[sh %in% c("S3", "Table S3")]
  if (length(sh) == 1) {
    s3 <- as.data.table(read_excel(SUPP, sheet = sh, range = readxl::cell_limits(c(4, 1), c(12, 3)), .name_repair = "minimal"))
    s3 <- s3[!is.na(item) & !is.na(genes)]
    exp <- c(thyroid = s3[item == "Thyroid-derived GWAS gene set", genes], liver = s3[item == "Liver-derived GWAS gene set", genes],
             shared = s3[item == "Shared (in both sets)", genes], union = s3[item == "Union of the two sets", genes])
    obs <- c(thyroid = set_counts$genes[1], liver = set_counts$genes[2], shared = set_counts$genes[3], union = set_counts$genes[4])
    cat("\nCheck against Table S3a (", basename(SUPP), ")\n", sep = "")
    print(data.table(set = names(obs), computed = as.integer(obs), table_S3a = as.integer(exp), difference = as.integer(obs - exp)), row.names = FALSE)
  }
}
## (b) the flags file used for the manuscript
if (file.exists(REF_FLAGS) && normalizePath(REF_FLAGS) != normalizePath(file.path(OUT_DIR, "gwas_gene_flags.tsv"), mustWork = FALSE)) {
  ref <- fread(REF_FLAGS)
  cat("\nCheck against ", REF_FLAGS, " (", nrow(ref), " genes)\n", sep = "")
  mm <- merge(by_sym[, c("Gene", fcols), with = FALSE], ref[, c("Gene", fcols), with = FALSE], by = "Gene", all = TRUE, suffixes = c("", ".ref"))
  only_new <- mm[is.na(GWAS_PDFF.ref), Gene]; only_ref <- mm[is.na(GWAS_PDFF), Gene]
  both <- mm[!is.na(GWAS_PDFF) & !is.na(GWAS_PDFF.ref)]
  diff_n <- sum(Reduce(`|`, lapply(fcols, function(cn) both[[cn]] != both[[paste0(cn, ".ref")]])))
  cat(sprintf("  genes only in the new run: %d, only in the reference: %d, genes with a different flag: %d\n", length(only_new), length(only_ref), diff_n))
  if (length(only_new)) cat("  only new: ", paste(head(only_new, 20), collapse = ", "), if (length(only_new) > 20) " ...", "\n", sep = "")
  if (length(only_ref)) cat("  only reference: ", paste(head(only_ref, 20), collapse = ", "), if (length(only_ref) > 20) " ...", "\n", sep = "")
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "00_gwas_gene_sets_sessionInfo.txt"))
message("Done. Output: ", normalizePath(OUT_DIR))
