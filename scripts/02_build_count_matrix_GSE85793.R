################################################################################
#  02_build_count_matrix_GSE85793.R
#  Builds the raw count matrix of the public mouse liver T3 dataset GSE85793 from the per-sample HTSeq count files
#  supplied by GEO (GSE85793_RAW), and the C/T subset used for the DESeq2 analysis.
#
#  NOTE  The original matrix (0_count_matrix.csv) was assembled before the analysis scripts were written and the assembling
#        step was not kept as a script. This script was re-written for the repository. It was checked against the original files:
#        the columns C1 and T1 of 0_count_matrix.csv are identical to the GEO files (same genes, same order, same counts), and
#        1_cleaned_count_matrix.csv is the C1-C4/T1-T4 subset of 0_count_matrix.csv (same rows, same order).
#
#  INPUT  (DATA_DIR/GSE85793_RAW_file_gunzip/)
#    GSM*_Solid5500_2015_02_05_htseq_<label>.csv : one file per sample (gunzipped GEO supplementary files); tab-separated,
#                                                  header "Symbol<TAB><label>", one row per gene. <label> = C1..C4 (vehicle), G*, T1..T4 (T3), GT*, CP*
#  OUTPUT (OUT_DIR, default = DATA_DIR because later steps read them from there)
#    0_count_matrix.csv          : all samples; gene x sample
#    1_cleaned_count_matrix.csv  : only the vehicle (C1-C4) and T3 (T1-T4) samples, used by 03_deseq2_mouse_liver_T3.R
#    1_cleaned_coldata.csv       : sample, condition (Control / T3)
#  Study design of GSE85793 (4 mice per group): vehicle (C), glucagon (G), T3 (T), glucagon + T3 (GT), chemical hybrid (CP).
#  Only vehicle vs T3 is analysed in the manuscript.
#
#  USAGE  DATA_DIR=/path/to/data Rscript 02_build_count_matrix_GSE85793.R
################################################################################
DATA_DIR <- Sys.getenv("DATA_DIR", "data")
RAW_DIR  <- file.path(DATA_DIR, "GSE85793_RAW_file_gunzip")
OUT_DIR  <- Sys.getenv("OUT_DIR", DATA_DIR)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

files <- list.files(RAW_DIR, pattern = "^GSM[0-9]+_.*_htseq_[A-Za-z0-9]+\\.csv$", full.names = TRUE)
if (!length(files)) stop("No per-sample count files found in ", RAW_DIR)
label <- sub("^.*_htseq_([A-Za-z0-9]+)\\.csv$", "\\1", basename(files))
read_one <- function(f) {
  x <- read.delim(f, header = TRUE, sep = "\t", quote = "", check.names = FALSE, colClasses = c("character", "integer"))
  setNames(x[[2]], x[[1]])
}
vecs <- lapply(files, read_one); names(vecs) <- label
genes <- names(vecs[[1]])
if (!all(vapply(vecs, function(v) identical(names(v), genes), TRUE))) stop("The gene lists of the sample files differ.")
if (anyDuplicated(genes)) stop("Duplicated gene symbols in the count files.")
mat <- do.call(cbind, vecs)
# column order of the original matrix: C1-4, G1-4, T1-4, GT1-4, CP1-4 (all labels present are kept)
ord <- c(paste0("C", 1:4), paste0("G", 1:4), paste0("T", 1:4), paste0("GT", 1:4), paste0("CP", 1:4))
mat <- mat[, c(intersect(ord, colnames(mat)), setdiff(colnames(mat), ord)), drop = FALSE]
write.csv(data.frame(gene = rownames(mat), mat, check.names = FALSE), file.path(OUT_DIR, "0_count_matrix.csv"), row.names = FALSE, quote = FALSE)

ct <- paste0(rep(c("C", "T"), each = 4), 1:4)
if (all(ct %in% colnames(mat))) {
  write.csv(data.frame(gene = rownames(mat), mat[, ct], check.names = FALSE), file.path(OUT_DIR, "1_cleaned_count_matrix.csv"), row.names = FALSE, quote = FALSE)
  write.csv(data.frame(sample = ct, condition = rep(c("Control", "T3"), each = 4)), file.path(OUT_DIR, "1_cleaned_coldata.csv"), row.names = FALSE, quote = FALSE)
} else message("C1-C4/T1-T4 are not all present; the cleaned files were not written.")
cat(sprintf("Samples: %s\nGenes: %d\nWritten to %s\n", paste(colnames(mat), collapse = ", "), nrow(mat), OUT_DIR))
