################################################################################
#  04_reortholog_synonym.R
#  Mouse -> human ortholog remapping of the liver-T3 DEGs, with optional resolution of former mouse gene names.
#  (Original: 261021_reortholog_synonym.R. The redox analysis is not part of this script.)
#
#  PURPOSE
#    Re-map the mouse genes of the liver-T3 differential-expression result (GSE85793) to human orthologues
#    with ONE Ensembl BioMart export, and record the reason why a gene could not be mapped.
#    With SYNONYM_FILE <- NA the script reproduces the primary analysis (448 human entries, Table S6a);
#    with the synonym table it reproduces the former-name sensitivity analysis (498 human entries).
#
#  INPUTS (all read from DATA_DIR, file names as they were used in the study)
#    ENSEMBL_FILE : Ensembl BioMart export (mouse dataset, GRCm39, release 116), columns:
#                   Gene stable ID, Gene stable ID version, Gene name, Human gene stable ID,
#                   Human gene name, Human homology type, Human orthology confidence
#    SYNONYM_FILE : (optional) BioMart table "Gene stable ID" + "Gene Synonym" (mouse). NA skips the name-resolution step.
#    DEG_FILE     : DESeq2 result of all genes (12,135 genes; log2FoldChange = apeglm-shrunken value)
#                   produced by 03_deseq2_mouse_liver_T3.R
#    RAW_COUNT    : count matrix before filtering (13,319 genes; first column = gene name)
#    OLD_HUMAN    : (optional) human responsive gene list of an earlier mapping, only for comparison
#
#  MAPPING RULES (these affect the interpretation and must be stated in the Methods)
#    - The mouse gene symbol is matched exactly to the Ensembl "Gene name".
#    - [name resolution] a gene whose name is not in Ensembl (M1) is looked up in the Ensembl synonym table (former names included).
#        accepted: the synonym maps to exactly one Ensembl gene AND the current name of that gene is not separately in the count matrix
#        rejected: no synonym / synonym maps to several genes / the current name is already in the count matrix (conflict)
#        if two or more DEG names resolve to the same current gene they are not treated as one-to-one (M3)
#    - Reason codes
#        M0  one-to-one: one human gene and homology type ortholog_one2one
#        M3  not one-to-one: orthologue exists but one2many / many2many or two or more human genes
#        M2  in Ensembl but no human orthologue assigned
#        M1  name not in Ensembl (former name, predicted name, ...)
#    - The homology type of the human side is the one computed by Ensembl (one2one = one-to-one in both directions).
#    - strict = one2one only / all = one2many and many2many included
#
#  OUTPUTS (OUT_DIR)
#    261021_ortholog_remap_synonym.xlsx (sheets described in README_QC), lists/*.txt,
#    liver_T3_human_responsive_genes_<n>.txt / liver_T3_human_universe_tested_<n>.txt (foreground and background of the
#    pathway analysis, including strict versions), settings.tsv, input_files_manifest.tsv (md5 of the inputs),
#    261021_remap_run_log.txt, sessionInfo.txt
#
#  USAGE   Rscript 04_reortholog_synonym.R        (DATA_DIR and OUT_DIR may be set as environment variables)
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

# ==============================================================================
# [SECTION 1]  Settings
# ==============================================================================
DATA_DIR     <- Sys.getenv("DATA_DIR", "data")                      # folder with the input files
ENSEMBL_FILE <- file.path(DATA_DIR, "261019_ensembl_GRCm39.txt")
DEG_FILE     <- file.path(DATA_DIR, "261020__DEG_all_genes_with_stat.csv")
RAW_COUNT    <- file.path(DATA_DIR, "1_raw_count_matrix_liver_T3.csv")
OLD_HUMAN    <- NA                                                  # optional: earlier human responsive gene list (e.g. 5_liver-T3_sig_merged.csv); NA = no comparison
PADJ_CUT <- 0.05
LFC_CUT  <- 1
SYNONYM_FILE <- file.path(DATA_DIR, "261021_mouse_synonym.txt")     # mouse former-name (synonym) table; NA skips name resolution (= primary 448-gene analysis)
ENSEMBL_RELEASE <- "116"          # Ensembl release of the downloaded tables (used in a status label; the download date is recorded in the Methods)
OUT_DIR  <- Sys.getenv("OUT_DIR", "results/04_ortholog_remap")

dir.create(file.path(OUT_DIR, "lists"), recursive = TRUE, showWarnings = FALSE)
OUT_XLSX <- file.path(OUT_DIR, "261021_ortholog_remap_synonym.xlsx")
STATUS_NO_ORTH <- sprintf("no human orthologue assigned in Ensembl release %s", ENSEMBL_RELEASE)

hdr_style <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", halign = "center", border = "Bottom")
add_sheet <- function(wb, name, df) {
  addWorksheet(wb, name); df <- as.data.frame(df)
  writeData(wb, name, df, headerStyle = hdr_style); freezePane(wb, name, firstRow = TRUE)
  if (ncol(df)) setColWidths(wb, name, cols = seq_len(ncol(df)),
                             widths = vapply(seq_len(ncol(df)), function(j) min(max(nchar(c(names(df)[j], as.character(head(df[[j]], 200)))), na.rm = TRUE) + 2, 45), numeric(1)))
}
qc <- data.table(item = character(), value = character())
qadd <- function(a, b) qc <<- rbind(qc, data.table(item = a, value = as.character(b)))

# ==============================================================================
# [SECTION 2]  Read inputs
# ==============================================================================
RUN_START <- Sys.time()
con_log <- file(file.path(OUT_DIR, "261021_remap_run_log.txt"), open = "wt", encoding = "UTF-8"); sink(con_log, split = TRUE)
status <- tryCatch({ withCallingHandlers({
cat(sprintf("Start: %s | working directory: %s\n", format(RUN_START, "%Y-%m-%d %H:%M:%S"), getwd()))
cat("\n[1] Reading inputs\n")
e <- data.table::fread(ENSEMBL_FILE, sep = "\t", quote = "", colClasses = "character", na.strings = c("", "NA"))
pick <- function(pat) { i <- grep(pat, names(e), ignore.case = TRUE)[1]; if (is.na(i)) stop("Column not found in the Ensembl file: ", pat); names(e)[i] }
setnames(e, c(pick("^Gene stable ID$"), pick("^Gene name$"), pick("^Human gene stable ID$"),
              pick("^Human gene name$"), pick("^Human homology type$"), pick("^Human orthology confidence")),
         c("mid", "mname", "hid", "hname", "htype", "hconf"))
e <- e[, .(mid, mname, hid, hname, htype, hconf)]
cat(sprintf("   Ensembl table: %d rows | %d mouse gene IDs | %d rows with a human orthologue\n", nrow(e), uniqueN(e$mid), sum(!is.na(e$hname))))
qadd("Ensembl table rows / mouse gene IDs / rows with a human orthologue", sprintf("%d / %d / %d", nrow(e), uniqueN(e$mid), sum(!is.na(e$hname))))
qadd("homology type distribution (rows)", paste(names(table(e$htype)), as.integer(table(e$htype)), sep = "=", collapse = ", "))

deg <- data.table::fread(DEG_FILE, data.table = TRUE)
need <- c("gene", "log2FoldChange", "padj")
if (!all(need %in% names(deg))) stop("Required columns missing in the DEG file: ", paste(setdiff(need, names(deg)), collapse = ", "))
deg[, responsive := !is.na(padj) & padj < PADJ_CUT & abs(log2FoldChange) >= LFC_CUT]
raw_genes <- unique(as.character(data.table::fread(RAW_COUNT, select = 1, data.table = TRUE)[[1]]))
cat(sprintf("   count matrix: %d genes | DEG results: %d genes | with padj: %d | responsive: %d\n", length(raw_genes), nrow(deg), sum(!is.na(deg$padj)), sum(deg$responsive)))
qadd("count matrix / DEG results / with padj (tested) / responsive genes", sprintf("%d / %d / %d / %d", length(raw_genes), nrow(deg), sum(!is.na(deg$padj)), sum(deg$responsive)))

# ==============================================================================
# [SECTION 3]  Reason code per mouse gene (all 13,319 genes)
# ==============================================================================
cat("\n[2] Reason code per mouse gene\n")
ep    <- e[!is.na(mname)]
pairs <- unique(ep[!is.na(hname), .(mname, mid, hname, hid, htype, hconf)])
ids   <- ep[, .(n_ensembl_ids = uniqueN(mid)), by = mname]
agg   <- pairs[, .(n_human = uniqueN(hname), human_genes = paste(sort(unique(hname)), collapse = ";"),
                   homology_types = paste(sort(unique(htype)), collapse = ";"),
                   low_confidence = any(hconf == "0")), by = mname]

ms <- data.table(gene = raw_genes)
ms[, map_name := fifelse(gene %in% ep$mname, gene, NA_character_)]
ms[, syn_status := fifelse(!is.na(map_name), "name match", NA_character_)]
ms[, syn_old_ids := NA_character_]
if (!is.na(SYNONYM_FILE) && file.exists(SYNONYM_FILE)) {
  sy <- data.table::fread(SYNONYM_FILE, sep = "\t", select = c("Gene stable ID", "Gene Synonym"), colClasses = "character", na.strings = c("", "NA"))
  setnames(sy, c("mid", "syn")); sy <- unique(sy[!is.na(syn)]); sy[, key := toupper(syn)]
  sy_by_key <- split(sy$mid, sy$key)
  cur_name <- unique(ep[, .(mid, mname)])[!duplicated(mid)]
  cm_keys <- toupper(raw_genes)
  for (i in which(is.na(ms$map_name))) {
    ii <- unique(sy_by_key[[toupper(ms$gene[i])]])
    if (is.null(ii) || length(ii) == 0) { set(ms, i, "syn_status", "no synonym"); next }
    if (length(ii) > 1) { set(ms, i, "syn_status", "synonym maps to several genes"); set(ms, i, "syn_old_ids", paste(ii, collapse = ";")); next }
    cn <- cur_name$mname[match(ii, cur_name$mid)]
    if (is.na(cn)) { set(ms, i, "syn_status", "no synonym"); next }
    if (toupper(cn) %in% cm_keys) { set(ms, i, "syn_status", "conflict (current name already in the count matrix)"); set(ms, i, "syn_old_ids", ii); next }
    set(ms, i, "map_name", cn); set(ms, i, "syn_status", "resolved by synonym (unique)"); set(ms, i, "syn_old_ids", ii)
  }
} else cat("   [note] SYNONYM_FILE not found; name-resolution step skipped (reproduces the primary 448-gene analysis).\n")
ms[, dup_target := !is.na(map_name) & .N > 1L, by = map_name]
n_res <- sum(ms$syn_status == "resolved by synonym (unique)", na.rm = TRUE); m1_all <- sum(ms$syn_status != "name match", na.rm = TRUE)
cat(sprintf("   Of %d genes whose name is not in Ensembl, %d were resolved by synonym\n", m1_all, n_res))
qadd("Name resolution: genes whose name is not in Ensembl / resolved (unique)", sprintf("%d / %d", m1_all, n_res))
qadd("Name-resolution status distribution (all)", paste(names(table(ms$syn_status)), as.integer(table(ms$syn_status)), sep = "=", collapse = ", "))
ms[, in_ensembl := !is.na(map_name)]
ms[, n_ensembl_ids := ids$n_ensembl_ids[match(map_name, ids$mname)]]
j <- match(ms$map_name, agg$mname)
ms[, `:=`(n_human = agg$n_human[j], human_genes = agg$human_genes[j], homology_types = agg$homology_types[j], low_confidence = agg$low_confidence[j])]
ms[, reason := fifelse(!in_ensembl, "M1 name not in Ensembl",
                       fifelse(is.na(n_human), "M2 no human orthologue",
                               fifelse(n_human == 1L & homology_types == "ortholog_one2one", "M0 one-to-one", "M3 not one-to-one")))]
ms[dup_target == TRUE & reason == "M0 one-to-one", reason := "M3 not one-to-one"]
k <- match(ms$gene, deg$gene)
ms[, `:=`(log2FC = deg$log2FoldChange[k], padj = deg$padj[k], pvalue = if ("pvalue" %in% names(deg)) deg$pvalue[k] else NA_real_)]
ms[, deg_status := fifelse(is.na(k), "excluded by pre-filter", fifelse(is.na(padj), "padj NA", "tested"))]
ms[, responsive := !is.na(padj) & padj < PADJ_CUT & abs(log2FC) >= LFC_CUT]
ms[, `:=`(in_ensembl = NULL)]

reason_tab <- rbind(
  ms[, .(scope = "count matrix (all)", .N), by = reason],
  ms[deg_status == "tested", .(scope = "tested (padj available)", .N), by = reason],
  ms[responsive == TRUE, .(scope = "responsive genes", .N), by = reason])
reason_wide <- dcast(reason_tab, reason ~ scope, value.var = "N", fill = 0L)
setcolorder(reason_wide, c("reason", "count matrix (all)", "tested (padj available)", "responsive genes"))
reason_wide <- reason_wide[order(reason)]
print(reason_wide)
for (i in seq_len(nrow(reason_wide)))
  qadd(paste0("reason code ", reason_wide$reason[i], " (all / tested / responsive)"),
       sprintf("%d / %d / %d", reason_wide$`count matrix (all)`[i], reason_wide$`tested (padj available)`[i], reason_wide$`responsive genes`[i]))
qadd("Mouse genes whose symbol has two or more Ensembl gene IDs (within the count matrix)", sum(ms$n_ensembl_ids > 1, na.rm = TRUE))

# Detail of M1 (name not found): predicted/clone-type names vs others (possible former names)
pred_pat <- "Rik$|^Gm[0-9]|^AK[0-9]|^BC[0-9]|^AI[0-9]|^LOC|^CT[0-9]|^A[0-9]{6}"
ms[, reason_detail := fifelse(reason == "M1 name not in Ensembl",
                              paste0(fifelse(grepl(pred_pat, gene), "predicted/clone-type name; ", ""), "synonym: ", syn_status),
                              fifelse(reason == "M2 no human orthologue",
                                      fifelse(grepl(pred_pat, gene), "predicted/clone-type name",
                                              fifelse(grepl("-ps[0-9]*$|^Snhg|^Mir|os$", gene), "pseudogene/non-coding-type name", "other (check whether the export lacks it)")),
                                      NA_character_))]

# ---- Status scheme used in the co-author's analysis (Table S7) ----
#   Order of evaluation: identifier conversion failed -> (no human orthologue assigned) -> ambiguous -> excluded for low expression -> tested normally
#   ("absent from the count matrix" is a status of the human-side audit (05_ortholog_audit.R) and does not occur here)
ms[, status := fifelse(reason == "M1 name not in Ensembl", "identifier conversion failed",
                       fifelse(reason == "M2 no human orthologue", STATUS_NO_ORTH,
                               fifelse(reason == "M3 not one-to-one", "ambiguous",
                                       fifelse(deg_status == "tested", "tested normally", "excluded for low expression"))))]
ms[, status_detail := fifelse(status %in% c("identifier conversion failed", STATUS_NO_ORTH), reason_detail,
                              fifelse(status == "ambiguous", paste0(homology_types, fifelse(n_human > 1L, paste0(" (", n_human, " human genes)"), "")),
                                      fifelse(status == "excluded for low expression",
                                              fifelse(deg_status == "excluded by pre-filter", "excluded by the pre-filter (count >= 10 in at least half of the samples)",
                                                      fifelse(is.na(pvalue), "padj NA (Cook's outlier)", "padj NA (independent filtering)")),
                                              fifelse(responsive, "responsive",
                                                      fifelse(!is.na(padj) & padj < PADJ_CUT, "below response threshold (padj significant only, |log2FC| < cutoff)", "non-responsive")))))]
ms[dup_target == TRUE & status == "ambiguous", status_detail := paste0(status_detail, "; several DEG names resolve to the same Ensembl gene")]
st_levels <- c("identifier conversion failed", STATUS_NO_ORTH, "ambiguous", "excluded for low expression", "tested normally")
st_long <- rbind(
  ms[, .(scope = "count matrix (all)", .N), by = status],
  ms[deg_status == "tested", .(scope = "tested (padj available)", .N), by = status],
  ms[responsive == TRUE, .(scope = "responsive genes", .N), by = status])
status_wide <- dcast(st_long, status ~ scope, value.var = "N", fill = 0L)
status_wide <- merge(data.table(status = st_levels), status_wide, by = "status", all.x = TRUE, sort = FALSE)
for (cc in c("count matrix (all)", "tested (padj available)", "responsive genes")) set(status_wide, which(is.na(status_wide[[cc]])), cc, 0L)
setcolorder(status_wide, c("status", "count matrix (all)", "tested (padj available)", "responsive genes"))
cat("\n[Status summary]\n"); print(status_wide)
for (i in seq_len(nrow(status_wide)))
  qadd(paste0("status ", status_wide$status[i], " (all / tested / responsive)"),
       sprintf("%d / %d / %d", status_wide$`count matrix (all)`[i], status_wide$`tested (padj available)`[i], status_wide$`responsive genes`[i]))
status_def <- data.table(status = c(st_levels, "(note) absent from the count matrix"),
                         meaning = c("the mouse name is not in the Ensembl table, so no human orthologue can be looked up (former name, predicted name, ...)",
                               sprintf("the name is in the Ensembl release %s table but no human orthologue is assigned. This does NOT mean that no human orthologue exists (other sources may supplement it)", ENSEMBL_RELEASE),
                               "the human orthologue is one2many or many2many, or there are two or more human genes",
                               "one-to-one orthologue but excluded from the DEG test (pre-filter, padj NA)",
                               "one-to-one orthologue and tested in the DEG analysis (responsive, below response threshold, non-responsive)",
                               "status of the human-side audit (05_ortholog_audit.R); it does not occur here because this table starts from the genes of the count matrix"))

# ==============================================================================
# [SECTION 4]  Responsive genes -> human genes (strict / all)
# ==============================================================================
cat("\n[3] Human counterparts of the responsive genes\n")
h_all <- pairs[, .(types_all = paste(sort(unique(htype)), collapse = ";"), n_mouse_in_ensembl = uniqueN(mname),
                   low_conf_any = any(hconf == "0")), by = hname]
h_all[, strict := types_all == "ortholog_one2one"]
resp_pairs <- merge(pairs[, .(mname, hname, htype, hconf)], ms[responsive == TRUE & !is.na(map_name), .(deg_gene = gene, mname = map_name, log2FC, padj, via_synonym = syn_status == "resolved by synonym (unique)", dup_target)], by = "mname")
resp_pairs <- merge(resp_pairs, h_all[, .(hname, types_all, strict)], by = "hname")
setorder(resp_pairs, hname, mname)

hum_resp <- resp_pairs[, .(mouse_genes = paste(sort(unique(deg_gene)), collapse = ";"), mouse_current_names = paste(sort(unique(mname)), collapse = ";"), n_mouse_responsive = uniqueN(deg_gene), via_synonym_any = any(via_synonym),
                           direction = fifelse(all(log2FC > 0), "up", fifelse(all(log2FC < 0), "down", "mixed")),
                           best_padj = min(padj), mean_log2FC = round(mean(log2FC), 3),
                           homology_types = paste(sort(unique(htype)), collapse = ";"),
                           strict_one2one = all(strict) & !any(dup_target)), by = hname]
setorder(hum_resp, hname)
cat(sprintf("   Responsive mouse genes with a human orthologue: %d -> human genes: %d (strict: %d)\n",
            uniqueN(resp_pairs$deg_gene), nrow(hum_resp), sum(hum_resp$strict_one2one)))
qadd("Responsive mouse genes with a human orthologue / corresponding human genes (all) / of which strict (one2one)",
     sprintf("%d / %d / %d", uniqueN(resp_pairs$deg_gene), nrow(hum_resp), sum(hum_resp$strict_one2one)))
qadd("Direction of the human responsive genes (all): up / down / mixed",
     sprintf("%d / %d / %d", sum(hum_resp$direction == "up"), sum(hum_resp$direction == "down"), sum(hum_resp$direction == "mixed")))

# Background set of the pathway analysis: human counterparts of the tested mouse genes
tested_pairs <- merge(pairs[, .(mname, hname)], ms[deg_status == "tested" & !is.na(map_name), .(mname = map_name)], by = "mname")
univ_all    <- sort(unique(tested_pairs$hname))
univ_strict <- sort(intersect(univ_all, h_all[strict == TRUE, hname]))
qadd("Pathway-analysis background (human counterparts of the tested mouse genes): all / strict", sprintf("%d / %d", length(univ_all), length(univ_strict)))
writeLines(sort(hum_resp$hname),                        file.path(OUT_DIR, "lists", "human_responsive_all.txt"))
writeLines(sort(hum_resp[strict_one2one == TRUE, hname]), file.path(OUT_DIR, "lists", "human_responsive_strict.txt"))
writeLines(univ_all,    file.path(OUT_DIR, "lists", "human_universe_all.txt"))
writeLines(univ_strict, file.path(OUT_DIR, "lists", "human_universe_strict.txt"))
# Also saved under the names (with gene counts) that the pathway analysis (06_pathway_ORA.R) reads directly
writeLines(sort(hum_resp$hname),                          file.path(OUT_DIR, sprintf("liver_T3_human_responsive_genes_%d.txt", nrow(hum_resp))))
writeLines(univ_all,                                      file.path(OUT_DIR, sprintf("liver_T3_human_universe_tested_%d.txt", length(univ_all))))
writeLines(sort(hum_resp[strict_one2one == TRUE, hname]), file.path(OUT_DIR, sprintf("liver_T3_human_responsive_genes_strict_%d.txt", sum(hum_resp$strict_one2one))))
writeLines(univ_strict,                                   file.path(OUT_DIR, sprintf("liver_T3_human_universe_tested_strict_%d.txt", length(univ_strict))))
writeLines(sort(hum_resp[via_synonym_any == TRUE, hname]), file.path(OUT_DIR, "lists", "human_responsive_added_by_synonym.txt"))   # responsive genes newly added by name resolution (for reference)

# Comparison with the earlier mapping (444 genes)
cmp_tab <- NULL
if (!is.na(OLD_HUMAN) && file.exists(OLD_HUMAN)) {
  # Human gene list of the earlier mapping. Any file works as long as the human-gene column has one of the names below
  # (Human gene name / human_gene etc.; case, spaces, dots and underscores are ignored). A value of "0" or empty means no mapping and is dropped.
  oh <- data.table::fread(OLD_HUMAN, data.table = TRUE)
  norm <- function(x) tolower(gsub("[ ._]", "", x))
  hc <- names(oh)[norm(names(oh)) %in% c("humangenename", "humangene", "humansymbol")][1]
  if (is.na(hc)) {
    cat("   [note] Human-gene column not found in the earlier list; comparison skipped. Column names:", paste(names(oh), collapse = ", "), "\n")
    qadd("Comparison with the earlier human responsive gene list", "skipped (human-gene column not found)")
  } else {
    old <- unique(as.character(oh[[hc]])); old <- old[!is.na(old) & old != "" & old != "0"]
    new <- hum_resp$hname
    only_old <- setdiff(old, new); only_new <- setdiff(new, old)
    qadd(sprintf("Comparison with the earlier list (%s, column %s): earlier / current / common / earlier only / current only", basename(OLD_HUMAN), hc),
         sprintf("%d / %d / %d / %d / %d", length(old), length(new), length(intersect(old, new)), length(only_old), length(only_new)))
    cmp_tab <- rbind(
      if (length(only_old)) data.table(human_gene = only_old, category = "only in the earlier list") else NULL,
      if (length(only_new)) data.table(human_gene = only_new, category = "only in the current list") else NULL)
    if (!is.null(cmp_tab)) {
      mm <- resp_pairs[, .(mouse_genes = paste(sort(unique(deg_gene)), collapse = ";")), by = hname]
      cmp_tab[, current_mouse_genes := mm$mouse_genes[match(human_gene, mm$hname)]]
    }
    cat(sprintf("   Versus the earlier list (%d): earlier only %d, current only %d\n", length(old), length(only_old), length(only_new)))
  }
}

# ==============================================================================
# [SECTION 5]  Save
# ==============================================================================
inp <- unique(c(ENSEMBL_FILE, DEG_FILE, RAW_COUNT, OLD_HUMAN, SYNONYM_FILE)); inp <- inp[!is.na(inp)]
man <- rbindlist(lapply(inp, function(f) if (file.exists(f)) data.table(file = f, size_bytes = file.size(f), modified = format(file.mtime(f), "%Y-%m-%d %H:%M:%S"), md5 = unname(tools::md5sum(f)))
                        else data.table(file = f, size_bytes = NA, modified = NA, md5 = "file not found")))
data.table::fwrite(man, file.path(OUT_DIR, "input_files_manifest.tsv"), sep = "\t")
data.table::fwrite(data.table(item = c("ENSEMBL_FILE", "DEG_FILE", "RAW_COUNT", "OLD_HUMAN", "SYNONYM_FILE", "PADJ_CUT", "LFC_CUT", "ENSEMBL_RELEASE", "OUT_DIR", "run start", "working directory"),
                               value = as.character(c(ENSEMBL_FILE, DEG_FILE, RAW_COUNT, OLD_HUMAN, SYNONYM_FILE, PADJ_CUT, LFC_CUT, ENSEMBL_RELEASE, OUT_DIR, format(RUN_START, "%Y-%m-%d %H:%M:%S"), getwd()))),
                   file.path(OUT_DIR, "settings.tsv"), sep = "\t")
cat("\n[5] Saving:", OUT_XLSX, "\n")
readme <- data.table(item = c("1_mouse_reason_codes", "2_responsive_pairs", "3_human_responsive", "4_reason_summary", "5_previous_list_comparison", "4b_status_summary", "0_status_definitions", "6_name_resolution"),
                     value = c("all 13,319 genes of the count matrix. M0 one-to-one, M3 not one-to-one, M2 no human orthologue, M1 name not in Ensembl",
                           "responsive mouse gene x human gene (homology type, confidence)", "direction and strict flag per human gene",
                           "number of genes per reason code by scope (all, tested, responsive)", "difference from the earlier human responsive gene list (444 genes)",
                           "counts per scope under the status scheme of the co-author's analysis (status and status_detail columns are in sheet 1_mouse_reason_codes)", "meaning of each status", "synonym-resolution result of the mouse genes whose name was not in Ensembl (resolved, none, several genes, conflict)"))
if (is.null(cmp_tab)) readme <- readme[!startsWith(readme[["item"]], "5_")]   # drop it from the README when there is no comparison sheet
wb <- createWorkbook()
add_sheet(wb, "README_QC", rbind(readme, data.table(item = "-----", value = "-----"), qc))
add_sheet(wb, "1_mouse_reason_codes", ms[order(reason, gene)])
add_sheet(wb, "2_responsive_pairs", resp_pairs)
add_sheet(wb, "3_human_responsive", hum_resp)
add_sheet(wb, "4_reason_summary", reason_wide)
add_sheet(wb, "4b_status_summary", status_wide)
add_sheet(wb, "0_status_definitions", status_def)
add_sheet(wb, "6_name_resolution", ms[syn_status != "name match", .(gene, syn_status, current_name = map_name, human_genes, homology_types, deg_status, responsive, log2FC, padj, reason, status)][order(syn_status, gene)])
if (!is.null(cmp_tab)) add_sheet(wb, "5_previous_list_comparison", cmp_tab)
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
data.table::fwrite(ms, file.path(OUT_DIR, "ortholog_mouse_status_all.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\n===== QC =====\n"); print(qc, nrows = 100)
cat("\nDone.\n")
cat(sprintf("\nEnd: %s (elapsed %.1f s)\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), as.numeric(difftime(Sys.time(), RUN_START, units = "secs"))))
}, warning = function(w) { cat("[warning]", conditionMessage(w), "\n"); invokeRestart("muffleWarning") },
   message = function(m) { cat("[message]", conditionMessage(m)); invokeRestart("muffleMessage") })
"OK"
}, error = function(e) { cat("\n[error]", conditionMessage(e), "\n"); paste("ERROR:", conditionMessage(e)) },
   finally = { sink(); close(con_log) })
cat("\nRun status:", status, "\nLog:", file.path(OUT_DIR, "261021_remap_run_log.txt"), "\n")
