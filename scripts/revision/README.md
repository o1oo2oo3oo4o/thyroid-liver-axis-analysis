# Revision analyses

Scripts of the additional analyses that were run for the revision of the manuscript (responses to the reviewers and the
editor). They are separate from the pipeline of the main analysis in `../` (steps 0 - 11, which produce the supplementary
tables of the submitted version). The scripts of this folder read their inputs from `data/` and write to `results/`;
**run every script from this folder** (`scripts/revision`), all paths are relative to it. `config.R` holds the file names,
the column names and the gene lists; the scripts read it with `source()`.

No raw data are part of the repository: the inputs are public GEO files, tables of the study or result files of the
earlier steps (`data/README.md`).

## Scripts

| Script | Question (reviewer comment) | Main output (`results/...`) |
|---|---|---|
| `00_prepare_GSE_inputs.R` | builds the count matrices and sample tables that `config.R` expects from the original GEO files (GSE85793, GSE130970, GSE135251, GSE270233); run once | `data/*_counts.csv`, `data/*_meta.csv`, `data/00_prepare_GSE_inputs_log.txt` |
| `01_KO_redox_interaction.R` | is the T3 response of the redox / glutathione genes attenuated by THRB loss in human hepatocytes? (R1 comment 5) | `01_KO_redox/` |
| `02_biopsy_stage_analysis.R` | are the fibrosis-stage associations of the candidates robust to fibrosis-stage heterogeneity? Models M1 - M7, meta-analysis with Cochran Q and I2 (R1 comment 6) | `02_biopsy_stage/` |
| `03_mouse_T3_robustness_v2-apeglm.R` | 4 versus 4 mice: estimated inter-gene correlation for cameraPR, leave-one-out subsets, limma-voom as a second engine (R1 comment 4) | `03_mouse_robustness_v2-apeglm/` |
| `04_human_final_tables_v2.R` | panel genes with 95% CI in the hepatocyte experiment, annotation of the 36 glutathione genes, maturity markers (R1 comments 2, 3, 5) | `04_human_final_tables_v2/` |
| `05_mouse_QC_leave_one_animal_out_v3-mainmap.R` | sample QC and eight leave-one-animal-out fits (R1 comment 4) | `05_mouse_QC_LOAO_v3-mainmap/` |
| `06_biopsy_F4_exclusion.R` | stage counts per cohort, all biopsy models with the hepatocyte score kept and F4 excluded (R1 comment 6, R2 major comment 4) | `06_biopsy_F4/` |
| `07_coloc_nonUKB_sensitivity.R` | colocalization windows without TSH dataset B, which contains UK Biobank participants (R2 major comment 1, Table S2) | `07_coloc_nonUKB/` |
| `10a_make_tss_hg19.R` | transcription start sites (hg19) of all genes; run once before `10` | `data/tss_hg19.csv` |
| `10_TRb_ChIP_annotation.R` | published HepG2 TR-beta BioChIP-seq peaks (Ayers et al. 2014) around the TSS of the prioritized genes and the redox sets (R2 major comment 3) | `10_TRb_ChIP/` |
| `11_glutathione_name_matched_concordance.R` | mouse-human direction agreement and human set tests for the glutathione genes that are tested in both species by gene symbol (20 of 36 genes); input of `../../figures/Fig6_glutathione_name_matched.R` | `11_glutathione_name_matched/` |
| `12_recompute_panel_stats_IVW.R` | recomputes the three statistics that depend on the number of genes in the prioritized panel (Fisher enrichment, Benjamini-Hochberg adjustment of the biopsy meta-analysis, TR-beta ChIP proportions) for the panel after the IVW-centred Mendelian randomization definition (22 genes) and for the 26-gene baseline | `12_panel_recompute_IVW_drop4/`, `12_panel_recompute_baseline_26/` |
| `13_recount_MR_support_IVW.R` | recounts, from Tables S7, S14, S15 and S16, every number that depends on the definition of genetic support (IVW FDR < 0.05): pairs and proteins, steatotic / fibrotic outcome groups, direction, MR-Egger rule, Steiger P, colocalization and MR-RAPS status, F statistics of the instruments; compares each value with the manuscript (definition of genetic support, Response to the Editor) | `13_MR_support_IVW/` |

The suffixes `v2-apeglm`, `v2` and `v3-mainmap` are kept in the file names so that they match the names of the result folders.
They are the versions used for the manuscript: a gene is called "responsive" with the apeglm-shrunken log2 fold change
(as in the DEG table of the manuscript), the glutathione / xenobiotic tests of `05` use the ortholog table of the main
analysis (448 foreground / 10,085 background genes), and `04` takes the DESeq2 size factors of the submitted run. Earlier
versions of these scripts (unshrunken log2 fold change, synonym-resolved ortholog table) are not part of the repository.

## Order

```bash
cd scripts/revision
Rscript 00_prepare_GSE_inputs.R                        # once: count matrices and sample tables from the GEO files
Rscript 01_KO_redox_interaction.R
Rscript 02_biopsy_stage_analysis.R
Rscript 03_mouse_T3_robustness_v2-apeglm.R
Rscript 04_human_final_tables_v2.R
Rscript 05_mouse_QC_leave_one_animal_out_v3-mainmap.R
Rscript 06_biopsy_F4_exclusion.R
Rscript 07_coloc_nonUKB_sensitivity.R
Rscript 10a_make_tss_hg19.R                            # once: needs Bioconductor annotation packages (internet)
Rscript 10_TRb_ChIP_annotation.R
Rscript 11_glutathione_name_matched_concordance.R
Rscript 12_recompute_panel_stats_IVW.R                 # after 04, 06 and 10
Rscript 13_recount_MR_support_IVW.R                  # needs the revised Supplementary Tables workbook (config.R: supp_revised_file, or SUPP_XLSX)
```

Script 12 reads the result files of `06` and `10` and `data/DEG_hepT3_raw.csv`. It uses the IVW-centred panel (22 genes) by
default; the baseline panel of the submitted analysis (26 genes) is run with an empty `STAGE_DROP`:

```bash
STAGE_DROP="" Rscript 12_recompute_panel_stats_IVW.R                     # baseline, 26 genes
STAGE_DROP="BTC,COL11A2,HSPA1A,LGALS4" Rscript 12_recompute_panel_stats_IVW.R   # same as the default
```

The baseline run must reproduce the values of the submitted manuscript; the script stops if the Benjamini-Hochberg values
stored by `06` are not reproduced.

## Software

R 4.3 or later. Packages: data.table, ggplot2, readxl, limma, DESeq2, apeglm (needed by `03` and `05`; they stop if it is
missing), dplyr and readr (`11`). Optional: clusterProfiler and org.Hs.eg.db (GO terms in `01` and `05`; the
steps are skipped with a message if they are missing). `00` additionally needs org.Hs.eg.db and AnnotationDbi, and `10a`
needs TxDb.Hsapiens.UCSC.hg19.knownGene, GenomicFeatures and org.Hs.eg.db (installed from Bioconductor if missing).
The session information of each run is written next to its results.

## What was checked

| Script | Check |
|---|---|
| `01`, `04`, `07`, `11`, `12` | re-run on the real result tables of the study from this folder layout (`docs/PROVENANCE.md`, section "Revision analyses"); the outputs agree with the result files that were used for the response letters (identical, or equal to the last printed digit) |
| `13` | run on the revised Supplementary Tables workbook (S7, S14, S15, S16): 669 IVW pairs with FDR < 0.05 (deCODE 411, ARIC 258) of 262 tested proteins; 224 proteins (deCODE 172, ARIC 100); steatotic outcomes 139, fibrotic 128, both 43 (one direction only: 136 and 122); 222 proteins if the MR-Egger intercept rule were added; Steiger P < 0.05 in 654 pairs, P >= 0.05 in 1, not computed in 14; MR-RAPS consistent in 643 of 644 estimable pairs (25 pairs with fewer than 3 variants); colocalization strong 51, suggestive 20; 9,356 unique instruments, F median 72.3 (IQR 44.2 - 155.6), none below 10; the 21 cells of Table S15 (overall, deCODE, ARIC) equal the recomputed values. ADH4, BTC, COL11A2, HSPA1A and LGALS4 have no IVW-significant pair (weighted-median-significant pairs: 3, 1, 2, 2, 2). All 32 values of `3_count_summary.csv` agree with the manuscript. The Mendelian randomization itself is not re-estimated |
| `00`, `02`, `03`, `05`, `06`, `10`, `10a` | run by the author on the real GEO / public files; the reproduction checks that the scripts print (for example the 541 responsive genes of `03`, the Table S12 values in `02`, the check of `00` against the count matrix of the submitted DEG analysis) agree. They were not re-run in the environment in which the files of this repository were cleaned |

## Not included

`08_instrument_F_statistics.R` (F statistics of the harmonized instruments) and `09_reverse_MR.R` (reverse-direction
Mendelian randomization) are not part of this release: they were only tested on simulated data, because the harmonized
instrument table, the LD-clumped liver-trait variants and the protein summary statistics (Mendelian randomization and
proteomics inputs of the collaborating group) are not available here.
