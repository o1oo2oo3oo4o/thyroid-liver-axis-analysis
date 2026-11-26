# Provenance of the Supplementary Tables and Figures

This file states, for every Supplementary Table of the manuscript, which input data and which script of this repository
produced it, and what was checked. "Not included" means the step was run outside the scripts of this repository (web
server, stand-alone program, or by a co-author) or its script has not been added to this release yet.

## Supplementary Tables

| Table | Content | Input | Script | Status in this repository |
|---|---|---|---|---|
| S1 | MAGMA competitive gene-set output, 12 tests | `gene_sets/gwas_gene_flags.tsv`; 12 GWAS summary statistics; MAGMA v1.10, `NCBI37.3.gene.loc`, `g1000_eur` | `00b_magma_gene_set_enrichment.R` | included (MAGMA, the reference files and the GWAS files are downloaded separately). Checked: all 12 rows (number of genes, beta, standard error, P) agree with the MAGMA output files (`gsa_*.gsa.out`) used for the manuscript. The script was checked with a stand-in for the MAGMA program: command lines, file handling (chromosome/position to rsID, rows without P or N dropped), resumption and the assembly of the 12-row table with its significance labels and trait-group columns reproduce Table S1 from MAGMA-format output. MAGMA itself was not re-run in the checking environment |
| S2a, S2b | Colocalization windows (458) and all dataset pairs (66 per window), prior sensitivity | public GWAS summary statistics (see `data/README.md`) for S2b; S2b for S2a | `01a_coloc_scan.R` (S2b), `01b_coloc_window_tiers.R` (S2a) | included. Checked: the 29,802 pair rows of S2b equal the output of the colocalization scan of the study; `01b` run on S2b gives the maximum cross-domain PP.H4, the tier, the supporting pairs and the prior-sensitivity values of all 458 windows of S2a without any difference (11 high, 16 moderate, 431 none). `01a` is the cleaned scan script (same parameters: windows lead +/- 500 kb, P < 1e-6 for lead variants, MAF >= 0.01, at least 50 shared variants, p1 = p2 = 1e-4, p12 = 1e-5); it was not re-run on the full GWAS files (about 16 GB of memory and 2 h), but its data handling and its output format were tested on small simulated files with a stand-in for `coloc.abf`, and its output is read by `01b` |
| S3a, S3b | Overlap statistics and the 199 shared genes | FUMA gene lists of both domains (`gene_sets/gwas_gene_flags.tsv`), S2a, `hg19_refgene.txt`, FUMA `genes.txt` | `00_gwas_gene_sets_from_FUMA.R` (gene sets), `01c_window_gene_mapping.R` (status of the shared genes) | included. Checked on the real FUMA files of the 12 jobs: `00` reproduces `gwas_gene_flags.tsv` gene by gene (0 genes only in one of the two, 0 different flags) and the Table S3a numbers thyroid 1,709 and shared 199; per trait PDFF 786, cT1 486, TSH 828, FT3 149, FT4 999. The liver-derived set has 1,120 symbols and the union 2,630 symbols; Table S3a prints 1,119 and 2,629 genes because ABHD16B and C20ORF135 are one gene (see "Gene counts" below). The 199 genes and their thyroid and liver traits are identical to the gene column and trait columns of Table S3b (the symbols of the flags file, i.e. FUMA symbols). Status of the shared genes (`01c` with the real FUMA `genes.txt` files): all 199 genes have a position; status, nearest window, supporting pairs and traits are identical to Table S3b (high 16, no 156, nominal 11, outside tested loci 16) |
| S4 | Gene to window assignment of the 27 supported windows | S2a, gene coordinates (GRCh37, UCSC refGene; FUMA `genes.txt` for symbols that are not in refGene) | `01c_window_gene_mapping.R` | included. Checked: with refGene alone 509 of the 526 gene-window rows are found; with the real FUMA `genes.txt` files all 526 rows are found and agree in coordinates, overlap, distance, biotype, number of windows, representative window and shared-gene flag (0 differences). The script marks four pairs of symbols in the same window as possible former / current names of one gene: SPX / C12orf39 and MYMK / TMEM8C (the two pairs of Table S4 and of the Results text; 60 % and 74 % overlap of the shorter gene) and, in addition, C19orf45 / TEX45 (70 %, same end position, same strand) and RP11-45M22.4 / FLCN (72 %, same strand). These two extra marks are candidates only and Table S4 does not make them (see open points). The goldenPath refGene layout was tested as well |
| S5a | Exploratory lead-variant cis-eQTL alignment | S2 lead variants, public eQTL resources | none | not included (not used in candidate selection) |
| S5b | Exploratory protein interaction edges | STRING, seed genes of the supported windows | none | not included (not used in candidate selection) |
| S6a | Orthologue assignment of the 448 human entries | `250912__DEG_all_genes.csv`, `0_count_matrix.csv`, Ensembl BioMart release 116 | `04_reortholog_synonym.R` | included. Checked: the 448 entries and their mapping fields agree with the table |
| S6b | Retrieval status of every evaluated positional gene | `250912__DEG_all_genes.csv`, `0_count_matrix.csv`, gene to window mapping, Ensembl BioMart | `05_ortholog_audit.R` | included. Checked: status labels agree with the audit output |
| S7 | Two-sample Mendelian randomization | pQTL summary statistics, liver MRI GWAS | co-author | not included. In the revised version a protein-outcome pair is called significant when its IVW FDR is < 0.05 (no condition on the other methods); `figures/Fig3_protein_candidates_supp.R` reads S7 in the layout of the revised workbook, and `scripts/revision/13_recount_MR_support_IVW.R` recounts the pairs, proteins, Steiger and instrument-strength numbers from S7, S14, S15 and S16 (see "Revision analyses") |
| S8 | Plasma proteomic associations | plasma proteomics cohort | co-author | not included |
| S9a | Human hepatocyte T3 experiment (GSE270233): wild-type effect, knockout effect, interaction, baseline contrast | GEO GSE270233 counts | `08_deseq2_human_hepatocyte_THRB_KO.R` | included. Checked: re-run in R 4.3.3 / DESeq2 1.42.0 on the count table gives the same genes, the same missing values and the same baseMean, log2 fold change, standard error, statistic, P and adjusted P as the original result files (largest absolute differences below 3e-5 for the statistic and 6e-7 for P values); the four contrast counts (799 / 35 / 207 / 719 responsive genes) agree with the table |
| S9b | 2 x 2 tables of the two candidate enrichment tests | `T3_effect_WT_KO_interaction_wide.csv` | `11_THRB_candidate_enrichment.R` | included. Checked: counts and both P values reproduced exactly |
| S9c | Mouse-human concordance of the glutathione-set genes | mouse and human DESeq2 results | `07_redox_cameraPR.R` (`S9_mouse_human_concordance_genes.csv`) | included |
| S10a | Gene lists of the over-representation analysis | `04_reortholog_synonym.R` output | `04_reortholog_synonym.R` | included |
| S10b | Full over-representation output (GO, KEGG, Reactome) | S10a lists | `06_pathway_ORA.R` | included. Checked: all 5,464 rows agree with the output workbook of the analysis |
| S11a-d | Redox gene sets: Fisher test, cameraPR, correlation and transferase sensitivity, membership, gene-level values | `261020__DEG_all_genes_with_stat.csv`, `DEG_WT_T3_effect.csv`, `DEG_KO_T3_effect.csv` | `07_redox_cameraPR.R` | included. Checked: S11a and S11d agree with the output files of the analysis |
| S12a, S12b | Fibrosis-stage results per gene and model; candidate meta-analysis | GSE130970, GSE135251 | `09_human_liver_fibrosis.R` | included. Checked: S12b agrees with the results workbook; the adjusted P values of S12b are Benjamini-Hochberg values within the candidate list (20 genes testable in both cohorts) |
| S13a-c | Set-level cameraPR results, set-score trends, cohort agreement | GSE130970, GSE135251, S11d set definitions | `09_human_liver_fibrosis.R` | included. Checked: S13a (77 model rows) agrees with the results workbook |
| S13d | Positive-control genes of the fibrosis analysis | primary-model gene results of both cohorts | `10_fibrosis_positive_controls.R` | included. Checked: the 20 genes, log2 fold changes, P values and flags agree with the table |

## Steps that were re-written for this repository

For four steps the code that made the original files was not kept in a form that can be run without interaction (console
history, interactive prompts, or no script). These scripts were written for the repository and compared with the original
result files:

| Script | Original situation | Check |
|---|---|---|
| `02_build_count_matrix_GSE85793.R` | the count matrix was assembled by hand before the analysis scripts were written | columns C1 and T1 are identical to the GEO files; `1_cleaned_count_matrix.csv` is the C1-C4/T1-T4 subset of `0_count_matrix.csv` (13,319 of 13,319 rows identical) |
| `03_deseq2_mouse_liver_T3.R` | the original analysis script (`2_run_DESeq2.R`) asks for its settings at the console; this script fixes the same settings (design ~ condition, reference Control, count >= 10 in half of the samples, apeglm shrinkage) | R 4.3 / DESeq2 1.42.0: 12,135 genes after filtering, identical baseMean and P values, the same genes with a missing adjusted P (11,819 with a value), adjusted P identical within 3e-7, Wald statistic identical. The shrunken log2 fold changes (apeglm) were not recomputed because apeglm was not available in the checking environment |
| `10_fibrosis_positive_controls.R` | the stage_check sheet was added to the results workbook after the main script had run | same 20 genes, log2 fold changes and P values to machine precision |
| `11_THRB_candidate_enrichment.R` | the code of Table S9b was not kept | 2 x 2 counts and both P values identical to Table S9b |

All other scripts are the analysis scripts of the study with English comments and messages, configurable paths and the
test hooks removed; the analysis code was not changed. They were checked for syntax and for identical code structure
against the original scripts, and the key outputs listed above were compared with the original result files. They were
not all re-run end to end on the full public data in the environment in which they were cleaned.

## Figures

Every figure script reads the Supplementary Tables workbook (Figure 6: the output files of `scripts/revision/11_*`) and, where stated, public files. Two workbook layouts are
accepted: one sheet per table (`Table S2a`) or one sheet per table number with the lettered tables side by side
(`Table S2`, as in the submitted workbook). Figures 3, 4 and 5 were run from start to end on the revised workbook (layout with the lettered tables side by side) and on a workbook with one sheet per table; the figure files of both runs are identical, and the numbers that the scripts print equal the numbers of the manuscript: 669 IVW pairs and 224 proteins; 7 decreased and 4 increased candidates in mouse liver; 6 measured proteins; 19 of the 22 panel genes tested in hepatocytes, 3 of 3 responsive candidates, conditional Fisher P = 0.014. Figure 6 was run on the output of `scripts/revision/11_*`; its image file is identical to the output of the script version that drew the revised figure (20 genes, 15 in the same direction). Figure 2 is not changed in this release. Figures 2 and 4 use `ggrastr` for the dense point layers when it is installed and draw them as ordinary points otherwise; the runs above were made without `ggrastr`, so the rasterised path was not tested. The gene labels placed by `ggrepel` (Figures 4 and 6) are fixed by `set.seed(1)`, so the layout is the same in every run.

| Figure | Script | Tables read | Further input |
|---|---|---|---|
| 2 | `Fig2_colocalization_and_T3_response_supp.R` | S2a, S2b, S4 | `GSE85793_mouse_liver_T3_DESeq2_all_genes.csv` (panel C); public GWAS files (optional, panel B) |
| 3 | `Fig3_protein_candidates_supp.R` | S7 (revised layout: one block of columns per outcome, rows per pQTL source, protein and method; IVW FDR < 0.05), S8 | `GSE85793_mouse_liver_T3_DESeq2_all_genes.csv` (panel B; read from sheet S6c if the workbook has one) |
| 4 | `Fig4_THRB_dependency_supp.R` | S9a | `T3_effect_WT_KO_interaction_wide.csv`, or the three files `DEG_WT_T3_effect.csv`, `DEG_KO_T3_effect.csv`, `DEG_Interaction.csv`, for the standard errors (panels B and C; a table S9d is used only if it has the columns `gene`, `WT_T3_lfcSE`, `KO_T3_lfcSE`, `Interaction_lfcSE`; the S9d of the submitted workbook is another table and is ignored) |
| 5 | `Fig5_pathway_and_redox_supp.R` | S10b, S11a, S11b, S11d | - |
| 6 | `Fig6_glutathione_name_matched.R` (replaces `Fig6_mouse_human_concordance_supp.R`, which used the earlier set of 18 genes) | none | `11b_concordance_genes_name_matched.csv` and `11d_concordance_summary.csv` (output of `scripts/revision/11_glutathione_name_matched_concordance.R`, 20 name-matched genes), `DEG_WT_T3_effect.csv` |
| S1 | `FigS1_MAGMA_gene_set_supp.R` | S1 | - |
| S2 | `FigS2_prior_sensitivity_supp.R` | S2a, S2b | - |
| S3 | `FigS3_window_layouts_supp.R` | S2a, S2b, S4 | public GWAS files (optional) |

Without the public GWAS summary statistics, the regional panels of Figures 2B and S3 are drawn as window and gene layouts
instead of regional association plots.

## Revision analyses (`scripts/revision/`)

The scripts of `scripts/revision/` produce the additional analyses of the responses to the reviewers (see
`scripts/revision/README.md`). Their inputs are listed in `scripts/revision/data/README.md`.

| Script | Reproduces / produces | Checked |
|---|---|---|
| `00_prepare_GSE_inputs.R` | count matrices and sample tables of GSE85793, GSE130970, GSE135251 and GSE270233 from the original GEO files | the 8-sample matrix of GSE85793 equals the count matrix of the submitted DEG analysis (genes, order, samples, values identical); sample numbers and stage counts agree with the series matrices (GSE130970: F0 - F4 = 25 / 28 / 9 / 14 / 2; GSE135251, 206 NAFLD samples: 38 / 47 / 53 / 54 / 14); the sample order of GSE270233 agrees with the series matrix |
| `01_KO_redox_interaction.R` | attenuation of the redox / glutathione T3 response in THRB-knockout hepatocytes | re-run on `DEG_hepT3_raw.csv`: the output tables agree with the result files used for the response (identical, except that one set-level column differs in the 15th significant digit) |
| `02_biopsy_stage_analysis.R` | models M1 - M7 in the two biopsy cohorts, heterogeneity, meta-analysis | run on the real counts; the reproduction check against the stored values of Table S12 gives a maximum difference of 0.0005 |
| `03_mouse_T3_robustness_v2-apeglm.R` | robustness of the mouse set tests | run on the real counts; reproduces 541 responsive genes, glutathione Fisher P = 2.07e-4 and cameraPR P = 1.8e-6 (rho = 0.01) of the manuscript |
| `04_human_final_tables_v2.R` | panel and glutathione tables, maturity markers | re-run on the result tables: T1 (panel) and T2 (glutathione annotation) identical to the result files; the T3 baseMean values use the size factors of the submitted DESeq2 run and equal Table S9a |
| `05_mouse_QC_leave_one_animal_out_v3-mainmap.R` | sample QC, eight leave-one-animal-out fits | run on the real counts, with the ortholog table of the main analysis (448 foreground / 10,085 background genes) for the GO term |
| `06_biopsy_F4_exclusion.R` | stage counts, all biopsy models with F4 excluded | run on the real counts; the Benjamini-Hochberg values recomputed from the stored P values agree with the stored values to 1e-15 (check in `12`) |
| `07_coloc_nonUKB_sensitivity.R` | colocalization windows without TSH dataset B | re-run on the pair table of Table S2: 27 windows (11 high, 16 moderate) with the rule of the manuscript; 28 windows (13 high, 15 moderate) without TSH dataset B; no window is lost; output identical to the result file |
| `10a_make_tss_hg19.R`, `10_TRb_ChIP_annotation.R` | hg19 TSS table; TR-beta peaks around the panel and the redox sets | run on the peaks of Ayers et al. 2014 and the UCSC knownGene TSS table |
| `11_glutathione_name_matched_concordance.R` | 20 name-matched glutathione genes (of 36): direction agreement and human set tests | re-run on the DESeq2 result tables: all output tables identical to the result files; the earlier 18-gene set is reproduced by the same code (14 of 18 concordant, P = 0.0154; cameraPR P = 0.337 at rho = 0.01) |
| `12_recompute_panel_stats_IVW.R` | Fisher enrichment, Benjamini-Hochberg adjustment and ChIP proportions for the panel size | re-run in both modes, identical to the result files. Baseline (26 genes, 23 expressed in hepatocytes): 3 of 23, expected 0.29, one-sided Fisher P = 0.0030. IVW-centred panel (22 genes, 19 expressed): 3 of 19, expected 0.24, P = 0.0017. The test restricted to wild-type responsive genes (3 of 3) is not panel-size dependent (P = 0.014) |
| `13_recount_MR_support_IVW.R` | the numbers that depend on the definition of genetic support (IVW FDR < 0.05), recounted from S7, S14, S15 and S16; the Mendelian randomization itself is not re-estimated | run on the revised workbook: 669 IVW pairs (deCODE 411, ARIC 258) of 262 tested proteins; 224 proteins (172 deCODE, 100 ARIC); steatotic 139, fibrotic 128, both 43 (one direction only 136 and 122); 222 proteins with the additional MR-Egger rule (ALDH2 and GCKR drop out); Steiger P < 0.05 in 654 pairs, P >= 0.05 in 1 (deCODE NCAN and Liver PDFFvat, P = 0.208), not computed in 14; colocalization strong 51 and suggestive 20; MR-RAPS consistent in 643 of 669 pairs and of 644 estimable pairs (25 pairs with fewer than 3 variants; Table S16 lists 24 pairs with fewer than 3 instruments, because deCODE CXCL10 and Liver cT1non has 3 instruments but 2 MR-RAPS variants); 9,356 unique instruments, F 29.7 - 58,031.9, median 72.3 (IQR 44.2 - 155.6), none below 10; the 21 cells of Table S15 (overall, deCODE, ARIC) equal the recomputed values. 35 of the 9,356 instruments have different F values in different outcome rows of Table S16; their median is used (the summary does not depend on this beyond the last decimal). ADH4, BTC, COL11A2, HSPA1A and LGALS4 have no IVW-significant pair (weighted-median-significant pairs 3, 1, 2, 2, 2) |

The panel after the IVW-centred definition has 22 genes: `BTC`, `COL11A2`, `HSPA1A` and `LGALS4`, which had no
IVW-significant pair, are removed from the 26 genes; `ADH4` also had no IVW-significant pair but stays because it is a
positional (colocalization) candidate.

`08_instrument_F_statistics.R` and `09_reverse_MR.R` (instrument F statistics and reverse-direction Mendelian
randomization) are not included in this release; they could only be tested on simulated data.

## Steps 0 and 1: what was checked, and how

| Script | Reproduces | Run on | Result |
|---|---|---|---|
| `00_gwas_gene_sets_from_FUMA.R` | gene sets | (a) the real FUMA files of the 12 jobs: `genes.txt`, `GenomicRiskLoci.txt` of all jobs and `eqtl.txt` of 11 jobs (the 670 MB `eqtl.txt` of TSH B was not copied to the checking environment; its result of the cis rule was taken from the step files that the original run wrote, `RESUME=1`); (b) simulated FUMA folders with decoy genes that fail each rule | (a) the number of genes after each rule equals the counts of the original run for all 11 recomputed jobs (cis PASS 337 / 473 / 515 / 446 / 332 / 344 / 350 / 372 / 149 / 1,003 / 891), and the gene sets equal `gwas_gene_flags.tsv` (2,630 genes, 0 differences); (b) exact recovery, decoys excluded |
| `00b_magma_gene_set_enrichment.R` | S1 | simulated reference panel and summary statistics with a stand-in for the MAGMA program | table layout, labels and trait-group columns identical to S1 (12 of 12 rows) |
| `00c_gene_overlap_hypergeometric.R` | overlap test | `gwas_gene_flags.tsv` | N = 19,000 / 20,000 / 25,000: expected 100.74 / 95.70 / 76.56, fold 1.975 / 2.079 / 2.599, P 1.03e-21 / 2.46e-24 / 7.26e-37, the values of the original result file (computed with the 1,120 symbols of the flags file; with 1,119 genes the expected values would be 100.65 / 95.62 / 76.49); permutation P < 1e-5 |
| `01a_coloc_scan.R` | S2b | small simulated GWAS files (one gzip file in a folder named by its accession), stand-in for `coloc.abf` | windows and pair table in the expected format, read by `01b` |
| `01b_coloc_window_tiers.R` | S2a | S2b | 0 differences in 458 windows |
| `01c_window_gene_mapping.R` | S4, S3b | S2a, UCSC refGene (hg19) in both file layouts, the 12 real FUMA `genes.txt` files | S4: 526 of 526 rows, 0 differences; S3b: 199 of 199 genes, 0 differences (see S3 and S4 above) |

## Gene counts: FUMA symbols and Table S3a

`gene_sets/gwas_gene_flags.tsv` has one row per gene symbol as FUMA (Ensembl v102) writes it: 2,630 rows, of which 1,709
thyroid-derived, 1,120 liver-derived and 199 shared. Table S3a, the Results text and the final integrated analysis file of the study
(version 11, 5 November 2026) give 1,119 liver-derived genes and a union of 2,629; thyroid (1,709) and shared (199) are the same.

The one-gene difference is a gene that FUMA lists twice. In the PDFF unadjusted run (GCST90446650) `genes.txt` has both ABHD16B
(ENSG00000183260, chr20:62,492,566-62,494,341, + strand) and C20ORF135 (ENSG00000268858, chr20:62,492,879-62,493,217, - strand, no
Entrez ID), an old name of ABHD16B that lies inside it. Both pass the cis rule, so `00` counts two symbols. The integrated file uses
the current name and counts the gene once, which gives 1,119 and 2,629. When the 16 old symbols below are replaced by the current symbols of the integrated file (C20ORF135
merged into ABHD16B), the flags file gives the GWAS-liver, GWAS-thyroid and shared gene sets of the integrated file without any
difference (1,119 / 1,709 / 199, union 2,629). Only the C20ORF135 / ABHD16B pair changes a count; the other 15 change the name only.
The console check of `00` against Table S3a therefore shows a difference of 1 in the liver and union rows and 0 in the other two.
The flags file and Table S3b keep the FUMA symbols.

| FUMA symbol (flags file, S3b) | Current symbol (integrated file) | Set |
|---|---|---|
| BRE | BABAM2 | liver |
| EMR1 | ADGRE1 | liver |
| FAM101A | RFLNA | liver |
| FAM64A | PIMREG | liver |
| SGK223 | PRAG1 | liver |
| C20ORF135 | ABHD16B (already a separate row: merged) | liver |
| FAM65B | RIPOR2 | liver and thyroid (shared) |
| HIST1H1E | H1-4 | liver and thyroid (shared) |
| TMEM8C | MYMK | liver and thyroid (shared) |
| ADSSL1 | ADSS1 | thyroid |
| C12ORF39 | SPX | thyroid |
| C17ORF53 | HROB | thyroid |
| CASC5 | KNL1 | thyroid |
| CCDC101 | SGF29 | thyroid |
| FAM109A | PHETA1 | thyroid |
| KIAA0101 | PCLAF | thyroid |

## Open points

- **Sample sizes.** The N used in the FUMA jobs and in the colocalization scan (FT3 57,107; FT4 117,166; TSH A 55,342; TSH B 124,358) differs from Table 1 of the manuscript
  (59,061; 119,120; 119,715; 247,107). The numbers for the eight liver datasets agree. S2b was produced with the first set of numbers.
- **Possible alias pairs in Table S4.** With the real FUMA files the script marks four symbol pairs in the same window with a partial
  overlap (SPX / C12orf39, MYMK / TMEM8C, C19orf45 / TEX45, RP11-45M22.4 / FLCN). Table S4 and the Results text count two of them as one gene
  (the first two) and list the other two as separate genes (501 distinct genes). If both further pairs were also one gene each, the number would be 499.
  The pairs are written to `fuma_possible_aliases.tsv` for the authors to decide; RP11-45M22.4 is a clone-based Ensembl name, TEX45 is a RefSeq symbol.
- **TSH pre-filtering.** The two TSH files were reduced to variants with P < 1e-5 before FUMA; no script for this step is in the
  repository (`docs/FUMA_settings.md`).
- `00c` is an exploratory analysis and is not part of the submitted tables.

## Not included in this release

The FUMA web analysis itself, the STRING, Mendelian randomization and plasma proteomics analyses, and the exploratory
annotations of Table S5 were run with other tools or by co-authors. Their settings are given in the Methods of the manuscript.
