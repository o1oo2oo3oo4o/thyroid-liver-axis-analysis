# Input data

Nothing in this folder is tracked by git except `gene_sets/gwas_gene_flags.tsv` (a result table of the study, see below).
Place the files listed below here (or point `DATA_DIR` to another folder). Files are found by their **name**, also in
sub-folders, so the folder layout does not matter.

## Which steps run from the Supplementary Tables workbook

The workbook `Antioxidant_ThyroidLiver_Supplementary_Tables_S1-S13.xlsx` (put it in this folder; it is found by the part
`Supplementary_Tables_S1-S13` of its name) is enough for the steps marked "workbook". For the other steps the public data of
the right-hand column have to be downloaded.

| Script | Reproduces | Runs from the workbook? | Otherwise download |
|---|---|---|---|
| `00_gwas_gene_sets_from_FUMA.R` | trait, thyroid-derived, liver-derived and shared gene sets | no (the workbook has the set sizes only; it is used for a check against S3a). The result of this step is provided: `gene_sets/gwas_gene_flags.tsv` | FUMA SNP2GENE output (`genes.txt`, `eqtl.txt`, `GenomicRiskLoci.txt`) of the 12 GWAS datasets, in `FUMA/<dataset>/` (or any folder whose path contains the dataset name or the GWAS Catalog ID); settings and job overview in `docs/FUMA_settings.md`. Large `eqtl.txt` files (TSH B: 670 MB) can be skipped on a second run with `RESUME=1` |
| `00b_magma_gene_set_enrichment.R` | S1 | no (S1 is its output; used for the check) | MAGMA v1.10, `NCBI37.3.gene.loc`, `g1000_eur.bed/.bim/.fam`, the 12 GWAS files, `gene_sets/gwas_gene_flags.tsv` |
| `00c_gene_overlap_hypergeometric.R` | exploratory overlap test (not in the submitted tables) | yes (set sizes of S3a), or from `gene_sets/gwas_gene_flags.tsv` | - |
| `01a_coloc_scan.R` | S2b (and the windows of S2a) | not needed: S2b is its output | the 12 GWAS files; R package `coloc` |
| `01b_coloc_window_tiers.R` | S2a | yes (from S2b), or from the output of `01a` | - |
| `01c_window_gene_mapping.R` | S4, S3b, gene-to-window table | yes (from S2a) plus `hg19_refgene.txt`; without the FUMA files 509 of the 526 rows of S4 and 151 of the 199 shared genes are found | `hg19_refgene.txt` (UCSC); FUMA `genes.txt` files (`FUMA_DIR`) for genes with former symbols or clone-based names (17 symbols of S4, 48 of the 199 shared genes); with them all rows are found |
| `02_*` ... `11_*` | S6 - S13 | see the table below | GEO and Ensembl files |
| `figures/*` | Figures 2 - 6, S1 - S3 | yes | optional public GWAS files for the regional panels |

## Public data

| Data | Source | Used by |
|---|---|---|
| Liver MRI GWAS: cT1 (GCST90446646-49), PDFF (GCST90446650-53); order unadjusted, BMI-, WFM-, VAT-adjusted | GWAS Catalog (https://www.ebi.ac.uk/gwas/); Nature Medicine 2024;30:3614-3623 | `00b_*`, `01a_*`, `figures/Fig2*`, `figures/FigS3*` |
| TSH GWAS: GCST010653 (dataset A; HUNT, Michigan Genomics Initiative and ThyroidOmics; manuscript references 40, 41) and GCST90296333 (dataset B; Williams AT et al., Nat Commun 2023;14:6713, doi:10.1038/s41467-023-42284-5; references 31, 42) | GWAS Catalog. For FUMA both files were reduced to variants with P < 1e-5 (`docs/FUMA_settings.md`); `00b_*` and `01a_*` are meant for the complete files (the colocalization scan needs all variants of a window) | `00b_*`, `01a_*`, `figures/Fig2*`, `figures/FigS3*` |
| FT3 / FT4 GWAS (`FT3_summary_data.txt.gz`, `FT4_summary_data.txt.gz`) | ThyroidOmics Consortium, "Thyroid function GWAS meta-analysis summary statistics (free thyroxine, FT4; free triiodothyronine, FT3)", 2024: https://genepi.med.uni-greifswald.de/thyroidomics/datasets/ . Study: Sterenborg RBTM et al., Multi-trait analysis characterizes the genetics of thyroid function and identifies causal associations with clinical implications, Nat Commun 2024;15:888, doi:10.1038/s41467-024-44701-9 (references 32 and 43 of the manuscript). Used under the terms of the consortium; the files are not redistributed here. Columns of the files as given to FUMA: `chr`, `position`, `rsid`, `Allele1`, `Allele2`, `Effect`, `StdErr`, `P.value` | `00b_*`, `01a_*`, `figures/Fig2*`, `figures/FigS3*` |
| `hg19_refgene.txt` | UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables): assembly hg19 (GRCh37), group "Genes and Gene Predictions", track "NCBI RefSeq", table `refGene`, output format "all fields from selected table" (or CSV with the fields name, chrom, strand, txStart, txEnd, name2); save as `hg19_refgene.txt`. The goldenPath file `refGene.txt.gz` of hg19, unpacked and renamed, is also accepted | `01c_*` |
| MAGMA v1.10 program, gene locations `NCBI37.3.gene.loc`, reference panel `g1000_eur.bed/.bim/.fam` (1000 Genomes Phase 3 EUR, build 37) | https://cncr.nl/research/magma/ (software and "auxiliary files") | `00b_*` |
| Mouse liver, vehicle vs T3 (GSE85793): per-sample HTSeq counts, `GSE85793_RAW` | GEO | `02_*`, `03_*` |
| Human hepatocytes, wild-type vs THRB knockout, vehicle vs T3 (GSE270233): count table `GSE270233_counts.csv` (12 samples, columns S1-S12 in the order of the GEO sample list) and the series matrix | GEO | `08_*` |
| Human liver biopsies, fibrosis stages: GSE130970 (series matrix and counts), GSE135251 (series matrix and `GSE135251_RAW.tar`) | GEO | `09_*` |
| Ensembl BioMart tables (release 116): mouse-to-human orthologues, human-to-mouse orthologues, all mouse genes (GRCm39) | Ensembl | `04_*`, `05_*`, `07_*` |

All GWAS files are on genome build GRCh37. `.gz` files can be used as they are (R package R.utils).

## Result tables kept in the repository

| File | Content | Made by | Read by |
|---|---|---|---|
| `gene_sets/gwas_gene_flags.tsv` | the 2,630 gene symbols of the GWAS gene sets, as FUMA writes them, with a 0/1 flag per trait (PDFF, cT1, TSH, FT3, FT4); thyroid-derived set = TSH, FT3 or FT4 (1,709 genes), liver-derived set = PDFF or cT1 (1,120 symbols = **1,119 genes**: ABHD16B and C20ORF135 are one gene that FUMA lists twice), shared = both (199 genes); union 2,630 symbols = **2,629 genes**, the numbers of Table S3a. See "Gene counts" in `docs/PROVENANCE.md` | `00_gwas_gene_sets_from_FUMA.R` (copy of the flags used for the manuscript) | `00b_*`, `00c_*`, `01c_*` |

## Files made by earlier steps (hand-over files)

| File | Made by | Read by |
|---|---|---|
| `0_count_matrix.csv`, `1_cleaned_count_matrix.csv`, `1_cleaned_coldata.csv` | `02_*` | `03_*`, `05_*` |
| `250912__DEG_all_genes.csv`, `261020__DEG_all_genes_with_stat.csv`, `GSE85793_mouse_liver_T3_DESeq2_all_genes.csv` | `03_*` | `04_*`, `05_*`, `07_*`, `figures/Fig2*` |
| `DEG_WT_T3_effect.csv`, `DEG_KO_T3_effect.csv`, `T3_effect_WT_KO_interaction_wide.csv` (written by `08_*` to `GSE270233_DESeq2/`) | `08_*` | `07_*`, `11_*` |
| `coloc_windows.tsv`, `coloc_pairs.tsv` | `01a_*` | `01b_*` |
| `window_tiers.tsv` | `01b_*` | `01c_*` |
| `261019_gene_locus_full_mapping.xlsx` (gene to window mapping) | `01c_*` | `05_*` |

## Supplementary Tables workbook

The figure scripts and the checks of steps 00 - 01 read `Antioxidant_ThyroidLiver_Supplementary_Tables_S1-S13.xlsx` (the
supplementary file of the manuscript). Put it in this folder; the scripts find it by the part
`Supplementary_Tables_S1-S13` of its name, or take its path from `SUPP_XLSX`.
