# Thyroid-liver axis analysis

Analysis code for the manuscript on genetic and transcriptomic links between thyroid hormone signalling and metabolic
dysfunction-associated steatotic liver disease (MASLD/MASH): shared thyroid-liver GWAS loci, hepatic T3 response in
mouse liver, THRB dependence in human hepatocytes, and redox gene-set and fibrosis-stage analyses in human liver.

The repository contains the R scripts that produced the supplementary tables and figures. It does not contain any
data; all inputs are public (see [`data/README.md`](data/README.md)) or are tables of the supplementary workbook.

## Layout

| Folder | Content |
|---|---|
| `scripts/` | analysis scripts, numbered in the order of the pipeline |
| `figures/` | one standalone script per figure; each reads the Supplementary Tables workbook (and, where stated, public files) |
| `docs/PROVENANCE.md` | which script and which input produced each supplementary table |
| `data/` | place the input files here (not tracked, except `gene_sets/gwas_gene_flags.tsv`); see `data/README.md` |
| `docs/FUMA_settings.md` | settings and job overview of the 12 FUMA web runs (from the `params.config` files of the jobs) |

## Pipeline

| Step | Script | What it does | Supplementary table |
|---|---|---|---|
| 0 | `00_gwas_gene_sets_from_FUMA.R` | FUMA output of the 12 GWAS datasets to cis-eQTL genes, trait gene sets, thyroid-derived / liver-derived / shared genes (FUMA itself is a web tool: `docs/FUMA_settings.md`); reproduces `data/gene_sets/gwas_gene_flags.tsv` | S3a (counts) |
| 0b | `00b_magma_gene_set_enrichment.R` | MAGMA gene-level analysis of the 12 GWAS datasets and competitive gene-set tests in both directions | S1 |
| 0c | `00c_gene_overlap_hypergeometric.R` | exploratory hypergeometric and permutation test of the shared-gene overlap (not part of the submitted tables) | - |
| 1a | `01a_coloc_scan.R` | colocalization scan (coloc.abf) of liver MRI and thyroid GWAS: 458 windows, 66 dataset pairs | S2b |
| 1b | `01b_coloc_window_tiers.R` | window tiers (high / moderate), supporting pairs, prior sensitivity | S2a |
| 1c | `01c_window_gene_mapping.R` | genes in the windows (GRCh37) and colocalization status of the 199 shared genes | S4, S3b |
| 2 | `02_build_count_matrix_GSE85793.R` | raw count matrix of the mouse liver T3 dataset (GSE85793) from the GEO per-sample files | - |
| 3 | `03_deseq2_mouse_liver_T3.R` | DESeq2, T3 vs vehicle, apeglm shrinkage | S6 (input), S10, S11 |
| 4 | `04_reortholog_synonym.R` | mouse to human orthologue remapping of the responsive genes (Ensembl release 116) | S6a |
| 5 | `05_ortholog_audit.R` | retrieval status of every positional gene in the mouse dataset | S6b |
| 6 | `06_pathway_ORA.R` | over-representation analysis (GO, KEGG, Reactome) | S10 |
| 7 | `07_redox_cameraPR.R` | redox gene sets: Fisher exact test and cameraPR in mouse and human hepatocytes | S11 |
| 8 | `08_deseq2_human_hepatocyte_THRB_KO.R` | DESeq2 (2 x 2 design) of the human hepatocyte THRB-knockout experiment (GSE270233) | S9a |
| 9 | `09_human_liver_fibrosis.R` | fibrosis-stage analysis in two human liver biopsy cohorts (GSE130970, GSE135251) | S12, S13 |
| 10 | `10_fibrosis_positive_controls.R` | positive-control genes of the fibrosis-stage analysis | S13d |
| 11 | `11_THRB_candidate_enrichment.R` | enrichment of the candidates among genes with a THRB-dependent T3 response | S9b |

Steps that were run outside R (FUMA web server, Mendelian randomization, plasma proteomics, STRING) are described in the
Methods of the manuscript; see `docs/PROVENANCE.md` for what is and is not included here. MAGMA is a stand-alone program that
`00b_magma_gene_set_enrichment.R` calls.

### Running from the Supplementary Tables

Steps 0c, 1b and 1c, and the checks of all steps of 0 - 1c, run from the Supplementary Tables workbook alone (step 1c also
needs the UCSC refGene table). Steps 0, 0b and 1a need the downloaded public data; `data/README.md` lists, for every script,
whether the workbook is enough and which files to download otherwise. All inputs are found by file name, no path has to be
edited:

```bash
DATA_DIR=data Rscript scripts/01b_coloc_window_tiers.R          # uses S2b of the workbook if no scan output is present
DATA_DIR=data Rscript scripts/01c_window_gene_mapping.R         # S4 and S3b; needs data/hg19_refgene.txt
DATA_DIR=data Rscript scripts/00c_gene_overlap_hypergeometric.R
```

## How to run

Each script is run from the repository root; the folder with the input files is given with `DATA_DIR` (default `data`),
and results are written to `results/<step>/` (or to `OUT_DIR`):

```bash
DATA_DIR=data Rscript scripts/03_deseq2_mouse_liver_T3.R
```

Figure scripts are run the same way; they find the Supplementary Tables workbook in `DATA_DIR` (or take its path from
`SUPP_XLSX`) and write the figure files to a `figures_out_*` folder in the current directory:

```bash
DATA_DIR=data Rscript figures/Fig5_pathway_and_redox_supp.R
```

The usage line at the top of every script lists its inputs, outputs and options. Hand-over files between steps keep the
names of the files used for the manuscript (for example `261020__DEG_all_genes_with_stat.csv`); the date prefixes are
file labels, not versions to be changed.

## Software

R 4.5.2 was used for the analyses. MAGMA v1.10 is called by `00b`. Main packages: coloc (01a), R.utils, readxl, DESeq2 1.50.2, apeglm, ashr, limma 3.66.0, clusterProfiler 4.18.4,
ReactomePA, org.Hs.eg.db, msigdbr, biomaRt, data.table, openxlsx, readr, dplyr. The figure scripts additionally use
ggplot2, patchwork, ggrepel, scales and readxl. The session information of each run is written next to its results.

## Reproduction checks

Where an original result file was available, the scripts were checked against it; these checks are listed in
`docs/PROVENANCE.md`. The scripts were cleaned (English comments and messages, configurable paths, removal of test hooks)
without changing the analysis code.

## License

MIT, see `LICENSE`.
