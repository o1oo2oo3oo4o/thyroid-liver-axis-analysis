
# thyroid-liver-axis-analysis

Code to draw the figures of the manuscript

> **Shared thyroid-liver genetics and hepatic triiodothyronine transcriptomics pinpoint thyroid hormone receptor beta-dependent candidate genes and a conserved xenobiotic-glutathione axis in steatotic liver disease**

The study combines thyroid-function and liver-imaging GWAS (colocalization, competitive gene-set analysis), hepatic T3 transcriptomics (mouse liver and human hepatocytes with and without *THRB*), plasma protein Mendelian randomization and liver biopsy cohorts. This repository contains the R scripts that draw Figures 2 to 6 and Supplementary Figures S1 to S3 from the Supplementary Tables of the manuscript.

Manuscript: *TBD (link or DOI after publication)* · Archived version of this code: *TBD (Zenodo DOI)*

## Repository contents

```
figures/                      R scripts, one per figure; each script runs on its own
data/                         input files (not tracked by git, see "Input data")
README.md  LICENSE  .gitignore
```

| Script | Figure | Main inputs |
|---|---|---|
| `Fig2_colocalization_and_T3_response_supp.R` | Figure 2: shared thyroid-liver signal and hormone responsiveness | Supplementary Tables S2a, S2b, S4; mouse DESeq2 results (panel C); public GWAS files (optional, panel B) |
| `Fig3_protein_candidates_supp.R` | Figure 3: protein-level candidate set | Supplementary Tables S6c, S7, S8 |
| `Fig4_THRB_dependency_supp.R` | Figure 4: effect of receptor beta loss in human hepatocytes | Supplementary Table S9a; standard errors from S9d or from `DEG_WT_T3_effect.csv`, `DEG_KO_T3_effect.csv`, `DEG_Interaction.csv` (panels B, C) |
| `Fig5_pathway_and_redox_supp.R` | Figure 5: functional composition and redox gene sets | Supplementary Tables S10b, S11a, S11b, S11d |
| `Fig6_mouse_human_concordance_supp.R` | Figure 6: mouse and human glutathione response | Supplementary Tables S9a, S9c, S11d (optional: `DEG_WT_T3_effect.csv`) |
| `FigS1_MAGMA_gene_set_supp.R` | Supplementary Figure S1: competitive gene-set results | Supplementary Table S1 |
| `FigS2_prior_sensitivity_supp.R` | Supplementary Figure S2: prior sensitivity of colocalization | Supplementary Tables S2a, S2b |
| `FigS3_window_layouts_supp.R` | Supplementary Figure S3: layouts of the 27 supported windows | Supplementary Tables S2a, S2b, S4; public GWAS files (optional) |

Every script writes its figures (`.pdf` and `.png`), a legend text file and `sessionInfo.txt` to its own folder `figures_out_supp_*` in the working directory.

## Requirements

R (tested with 4.3.3) and the packages `ggplot2`, `data.table`, `readxl`, `patchwork`, `ggrepel`, `scales`, `ragg`, plus `stringr` (Figure 5) and `cowplot` (Supplementary Figure S3):

```r
install.packages(c("ggplot2", "data.table", "readxl", "patchwork", "ggrepel", "scales", "ragg", "stringr", "cowplot"))
```

## Input data

Put the input files in `data/` (or set the environment variable `DATA_DIR` to the folder that holds them). Files are found by name, also inside sub-folders.

| File | Needed for | Source |
|---|---|---|
| Supplementary Tables workbook (file name containing `Supplementary_Tables_S1-S13`, `.xlsx`) | all figures | Supplementary material of the manuscript (set `SUPP_XLSX` to give the full path instead) |
| `DEG_WT_T3_effect.csv`, `DEG_KO_T3_effect.csv`, `DEG_Interaction.csv` (columns `gene`, `lfcSE`) | Figure 4B, C (if the workbook has no sheet S9d); Figure 6C (optional) | DESeq2 results of GEO [GSE270233](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270233) |
| `GSE85793_mouse_liver_T3_DESeq2_all_genes.csv` (columns `gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `pvalue`, `padj`) | Figure 2C | DESeq2 results of GEO [GSE85793](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85793) (Methods, section 2.4) |
| GWAS summary statistics, e.g. `GCST90446646.tsv` … `GCST90446653.tsv`, `GCST010653.txt`, `GCST90296333.tsv`, `FT3_summary_data.txt`, `FT4_summary_data.txt` | regional association plots of Figure 2B and Supplementary Figure S3 (optional) | see below |

Without the optional inputs the scripts still run: the regional association plots are replaced by window and gene layouts, and panels that need missing standard errors or the mouse results file are replaced by a note.

### Public data used in the study

| Data | Accession | Where to get it |
|---|---|---|
| TSH GWAS (dataset A) | GCST010653 | [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/) |
| TSH GWAS (dataset B) | GCST90296333 | NHGRI-EBI GWAS Catalog |
| Liver PDFF and cT1 GWAS | GCST90446646 to GCST90446653 | NHGRI-EBI GWAS Catalog |
| FT4 and FT3 GWAS | – | [ThyroidOmics Consortium](https://genepi.med.uni-greifswald.de/thyroidomics/datasets/) |
| Mouse liver RNA-seq, T3 vs vehicle | GSE85793 | NCBI GEO |
| Human hepatocytes, *THRB* knockout, T3 vs vehicle | GSE270233 | NCBI GEO |
| Liver biopsy cohorts | GSE135251, GSE130970 | NCBI GEO |

The ThyroidOmics files must not be redistributed or uploaded to a public repository (terms of the data provider), and the GWAS files are too large for this repository. Download them from the sources above.

## How to run

```bash
mkdir -p data                      # copy the input files into data/
Rscript figures/Fig6_mouse_human_concordance_supp.R
```

Each script prints the numbers that appear in the manuscript, for example:

| Script | Printed check |
|---|---|
| `Fig2_…` | 458 windows, 451 with a cross-domain result, 11 high and 16 moderate |
| `Fig4_…` | 15,297 tested genes, 799 wild-type responsive, 35 knockout responsive, 194 altered responses |
| `Fig6_…` | 14 of 18 genes concordant, one-sided binomial P = 0.015 |

## Scope

This repository contains the code that draws the figures from the Supplementary Tables. The upstream analyses (genome-wide association processing with FUMA, MAGMA gene-set analysis, colocalization with `coloc.abf`, Mendelian randomization, DESeq2 and the biopsy cohort analyses) are described in the Methods of the manuscript; *TBD: state here whether their scripts are added to this repository.*

## Data availability

The Supplementary Tables S1 to S13 contain the full results and gene-set membership. All public datasets are listed above with their accession numbers.

## License and citation

Code: MIT License (see `LICENSE`). If you use this code, please cite the manuscript: *TBD*.
