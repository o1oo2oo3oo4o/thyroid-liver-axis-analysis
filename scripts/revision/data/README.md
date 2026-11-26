# Input files of the revision analyses

Put the files of the table below in this folder (`scripts/revision/data`). Only this README and the three small gene-set
tables are tracked by git (see `.gitignore`); counts, result tables and the supplementary workbook are not. Output goes to
`scripts/revision/results/` (not tracked).

## Files that are part of the repository

| File | Content | Source |
|---|---|---|
| `redox_gene_sets_human.csv` | columns `set`, `human_symbol`: the six redox gene sets (human symbols) | Supplementary Table S11 (membership) |
| `redox_gene_sets_mouse.csv` | columns `set`, `mouse_symbol`: the same sets with mouse symbols | Supplementary Table S11 (membership) |
| `glutathione_sensitivity_conditions_mouse.csv` | columns `condition`, `gene`, `included`, `tested`: membership of the glutathione set under the sensitivity conditions | gene-set membership of the supplementary tables (glutathione sensitivity conditions) |

## Files to add

| File | Content | How to obtain | Read by |
|---|---|---|---|
| `GSE85793_counts.csv`, `GSE85793_meta.csv` | raw counts (gene symbol x 8 samples: C1 - C4 vehicle, T1 - T4 T3) and sample table | `00_prepare_GSE_inputs.R` from GEO GSE85793 (`GSE85793_RAW.tar`, series matrix) | `03`, `05` |
| `GSE130970_counts.csv`, `GSE130970_meta.csv`, `GSE135251_counts.csv`, `GSE135251_meta.csv` | raw counts (HGNC symbol x samples) and sample tables (fibrosis stage, NAS, age and sex where deposited) | `00_prepare_GSE_inputs.R` from GEO GSE130970 and GSE135251 | `02`, `06` |
| `GSE270233_counts.csv`, `GSE270233_meta.csv` | raw counts of the human hepatocyte experiment and sample table | `00_prepare_GSE_inputs.R` from GEO GSE270233 | `04` (maturity table, optional) |
| `DEG_hepT3_raw.csv` | one row per gene (15,915): `log2FoldChange`, `lfcSE`, `pvalue`, `padj` for the wild-type (`_WT`), knockout (`_KO`) and interaction (`_int`) contrasts, `lfc_shrunk_WT`, `lfc_shrunk_KO`, `log2FoldChange_KOvsWT_basal`, `padj_KOvsWT_basal` | DESeq2 result tables written by `../08_deseq2_human_hepatocyte_THRB_KO.R` (`DEG_WT_T3_effect.csv`, `DEG_KO_T3_effect.csv`, `DEG_Interaction.csv`, `DEG_KO_vs_WT_basal.csv`) joined by gene; the three contrasts without the shrunken and baseline columns are also in `T3_effect_WT_KO_interaction_wide.csv` | `01`, `04`, `12` |
| `mouse_T3_results.csv` | mouse T3 result with human orthologues: `gene`, `human_genes`, `homology_types`, `n_human`, `log2FC`, `padj`, `pvalue`, `responsive` | table of the study (DESeq2 result of `../03_deseq2_mouse_liver_T3.R` joined with the orthologue table of `../04_reortholog_synonym.R`). **The `human_genes` column of the file used for the response letters is the synonym-resolved sensitivity mapping (498 foreground / 11,290 background genes), not the main mapping of the manuscript (448 / 10,085 genes).** This affects the microsomal P450 set in `01` (section B4, one gene) and three genes of the glutathione annotation in `04` (Ggt1, Gpx6, Gstp2); `03` and `05` do not use the column | `01` (optional, B4), `03` (optional comparison), `04` |
| `ortholog_mouse_status_all.tsv` | main-analysis ortholog table (448 foreground / 10,085 background genes) | written by `../04_reortholog_synonym.R` | `05` (required), `11` (optional) |
| `Antioxidant_ThyroidLiver_Supplementary_Tables_S1-S13.xlsx` | Supplementary Tables workbook, sheet `Table S2` (window summary on the left, pair table on the right) | supplementary file of the manuscript | `07` |
| `Antioxidant_ThyroidLiver_Supplementary_Tables_revised.xlsx` | revised Supplementary Tables workbook with the sheets `Table S7`, `Table S14`, `Table S15`, `Table S16` (any file name: set `supp_revised_file` in `config.R` or the environment variable `SUPP_XLSX`) | supplementary file of the revised manuscript | `13` |
| `Ayers2014_TRb_peaks.csv` | columns `chr`, `start`, `end` (hg19): HepG2 TR-beta BioChIP-seq peaks | Supplementary Table S5 of Ayers et al., PLoS ONE 2014, doi:10.1371/journal.pone.0081186, saved as CSV | `10` |
| `tss_hg19.csv` | columns `gene`, `chr`, `tss`, `strand` for all genes | `10a_make_tss_hg19.R` | `10` |
| `250912__DEG_all_genes.csv` | mouse T3 DESeq2 result (all 12,135 tested genes; `gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `pvalue`, `padj`) | written by `../03_deseq2_mouse_liver_T3.R` | `11` |
| `DEG_WT_T3_effect.csv`, `T3_effect_WT_KO_interaction_wide.csv` | human hepatocyte DESeq2 results (wild-type T3 effect; the three contrasts side by side) | written by `../08_deseq2_human_hepatocyte_THRB_KO.R` | `11` (`T3_effect_WT_KO_interaction_wide.csv` optional) |

Count-matrix format: first column `gene` (symbol), the other columns are samples; the column names must match the `sample`
column of the sample table. Column names of the sample tables can be changed in `config.R` (`meta_cols`).

`10_TRb_ChIP_annotation.R` and `12_recompute_panel_stats_IVW.R` read the result files of the earlier scripts in
`results/` (`results/06_biopsy_F4/`, `results/10_TRb_ChIP/`); run `06` and `10` first.
