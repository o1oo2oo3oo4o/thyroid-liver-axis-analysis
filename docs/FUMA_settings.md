# FUMA settings (SNP2GENE)

FUMA is a web tool (https://fuma.ctglab.nl); there is no script for it. The result files of each run (`genes.txt`,
`eqtl.txt`, `GenomicRiskLoci.txt`) are read by `scripts/00_gwas_gene_sets_from_FUMA.R`. This page records how the 12 FUMA
jobs were run. Every value below was read from the `params.config` file of the job (one per job, from the FUMA result
download) and compared between the 12 jobs; the settings of the 12 jobs are identical except for the input file, its
column names and the sample size.

## Jobs

One SNP2GENE job per GWAS dataset, run in July - August 2025. The first column is the folder name that
`00_gwas_gene_sets_from_FUMA.R` expects (a folder whose path contains the GWAS Catalog ID is accepted as well).

| Folder (dataset) | GWAS | File given to FUMA | N entered in FUMA | FUMA run | Genomic risk loci | Genes in `genes.txt` | Genes passing the cis rule |
|---|---|---|---|---|---|---|---|
| `PDFF_650` | PDFF unadjusted | `GCST90446650.tsv.gz` | 36,394 | 2025-07-15 | 27 | 494 | 337 |
| `PDFF_651` | PDFF BMI-adjusted | `GCST90446651.tsv.gz` | 35,146 | 2025-07-15 | 38 | 701 | 473 |
| `PDFF_652` | PDFF WFM-adjusted | `GCST90446652.tsv.gz` | 34,367 | 2025-07-15 | 37 | 754 | 515 |
| `PDFF_653` | PDFF VAT-adjusted | `GCST90446653.tsv.gz` | 33,540 | 2025-07-16 | 31 | 626 | 446 |
| `cT1_646` | cT1 unadjusted | `GCST90446646.tsv.gz` | 30,481 | 2025-07-15 | 25 | 481 | 332 |
| `cT1_647` | cT1 BMI-adjusted | `GCST90446647.tsv.gz` | 29,312 | 2025-07-15 | 21 | 488 | 344 |
| `cT1_648` | cT1 WFM-adjusted | `GCST90446648.tsv.gz` | 28,638 | 2025-07-15 | 24 | 490 | 350 |
| `cT1_649` | cT1 VAT-adjusted | `GCST90446649.tsv.gz` | 27,782 | 2025-07-15 | 30 | 538 | 372 |
| `FT3` | free T3 | `FT3_summary_data.txt.gz` | 57,107 | 2025-08-04 | 18 | 192 | 149 |
| `FT4` | free T4 | `FT4_summary_data.txt.gz` | 117,166 | 2025-08-04 | 99 | 1,475 | 1,003 |
| `TSH_A` | TSH, dataset A (GCST010653) | `GCST010653_filtered_p_1e-5.txt` | 55,342 | 2025-07-21 | 101 | 1,842 | 891 |
| `TSH_B` | TSH, dataset B (GCST90296333) | `GCST90296333_filtered_p_1e-5.txt` | 124,358 | 2025-07-22 | 224 | 3,690 | 2,322 |

"Genes passing the cis rule" is the result of the three filters described at the end of this page; it equals the
`FINAL cis PASS` counts that were written when the analysis was run.

The N entered in FUMA is the sample size used for the FUMA run and for the colocalization scan (`01a_coloc_scan.R` uses the
same numbers). For FT3, FT4 and the two TSH datasets it differs from the sample size printed in Table 1 of the manuscript
(59,061 / 119,120 / 119,715 / 247,107); the numbers for the liver datasets are the same in both places. FUMA uses N only for
the optional MAGMA gene analysis inside FUMA, which is not used here, so the FUMA gene lists do not depend on it.

## Input files

| Item | Value |
|---|---|
| Genome build of the summary statistics | GRCh37 (hg19) for all 12 files (`GRCh38 = 0`) |
| Liver files (`GCST904466xx.tsv.gz`) | GWAS Catalog files; columns `chromosome`, `base_pair_location`, `ID` (rsID), `effect_allele`, `other_allele`, `beta`, `standard_error`, `p_value` |
| FT3, FT4 (`FT3_summary_data.txt.gz`, `FT4_summary_data.txt.gz`) | ThyroidOmics Consortium files; columns `chr`, `position`, `rsid`, `Allele1` (effect), `Allele2`, `Effect`, `StdErr`, `P.value` |
| TSH A | `GCST010653_filtered_p_1e-5.txt`; columns `Chr`, `Pos`, `Allele1` (effect), `Allele2`, `Effect`, `SE`, `P-value`; no rsID column |
| TSH B | `GCST90296333_filtered_p_1e-5.txt`; columns `chromosome`, `base_pair_location`, `effect_allele`, `other_allele`, `beta`, `standard_error`, `p_value`; no rsID column |
| Pre-filtering | the two TSH files were reduced to variants with P < 1e-5 before the upload (as the file names say); the other ten files carry their original names. No script for this step is in the repository: to repeat it, keep the rows of the two public TSH files with P < 1e-5. The lead variants (P < 1e-6) are not affected, but FUMA then only sees variants with P < 1e-5 as GWAS-tagged candidate variants (`gwasP = 0.05` has no effect below that), so the TSH gene lists would differ if the unfiltered files were uploaded |

## SNP2GENE parameters (`params.config`, identical in the 12 jobs)

| Parameter (`params.config` key) | Value |
|---|---|
| FUMA version | v1.5.2 |
| MAGMA version inside FUMA | v1.08 (switched on, `magma = 1`; its output is not used) |
| GWAS Catalog / ANNOVAR versions | e0_r2022-11-29 / 2017-07-17 |
| Ensembl annotation (`ensembl`) | v102 |
| Gene types (`genetype`) | protein_coding |
| LD reference panel (`refpanel`, `pop`) | 1000 Genomes Phase 3, EUR |
| Minimum MAF (`MAF`) | 0.01 |
| Maximum P of lead / independent significant variants (`leadP`) | 1e-6 |
| P-value threshold of GWAS-tagged variants (`gwasP`) | 0.05 |
| r2 of independent significant variants (`r2`) | 0.2 |
| r2 of lead variants (`r2_2`) | 0.1 |
| Merge distance of LD blocks into one genomic risk locus (`mergeDist`) | 250 kb |
| Reference-panel SNPs added to the annotation (`refSNPs`) | 1 (yes) |
| MHC region (`exMHC`, `MHCopt`, `extMHC`) | `exMHC = 1`, `MHCopt = annot`, no extended MHC region. The manuscript Methods do not mention this setting. In the downloaded files `genes.txt` still lists genes that overlap chr6:28,477,797-33,448,354 (the FUMA default MHC coordinates) in 9 of the 12 jobs (1 - 25 genes), so MHC genes were not removed from the FUMA output; the exact meaning of the option `annot` is given in the FUMA tutorial ("Exclude MHC region") |
| Positional mapping (`posMap`) | on, window 10 kb (`posMapWindowSize = 10`), no CADD / RegulomeDB / chromatin-state filter |
| eQTL mapping (`eqtlMap`) | on; 237 eQTL data sets selected (GTEx v6/v7/v8, eQTLGen cis and trans, eQTLcatalogue, BIOS, BloodeQTL, CMC, PsychENCODE, BRAINEAC, DICE, scRNA-eQTLs, MuTHER, xQTLServer, InsPIRE, TIGER and others; the manuscript calls this "all eQTL resources integrated in FUMA"); `eqtlMapSig = 1` (FUMA option: use only significant SNP-gene pairs; the FDR is the column `FDR` of `eqtl.txt`), `eqtlMapP = 1`, no CADD / RegulomeDB / chromatin-state filter. The list of data sets is the same in all 12 jobs (identical `eqtlMaptss` string) |
| Chromatin interaction mapping (`ciMap`) | off |
| Genome build (`GRCh38`) | 0 (GRCh37) |

## Agreement with the manuscript

The Methods of the manuscript state: 1000 Genomes Phase 3 EUR, Ensembl v102, independent significant variants at P < 1e-6
with r2 < 0.2, lead variants at r2 < 0.1, MAF >= 0.01, positional mapping within 10 kb and eQTL mapping against all eQTL
resources integrated in FUMA. All of this agrees with `params.config`. Not stated in the manuscript, but recorded in
`params.config`: FUMA v1.5.2, gene type protein_coding, merge distance 250 kb, MHC option (`annot`) and that chromatin
interaction mapping was off.

## How the FUMA results are used

`00_gwas_gene_sets_from_FUMA.R` applies three filters to the files downloaded from FUMA. They are not FUMA settings:

1. `genes.txt`: keep genes with `eqtlMapminQ` < 0.05 (minimum eQTL FDR over all eQTL pairs of the gene);
2. the gene lies on the chromosome of its genomic risk locus (`GenomicRiskLoci.txt`);
3. every eQTL pair of the gene in `eqtl.txt` is on the chromosome of the gene (cis only).

The trait gene set is the union over the datasets of the trait (PDFF: 4, cT1: 4, FT3: 1, FT4: 1); TSH is the intersection
of dataset A and dataset B. Genes are counted by upper-case symbol (FUMA lists `C6orf165` and `C6ORF165` as two Ensembl
genes; PIK3R2, SOX7, ALG9, CLN3, NDUFB8 and ZNF286A each have two Ensembl IDs with the same symbol).

Result with the FUMA files of the manuscript (checked on the real files, see `docs/PROVENANCE.md`): PDFF 786, cT1 486, TSH 828,
FT3 149 and FT4 999 genes; thyroid-derived set 1,709, shared 199, liver-derived set 1,120 symbols and union 2,630 symbols. FUMA
lists ABHD16B and its old name C20ORF135 as two genes in the PDFF runs; counted as one gene they give the 1,119 liver genes and the
union of 2,629 of Table S3a (`docs/PROVENANCE.md`, "Gene counts").
