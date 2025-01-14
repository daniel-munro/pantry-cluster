# Pantry Data

## Data descriptions

### General info

`info/` contains the following metadata:
- `tissues.phenos.txt`: A list of all tissues with RNA phenotypes
- `tisssues.pheast.txt`: A list of all tissues with cis-QTL and xTWAS results
- `phenotypes_per_tissue.tsv.gz`: A table of all RNA phenotypes per tissue used for genetic analyses
- `traits.txt`: A list of GWAS traits tested in xTWAS
- `gwas_metadata.txt`: Metadata originally included with the GWAS traits summary stats

### RNA phenotype tables

RNA phenotypes were quantified for Geuvadis and all 54 GTEx tissues, including the five tissues excluded from downstream genetic analyses due to small sample size. The tables are in BED format, with the first three columns specifying the TSS coordinate of each phenotype's gene, i.e., the input format for `tensorQTL`. These are the raw phenotypes, i.e. the actual ratios or quantities. Processed tables were prepared for QTL mapping and xTWAS, but are not included here. They are the result of quantile-normalizing samples to the average empirical distribution, followed by rank-based inverse-normal transformation per phenotype. For GTEx tissues, we included only samples with genotypes, converted sample IDs to individual IDs to match the genotypes, and omitted the five tissues originally excluded from the GTEx QTL mapping.

### Covariates

For each tissue-modality RNA phenotype table, the covariates include the first 20 principal components (PCs) of the phenotype table and the first 5 PCs of the genotype alternative allele count matrix. The covariates for the `combined` pseudo-modality are computed in the same way from the corresponding `combined` phenotype table, and are not concatenations of each modality's covariates. We include two formats for each covariate table, one compatible with tensorQTL (`*.covar.tsv`), and one in plink format for compatibility with FUSION (`*.covar.plink.tsv`).

We include covariates and cis-QTL results for a pseudo-modality called `combined`, which is a concatenation of phenotypes from all six modalities. We used these for 'combined-modality mapping' of cis-QTLs, but not for xTWAS.

### cis-QTLs

`QTLs/` includes output files from tensorQTL run in `cis` mode (`*.cis_qtl.txt.gz`, top association per gene, even if not significant) and `cis_independent` mode (`*.cis_independent_qtl.txt.gz`, all conditionally-independent cis-QTLs). Files for the pseudo-modality `combined` are the result of running tensorQTL on combined-modality phenotype tables and group files, so that cis-QTLs are conditionally independent across modalities per gene.

### xTWAS

A set of TWAS transcriptomic models was generated for each tissue-modality pair using FUSION. The weights and summary files are provided in compressed archives in `TWAS_weights/`, and the outputs of the FUSION script `FUSION.profile_wgt.R`, i.e. tables giving cis-heritability and model stats per phenotype, are also included separately (`*.twas_weights.profile`).

TWAS associations were calculated for each set of models against summary stats for each of [114 GWAS traits previously harmonized to GTEx v8 variant referece](https://zenodo.org/records/3629742#.XjCh9OF7m90). Outputs of `FUSION.assoc_test.R` were concatenated across chromosomes to get one table per trait-tissue-modality combination. Tables are provided in compressed archives per trait in `TWAS_associations/`.

### Aggregated results

For convenience, some files with aggregated, filtered, and otherwise processed results are provided in `processed/`:

- `genes_pcg_lncrna.txt`: Ensembl gene IDs for the subset of genes that are protein-coding or lncRNA. We include all genes in the Ensembl annotations in phenotyping and genetic analyses, and filter to this set of genes for all files in the `processed/` directory.
- `geuvadis.sep.assoc.tsv.gz`: Concatenation of the top associations per gene (`*.cis_qtl.txt.gz`) for the six modalities for Geuvadis, with a subset of the original columns.
- `geuvadis.sep.qtls.tsv.gz`: Concatenation of the conditionally independent cis-QTLs (`*.cis_independent_qtl.txt.gz`) for the six modalities for Geuvadis, with a subset of the original columns.
- `geuvadis.comb.qtls.tsv.gz`: The `GEUVADIS.combined.cis_independent_qtl.txt.gz` file processed to match the format of `geuvadis.sep.qtls.tsv.gz`, for easier comparison of cis-QTLs from separate vs. combined-modality mapping.
- `geuvadis.hsq.tsv.gz`: Concatenation and processing of the `*.twas_weights.profile` TWAS transcriptomic model summary tables for all modalities for Geuvadis.
- `geuvadis.twas.tsv.gz`: Concatenation and processing of all TWAS associations (even if not significant) for all modalities and traits for Geuvadis.
- `gtex.sep.assoc.tsv.gz`: Concatenation of the top associations per gene (`*.cis_qtl.txt.gz`) for the six modalities for all GTEx tissues, with a subset of the original columns.
- `gtex.sep.qtls.tsv.gz`: Concatenation of the conditionally independent cis-QTLs (`*.cis_independent_qtl.txt.gz`) for the six modalities for all GTEx tissues, with a subset of the original columns.
- `gtex.comb.qtls.tsv.gz`: Concatenation of the `{tissue}.combined.cis_independent_qtl.txt.gz` files for all GTEx tissues processed to match the format of `gtex.sep.qtls.tsv.gz`, for easier comparison of cis-QTLs from separate vs. combined-modality mapping.
- `gtex.hsq.tsv.gz`: Concatenation and processing of the `*.twas_weights.profile` TWAS transcriptomic model summary tables for the six modalities for all GTEx tissues.
- `gtex.twas_hits.tsv.gz`: Concatenation and processing of all TWAS hits for all modalities and traits for all GTEx tissues. Unlike `geuvadis.twas.tsv.gz`, this is filtered to hits passing a genome-wide P-value threshold, adjusted for number of RNA modalities, of 8.33e-9 (5e-8 / 6), to reduce file size.

## Files

This is a summarized outline of the data files. Generally, files exist for every combination of values for the variables in curly braces, with exceptions, for example not all phenotyped tissues were used for genetic analyses. `GEUVADIS` is included as a tissue alongside the GTEx tissues. See `file_tree.txt` for the complete tree, or see `check_repo.py` for code that generates the exact list of files, and run it to check for their existence.

```
├── info
│   ├── tissues.phenos.txt
│   ├── tissues.pheast.txt
│   ├── phenotypes_per_tissue.tsv.gz
│   ├── traits.txt
│   └── gwas_metadata.txt
├── RNA_phenotypes
│   └── {tissue}.{modality}.unnorm.bed.gz
├── covariates
│   ├── {tissue}.{modality}.covar.tsv (compressed into covariates.tar.bz2)
│   └── {tissue}.{modality}.covar.plink.tsv (compressed into covariates.plink.tar.bz2)
├── QTLs
│   ├── {tissue}.{modality}.cis_qtl.txt.gz
│   └── {tissue}.{modality}.cis_independent_qtl.txt.gz
├── TWAS_weights
│   ├── {tissue}.{modality}.twas_weights.tar.bz2
│   └── {tissue}.{modality}.twas_weights.profile
├── TWAS_associations
│   └── {trait}/fusion.{tissue}.{modality}.{trait}.tsv (compressed into {trait}.tar.bz2)
└── processed
    ├── genes_pcg_lncrna.txt
    ├── {geuvadis,gtex}.sep.assoc.tsv.gz
    ├── {geuvadis,gtex}.{sep,comb}.qtls.tsv.gz
    ├── {geuvadis,gtex}.hsq.tsv.gz
    ├── geuvadis.twas.tsv.gz
    └── gtex.twas_hits.tsv.gz
```
