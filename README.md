# GeneDepresion

RNA-seq preprocessing and normalization workflow for GSE102556.

The current Step 1 script prepares region-specific expression matrices using one normalization path:

```text
filtered raw counts -> GeTMM -> log2(GeTMM + 1)
```

## Main Scripts

- `run_step1.R` - easiest entry point; runs the current recommended Step 1 workflow.
- `code/Step1_GSE102556_GeTMM_brain_area_groups.R` - recommended Step 1; groups exact GEO tissue labels into broader editable brain areas before GeTMM-log2 normalization.
- `code/Step2_Covariate_Exploration_FrontalCortex.R` - PCA/QC and covariate exploration for Frontal Cortex.
- `code/Step3_LMM_Frontal_Cortex_simple.R` - simple Frontal Cortex differential expression script using one mixed model per gene.
- `code/Step3_LM_Small_Groups_simple.R` - simple differential expression script for Insula, Ventral Striatum, and Hippocampal Formation.

## Brain-Area Mapping

The broad grouping is controlled by:

- `config/GSE102556_brain_area_mapping_template.csv`

Default grouping:

- OFC, dlPFC, Cg25 -> `Frontal_Cortex`
- aINS -> `Insula`
- NAc -> `Ventral_Striatum`
- Subiculum -> `Hippocampal_Formation`

Edit the mapping CSV and rerun the brain-area script if a different anatomical grouping is needed.

## Documentation

- `docs/Step1_GSE102556_GeTMM_log2_Notion.md`
- `docs/Data_Treatment_GSE102556_Notion.md`
- `docs/Full_Process_Simple_Explanation.md`

## Current Analysis Summary

The current workflow is:

```text
raw counts
-> metadata cleaning
-> brain-area grouping
-> GeTMM normalization
-> log2(GeTMM + 1)
-> PCA/QC inspection
-> regression-based differential expression
```

The current differential expression step tested:

- 24,287 genes in `Frontal_Cortex`
- 23,375 genes in `Insula`
- 24,275 genes in `Ventral_Striatum`
- 24,700 genes in `Hippocampal_Formation`

Using FDR < 0.05, no genes were significant in the current models.

Using FDR < 0.10, one Frontal Cortex gene passed the exploratory threshold.

## Inputs

The script expects local input files under:

```text
D:/genes/genedata
```

Required:

- `GSE102556_raw_counts_GRCh38.p13_NCBI.tsv`
- `Human.GRCh38.p13.annot.tsv`

Optional:

- `GSE102556_metadata_clean.csv`

If the optional metadata file is missing, the script fetches GEO metadata using:

```r
GEOquery::getGEO("GSE102556", GSEMatrix = TRUE)
```

## Outputs

The brain-area grouped Step 1 writes outputs to:

```text
D:/genes/genedata/Step1_GSE102556_GeTMM_brain_area_groups
```

Output folders are organized by file type:

```text
00_config/                 mapping used for the run
01_metadata/               matched metadata and per-group metadata
02_normalized_expression/  final log2(GeTMM + 1) matrices
03_filtered_raw_counts/    raw counts after gene filtering
04_normalization_factors/  TMM and effective-library-size tables
05_gene_annotation/        annotation for retained genes
06_repetition_reports/     repeated donor / gene-entry reports
07_qc/                     compact QC and model-decision summaries
08_rds/                    compact R object for later steps
```

Large input data and generated CSV outputs are intentionally ignored by Git.
