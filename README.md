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

## Simple Step-by-Step Explanation

This section explains the whole data process from normalization until the current
differential expression analysis.

### 1. Start With Raw RNA-seq Counts

The starting data are raw RNA-seq gene counts from GSE102556.

The raw count table has:

- genes in rows
- samples in columns
- count values showing how many sequencing reads mapped to each gene

Raw counts cannot be compared directly across samples because samples can have
different sequencing depths, and genes have different lengths.

### 2. Clean and Standardize Metadata

The workflow loads sample metadata from GEO.

Important metadata columns are cleaned into standard names:

- `SampleID`
- `Diagnosis`
- `Age`
- `Sex`
- `RIN`
- `OriginalRegion`
- `BrainID`

`BrainID` identifies the donor/person. This matters because the same donor can
appear more than once if multiple brain regions were sampled.

### 3. Group Brain Regions

The original tissue labels are grouped into broader analysis groups.

| Original region | Analysis group |
|---|---|
| dlPFC | Frontal_Cortex |
| OFC | Frontal_Cortex |
| Cg25 | Frontal_Cortex |
| Anterior insula | Insula |
| Nucleus accumbens | Ventral_Striatum |
| Subiculum | Hippocampal_Formation |

This mapping is controlled by:

```text
config/GSE102556_brain_area_mapping_template.csv
```

### 4. Match Expression Samples to Metadata

The workflow keeps only samples that are present in both:

- the expression count matrix
- the metadata table

This avoids mismatches between expression columns and metadata rows.

### 5. Filter and Annotate Genes

Genes are matched to the annotation file so gene length can be used.

Very short or unusable genes are removed before normalization.

### 6. Normalize With GeTMM-log2

The project uses one normalization path:

```text
raw counts
-> gene-length adjustment
-> TMM library-size adjustment
-> GeTMM
-> log2(GeTMM + 1)
```

GeTMM adjusts for:

- gene length
- sequencing-depth differences between samples

The `log2(GeTMM + 1)` transform makes expression values easier to model and
keeps zero values valid.

### 7. Save Normalized Expression and Metadata

Step 1 writes the main downstream files:

```text
normalized_GeTMM_log2_<group>.csv
metadata_<group>.csv
```

The expression files have:

- genes in rows
- samples in columns

The metadata files have:

- one row per sample

### 8. Run QC and PCA

PCA/QC was used to inspect sample structure before differential expression.

The important finding was that `OriginalRegion` strongly matters inside
`Frontal_Cortex`, because that group combines dlPFC, OFC, and Cg25.

Therefore, `OriginalRegion` is included as a covariate in the Frontal Cortex
model.

### 9. Choose the Differential Expression Models

Frontal Cortex uses a mixed model because it can contain repeated donors:

```r
expression_value ~ Diagnosis + Age + Sex + RIN + OriginalRegion + (1 | BrainID)
```

The other groups use simple linear models:

```r
expression_value ~ Diagnosis + Age + Sex + RIN
```

The model is fit once per gene.

### 10. Correct for Testing Thousands of Genes

Each gene gets a diagnosis p-value.

Because thousands of genes are tested, the workflow corrects p-values using
Benjamini-Hochberg FDR.

FDR means false discovery rate.

Simple interpretation:

- p-value: evidence for one gene by itself
- FDR: corrected evidence after testing thousands of genes
- FDR < 0.05: strict significance threshold
- FDR < 0.10: exploratory threshold

### 11. Current Result

The current models tested:

| Group | Samples used | Genes tested | FDR < 0.05 | FDR < 0.10 |
|---|---:|---:|---:|---:|
| Frontal_Cortex | 122 | 24,287 | 0 | 1 |
| Insula | 48 | 23,375 | 0 | 0 |
| Ventral_Striatum | 48 | 24,275 | 0 | 0 |
| Hippocampal_Formation | 43 | 24,700 | 0 | 0 |

Using the strict FDR < 0.05 threshold, no genes were significant.

Using the exploratory FDR < 0.10 threshold, one Frontal Cortex gene passed.

### 12. Important Caution

The Frontal Cortex mixed model completed, but it printed many singular-fit
messages and two convergence warnings.

This means the output table was created, but the Frontal Cortex mixed-model
results should be interpreted carefully.

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
