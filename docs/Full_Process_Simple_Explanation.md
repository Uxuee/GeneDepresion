# Full Process Simple Explanation

This document explains, in simple terms, what has happened to the GSE102556 data
from normalization until the current differential expression step.

The goal of the project is to compare gene expression between depression
samples and control samples while accounting for important sample information
such as age, sex, RNA quality, brain region, and repeated donors.

## 1. Starting point

The project starts from RNA-seq raw count data.

The raw count file has:

- one row per gene entry
- one column per sample
- count values showing how many RNA-seq reads were assigned to each gene

The project also uses a gene annotation file with gene lengths.

Gene length is needed because longer genes naturally collect more reads than
shorter genes.

## 2. Metadata cleaning

The sample metadata was loaded from GEO.

The metadata contains sample information such as:

- sample ID
- diagnosis
- age
- sex
- RIN
- pH
- PMI
- original brain region

The project created clean standard columns such as:

```text
SampleID
Diagnosis
Age
Sex
RIN
OriginalRegion
BrainID
```

`BrainID` identifies the donor/person.

This matters because some donors can appear more than once if multiple brain
regions were sampled from the same person.

## 3. Brain-region grouping

The original GEO tissue labels were mapped into broader analysis groups.

The broad groups are:

| Original region | Analysis group |
|---|---|
| Dorsolateral prefrontal cortex | Frontal_Cortex |
| Orbitofrontal cortex | Frontal_Cortex |
| Cingulate gyrus 25 | Frontal_Cortex |
| Anterior insula | Insula |
| Nucleus accumbens | Ventral_Striatum |
| Subiculum | Hippocampal_Formation |

This grouping is controlled by:

```text
config/GSE102556_brain_area_mapping_template.csv
```

The Frontal Cortex group is special because it combines three original regions.
For that reason, `OriginalRegion` is kept as a model covariate later.

## 4. Matching samples

The code kept only samples found in both:

- the raw count matrix
- the cleaned metadata

This prevents expression columns and metadata rows from becoming mismatched.

The scripts always check that sample IDs line up before modeling.

## 5. Gene filtering and annotation

The workflow matched genes to the annotation file and kept genes with usable
gene-length information.

Very short genes were removed.

This step helps avoid unstable normalization for genes where gene length is not
useful or reliable.

## 6. Normalization

The project uses one normalization path:

```text
raw counts
-> gene-length adjustment
-> TMM library-size adjustment
-> GeTMM
-> log2(GeTMM + 1)
```

There are no parallel CPM, TPM, or alternative normalization branches.

This is important because the project stays focused on one consistent expression
scale.

## 7. What GeTMM-log2 means

GeTMM combines two ideas:

1. adjust for gene length
2. adjust for differences in sequencing depth between samples

After GeTMM is calculated, the values are transformed with:

```r
log2(GeTMM + 1)
```

The `+ 1` prevents problems with zeros.

The log transform makes the data easier to model because very large expression
values become less extreme.

## 8. Output files from normalization

The Step 1 normalization writes outputs here:

```text
D:/genes/genedata/Step1_GSE102556_GeTMM_brain_area_groups
```

Important folders:

```text
01_metadata
02_normalized_expression
06_repetition_reports
07_qc
```

The most important files for the next analysis step are:

```text
normalized_GeTMM_log2_<group>.csv
metadata_<group>.csv
```

The expression files contain genes in rows and samples in columns.

The metadata files contain one row per sample.

## 9. Repeated-donor checks

The workflow checked whether the same `BrainID` appeared more than once.

This matters because repeated samples from the same donor are not fully
independent.

The important conclusion was:

- Frontal_Cortex can have repeated donors because it combines dlPFC, OFC, and Cg25.
- Insula, Ventral_Striatum, and Hippocampal_Formation do not have the same issue.

This is why Frontal Cortex uses a mixed model later.

## 10. PCA and QC

PCA/QC was already completed before the current differential expression step.

PCA was used to inspect sample structure.

The key lesson from PCA/QC was:

- `OriginalRegion` strongly matters inside Frontal_Cortex.
- Therefore, `OriginalRegion` should be included in the Frontal Cortex model.

PCA is not repeated in the current Step 3 regression scripts.

## 11. Current step: transcriptomics differential expression

The current step tests whether each gene is associated with depression diagnosis.

The first simple LM/LMM scripts were useful for checking the idea, but they are
now superseded by transcriptomics-specific Step 3 scripts.

The current recommended approach uses:

- limma-voom for the small groups
- dream / variancePartition for repeated-donor analyses
- p-value histograms to inspect whether there is left-skewed signal

## 12. Why use filtered raw counts for Step 3?

The Step 1 GeTMM-log2 files are useful normalized expression matrices.

However, limma-voom and dream are designed to start from raw count data so they
can estimate count-based mean-variance weights.

Therefore the new Step 3 scripts read:

```text
filtered_raw_counts_<group>.csv
metadata_<group>.csv
```

This does not redo the earlier GeTMM-log2 preprocessing. It uses the already
filtered raw counts in a transcriptomics-appropriate differential expression
model.

## 13. Small group models

Insula, Ventral Striatum, and Hippocampal Formation use limma-voom:

```r
~ Diagnosis + Age + Sex + RIN
```

The script is:

```text
code/Step3_small_groups_limma_voom_simple.R
```

limma-voom:

- starts from filtered raw counts
- estimates voom precision weights
- fits the linear model
- uses empirical Bayes moderation
- saves one result CSV per group
- saves one diagnosis p-value histogram per group

## 14. Frontal Cortex dream model

Frontal Cortex uses dream because it can contain repeated donors:

```r
~ Diagnosis + Age + Sex + RIN + OriginalRegion + (1 | BrainID)
```

This means:

- `Diagnosis` tests the MDD versus control expression difference.
- `Age`, `Sex`, and `RIN` adjust for sample differences.
- `OriginalRegion` adjusts for dlPFC/OFC/Cg25 differences.
- `(1 | BrainID)` accounts for repeated donors.

The script is:

```text
code/Step3_frontal_dream_simple.R
```

## 15. Combined all-region model

The combined model tests whether a shared depression effect across regions is
easier to detect than separate regional analyses.

It uses dream because donors can repeat across regions:

```r
~ Diagnosis + BrainArea + Age + Sex + RIN + (1 | BrainID)
```

The script is:

```text
code/Step3_all_regions_combined_simple.R
```

`BrainArea` is included so the diagnosis effect is tested while adjusting for
broad region differences.

## 16. How many genes are tested?

The number of tested genes can differ by model because each filtered count
matrix has its own retained genes.

| Group | Samples used | Genes tested |
|---|---:|---:|
| Frontal_Cortex | 122 | 24,287 |
| Insula | 48 | 23,375 |
| Ventral_Striatum | 48 | 24,275 |
| Hippocampal_Formation | 43 | 24,700 |

For the combined all-region model, the script keeps only genes present in all
four filtered count matrices.

## 17. What FDR means

FDR means false discovery rate.

Because thousands of genes are tested, some raw p-values can look small by
chance.

The scripts correct the p-values using Benjamini-Hochberg FDR:

```r
p.adjust(p_values, method = "BH")
```

Simple interpretation:

- p-value: evidence for one gene by itself
- FDR: corrected evidence after testing thousands of genes
- FDR < 0.05: strict significance threshold
- FDR < 0.10: more exploratory threshold

## 18. Why make p-value histograms?

Each new Step 3 script saves a simple diagnosis p-value histogram.

The histogram helps answer:

- are p-values almost uniform?
- are many p-values close to zero?
- is there left-skewed signal that FDR correction may be reducing?

## 19. Where the new result files are saved

The regression result CSV files are saved locally in:

```text
outputs/step3_limma_voom_dream
```

Expected output CSV files:

```text
limma_voom_Insula_results.csv
limma_voom_Ventral_Striatum_results.csv
limma_voom_Hippocampal_Formation_results.csv
dream_Frontal_Cortex_results.csv
dream_all_regions_combined_results.csv
```

Generated output files are intentionally not committed to GitHub.
