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

## 11. Current step: differential expression

The current step tests whether each gene is associated with depression diagnosis.

This is done one gene at a time.

For each gene:

1. take that gene's normalized expression values across samples
2. combine them with the aligned metadata
3. fit a regression model
4. extract the diagnosis effect
5. save the coefficient and p-value

## 12. Frontal Cortex model

Frontal Cortex uses a linear mixed model:

```r
expression_value ~ Diagnosis + Age + Sex + RIN + OriginalRegion + (1 | BrainID)
```

This means:

- `Diagnosis` tests the MDD versus control expression difference.
- `Age`, `Sex`, and `RIN` adjust for sample differences.
- `OriginalRegion` adjusts for dlPFC/OFC/Cg25 differences.
- `(1 | BrainID)` accounts for repeated donors.

The script is:

```text
code/Step3_LMM_Frontal_Cortex_simple.R
```

## 13. Other group models

Insula, Ventral Striatum, and Hippocampal Formation use simple linear models:

```r
expression_value ~ Diagnosis + Age + Sex + RIN
```

The script is:

```text
code/Step3_LM_Small_Groups_simple.R
```

These groups do not include `(1 | BrainID)` because they do not have the same
within-group repeated-donor issue.

## 14. How many genes were tested?

The number of tested genes differs by group because each group has its own
filtered normalized expression matrix.

| Group | Samples used | Genes tested |
|---|---:|---:|
| Frontal_Cortex | 122 | 24,287 |
| Insula | 48 | 23,375 |
| Ventral_Striatum | 48 | 24,275 |
| Hippocampal_Formation | 43 | 24,700 |

## 15. What FDR means

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

## 16. Current result summary

All models completed for all tested genes.

| Group | Genes tested | FDR < 0.05 | FDR < 0.10 | Best FDR | Top GeneID |
|---|---:|---:|---:|---:|---:|
| Frontal_Cortex | 24,287 | 0 | 1 | 0.06553346 | 5147 |
| Insula | 23,375 | 0 | 0 | 0.9344865 | 10090 |
| Ventral_Striatum | 24,275 | 0 | 0 | 0.1279001 | 80705 |
| Hippocampal_Formation | 24,700 | 0 | 0 | 0.9983605 | 65083 |

Using FDR < 0.05, no genes are significant.

Using FDR < 0.10, one Frontal Cortex gene passes the exploratory threshold.

## 17. Important caution

The Frontal Cortex mixed model completed, but it printed many singular-fit
messages and two convergence warnings.

This means the output table exists, but the Frontal Cortex mixed-model result
should be interpreted carefully.

It does not mean the script failed. It means the random donor effect may be weak
or hard to estimate for many genes.

## 18. Where the current result files are saved

The regression result CSV files are saved locally in:

```text
outputs/differential_expression
```

The output files are:

```text
LMM_Frontal_Cortex_results.csv
LM_Insula_results.csv
LM_Ventral_Striatum_results.csv
LM_Hippocampal_Formation_results.csv
```

Generated output files are intentionally not committed to GitHub.
