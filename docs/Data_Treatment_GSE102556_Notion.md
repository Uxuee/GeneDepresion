# Data treatment for GSE102556

## Short summary

This project prepares GSE102556 human postmortem brain RNA-seq data for downstream depression biomarker analysis.

The current recommended Step 1 workflow is:

```text
raw counts
-> metadata matching
-> donor and brain-area harmonization
-> gene annotation and length filtering
-> broad brain-area grouping
-> region/group-specific low-expression filtering
-> GeTMM normalization
-> log2(GeTMM + 1)
```

The workflow does not run differential expression, PCA, WGCNA, or machine learning. It only creates clean expression matrices, metadata tables, normalization factors, QC summaries, and repetition reports.

## Main script

Recommended script:

```text
code/Step1_GSE102556_GeTMM_brain_area_groups.R
```

Earlier exact-region script:

```text
code/Step1.R
```

The recommended script groups exact GEO tissue labels into broader brain areas before normalization.

## Input files

The script expects the main data files under:

```text
D:/genes/genedata
```

Required input files:

```text
GSE102556_raw_counts_GRCh38.p13_NCBI.tsv
Human.GRCh38.p13.annot.tsv
```

Optional metadata file:

```text
GSE102556_metadata_clean.csv
```

If the optional metadata file is absent, metadata is fetched from GEO using:

```r
GEOquery::getGEO("GSE102556", GSEMatrix = TRUE)
```

During the run, GEOquery found:

```text
GSE102556-GPL11154_series_matrix.txt.gz
GSE102556-GPL13112_series_matrix.txt.gz
```

These are GEO series-matrix metadata files for platform-specific parts of the GSE102556 study.

## Output folder

The broad brain-area workflow writes outputs to:

```text
D:/genes/genedata/Step1_GSE102556_GeTMM_brain_area_groups
```

Generated outputs are not committed to Git because they are large and reproducible.

## Raw count treatment

The raw count matrix is read from:

```text
GSE102556_raw_counts_GRCh38.p13_NCBI.tsv
```

The first column is treated as `GeneID`.

All other columns are treated as sample-level raw count columns.

The script checks for duplicated `GeneID` rows.

If duplicated `GeneID` entries exist, they are collapsed by summing raw counts across duplicate rows.

The duplicate report is exported as:

```text
duplicate_gene_ids_in_raw_counts.csv
```

This is important because all downstream matrices must have one row per gene.

## Gene annotation treatment

The gene annotation table is read from:

```text
Human.GRCh38.p13.annot.tsv
```

The script requires:

```text
GeneID
Length
```

Genes are matched between the raw count matrix and the annotation table by `GeneID`.

Gene length is used for GeTMM normalization.

Genes with missing, invalid, or non-positive length are removed.

Genes shorter than or equal to 20 bp are also removed:

```r
Length > 20
```

In the successful broad brain-area run, the raw count matrix after annotation matching contained:

```text
39,376 genes x 261 samples
```

## Metadata treatment

Metadata is loaded in this order:

1. Use `GSE102556_metadata_clean.csv` if present.
2. Use `metadata_clean_matched_all_samples.csv` if present.
3. Use a previous Step 1 matched metadata file if present.
4. Otherwise fetch metadata from GEO with `GEOquery::getGEO("GSE102556", GSEMatrix = TRUE)`.

For the run documented here, no local clean metadata file was used, so GEOquery fetched the metadata.

The script standardizes useful metadata columns, including:

```text
SampleID
OriginalRegion
RegionID
BrainID
RIN
Diagnosis
Sex
Age
PMI
pH
```

## OriginalRegion

`OriginalRegion` stores the exact tissue label from GEO.

Examples:

```text
Orbitofrontal (OFC; BA11)
Dorsolateral prefrontal cortex (dlPFC; BA8/9)
Cingulate gyrus 25 (Cg25)
Anterior Insula (aINS)
Nucleus Accumbens (Nac)
Subiculum (Sub)
```

This field is preserved so the exact anatomical source is never lost.

## BrainID

GSE102556 metadata does not expose a simple explicit donor/person column under a stable name.

The script therefore derives `BrainID` from the donor/person prefix in the GEO `title` field.

Example:

```text
14: Orbitofrontal (OFC; BA11)
```

becomes:

```text
BrainID = 14
OriginalRegion = Orbitofrontal (OFC; BA11)
```

This lets the workflow detect repeated donors across brain regions.

## Broad brain-area grouping

The exact GEO tissue labels are mapped to broader analysis groups using:

```text
GSE102556_brain_area_mapping_template.csv
```

This file is intentionally editable.

Default mapping:

| OriginalRegion | BrainArea | MajorBrainPart |
|---|---|---|
| Orbitofrontal (OFC; BA11) | Frontal_Cortex | Cerebral_Cortex |
| Dorsolateral prefrontal cortex (dlPFC; BA8/9) | Frontal_Cortex | Cerebral_Cortex |
| Cingulate gyrus 25 (Cg25) | Frontal_Cortex | Cerebral_Cortex |
| Anterior Insula (aINS) | Insula | Cerebral_Cortex |
| Nucleus Accumbens (Nac) | Ventral_Striatum | Subcortical_Striatum |
| Subiculum (Sub) | Hippocampal_Formation | Temporal_Limbic_Hippocampal |

The current default grouping column is:

```r
analysis_group_column <- "BrainArea"
```

The script creates:

```text
OriginalRegion
BrainArea
MajorBrainPart
AnalysisGroup
```

where `AnalysisGroup` is copied from the selected grouping column.

## Why group into BrainArea?

The exact GEO labels are anatomically specific. Treating each exact label as a fully independent analysis can make sample sizes small.

The broad grouping allows related anatomical regions to be normalized and modeled together when scientifically appropriate.

For example:

```text
OFC + dlPFC + Cg25 -> Frontal_Cortex
```

The exact region is still preserved as `OriginalRegion`, so later models can still include or inspect the fine anatomical label.

## Sample-handling modes

The script produces two versions of outputs.

## Mode 1: keep_all_samples

This mode keeps all matched samples.

It preserves repeated donors across brain areas and within broad brain areas.

This is useful if the downstream model will explicitly handle repeated donor measurements, for example with:

```r
(1 | BrainID)
```

## Mode 2: keep_highest_RIN_within_analysis_group

This mode keeps only one sample per `BrainID` inside each broad `AnalysisGroup`.

If the same donor appears multiple times inside the same broad brain area, the sample with the highest RIN is kept.

This reduces repeated-measure structure inside a broad brain-area analysis.

It does not remove the same donor across unrelated analysis groups.

## Low-expression filtering

Low-expression filtering is done separately inside each broad `AnalysisGroup`.

The rule is:

```r
quantile(raw_counts, 0.80) > 10
```

In words:

```text
keep genes whose 80th percentile raw count is greater than 10 inside that analysis group
```

This follows the HTSeq-style rule:

```text
exclude genes with fewer than 10 counts in at least 80% of samples
```

Filtering is applied before normalization.

The workflow does not apply CPM before GeTMM.

## GeTMM normalization

GeTMM is computed directly from filtered raw counts.

## Step 1: gene length correction

Gene length is converted from base pairs to kilobases:

```r
length_kb = Length / 1000
```

Reads per kilobase are computed:

```r
RPK = counts / length_kb
```

## Step 2: RNA composition correction

TMM normalization factors are calculated on the RPK matrix:

```r
edgeR::calcNormFactors(method = "TMM")
```

This adjusts for RNA composition differences between samples.

## Step 3: effective library size

For each sample:

```r
effective_library_size = colSums(RPK) * TMM_factor
```

## Step 4: GeTMM value

For each gene and sample:

```r
GeTMM = 1e6 * RPK / effective_library_size
```

## Step 5: final log transform

The final normalized matrix is:

```r
log2(GeTMM + 1)
```

This is the main normalized expression output.

## Main exported files

At the top level of the output folder, the script writes:

```text
brain_area_mapping_used.csv
metadata_clean_matched_all_samples_with_brain_area.csv
analysis_group_composition_report.csv
model_decision_repetition_summary.csv
Step1_QC_summary_by_analysis_group_all_modes.csv
Step1_metadata_QC_and_repetition_reports_by_brain_area.rds
```

It also writes donor and gene repetition reports before group-level deduplication:

```text
all_matched_samples_before_group_dedup_donor_repetition_overall.csv
all_matched_samples_before_group_dedup_donor_repetition_by_analysis_group.csv
all_matched_samples_before_group_dedup_donor_repetition_by_original_region.csv
all_matched_samples_before_group_dedup_gene_entry_repetition_summary.csv
all_matched_samples_before_group_dedup_gene_entry_repetition_summary_compact.csv
```

For each sample-handling mode and analysis group, the script writes:

```text
normalized_GeTMM_log2_<analysis_group>.csv
filtered_raw_counts_<analysis_group>.csv
metadata_<analysis_group>.csv
gene_annotation_<analysis_group>.csv
GeTMM_normalization_factors_<analysis_group>.csv
```

It also writes repetition reports for each group.

## Brain-area composition in the successful run

| AnalysisGroup | Samples | Unique BrainID | Original regions |
|---|---:|---:|---|
| Frontal_Cortex | 122 | 48 | dlPFC, OFC, Cg25 |
| Hippocampal_Formation | 43 | 43 | Subiculum |
| Insula | 48 | 48 | Anterior Insula |
| Ventral_Striatum | 48 | 48 | Nucleus Accumbens |

Detailed composition:

| AnalysisGroup | Original region counts |
|---|---|
| Frontal_Cortex | dlPFC = 48, OFC = 47, Cg25 = 27 |
| Hippocampal_Formation | Subiculum = 43 |
| Insula | Anterior Insula = 48 |
| Ventral_Striatum | Nucleus Accumbens = 48 |

## QC summary from the successful run

## keep_all_samples

| AnalysisGroup | Samples | Unique BrainID | Genes retained | CTRL | MDD |
|---|---:|---:|---:|---:|---:|
| Frontal_Cortex | 122 | 48 | 24,287 | 59 | 63 |
| Hippocampal_Formation | 43 | 43 | 24,700 | 19 | 24 |
| Insula | 48 | 48 | 23,375 | 22 | 26 |
| Ventral_Striatum | 48 | 48 | 24,275 | 22 | 26 |

## keep_highest_RIN_within_analysis_group

| AnalysisGroup | Samples before | Samples after | Dropped | Genes retained | CTRL | MDD |
|---|---:|---:|---:|---:|---:|---:|
| Frontal_Cortex | 122 | 48 | 74 | 24,147 | 22 | 26 |
| Hippocampal_Formation | 43 | 43 | 0 | 24,700 | 19 | 24 |
| Insula | 48 | 48 | 0 | 23,375 | 22 | 26 |
| Ventral_Striatum | 48 | 48 | 0 | 24,275 | 22 | 26 |

The largest change is in `Frontal_Cortex`, because this group combines multiple frontal regions from the same donors.

## Repeated-measure interpretation

The model-decision report checks donor repetition at three levels.

| Scope | Repetition result | Suggested model |
|---|---|---|
| Original-region-wise | No repeated BrainID within the same exact region | Linear model is probably OK inside each original fine region |
| Brain-area-wise | Repeated BrainID appears within broad BrainArea groups | For keep-all samples, consider LMM/dream with `(1 | BrainID)`, or use highest-RIN mode |
| All groups combined | Donors repeat across all samples/groups | LMM/dream with `(1 | BrainID)` is advisable |

The key issue is `Frontal_Cortex`.

Because OFC, dlPFC, and Cg25 are grouped together, many donors have multiple frontal-cortex samples.

Therefore:

- If using `keep_all_samples`, the downstream model should handle repeated donor measurements.
- If using `keep_highest_RIN_within_analysis_group`, frontal cortex is reduced to one sample per donor.

## Recommended use for downstream analysis

For a first broad brain-area differential expression analysis:

1. Use the output folder for the desired mode.
2. Use `normalized_GeTMM_log2_<analysis_group>.csv` as expression input.
3. Use `metadata_<analysis_group>.csv` as metadata input.
4. For `keep_all_samples`, use mixed models when repeated `BrainID` is present.
5. For `keep_highest_RIN_within_analysis_group`, a linear model may be acceptable from the donor-repetition perspective.

For all-group or multi-region modeling, include a donor random intercept:

```r
(1 | BrainID)
```

## What was deliberately not done

Step 1 does not perform:

- differential expression
- PCA
- WGCNA
- machine learning
- biomarker selection
- pathway enrichment

Those analyses belong in later steps after this preprocessing and normalization layer is stable.

## How to change the anatomical grouping

Edit:

```text
GSE102556_brain_area_mapping_template.csv
```

Then rerun:

```r
code/Step1_GSE102556_GeTMM_brain_area_groups.R
```

The script will use the edited mapping and regenerate all grouped outputs.

## Final interpretation

The data treatment creates reproducible, broad brain-area expression matrices from GSE102556.

It preserves exact GEO tissue labels, adds editable anatomical groupings, detects repeated donors, and normalizes expression using one coherent GeTMM-log2 path.

The most important downstream decision is whether to use all repeated frontal-cortex samples with a mixed model, or to use the highest-RIN version to keep one frontal-cortex sample per donor.
