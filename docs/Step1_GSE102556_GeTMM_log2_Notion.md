# Step 1: GSE102556 GeTMM-log2 normalization

## Purpose

This Step 1 script prepares GSE102556 RNA-seq count data for downstream differential expression or region-wise modeling.

It does not run differential expression, PCA, WGCNA, or machine learning. Its only job is to produce clean, region-specific expression matrices and metadata reports.

## Inputs

The script reads data from `D:/genes/genedata`.

Required files:

- `GSE102556_raw_counts_GRCh38.p13_NCBI.tsv`
- `Human.GRCh38.p13.annot.tsv`

Optional file:

- `GSE102556_metadata_clean.csv`

If the optional metadata file is not present, the script fetches metadata with:

```r
GEOquery::getGEO("GSE102556", GSEMatrix = TRUE)
```

In the successful run, the local metadata file was not present, so GEOquery downloaded the GEO Series metadata and parsed it.

## Main Design Decision

The script uses one coherent normalization path:

```text
raw counts
-> gene/sample filtering
-> region split
-> GeTMM
-> log2(GeTMM + 1)
```

There are no parallel CPM, FPKM, TPM, or alternative normalization branches.

## Step-by-step workflow

## 1. Read raw counts

The raw count matrix is loaded from:

```text
GSE102556_raw_counts_GRCh38.p13_NCBI.tsv
```

The file is expected to have:

- one `GeneID` column
- one column per GSM sample

Rows with missing or blank `GeneID` are removed.

## 2. Report and collapse duplicated GeneIDs

The script checks whether the raw count table contains multiple rows for the same `GeneID`.

It writes:

```text
duplicate_gene_ids_in_raw_counts.csv
```

Then duplicated `GeneID` rows are collapsed by summing raw counts across duplicate rows.

In the current run, there were no duplicated GeneID rows.

## 3. Read gene annotation

The script reads:

```text
Human.GRCh38.p13.annot.tsv
```

It uses the annotation columns:

- `GeneID`
- `Length`

Genes are retained only if they match the annotation and have valid positive length.

## 4. Remove very short genes

The script keeps genes with:

```r
Length > 20
```

This follows the short-gene filter used in the previous Step1.R logic.

## 5. Fetch or load metadata

If `GSE102556_metadata_clean.csv` exists, the script reads it.

If it does not exist, the script fetches GEO metadata using GEOquery.

For GSE102556, GEOquery finds two platform-specific series matrix files:

```text
GSE102556-GPL11154_series_matrix.txt.gz
GSE102556-GPL13112_series_matrix.txt.gz
```

The script parses the GEO metadata fields such as:

- sample accession
- tissue or brain region
- phenotype
- RIN
- age
- PMI
- pH

## 6. Create standard metadata columns

The script standardizes metadata into these columns:

- `SampleID`
- `RegionID`
- `BrainID`
- `RIN`

For this dataset, no explicit donor column was found. The script therefore uses the numeric prefix in the GEO title field as `BrainID`.

Example:

```text
14: Orbitofrontal (OFC; BA11)
```

becomes:

```text
BrainID = 14
RegionID = Orbitofrontal (OFC; BA11)
```

## 7. Match metadata to count samples

Only samples present in both the raw count matrix and metadata are kept.

The script writes:

```text
metadata_clean_matched_all_samples.csv
```

In the successful run:

- 263 GEO metadata rows were loaded
- 261 samples matched the count matrix and had a brain region

## 8. Write donor and gene repetition reports before deduplication

Before any within-region duplicate handling, the script writes donor repetition reports:

```text
all_matched_samples_before_region_dedup_donor_repetition_overall.csv
all_matched_samples_before_region_dedup_donor_repetition_by_region.csv
```

It also writes a gene-entry repetition report:

```text
all_matched_samples_before_region_dedup_gene_entry_repetition_summary.csv
```

This report records, for each gene:

- how many expression entries exist
- how many unique BrainID/person IDs are present
- the maximum number of entries from the same person

## 9. Write model-design recommendation

The script writes:

```text
model_decision_repetition_summary.csv
```

For the current run, the recommendation was:

```text
Repeated BrainID appears across all regions; an all-region model should use LMM/dream with (1 | BrainID).
Region-wise linear models are probably OK if each region has no repeated BrainID.
```

This is because donors repeat across brain regions, but not within a single region.

## 10. Run two sample-handling modes

The script produces two versions of every region-level output.

### Mode A: keep_all_samples

Keeps every matched sample.

This preserves repeated BrainID entries across different brain regions.

### Mode B: keep_highest_RIN_within_region

Within each brain region, if the same `BrainID` appears more than once, the script keeps the sample with the highest RIN.

The same person is not removed across different brain regions.

In the current GSE102556 run, no within-region repeated BrainID values were detected, so the sample counts were the same in both modes.

## 11. Split by brain region before normalization

Samples are split by `RegionID` before low-expression filtering and normalization.

The regions found in the current run were:

- Anterior Insula (aINS)
- Cingulate gyrus 25 (Cg25)
- Dorsolateral prefrontal cortex (dlPFC; BA8/9)
- Nucleus Accumbens (Nac)
- Orbitofrontal (OFC; BA11)
- Subiculum (Sub)

## 12. Apply region-wise low-expression filtering

Within each region, genes are kept if:

```r
quantile(raw_counts, 0.80) > 10
```

This means a gene must have enough raw-count evidence in that region before normalization.

## 13. Run GeTMM normalization

For each region and mode, the script computes GeTMM from filtered raw counts.

### 13.1 Gene length correction

Gene length is converted to kilobases:

```r
length_kb = Length / 1000
```

Reads per kilobase are computed:

```r
RPK = counts / length_kb
```

### 13.2 RNA composition correction

TMM normalization factors are computed with edgeR:

```r
edgeR::calcNormFactors(method = "TMM")
```

The TMM factors are calculated on the RPK matrix.

### 13.3 Effective library size

For each sample:

```r
effective_library_size = colSums(RPK) * TMM_factor
```

### 13.4 GeTMM

For each gene and sample:

```r
GeTMM = 1e6 * RPK / effective_library_size
```

### 13.5 Final log transform

The final normalized matrix is:

```r
log2(GeTMM + 1)
```

## 14. Per-region outputs

For each mode and region, the script writes:

```text
normalized_GeTMM_log2_<region>.csv
filtered_raw_counts_<region>.csv
metadata_<region>.csv
gene_annotation_<region>.csv
GeTMM_normalization_factors_<region>.csv
region_<region>_<mode>_donor_repetition_overall.csv
region_<region>_<mode>_donor_repetition_by_region.csv
region_<region>_<mode>_gene_entry_repetition_summary.csv
```

## 15. QC outputs

The script writes one QC summary per mode:

```text
Step1_QC_summary_by_region_keep_all_samples.csv
Step1_QC_summary_by_region_keep_highest_RIN_within_region.csv
```

It also writes one combined QC file:

```text
Step1_QC_summary_by_region_all_modes.csv
```

## Current run summary

The script completed successfully.

Output folder:

```text
D:/genes/genedata/Step1_GeTMM_log2_outputs
```

The run produced 105 files.

Matched samples:

```text
261
```

Genes after annotation matching and `Length > 20`:

```text
39,356
```

Region-level retained genes after expression filtering:

| Region | Samples | Genes retained |
|---|---:|---:|
| Anterior Insula (aINS) | 48 | 23,375 |
| Cingulate gyrus 25 (Cg25) | 27 | 24,663 |
| Dorsolateral prefrontal cortex (dlPFC; BA8/9) | 48 | 24,195 |
| Nucleus Accumbens (Nac) | 48 | 24,275 |
| Orbitofrontal (OFC; BA11) | 47 | 23,938 |
| Subiculum (Sub) | 43 | 24,700 |

## Recommended next step

Use the region-specific `normalized_GeTMM_log2_<region>.csv` files for downstream exploratory analysis or differential expression.

For region-wise differential expression, a linear model may be acceptable from the donor-repetition perspective.

For an all-region model, use a mixed model with a donor random intercept:

```r
(1 | BrainID)
```
