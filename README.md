# GeneDepresion

RNA-seq preprocessing and normalization workflow for GSE102556.

The current Step 1 script prepares region-specific expression matrices using one normalization path:

```text
filtered raw counts -> GeTMM -> log2(GeTMM + 1)
```

## Main Scripts

- `code/Step1_GSE102556_GeTMM_brain_area_groups.R` - recommended Step 1; groups exact GEO tissue labels into broader editable brain areas before GeTMM-log2 normalization.
- `code/Step1.R` - earlier Step 1 version that normalizes each exact GEO tissue label separately.

## Brain-Area Mapping

The broad grouping is controlled by:

- `GSE102556_brain_area_mapping_template.csv`

Default grouping:

- OFC, dlPFC, Cg25 -> `Frontal_Cortex`
- aINS -> `Insula`
- NAc -> `Ventral_Striatum`
- Subiculum -> `Hippocampal_Formation`

Edit the mapping CSV and rerun the brain-area script if a different anatomical grouping is needed.

## Documentation

- `docs/Step1_GSE102556_GeTMM_log2_Notion.md`

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

The earlier exact-region Step 1 writes outputs to:

```text
D:/genes/genedata/Step1_GeTMM_log2_outputs
```

Large input data and generated CSV outputs are intentionally ignored by Git.
