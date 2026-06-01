# GeneDepresion

RNA-seq preprocessing and normalization workflow for GSE102556.

The current Step 1 script prepares region-specific expression matrices using one normalization path:

```text
filtered raw counts -> GeTMM -> log2(GeTMM + 1)
```

## Main Script

- `code/Step1.R`

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

Generated outputs are written to:

```text
D:/genes/genedata/Step1_GeTMM_log2_outputs
```

Large input data and generated CSV outputs are intentionally ignored by Git.
