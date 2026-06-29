# AHBA Module Analysis Stage

This folder adds a post-hoc AHBA gene-module analysis stage to the clean
Length > 100 MDD pipeline. It does not rerun LMM, ML, or enrichment.

Run from the project root:

```bash
python code/ahba_modules/00_fetch_ahba_data.py
python code/ahba_modules/01_preprocess_ahba_microarray.py
Rscript code/ahba_modules/02_build_ahba_modules.R
python code/ahba_modules/03_test_ahba_module_enrichment.py
```

By default, the WGCNA step is memory-aware:

- it restricts module construction to genes tested in the Length > 100 MDD pipeline and present in AHBA;
- if more than 10,000 genes remain, it keeps the top variable genes before WGCNA.

Optional overrides:

```bash
# Use all processed AHBA genes that pass missingness/variance filters.
AHBA_RESTRICT_TO_MDD_TESTED=false Rscript code/ahba_modules/02_build_ahba_modules.R

# Change the top-variable-gene cap; use 0 for no cap.
AHBA_MAX_GENES_FOR_WGCNA=8000 Rscript code/ahba_modules/02_build_ahba_modules.R
```

The enrichment step uses unique gene symbols for candidate sets. Region-specific
and direction-specific sets are preserved, but duplicated gene-region rows are
deduplicated within each tested set before hypergeometric testing.

Outputs are written to:

```text
outputs/length_gt_100_clean_pipeline/ahba_modules/
```

Interpretation: AHBA module analysis provides healthy-brain coexpression context
for ML-prioritized MDD candidate genes. It does not validate them as disease
biomarkers.


## Published AHBA Canonical Module Branch

The de novo WGCNA sensitivity analysis produced one giant module across the tested parameter grid, so biological enrichment should use the published AHBA canonical module branch instead of the one-module de novo result.

Place Hawrylycz et al. 2015 Supplementary Table 6 in:

```text
data/ahba/published_modules/
```

Then run:

```bash
python code/ahba_modules/04_fetch_or_load_published_ahba_modules.py
python code/ahba_modules/05_test_published_ahba_module_enrichment.py
```

The published-module enrichment step deduplicates gene-region pairs to unique gene symbols within each candidate set. Its background is the intersection of the clean MDD tested-gene universe and genes present in the published AHBA module assignment table.
