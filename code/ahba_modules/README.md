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

Outputs are written to:

```text
outputs/length_gt_100_clean_pipeline/ahba_modules/
```

Interpretation: AHBA module analysis provides healthy-brain coexpression context
for ML-prioritized MDD candidate genes. It does not validate them as disease
biomarkers.
