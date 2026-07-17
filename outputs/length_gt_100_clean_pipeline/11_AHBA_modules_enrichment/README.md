# AHBA Module Analysis Outputs

This folder is reserved for the AHBA module-analysis stage added after the
extended ML-prioritized Tier 2 extraction.

Run from the project root:

```bash
python code/ahba_modules/00_fetch_ahba_data.py
python code/ahba_modules/01_preprocess_ahba_microarray.py
Rscript code/ahba_modules/02_build_ahba_modules.R
python code/ahba_modules/03_test_ahba_module_enrichment.py
```

Expected output subfolders:

- `processed/`: cleaned AHBA gene-by-sample matrix and probe/gene mapping report.
- `modules/`: AHBA gene-to-module assignments, module sizes, eigengenes, and parameters.
- `enrichment/`: module-overlap/enrichment tables using the correct MDD/AHBA background.
- `figures/`: enrichment dot plots and region-by-module heatmaps.

Interpretation: AHBA module analysis provides healthy-brain coexpression context
for ML-prioritized MDD candidate genes. It does not validate them as disease
biomarkers.
