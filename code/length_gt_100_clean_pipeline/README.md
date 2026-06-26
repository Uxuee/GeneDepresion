# Length > 100 Clean Tier2-Only Pipeline

This folder provides a reproducible code path from completed Length > 100 raw/normalized gene-expression outputs through Tier 2 candidates, enrichment, ML v1 using Tier2_fixed features only, and final Tier2FixedOnly ML-prioritized genes.

The pipeline reuses completed outputs from `outputs/length_gt_100_elisa_style` and writes clean organized outputs to `outputs/length_gt_100_clean_pipeline`.

## Run Commands

Run the full clean pipeline:

```bash
python code/length_gt_100_clean_pipeline/run_length_gt_100_clean_pipeline.py --skip-existing
```

Run only the Tier2-only ML prioritization branch:

```bash
python code/length_gt_100_clean_pipeline/run_length_gt_100_clean_pipeline.py --from-step 9 --to-step 10 --skip-existing
```

Smoke-test without writing:

```bash
python code/length_gt_100_clean_pipeline/run_length_gt_100_clean_pipeline.py --dry-run
```

## Numbered Scripts

- `00_pipeline_config.py`: central paths and configuration.
- `01_prepare_length_gt_100_inputs.py`: copies/filter summaries and metadata; creates tested-gene universe.
- `02_run_getmm_log2_normalization.py`: reuses log2(GeTMM + 1) normalized expression matrices.
- `03_run_broad_region_LMM.py`: reuses completed broad-region LMM outputs.
- `04_compute_standardized_effect_sizes.py`: reuses standardized effect sizes and writes a sign-check summary.
- `05_create_evidence_tiers.py`: creates clean Tier 1, Tier 2, and Tier 3 outputs. Tier 4 is excluded.
- `06_run_frontal_subregion_sensitivity.py`: reuses fixed-effect frontal subregion sensitivity outputs.
- `07_package_Tier2_candidate_gene_sets.py`: packages final Tier 2 candidate pools.
- `08_run_Tier2_enrichment.py`: copies completed Tier 2 enrichment tables and selected report-ready plots.
- `09_run_ML_v1_Tier2_fixed_only.py`: filters existing ML v1 outputs to `feature_set == "Tier2_fixed"`.
- `10_extract_Tier2FixedOnly_ML_prioritized_genes.py`: creates final Tier2FixedOnly ML-prioritized genes.
- `11_collect_report_images.py`: collects report-ready images and writes an image inventory.
- `12_write_integrated_clean_report.py`: writes the integrated report and reproducibility manifest.
- `run_length_gt_100_clean_pipeline.py`: ordered runner with `--from-step`, `--to-step`, `--skip-existing`, and `--dry-run`.

## Expected Outputs

Primary outputs are written under:

`outputs/length_gt_100_clean_pipeline`

The integrated report is written to:

`docs/MDD_Length_gt_100_clean_Tier2_to_ML_prioritized_Report.md`

## Interpretation Warnings

Tier 2 genes are exploratory candidate genes, not FDR-significant discoveries.

ML was used as a Tier2-guided gene-prioritization tool. Classification performance is reported only descriptively and should not be interpreted as unbiased diagnostic prediction.

The final prioritized genes should be interpreted as internally supported candidate genes, not validated biomarkers.

Tier 4 was explored previously but is excluded from this clean Tier2-only pipeline and report story.

