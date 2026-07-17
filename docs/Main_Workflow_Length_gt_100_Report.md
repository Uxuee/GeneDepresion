# Main Workflow Report: Length > 100 Clean MDD Transcriptomics Pipeline

Generated for the current main workflow in:

`D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline`

## Start Here

This folder is the clean, report-ready workflow for the MDD/control transcriptomics analysis. It starts from Length > 100 filtered gene-expression inputs, keeps one coherent GeTMM-log2 normalization path, runs statistical models, defines exploratory Tier 2 candidate genes, performs enrichment, prioritizes Tier 2 genes with ML support, and now contains the first MDD WGCNA validation branch.

The main interpretation rule is important: Tier 2 and ML-prioritized genes are exploratory candidate signatures. They are not externally validated biomarkers.

## Main Locations

| Purpose | Location |
| --- | --- |
| Main workflow output root | `outputs/length_gt_100_clean_pipeline` |
| Pipeline manifest | `outputs/length_gt_100_clean_pipeline/pipeline_manifest.json` |
| Clean technical report | `outputs/length_gt_100_clean_pipeline/MDD_Length_gt_100_clean_Tier2_to_ML_prioritized_Report.md` |
| GitHub docs copy of clean report | `docs/MDD_Length_gt_100_clean_Tier2_to_ML_prioritized_Report.md` |
| Step-by-step narrative report | `docs/Step_by_step_MDD_Length_gt_100_clean_report.md` |
| Main workflow notebook | `outputs/length_gt_100_clean_pipeline/notebooks/01_clear_length_gt_100_workflow_notebook.md` |
| Main workflow scripts | `code/length_gt_100_clean_pipeline` |
| Report-ready figures | `outputs/length_gt_100_clean_pipeline/images` |

## Inputs

The clean pipeline reuses completed Length > 100 analysis inputs from:

`D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style`

The key reused inputs recorded in the manifest are:

- `LMM_all_genes_length_gt_100.csv`
- `standardized_effect_size_4groups\LMM_all_genes_with_standardized_effect_size_4groups.csv`
- `frontal_subregion_sensitivity\frontal_subregion_LMM_all_genes_with_standardized_effect_size.csv`
- `ML_pilot_Tier2\feature_importance\feature_importance_summary.csv`
- `enrichment_Tier2`

The normalized expression matrices used by downstream work live here:

`outputs/length_gt_100_clean_pipeline/02_normalized_expression`

They are GeTMM normalized and transformed as:

```text
log2(GeTMM + 1)
```

## Methods in One Line

```text
Length > 100 filtered expression
-> GeTMM-log2 normalization
-> broad-region LMM
-> standardized diagnosis effect size
-> evidence tiers
-> frontal subregion sensitivity
-> final Tier 2 candidate sets
-> Tier 2 enrichment
-> extended ML v1 prioritization
-> ML-prioritized Tier 2 genes
-> AHBA and MDD WGCNA module follow-up
```

## Script Order

The clean workflow scripts are in:

`code/length_gt_100_clean_pipeline`

Run order:

1. `01_prepare_length_gt_100_inputs.py`
2. `02_run_getmm_log2_normalization.py`
3. `03_run_broad_region_LMM.py`
4. `04_compute_standardized_effect_sizes.py`
5. `05_create_evidence_tiers.py`
6. `06_run_frontal_subregion_sensitivity.py`
7. `07_package_Tier2_candidate_gene_sets.py`
8. `08_run_Tier2_enrichment.py`
9. `09_run_ML_v1_extended_prioritization.py`
10. `10_extract_extended_ML_prioritized_Tier2_genes.py`
11. `11_collect_report_images.py`
12. `12_write_integrated_clean_report.py`

The wrapper is:

`code/length_gt_100_clean_pipeline/run_length_gt_100_clean_pipeline.py`

## Main Result Counts

### Broad-region evidence tiers

| Brain region | Genes tested | Tier 1 | Tier 2 | Tier 3 |
| --- | ---: | ---: | ---: | ---: |
| Frontal_Cortex | 23,941 | 0 | 1 | 99 |
| Hippocampal_Formation | 24,371 | 0 | 127 | 70 |
| Insula | 23,090 | 0 | 70 | 92 |
| Ventral_Striatum | 23,917 | 4 | 146 | 167 |

### Final Tier 2 candidate pools

| Region or subregion | Model source | Tier 2 genes | Up in MDD | Down in MDD |
| --- | --- | ---: | ---: | ---: |
| Cingulate gyrus 25 | frontal_subregion_fixed_effect_sensitivity | 186 | 137 | 49 |
| Hippocampal Formation | grouped_region_LMM | 127 | 80 | 47 |
| Insula | grouped_region_LMM | 70 | 37 | 33 |
| OFC BA11 | frontal_subregion_fixed_effect_sensitivity | 70 | 38 | 32 |
| Ventral Striatum | grouped_region_LMM | 146 | 76 | 70 |
| dlPFC BA8/9 | frontal_subregion_fixed_effect_sensitivity | 72 | 24 | 48 |

### Final ML-prioritized Tier 2 output

The final extended ML v1 ranking contains:

- 671 final ranked Tier 2 genes
- 35 high-confidence ML-prioritized Tier 2 genes
- 130 moderate-confidence ML-prioritized Tier 2 genes

Main tables:

- `10_ML_prioritized_Tier2_extended/final_ranked_ML_prioritized_Tier2_genes.csv`
- `10_ML_prioritized_Tier2_extended/final_high_confidence_ML_prioritized_Tier2_genes.csv`
- `10_ML_prioritized_Tier2_extended/final_moderate_confidence_ML_prioritized_Tier2_genes.csv`
- `10_ML_prioritized_Tier2_extended/final_top_genes_by_region.csv`

## Current MDD WGCNA Branch

The MDD WGCNA module branch is in:

`outputs/length_gt_100_clean_pipeline/mdd_wgcna_modules`

Script:

`code/Step6_MDD_WGCNA_modules.R`

This branch builds modules from normalized MDD/control expression matrices, not DEG-only or Tier 2-only genes. The current saved output is a limited validation run for Frontal Cortex, clearly labeled in:

`outputs/length_gt_100_clean_pipeline/mdd_wgcna_modules/MDD_WGCNA_Module_Report.md`

The full WGCNA grid is intentionally separate because it is longer-running.

## What to Push to GitHub

For a clean GitHub update, prioritize these human-readable workflow files:

- `README.md`
- `docs/Step_by_step_MDD_Length_gt_100_clean_report.md`
- `docs/MDD_Length_gt_100_clean_Tier2_to_ML_prioritized_Report.md`
- `docs/Main_Workflow_Length_gt_100_Report.md`
- `outputs/length_gt_100_clean_pipeline/MAIN_WORKFLOW_REPORT.md`
- `outputs/length_gt_100_clean_pipeline/notebooks/01_clear_length_gt_100_workflow_notebook.md`
- `outputs/length_gt_100_clean_pipeline/pipeline_manifest.json`

Avoid pushing large raw/intermediate CSV outputs unless they are intentionally selected for GitHub.

## Limitations

- Tier 2 genes are exploratory candidates, not final FDR-significant discoveries.
- ML v1 is internal prioritization, not validated diagnostic classification.
- No external validation cohort is included in this clean workflow.
- The MDD WGCNA branch currently has a saved validation run, not the full completed parameter grid.
