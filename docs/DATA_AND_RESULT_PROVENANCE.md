> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Data and Result Provenance

## Canonical result

The authoritative ranking is `outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv`.

- Population: 671 Tier 2 gene-region pairs
- Classes: 30 high, 119 moderate, 522 exploratory/low
- Ranking SHA-256: `6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1`
- Scorer: `code/extract_stable_ML_prioritized_Tier2_genes.py`
- Wrapper: `code/run_corrected_ml_priority_scoring.py`
- Reconstruction: `code/analyze_corrected_ml_priority_results.py`
- Tests: `tests/test_ml_priority_scoring.py`, `tests/test_canonical_ml_priority_result.py`

## Raw event provenance

The large event files remain in the existing Elisa-style tree and are not duplicated.

| Evidence | File | Git blob | SHA-256 |
|---|---|---|---|
| v1 fold events | `outputs/length_gt_100_elisa_style/ML_pilot_Tier2/feature_importance/fold_level_feature_importance.csv` | `8e753431b47173c37736d84c6c1fb4160f2251fc` | `406ba687a9f3f4244de0d296550e1b2fcfb043f504c175c56328b447133cde1e` |
| v2 fold/model selections | `outputs/length_gt_100_elisa_style/ML_strict_v2/feature_stability/selected_genes_per_fold_all_models.csv` | `cd320246c2c17006992736c7e3710155cadb8888` | `7b33e7d0a6eb2ff2efaffdc223aced7f4722f18236a2902cc9102257b2215ef8` |

Six Tier 2 tables in `outputs/length_gt_100_elisa_style/DE_candidate_gene_sets/Tier2_working_exploratory/` define the population.

v1 uses 125 events per region and recorded repeat seeds 2026–2030. Strict v2 uses 150 non-baseline events per region and recorded outer seed 31001. Corrected scoring is deterministic post-processing and fits no models.

## Legacy provenance

`10_ML_prioritized_Tier2_extended/` is pre-correction: 35 high, 130 moderate, 506 low. It exists only for auditability and must not support current figures, tables, manuscripts, or biological conclusions.

The correction produced 5 high→moderate, 16 moderate→low, and 650 unchanged pairs. Majority-baseline removal caused the largest category-level effect; corrected weights caused one independent category change; RF sign correction and the missing-sign gate caused no independent category changes here.

## Git provenance

Audit starting point: branch `codex/clean-work-path-review-images`, HEAD `f4f60c2d12cae1e9ee600fc3ef0e6c74bfb450ef`. This cleanup is not committed. Record the approved publication commit in `canonical_manifest.json` later.
