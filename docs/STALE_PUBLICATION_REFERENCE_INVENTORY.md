> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Stale publication reference inventory

This inventory separates current canonical references from historical material.
Phase 3A deliberately did not rewrite the pre-existing dirty reports or the
preliminary manuscript package.

| Reference class | Current status | Required handling |
|---|---|---|
| `10_ML_prioritized_Tier2_extended` and the old 35 high / 130 moderate list | superseded | Replace with corrected 30 high / 119 moderate / 522 low; do not restore as current |
| `14_Preliminary_Manuscript_Package` | preliminary and stale | Exclude from Phase 3A staging; refresh only in a later controlled publication phase |
| Step 13 MDD module-level ML2 outputs | generator missing | Exclude from publication |
| Tier4 GSEA claims | method/results missing | Exclude from publication |
| AHBA FDR005/FDR010/M26-only module-score models | historical post-selection sensitivity | Keep for provenance; only all-32 is active |
| Broad pooled Frontal Cortex WGCNA | supplementary | Do not present as a main original region/subregion result |
| WGCNA biological annotations | incomplete mapping | Use cautiously as supportive interpretation pending controlled reannotation |

The old pipeline runner is retained for historical reproducibility but now
requires `--allow-legacy` for steps 9–12. The new canonical runner never invokes
superseded or excluded stages.

## File-level candidate classification

The Phase 3A content scan searched for old counts, denominator 180, the legacy
ranking path, and legacy M26/MEdarkred claims. Candidate files were classified as
follows; a match is not automatically an error because many files are frozen
historical evidence.

### current-stale — refresh later, deliberately not edited in Phase 3A

- `docs/Step_by_step_MDD_Length_gt_100_clean_report.md`
- `docs/Main_Workflow_Length_gt_100_Report.md`
- `docs/MDD_Length_gt_100_clean_Tier2_to_ML_prioritized_Report.md`
- `docs/MDD_WGCNA_Module_Report.md`
- `docs/ahba_modules/AHBA_gene_module_analysis_report.md`
- `outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/MDD_WGCNA_Module_Report.md`
- `outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/MDD_WGCNA_Module_Report_FULL.md`
- `outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/MDD_WGCNA_Top_Module_Interpretation.md`

### manuscript-stale — excluded from Phase 3A staging

- `docs/MDD_Transcriptomic_Module_Analysis_Report_v2_NOTION_READY.md`
- `docs/MDD_Transcriptomic_Module_Analysis_Report_v2_clean_supervisor.md`
- `docs/Supervisor_Explanation_PreExternalValidation.md`
- `docs/WGCNA_Module_Analysis_Supervisor_Deep_Dive.md`
- `docs/WGCNA_Module_Analysis_Pseudocode_Walkthrough.md`
- `outputs/length_gt_100_clean_pipeline/14_Preliminary_Manuscript_Package/**`

### historical-valid — retain as evidence, not current navigation

- `docs/ML_PRIORITY_SCORE_CODE_AUDIT.md`
- `docs/ML_PRIORITY_SCORE_FINAL.md`
- `docs/ML_PRIORITY_SCORING_CORRECTNESS_FIX_CHANGELOG.md`
- `docs/CANONICALIZATION_AUDIT.md`
- `docs/DOWNSTREAM_CORRECTION_IMPACT.md`
- `docs/DOWNSTREAM_CORRECTED_RESULT_COMPARISON.md`
- `docs/PUBLICATION_AUDIT_NUMERICAL_CHECKS.md`
- `docs/PUBLICATION_CLAIM_SUPPORT_MATRIX.md`
- `docs/PUBLICATION_END_TO_END_AUDIT.md`
- `docs/PUBLICATION_REPRODUCIBILITY_GAPS.md`
- `outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/audit_all_regions_wgcna_algorithm_and_coverage/**`
- `outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/covariate_justification_audit/**`
- `outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/frontal_subregion_wgcna/**`
- `outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/module_level_ml2/**`
- `outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/biological_annotation/**`

### comparison-only code or historical runner

- `code/analyze_corrected_ml_priority_results.py` and
  `code/summarize_corrected_ml_downstream.py`: legacy values are used only for
  before/after comparison, not as canonical inputs.
- `code/length_gt_100_clean_pipeline/pipeline_lib.py` and its README: historical
  runner; steps 9–12 now require `--allow-legacy`.
- `outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/report_ready_original_groups_only_figures/rebuild_original_groups_only_figures.R`:
  retained historical figure generator, not the canonical pipeline entry point.

### root README

Legacy matches occur only below `<!-- HISTORICAL_CONTENT_START -->`. The active
top section carries the corrected 30/119/522 counts and canonical path. The
validator scans only the current-facing section and fails on old counts there.

No active corrected downstream script consumes the legacy ranking after Phase
3A. `code/canonical_downstream_paths.py` no longer defines legacy ranking paths.
