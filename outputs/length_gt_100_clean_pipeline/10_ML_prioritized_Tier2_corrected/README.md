# Corrected ML-prioritised Tier 2 results

This directory contains the September 2026 correctness regeneration of the exploratory ML-prioritised Tier 2 ranking. The frozen legacy result remains unchanged in `10_ML_prioritized_Tier2_extended`.

`canonical_manifest.json` is the machine-readable publication manifest. The authoritative human-readable method is `docs/ML_PRIORITY_SCORE_FINAL.md`.

## Inputs and scope

- Candidate pool: existing Tier 2 gene-region pairs only.
- v1 evidence: `Tier2_fixed` only. Historical v1 Tier4 top-k runs do not contribute to this final score.
- v1 event denominator: 125 per region, from five feature-using models × five repeats × five folds.
- v2 evidence: nested train-only ANOVA-F and standardised-effect selection at top 25, 50, 100, 250, and 500 across three outer folds.
- v2 denominator: the actual distinct non-baseline event grid. It is 150 per region in this run.
- No GeTMM, LMM, tier, enrichment, v1 training, or v2 training stage was rerun.

## Correctness rules

- The majority-class baseline supplies no gene-level selection evidence.
- Random forest contributes nondirectional selection, importance/rank, and model-support evidence only.
- Sign consistency and ML direction use signed coefficients from L1 logistic, L2 logistic, elastic-net logistic, and linear SVM only.
- Missing sign consistency is scored neutrally as 0.5 inside the weighted score but cannot satisfy the high-confidence directional gate.
- ROC-AUC, balanced accuracy, and other predictive-performance metrics do not enter the gene-priority score.

## Score

~~~text
ML_priority_score =
    0.40 * combined_selection_frequency
  + 0.15 * model_support_scaled
  + 0.10 * sparse_model_support_binary
  + 0.15 * importance_rank_score
  + 0.10 * sign_consistency_score
  + 0.05 * lmm_evidence_scaled
  + 0.05 * v2_support_score
~~~

The frozen legacy CSV is reproducible with the historical first-three weights 0.30/0.20/0.15. This corrected export uses the confirmed 0.40/0.15/0.10 weights above.

## Main audit outputs

- `ML_prioritized_Tier2_genes_all_ranked.csv`: complete corrected ranking and score components.
- `before_after_all_genes.csv`: old-versus-corrected comparison for all 671 rows.
- `confidence_transition_summary.csv`: confidence-label transition counts.
- `isolated_correction_confidence_counts.csv`: counts under isolated corrections.
- `isolated_correction_effect_summary.csv`: score and label impact by correction.
- `top20_absolute_score_changes.csv`: largest score changes.
- `top20_absolute_rank_changes.csv`: largest within-region rank changes.
- `top20_RF_sign_consistency_changes.csv`: largest directional-consistency changes.
- `top20_baseline_removal_v2_frequency_changes.csv`: largest v2-frequency changes.
- `real_gene_score_reconstruction.csv`: exact arithmetic for representative genes.
- `corrected_regeneration_report.md`: concise audit report.

## Interpretation

These are exploratory within-study gene-prioritisation results, not validated biomarkers or externally validated diagnostic predictors.
