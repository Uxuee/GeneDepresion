# ML-prioritisation scoring correctness fix

> **Correction audit trail, not the canonical method.** Early sections describe the retained pre-correction 35/130 result; the authoritative corrected result is 30/119/522 and is documented in `ML_PRIORITY_SCORE_FINAL.md`.

Date: 2026-09-24

## Scope

This change makes three narrowly scoped correctness fixes:

1. unsigned random-forest importance no longer enters sign consistency or ML direction;
2. majority-class baseline events no longer enter strict-v2 gene-level evidence or its denominator; and
3. missing sign consistency cannot satisfy the high-confidence sign gate.

It does not change the Tier2 candidate pool, v1 Tier2_fixed restriction, nested-v2 top-k design, selection rules, model/sparse support definitions, rank transformation, v1/v2 frequency weighting, score-component weights, or confidence thresholds.

## Pre-edit impact analysis

The retained ranking contains 671 gene-region rows, including 35 high-confidence and 130 moderate-confidence rows.

| Question | Retained-data result |
|---|---|
| Genes exposed to RF sign contamination | 626 rows have RF support; 539 have both RF and linear support. The exact number whose numeric sign consistency changes cannot be recovered without the missing event-level coefficients. |
| Gene-event selections from majority baseline | The executed raw grid implies 5,550 baseline feature rows per region and 33,300 across six regions: (25+50+100+250+500) × 2 selectors × 3 folds × 6 regions. The exact subset merging to the 671 Tier2 gene-region rows cannot be recovered without the missing v2 event table. |
| Current high-confidence genes with missing sign consistency | 0 |
| Current moderate-confidence genes with missing sign consistency | 0 |
| All retained rows with missing sign consistency | 57 |

The 626 RF-supported rows are the number exposed to the faulty calculation, not a claim that all 626 numeric values necessarily change. A value is unchanged when the signed coefficients and the added positive RF values already have compatible proportions.

## Denominator effect

The retained executed v2 frequencies reproduce an old denominator of 180 events per region:

~~~text
2 selectors × 5 top-k values × 6 models × 3 outer folds = 180
~~~

Removing the majority baseline leaves an expected executed grid of:

~~~text
2 selectors × 5 top-k values × 5 feature-using models × 3 outer folds = 150
~~~

The implementation does not hard-code 150. It removes baseline rows and then counts the distinct non-baseline event combinations actually present for each region. The regression test likewise constructs an event grid and checks the derived denominator.

## Weighted-score missing-sign handling

The weighted score continues to replace missing sign consistency with the neutral value 0.5. This was intentionally not changed because the requested correctness issue concerns the high-confidence gate, and changing the score imputation would be a separate methodological decision.

The high-confidence gate now requires both:

~~~text
sign consistency is present
sign consistency >= 0.70
~~~

Missing values therefore fail the gate rather than being silently treated as 1.

## Score formula retained

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

The tracked scorer had stale 0.30/0.20/0.15 literals even though the retained exported score and confirmed audit use the formula above. The literals were replaced with a single named SCORE_WEIGHTS mapping so the implementation and regression test now enforce the confirmed formula.

## Regeneration status

The exact original event tables were subsequently recovered from the local Git `main` branch and verified by SHA-256. No ML training was rerun. Corrected results were written separately to `outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/`, leaving the legacy outputs intact.

- Corrected v2 denominator: 150 events per region.
- Legacy high/moderate/low: 35/130/506.
- Corrected high/moderate/low: 30/119/522.
- Transitions: five high-to-moderate, sixteen moderate-to-low, and no other cross-category entries.
- Full tables, isolated-effect analyses, and exact real-gene arithmetic are stored in the corrected-results directory.

The regeneration also confirmed that the frozen CSV used historical first-three score weights 0.30/0.20/0.15. The corrected export uses the requested and documented 0.40/0.15/0.10, with the other four weights unchanged.

## Regression tests added

tests/test_ml_priority_scoring.py verifies:

1. random forest never enters sign consistency or direction;
2. majority baseline never enters v2 gene-selection counts or frequencies;
3. missing sign consistency cannot pass the high-confidence gate;
4. score weights are exactly 0.40/0.15/0.10/0.15/0.10/0.05/0.05; and
5. the v2 denominator equals the distinct non-baseline event grid.

The bundled Python runtime's original location remained blocked by an executable/DLL ACL. An ACL-clean temporary copy was used without changing the tests. All five tests passed.

## Files modified

- code/extract_stable_ML_prioritized_Tier2_genes.py
  - filters majority baseline before v2 ranking, selection, and denominator calculation;
  - creates a signed-only importance field for direction calculations;
  - makes missing signs fail the high-confidence gate;
  - centralises the confirmed score weights;
  - exports denominators and score components for future reproducibility;
  - validates required event inputs before clearing an existing output folder.
- tests/test_ml_priority_scoring.py
  - adds five focused regression tests.
- code/length_gt_100_clean_pipeline/pipeline_lib.py
  - corrects the generated manifest’s source-strategy wording.
- code/length_gt_100_clean_pipeline/README.md
  - corrects final v1/v2 method wording.
- outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended/README.md
  - corrects the retained output-package description.
- outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended/final_gene_ranking_method.md
  - documents the corrected sources, direction rules, missing-sign rule, performance exclusion, and score weights.
- docs/ML_PRIORITY_SCORING_CORRECTNESS_FIX_CHANGELOG.md
  - records this impact analysis and regeneration limitation.

## Interpretation

The corrected code remains an exploratory gene-prioritisation procedure. ROC-AUC and other predictive-performance measures do not enter the priority score, and the result must not be described as external validation or a validated biomarker set.
