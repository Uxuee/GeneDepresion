> **Method scope, 2026-11-21:** This unchanged scoring specification describes the historical canonical corrected weighted **v1 + strict-v2** benchmark, not v2 alone. Tier2-RRA is now the current primary exploratory ranking; the weighted confidence labels do not transfer. See [PIPELINE.md](PIPELINE.md).

# Final ML-Priority Score Specification

## Status

This is the authoritative specification for the corrected exploratory ML-prioritisation result in `outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/`. It supersedes the pre-correction numerical ranking retained in `10_ML_prioritized_Tier2_extended/` for provenance.

The candidate population contains 671 Tier 2 gene-region pairs: 30 high-confidence, 119 moderate-confidence, and 522 exploratory/low-confidence pairs. These labels describe internal prioritisation strength; they are not biomarker-validation claims.

## v1 evidence

Only `Tier2_fixed` contributes to the final v1 score. Historical Tier4 top-k and strict-univariate v1 runs are sensitivity/performance analyses and do not contribute.

One v1 event is one region × `Tier2_fixed` × feature-using model × repeat × fold combination. The five models are L1 logistic, L2 logistic, elastic-net logistic, linear SVM, and exploratory random forest. Thus `5 models × 5 repeats × 5 folds = 125 events` per region.

For L1 and elastic-net, a gene is selected when its coefficient is nonzero above tolerance. For L2 logistic, linear SVM, and random forest, selection is conservative top-20%-within-event importance support. Random-forest importance is unsigned.

## Strict v2 evidence

v2 performs feature selection in each outer training fold only. The selectors are train-only ANOVA F and train-only standardised effect, each at top-k 25, 50, 100, 250, and 500.

One v2 event is one region × selector × top-k × non-baseline model × outer fold. The denominator is dynamically obtained from the non-baseline grid: `2 selectors × 5 top-k values × 5 models × 3 folds = 150 events` per region.

This is a global-grid frequency: nested top-k settings can repeatedly select the same high-ranked gene. It is an internal stability/support measure, not a selection probability for an arbitrary model. The majority baseline uses no expression features and contributes no gene evidence or denominator event.

## Score components

- `selection_frequency_v1`: selected v1 events divided by 125.
- `selection_frequency_v2`: selected v2 events divided by 150 here.
- `combined_selection_frequency = 0.60 × v2_selection_frequency + 0.40 × v1_selection_frequency`.
- `model_support_scaled`: distinct supported model families, counted once across v1/v2, divided by the regional maximum.
- `sparse_model_support_binary`: 1 when L1 or elastic-net supports the gene, otherwise 0.
- `importance_rank_score = 1 - ((mean_rank - 1) / regional_max_mean_rank)`, clipped to 0–1.
- `sign_consistency`: the larger positive/negative coefficient fraction among L1, L2, elastic-net, and linear SVM. Random forest is excluded.
- `sign_consistency_score`: observed consistency, or neutral 0.5 when missing. Missing evidence cannot pass the high gate.
- `lmm_evidence_scaled`: within-region min-max scaling of `-log10(Diagnosis p) × abs(standardised Diagnosis effect)`.
- `v2_support_score`: 1 for any v2 support, otherwise 0.

```text
ML_priority_score =
    0.40 × combined_selection_frequency
  + 0.15 × model_support_scaled
  + 0.10 × sparse_model_support_binary
  + 0.15 × importance_rank_score
  + 0.10 × sign_consistency_score
  + 0.05 × lmm_evidence_scaled
  + 0.05 × v2_support_score
```

ROC-AUC, balanced accuracy, permutation performance, and other predictive metrics do not enter the gene score.

## Confidence gates

High confidence requires score ≥0.70, combined frequency ≥0.50, at least two model families, sparse support, observed sign consistency ≥0.70, and no observed disagreement with the gene-wise Diagnosis direction.

Moderate confidence requires score ≥0.45, combined frequency ≥0.25, at least one model family, and no observed direction disagreement. Remaining pairs are exploratory/low confidence.

## Worked examples

### LOC107985795, CG25

- v1: `108 / 125 = 0.864`
- v2: `70 / 150 = 0.4666666667`
- combined frequency: `0.6256`
- model support scaled: 1
- sparse support: 1
- importance-rank score: 0.7992512728361785
- corrected sign consistency: 1
- LMM evidence scaled: 0.493570482665677
- v2 support: 1
- corrected score: `0.7948062150587105`
- class: high confidence

### PRELID2, Insula

- v1: `34 / 125 = 0.272`
- v2: `30 / 150 = 0.20`
- combined frequency: 0.2288
- corrected score: `0.5626053853815516`
- class: exploratory/low confidence

PRELID2 shows why the gates matter: its score exceeds 0.45, but its combined frequency is below the moderate threshold of 0.25.
