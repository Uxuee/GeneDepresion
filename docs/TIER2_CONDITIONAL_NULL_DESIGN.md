# Tier2 conditional label-permutation/null design

## Status and question

This design was frozen before inspecting any permuted-label result. It starts from benchmark commit `083c917170f577d06514e04fc51e34db7dd87a35` and does not alter the observed Tier2-RRA, LASSO benchmark, Stability Selection benchmark, all-gene RRA, Tier2 identities, or upstream analyses.

The test is conditional on the 671 Tier2 region–gene pairs that were selected using the real Diagnosis labels. It asks whether, once this enriched candidate universe is fixed, the true labels yield more coherent downstream ML prioritisation than labels shuffled within region. It is not a global-null test of expression → gene-wise model → Tier2 → ML. Such a test would have to regenerate Tier2 after every shuffle and is outside this phase.

## Frozen candidate universe

The sole candidate source is `outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv`, audited SHA-256 `1672ddf95a4f4c9e44ae42e7accc208c311e2e1da8533e1808dcf91bbf7b67a1`.

| Analysis unit | Tier2 pairs |
|---|---:|
| Ventral Striatum | 146 |
| Insula | 70 |
| Hippocampal Formation | 127 |
| CG25 | 186 |
| dlPFC BA8/9 | 72 |
| OFC BA11 | 70 |
| **Total** | **671** |

Every permutation retains the same expression matrix, sample identifiers, gene identifiers, and region-specific class counts. Only the Diagnosis vector is permuted, separately within each region. Seeds are `41001` through `41020`; a stable hash of region and permutation seed produces each region-specific shuffle. The complete sample-to-permuted-label mapping is saved.

## Frozen Tier2-RRA implementation

The null run reproduces `code/tier2_ml_rra/01_run_tier2_repeated_ml.py` and `02_run_tier2_rra.R`:

- five outer folds and five repeats, with outer split seeds 2026–2030 regenerated from each permuted label vector;
- three-fold stratified inner CV (two only if a training class has fewer than three samples), tuned on balanced accuracy;
- training-only median imputation and, for non-random-forest models, training-only standardisation;
- no feature selector or top-k screen inside the frozen Tier2 analysis;
- L1 logistic: `C = 0.01, 0.1, 1`, liblinear, balanced class weights;
- Elastic Net logistic: the same C grid, `l1_ratio = 0.5`, saga, balanced class weights;
- L2 logistic: the same C grid, liblinear, balanced class weights;
- linear SVM: the same C grid and balanced class weights;
- Random Forest: 200 trees, `max_depth = 3` or unlimited, `min_samples_leaf = 2`, balanced-subsample weights;
- deterministic inner/model seeds derived by the existing stable SHA-256 seed function;
- sparse L1/Elastic-Net lists retain coefficients with `abs(coef) > 1e-12`; L2 and SVM rank all features; Random Forest retains strictly positive importance;
- each non-empty fold list is ranked by decreasing absolute importance, then stable gene ID. Empty sparse lists remain empty evidence;
- level 1 RRA aggregates the 25 outer-fit lists independently for each region/model with `RobustRankAggreg::aggregateRanks`, `N =` the fixed regional universe, `method = "RRA"`, `full = FALSE`, `exact = FALSE`;
- level 2 RRA aggregates the five model-level lists using the same function and universe;
- the most-frequent DummyClassifier is prediction context only and cannot enter gene rankings or RRA.

Each permutation therefore contains 750 feature-model fits (6 regions × 25 splits × 5 models) plus 150 dummy references. No observed preprocessing, parameters, coefficients, importances, rankings, or RRA outputs are reused.

## Frozen LASSO-alone null benchmark

The LASSO null uses the `logistic_l1` rows generated within the complete permuted Tier2-RRA run; it does not fit a different model:

- the same 25 permuted-label outer splits per region;
- the same training-only preprocessing and three-fold tuning over `C = 0.01, 0.1, 1`;
- selection is `abs(coef) > 1e-12`;
- primary evidence is selection count divided by exactly 25;
- deterministic ordering is decreasing selection frequency, lower mean selected rank, higher mean absolute coefficient, then gene ID;
- absent genes receive zero selection frequency and a shared worst-rank block for correlation summaries.

Thus each complete permutation contains 150 LASSO outer fits as a strict subset of the 750 feature-model fits; no extra LASSO retraining is required.

## Pilot protocol and resumability

One complete permutation is run first, validated, timed, and rerun with `--resume` to prove no completed fit is recomputed. Its raw-output checksums must remain identical. Only then may the run extend to the prospectively fixed 20-permutation pilot. A later 99-permutation extension can append seeds without recomputing completed directories.

The pilot empirical tail probability is `(1 + count(null >= observed)) / (1 + 20)` for statistics where larger means stronger structure. For RRA-score concentration, where smaller means stronger structure, the inequality is reversed. The minimum attainable pilot value is `1/21 = 0.0476`; all probabilities are labelled pilot estimates.

## Predeclared endpoints

No composite score is permitted. Separate endpoints are:

1. RRA resample stability: pairwise top-10, top-25, and top-20% Jaccard plus rank correlations.
2. Model-to-model agreement and recurrence/support across models.
3. RRA consensus strength: minimum score, median score at each top-set definition, and counts supported by at least 3, 4, or all 5 model-level lists.
4. Observed top-10/top-25 gene rank and RRA-score behavior under the null.
5. LASSO support and stability: maximum and top-set mean selection frequency, threshold counts, and resample Jaccard.
6. Secondary prediction diagnostics for each model and the region-wise best-of-five balanced accuracy, compared with the corresponding null maximum.

Prediction is not the primary endpoint. Because Tier2 was created from the full cohort using real labels, observed prediction separation is expected to be optimistic and cannot support a clinical classifier claim.

## Prospective interpretation

- **A:** RRA stability/concentration exceeds null and LASSO agrees → conditional label-dependent structure supports exploratory prioritisation.
- **B:** RRA exceeds null but LASSO does not → ensemble-only structure, interpreted cautiously.
- **C:** LASSO exceeds null but RRA does not → the sparse single model may be more defensible.
- **D:** observed ranking stability resembles null → downgrade or remove ML prioritisation as a major result.
- **E:** prediction exceeds null but ranking stability does not → prediction alone does not justify gene prioritisation.

The interpretation rules, model grids, seeds, summaries, and tails are fixed before seeing permutation results.
