# ML–RRA redesign: pre-implementation design note

Date: 2026-10-05  
Source commit: `11d0c63167090aef2c6aa014859bca8e27128eed`  
Experimental branch: `codex/rra-prioritization-redesign`

## Boundary

This is a new experimental prioritisation method. It does not edit, replace, relabel, benchmark, or feed back into the frozen corrected ML score. The canonical `F/M/S/R/D/L/V` scorer and its 671 gene-region rows remain unchanged.

Tier2, Tier4, the legacy 671 pairs, and old high/moderate labels are excluded from feature selection, fitting, tuning, and RRA. They may only be joined later as descriptive annotations, and this implementation does not perform that join.

## Exact input units and preserved input paths

The preserved strict-v2 loader defines six independent analysis units:

| Unit | Normalized-expression input |
|---|---|
| Ventral Striatum | `outputs/length_gt_100_elisa_style/normalized_expression_length_gt_100_by_region/normalized_GeTMM_log2_ventral_striatum.csv` |
| Insula | `outputs/length_gt_100_elisa_style/normalized_expression_length_gt_100_by_region/normalized_GeTMM_log2_insula.csv` |
| Hippocampal Formation | `outputs/length_gt_100_elisa_style/normalized_expression_length_gt_100_by_region/normalized_GeTMM_log2_hippocampal_formation.csv` |
| CG25 | `outputs/length_gt_100_elisa_style/frontal_subregion_sensitivity/subregion_normalized_expression_matrices/Cingulate_gyrus_25_Cg25_normalized_GeTMM_log2_length_gt_100.csv` |
| dlPFC BA8/9 | `outputs/length_gt_100_elisa_style/frontal_subregion_sensitivity/subregion_normalized_expression_matrices/Dorsolateral_prefrontal_cortex_dlPFC_BA8_9_normalized_GeTMM_log2_length_gt_100.csv` |
| OFC BA11 | `outputs/length_gt_100_elisa_style/frontal_subregion_sensitivity/subregion_normalized_expression_matrices/Orbitofrontal_OFC_BA11_normalized_GeTMM_log2_length_gt_100.csv` |

Metadata comes from `step1_from_raw/Length_gt_100/01_metadata/keep_all_samples/`: the three region-specific metadata files and `metadata_frontal_cortex.csv`, filtered by `OriginalRegion` for each frontal subregion.

The clean release clone intentionally omits these local canonical inputs. For this experiment they will be restored byte-for-byte from the preserved unpublished local-main Git objects into ignored paths in this isolated clone. Their blob IDs and SHA-256 values will be recorded in the run manifest. This is an input-availability constraint, not a change in scientific data.

## Labels

The preserved encoding is retained exactly:

- `MDD` → `1`
- `CTRL` or `CONTROL` → `0`
- any other or missing diagnosis → excluded

Expression is loaded as genes × samples and transposed to samples × genes after metadata alignment. Pooled Frontal Cortex is excluded.

## Audited preserved implementation

The preserved v1 pilot uses fixed Tier2 features and seeds 2026–2030. Its performance is optimistic because the Tier2 pool was defined on the full cohort. The corrected final scorer uses only v1 `Tier2_fixed` events, denominator 125.

The preserved strict-v2 script recomputes ANOVA-F or absolute standardized-effect rankings in each outer training split and uses median imputation, standard scaling for linear models, balanced class weights, and balanced accuracy for tuning. Its historical export counts selector × top-k combinations as separate evidence; the new design will not do that. Selector and k become hyperparameters chosen inside inner CV.

The frozen corrected score is:

`0.40F + 0.15M + 0.10S + 0.15R + 0.10D + 0.05L + 0.05V`

That scorer remains untouched.

## New repeated nested-CV design

- Outer CV: five stratified folds × five repeats.
- Outer seeds: 2026, 2027, 2028, 2029, 2030.
- Inner CV: three stratified folds; reduced only if an outer-training class has fewer than three samples, with the change recorded.
- Tune metric: balanced accuracy.
- Selector candidates: train-only ANOVA F and train-only absolute standardized mean difference.
- Candidate k: 25, 50, 100, 250, 500, capped at the number of genes.
- The selector is fitted separately inside every inner-training fold. The inner validation fold and outer test fold never fit imputation, scaling, selection, or model parameters.

## Five model families and grids

The runtime-bounded grids are the preserved strict-v2 compact grids:

| Model | Implementation | Grid |
|---|---|---|
| L1 logistic | `LogisticRegression(penalty="l1", solver="liblinear", class_weight="balanced", max_iter=4000)` | `C ∈ {0.01, 0.1, 1}` |
| Elastic Net logistic | `LogisticRegression(penalty="elasticnet", solver="saga", class_weight="balanced", max_iter=5000)` | `C ∈ {0.01, 0.1, 1}`, `l1_ratio=0.5` |
| L2 logistic | `LogisticRegression(penalty="l2", solver="liblinear", class_weight="balanced", max_iter=4000)` | `C ∈ {0.01, 0.1, 1}` |
| Linear SVM | `LinearSVC(class_weight="balanced", max_iter=6000)` | `C ∈ {0.01, 0.1, 1}` |
| Random Forest | `RandomForestClassifier(n_estimators=200, class_weight="balanced_subsample")` | `max_depth ∈ {3, None}`, `min_samples_leaf=2` |

Every hyperparameter candidate is crossed with selector method and k inside inner CV. Model-specific seeds are deterministically derived from region, outer repeat, outer fold, and model. No held-out metric weights gene rankings.

Median imputation is used for all models. Standard scaling follows imputation for the four linear models. Random Forest is not scaled.

## Fold ranking and two RRA levels

After inner tuning, the best complete training pipeline is refitted on the full outer-training fold. L1/Elastic Net omit coefficients with absolute value ≤ `1e-12`; L2/SVM rank all selected genes by absolute coefficient. Random Forest omits exact-zero importance and ranks the remaining selected genes. Linear signs are exported separately and never enter RRA.

Level 1 aggregates the 25 outer lists for each region/model. Level 2 aggregates the five model-level lists for each region. Both levels use the canonical R package `RobustRankAggreg::aggregateRanks(method="RRA")`. RRA scores are ranking/aggregation scores, not MDD association p-values.

## Environment audit

- R: 4.6.1
- `RobustRankAggreg`: 1.2.1, installed in an ignored project-local library
- Deterministic toy RRA wrapper test: passed
- Python runtime: 3.12.14, isolated copied runtime because the bundled source runtime has a machine ACL failure
- Python package pins are recorded in `environment/ml_rra/requirements.txt`

## Reusable code and blockers

Reusable elements are the six-unit input specification, diagnosis mapping, train-only ANOVA and standardized-effect selectors, median/scaling preprocessing, compact model grids, and existing canonical integrity tests.

The main reproducibility blocker is that normalized expression and metadata are intentionally absent from a clean release clone. The run must therefore use explicitly restored local canonical inputs and record their hashes. No other upstream analysis will be rerun.
