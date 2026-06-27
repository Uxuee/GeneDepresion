# ML v1 Exact Model Settings

Source inspected: `D:/genes/GeneDepresion/code/run_ML_pilot_Tier2_length_gt_100.py` and clean Step 9 reuse wrapper `code/length_gt_100_clean_pipeline/09_run_ML_v1_extended_prioritization.py`.

This is a documentation/audit output only. ML, LMM, and enrichment were not rerun.

## Cross-validation

- CV was performed separately for each analyzed region/subregion.
- Diagnosis labels were encoded as MDD=1 and Control=0, and splits were stratified by Diagnosis.
- The code selected 5 folds when the smaller class had at least 5 samples; otherwise 3 folds when the smaller class had at least 3 samples.
- In the completed ML v1 outputs, all rows used `RepeatedStratifiedKFold_5_fold_5_repeats_runtime_bounded`.
- Repeats used `RANDOM_SEEDS = [2026, 2027, 2028, 2029, 2030]`, one `StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)` per seed.
- Subject/group-aware splitting was available in code using `StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=2026)` if `BrainID` existed, had no missing values, and duplicated subjects were present.
- In the completed clean ML v1 performance table, `group_splitting_used=False` for all 252 rows, so the actual run did not use grouped splits.

## Preprocessing

- All models were wrapped in an sklearn `Pipeline`.
- Every pipeline started with `SimpleImputer(strategy="median")`.
- Linear/scaled models then used `StandardScaler()` inside the same Pipeline.
- Models with scaling: `logistic_l2`, `logistic_l1`, `elastic_net_logistic`, and `linear_svm`.
- Models without scaling: `majority_class_baseline`, `random_forest_exploratory`, and optional `xgboost_exploratory`.
- Because imputation and scaling were inside the sklearn Pipeline, preprocessing was fit on each training fold only and then applied to the corresponding held-out fold.
- No variance filter was applied in the fixed feature-set modes.
- For `strict_univariate_top100`, missing values were median-imputed within the current training fold before `f_classif`; feature selection was recomputed inside each fold.
- Missing expression values were handled by median imputation.

## Feature Sets

| feature_set | Gene selection | Selection timing/leakage status | Report role |
|---|---|---|---|
| `Tier2_fixed` | Tier 2 exploratory candidate genes: p_value_diagnosis < 0.01 and abs_standardized_effect_size_diagnosis > 0.8, available in the expression matrix. | Fixed/full-dataset gene set; exploratory and potentially optimistic for classifier performance. | Main biology-driven ML candidate pool and main downstream prioritization source. |
| `Tier4_top25` | Top 25 genes from the precomputed Tier 4 ranked signature. Genes are ordered by abs(signed_standardized_logp_score). | Fixed/full-dataset ranking; preliminary support/sensitivity feature set. | Support/sensitivity only. |
| `Tier4_top50` | Top 50 genes from the precomputed Tier 4 ranked signature by abs(signed_standardized_logp_score). | Fixed/full-dataset ranking; preliminary support/sensitivity feature set. | Support/sensitivity only. |
| `Tier4_top100` | Top 100 genes from the precomputed Tier 4 ranked signature by abs(signed_standardized_logp_score). | Fixed/full-dataset ranking; preliminary support/sensitivity feature set. | Support/sensitivity only. |
| `Tier4_top250` | Top 250 genes from the precomputed Tier 4 ranked signature by abs(signed_standardized_logp_score). | Fixed/full-dataset ranking; preliminary support/sensitivity feature set. | Support/sensitivity only. |
| `Tier4_top500` | Top 500 genes from the precomputed Tier 4 ranked signature by abs(signed_standardized_logp_score). | Fixed/full-dataset ranking; preliminary support/sensitivity feature set. | Support/sensitivity only. |
| `strict_univariate_top100` | Candidate universe is ranked-signature genes available in expression; within each training fold, SimpleImputer(strategy="median") is fit on training data and f_classif is used to select the top 100 genes. | Fold-internal feature selection; leakage-controlled pilot reference, still without external validation. | Reference/sanity branch, not the main candidate source. |


Tier 4 ranking used the signed standardized log-p score for direction and absolute value for top-k selection: `signed_standardized_logp_score = sign(standardized_effect_size_diagnosis) * -log10(p_value_diagnosis)`, and top-k genes were ordered by `abs(signed_standardized_logp_score)`. This allowed both up-in-MDD and down-in-MDD genes into Tier4_topK feature sets.

## Models and Hyperparameters

| model_name | sklearn_class | key settings | feature importance | notes |
|---|---|---|---|---|
| `majority_class_baseline` | `sklearn.dummy.DummyClassifier` | strategy="most_frequent". Other parameters use sklearn defaults. | None; no coefficients or feature_importances_ extracted. | Included to compare ML models against a trivial majority-class baseline. |
| `logistic_l2` | `sklearn.linear_model.LogisticRegression` | penalty="l2"; solver="liblinear"; class_weight="balanced"; max_iter=5000; random_state=CV seed; C not specified, so sklearn default C=1.0. Other parameters use sklearn defaults. | model.coef_; nonzero coefficients where abs(coef) > 1e-12 saved with signed direction. | Preprocessing and model are fit inside each training fold only. |
| `logistic_l1` | `sklearn.linear_model.LogisticRegression` | penalty="l1"; solver="liblinear"; class_weight="balanced"; max_iter=5000; random_state=CV seed; C not specified, so sklearn default C=1.0. Other parameters use sklearn defaults. | model.coef_; nonzero coefficients where abs(coef) > 1e-12 saved with signed direction. | Used for sparse-model support in downstream prioritization. |
| `elastic_net_logistic` | `sklearn.linear_model.LogisticRegression` | penalty="elasticnet"; solver="saga"; l1_ratio=0.5; class_weight="balanced"; max_iter=8000; random_state=CV seed; C not specified, so sklearn default C=1.0. Other parameters use sklearn defaults. | model.coef_; nonzero coefficients where abs(coef) > 1e-12 saved with signed direction. | Used for sparse-model support in downstream prioritization. |
| `linear_svm` | `sklearn.svm.LinearSVC` | class_weight="balanced"; max_iter=8000; random_state=CV seed; C not specified, so sklearn default C=1.0. probability is not available for LinearSVC; decision_function was used for ROC/PR scoring. Other parameters use sklearn defaults. | model.coef_; nonzero coefficients where abs(coef) > 1e-12 saved with signed direction. | Best model in several regions by ROC-AUC in the completed run. |
| `random_forest_exploratory` | `sklearn.ensemble.RandomForestClassifier` | n_estimators=100; max_features="sqrt"; class_weight="balanced_subsample"; min_samples_leaf=2; random_state=CV seed; n_jobs=-1. max_depth, min_samples_split, criterion, bootstrap, and other parameters were not specified and therefore use sklearn defaults. | model.feature_importances_; positive importances saved as unsigned tree_importance_exploratory with direction="not_signed". | Labeled exploratory only; tree importances are not signed effect directions. |
| `xgboost_exploratory` | `xgboost.XGBClassifier` | Code path if installed: n_estimators=150; max_depth=2; learning_rate=0.05; subsample=0.8; colsample_bytree=0.8; eval_metric="logloss"; random_state=CV seed; n_jobs=1. Actual completed ML v1 run skipped this model because XGBoost was not installed/available. | If run, feature_importances_ would be saved as unsigned tree_importance_exploratory with direction="not_signed". | Not present in completed performance outputs. |

## Metrics

The fold-level metrics computed were:

- accuracy
- balanced accuracy
- ROC-AUC
- PR-AUC, via `average_precision_score`
- F1 score with `zero_division=0`
- sensitivity/recall for MDD: TP / (TP + FN)
- specificity for Control: TN / (TN + FP)
- confusion matrix totals: TN, FP, FN, TP

Per-model summary tables report mean and sample standard deviation across folds/repeats for accuracy, balanced accuracy, ROC-AUC, PR-AUC, F1, sensitivity, and specificity.

## Feature Importance and Direction

For linear models, feature importance was extracted from `model.coef_` after fitting the fold-specific pipeline. Coefficients with `abs(coef) > 1e-12` were saved. The sign was preserved as `positive_MDD` when the coefficient was greater than zero and `negative_MDD` when less than zero.

For Random Forest, feature importance was extracted from `model.feature_importances_`; values greater than zero were saved as `tree_importance_exploratory` with `direction="not_signed"`.

The output `feature_importance/fold_level_feature_importance.csv` stores per-fold importance events. The summary file aggregates by analysis level, region/subregion, feature set, model, and gene using:

- `selection_frequency`
- `mean_importance`
- `mean_abs_importance`
- `median_abs_importance`
- `n_positive_direction`
- `n_negative_direction`
- `direction_consistency`
- `dominant_direction`

For signed models, `direction_consistency` is the larger of positive or negative counts divided by the number of times the gene was selected. Tree importances are unsigned and recorded as `not_signed`.

## Permutation Testing

The script selected the best non-baseline model per region/subregion by sorting on `roc_auc_mean` and then `balanced_accuracy_mean`. A preliminary label-permutation test was run with `PERMUTATION_SEEDS = [9001, ..., 9020]`, giving 20 permutations. For each permutation, labels were shuffled, CV splits were rebuilt, and only the first three CV splits were used to bound runtime.

| region/subregion | feature_set | model | permutations | CV splits per permutation | empirical p ROC-AUC | empirical p balanced accuracy |
|---|---|---|---:|---:|---:|---:|
| Ventral Striatum | `Tier4_top500` | `linear_svm` | 20 | 3 | 0.047619 | 0.047619 |
| Insula | `Tier2_fixed` | `linear_svm` | 20 | 3 | 0.047619 | 0.047619 |
| Hippocampal Formation | `Tier4_top500` | `logistic_l2` | 20 | 3 | 0.047619 | 0.047619 |
| Cingulate gyrus 25 | `Tier4_top500` | `linear_svm` | 20 | 3 | 0.047619 | 0.047619 |
| dlPFC BA8/9 | `Tier4_top250` | `linear_svm` | 20 | 3 | 0.047619 | 0.047619 |
| OFC BA11 | `Tier2_fixed` | `linear_svm` | 20 | 3 | 0.047619 | 0.047619 |


Output file: `D:/genes/GeneDepresion/outputs/length_gt_100_elisa_style/ML_pilot_Tier2/permutation_test_results.csv`.

## Interpretation Caveats

ML v1 was used for exploratory gene prioritization, not validated diagnostic prediction. Fixed Tier 2 and Tier 4 feature sets were derived from the full dataset, so classifier performance from those feature sets may be optimistic. The `strict_univariate_top100` branch is more leakage-controlled because feature selection is recomputed inside each training fold, but it remains an internal pilot without external validation.
