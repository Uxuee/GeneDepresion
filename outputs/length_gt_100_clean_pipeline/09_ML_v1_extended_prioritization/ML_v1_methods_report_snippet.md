### ML v1 model settings and preprocessing

ML v1 was used for exploratory gene prioritization, not validated diagnostic prediction. The implementation used classical sklearn models only: a majority-class baseline, L2 logistic regression, L1 logistic regression, Elastic Net logistic regression, linear SVM, and exploratory Random Forest. XGBoost was coded as optional but was skipped in the completed run because `xgboost_available=false`.

Cross-validation was performed separately per region/subregion. The completed run used `RepeatedStratifiedKFold_5_fold_5_repeats_runtime_bounded`, with `StratifiedKFold(n_splits=5, shuffle=True)` repeated across random seeds 2026-2030. Splits were stratified by Diagnosis. Group-aware splitting via `StratifiedGroupKFold` was available in code when duplicated `BrainID` groups were present, but the completed ML v1 output reports `group_splitting_used=False` for all rows.

Preprocessing steps were fit within each cross-validation training fold and then applied to the corresponding held-out fold, preventing direct leakage through imputation or scaling. All models used `SimpleImputer(strategy="median")` inside an sklearn Pipeline. Logistic regression and linear SVM models also used `StandardScaler()` inside the Pipeline; Random Forest and the majority baseline were not scaled.

Feature sets were `Tier2_fixed`, `Tier4_top25`, `Tier4_top50`, `Tier4_top100`, `Tier4_top250`, `Tier4_top500`, and `strict_univariate_top100`. `Tier2_fixed` and Tier4 top-k sets are exploratory because they use full-dataset candidate/ranking information. Tier4 top-k sets were ordered by `abs(signed_standardized_logp_score)`, while the signed score preserved up-in-MDD versus down-in-MDD direction. `strict_univariate_top100` recomputed univariate feature selection inside each training fold and was retained as a leakage-controlled reference/sanity branch.

| model | exact settings |
|---|---|
| majority_class_baseline | `DummyClassifier(strategy="most_frequent")` |
| logistic_l2 | `LogisticRegression(penalty="l2", solver="liblinear", class_weight="balanced", max_iter=5000, random_state=seed)`; `C=1.0` by sklearn default |
| logistic_l1 | `LogisticRegression(penalty="l1", solver="liblinear", class_weight="balanced", max_iter=5000, random_state=seed)`; `C=1.0` by sklearn default |
| elastic_net_logistic | `LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, class_weight="balanced", max_iter=8000, random_state=seed)`; `C=1.0` by sklearn default |
| linear_svm | `LinearSVC(class_weight="balanced", max_iter=8000, random_state=seed)`; `C=1.0` by sklearn default; `decision_function` used for ROC/PR scoring |
| random_forest_exploratory | `RandomForestClassifier(n_estimators=100, max_features="sqrt", class_weight="balanced_subsample", min_samples_leaf=2, random_state=seed, n_jobs=-1)`; other unspecified parameters used sklearn defaults |
