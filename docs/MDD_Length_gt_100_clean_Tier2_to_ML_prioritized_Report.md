> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Length > 100 MDD transcriptomics: Tier 2 candidate genes, enrichment, and Tier2-only ML prioritization

> **Superseded ML section:** this report uses the legacy 35/130 ranking and must be regenerated before publication. Current method/counts are in `ML_PRIORITY_SCORE_FINAL.md`.

## 1. Executive summary
This clean pipeline organizes the completed Length > 100 analysis into a reproducible path from filtered/normalized gene expression through LMM-derived Tier 2 candidate genes, enrichment, extended ML v1 prioritization, and stable ML-prioritized Tier 2 genes.

Tier 2 genes are exploratory candidate genes, not FDR-significant discoveries. Tier 2 genes define the main exploratory candidate pool. Tier 4 ranked signatures are used only as an additional ML-support layer to test whether candidate-gene prioritization is stable across feature-set definitions.

ML v1 is interpreted as gene prioritization, not validated diagnostic prediction. Classification performance is descriptive because feature sets were derived from the same dataset. The final prioritized genes should be interpreted as internally supported candidate genes, not validated biomarkers.

## 2. Clean pipeline overview
The clean path is: raw counts/expression inputs -> Length > 100 gene filtering -> GeTMM -> log2(GeTMM + 1) -> broad-region LMM -> standardized effect sizes -> evidence tiers -> frontal subregion fixed-effect sensitivity -> final Tier 2 candidate pools -> Tier 2 enrichment -> extended ML v1 using Tier2_fixed plus Tier4_top-k feature sets -> stable ML-prioritized Tier 2 genes -> report-ready images.

## 3. Input data and normalization
The clean pipeline reuses completed Length > 100 outputs from `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style`. Normalized expression uses filtered raw counts, GeTMM normalization, and log2(GeTMM + 1).

## 4. Broad-region LMM
Broad-region Tier 2 candidates come from the grouped-region LMM for Ventral Striatum, Insula, Hippocampal Formation, and broad Frontal Cortex.

## 5. Standardized effect sizes
Standardized effect size was defined as `beta_diagnosis / SD(expression of gene)`. The standardized-effect branch was used to define Tier 2 genes.

## 6. Evidence tiers
Tier 1: Storey q < 0.10 and abs(standardized effect size) > 0.8.
Tier 2: p_value_diagnosis < 0.01 and abs(standardized effect size) > 0.8.
Tier 3: p_value_diagnosis < 0.01 and abs(beta_diagnosis) > 0.2.
Tier 4 is a ranking of all tested genes, not a separate selected gene list or biological conclusion. It helps avoid relying on one hard cutoff by providing ranked-feature support in the ML prioritization stage.

## 7. Frontal subregion fixed-effect sensitivity analysis
Frontal subregion Tier 2 candidates come from the frontal fixed-effect sensitivity model because the random effect was not estimable after splitting. The documented model is Diagnosis + Age + Sex + RIN.

## 8. Final Tier 2 candidate gene pools
| region_or_subregion   | model_source                               |   Tier2_genes |   up_in_MDD |   down_in_MDD |
|:----------------------|:-------------------------------------------|--------------:|------------:|--------------:|
| Cingulate gyrus 25    | frontal_subregion_fixed_effect_sensitivity |           186 |         137 |            49 |
| Hippocampal Formation | grouped_region_LMM                         |           127 |          80 |            47 |
| Insula                | grouped_region_LMM                         |            70 |          37 |            33 |
| OFC BA11              | frontal_subregion_fixed_effect_sensitivity |            70 |          38 |            32 |
| Ventral Striatum      | grouped_region_LMM                         |           146 |          76 |            70 |
| dlPFC BA8/9           | frontal_subregion_fixed_effect_sensitivity |            72 |          24 |            48 |

## 9. Tier 2 enrichment
Tier 2 enrichment was reused from completed outputs. It used Tier 2 gene lists and region/subregion-specific tested-gene backgrounds. Enrichment is exploratory biological interpretation.

## 10. ML v1 extended prioritization
In ML, samples are rows and genes are features. ML v1 used Tier2_fixed as the main candidate set and Tier4_top-k ranked signatures as additional support/sensitivity. Tier 4 is a ranking of all tested genes, not a separate selected gene list. Training/testing was by cross-validation within the same dataset. Every sample can be held out as test in some fold, but the test data are not external.

Tier 2 defines the main exploratory candidate-gene pool using nominal diagnosis evidence and standardized effect size. Tier 4 is a ranked signature of all tested genes and is used only as an additional ML-support layer. In the extended ML v1 prioritization, models are trained using Tier2_fixed and Tier4_top-k feature sets. The final prioritization remains Tier 2-centered: genes are reported as ML-prioritized Tier 2 candidates when they belong to the Tier 2 pool and receive repeated support from the ML feature-importance outputs.

Tier2FixedOnly was tried as a simplified sensitivity analysis but produced less informative prioritization plots, so the extended ML v1 strategy is used as the main workflow.

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


## 11. Stable ML-prioritized Tier 2 genes
The clean prioritization reuses the previous extended ML v1 strategy. The final candidates are ML-prioritized Tier 2 genes. Tier 4 only provides additional ranked-feature support.

### ML-priority score and confidence-class thresholds

Genes were assigned to high- or moderate-confidence classes using predefined thresholds on the composite ML-priority score and supporting stability metrics. The exact thresholds were extracted from the prioritization script `D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py`.

The composite score was:

`ML_priority_score = 0.30 * selection_frequency_combined + 0.20 * model_support_scaled + 0.15 * sparse_model_support_binary + 0.15 * importance_rank_score + 0.10 * sign_consistency_score + 0.05 * lmm_evidence_scaled + 0.05 * v2_support_score`, clipped to 0-1.

High-confidence required `ML_priority_score >= 0.70`, `selection_frequency_combined >= 0.50`, `model_support_count >= 2`, sparse-model support, `sign_consistency >= 0.70`, and no explicit statistical/ML sign disagreement. Moderate-confidence required `ML_priority_score >= 0.45`, `selection_frequency_combined >= 0.25`, `model_support_count >= 1`, and no explicit statistical/ML sign disagreement. Labels are mutually exclusive because high-confidence is assigned first with `np.select([high, moderate], ...)`; moderate-confidence therefore excludes high-confidence rows.

Current final counts are 671 ranked Tier 2 gene-region pairs, 35 high-confidence, 130 additional moderate-confidence, and 165 high/moderate candidates total.


### Final ML-prioritized Tier 2 gene ranking

The final ranked table consists of Tier 2 candidate genes ordered by the composite `ML_priority_score`. The score integrates statistical diagnosis evidence with repeated ML feature-importance support across `Tier2_fixed` and `Tier4_top-k` feature-set definitions, model-family support, rank stability/importance, sign consistency, and strict-v2 support. Tier 4 is used only as an additional ranked-support layer; the reported final genes remain Tier 2 candidate genes.

These are internally supported exploratory candidates, not validated biomarkers. The full ranked table is available at `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended\final_ranked_ML_prioritized_Tier2_genes.csv`.

Extended ML prioritization identified 35 high-confidence and 130 additional moderate-confidence ML-prioritized Tier 2 gene-region pairs, for a total of 165 high/moderate candidates. Here, moderate-confidence is exclusive of high-confidence; high-confidence rows are not duplicated in the moderate-confidence file.

High/moderate counts by region:

| region_or_subregion   |   total_Tier2_genes_ranked |   high_confidence_count |   moderate_confidence_count |   exploratory_low_priority_count |   median_ML_priority_score |   max_ML_priority_score |   number_up_in_MDD |   number_down_in_MDD |   sign_agreement_rate |
|:----------------------|---------------------------:|------------------------:|----------------------------:|---------------------------------:|---------------------------:|------------------------:|-------------------:|---------------------:|----------------------:|
| Cingulate gyrus 25    |                        186 |                       8 |                          29 |                              149 |                   0.557238 |                0.831105 |                137 |                   49 |              0.901961 |
| Hippocampal Formation |                        127 |                       3 |                          25 |                               99 |                   0.535342 |                0.846446 |                 80 |                   47 |              0.827273 |
| Insula                |                         70 |                       6 |                          18 |                               46 |                   0.611851 |                0.889593 |                 37 |                   33 |              0.969231 |
| OFC BA11              |                         70 |                       3 |                          28 |                               39 |                   0.638499 |                0.83251  |                 38 |                   32 |              0.985294 |
| Ventral Striatum      |                        146 |                       8 |                          18 |                              120 |                   0.531131 |                0.871387 |                 76 |                   70 |              0.741667 |
| dlPFC BA8/9           |                         72 |                       7 |                          12 |                               53 |                   0.598543 |                0.831354 |                 24 |                   48 |              0.761194 |

Top 5 high/moderate-confidence genes per region:

| Region/subregion      | Gene         | Direction   | ML priority class                         |   ML priority score | Statistical evidence          | ML support summary                                        |
|:----------------------|:-------------|:------------|:------------------------------------------|--------------------:|:------------------------------|:----------------------------------------------------------|
| Cingulate gyrus 25    | LOC107985795 | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.831 | p=0.000963; std_effect=1.281  | freq=0.632; models=5; sparse=True; sign_consistency=1.000 |
| Cingulate gyrus 25    | LOC105377644 | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.823 | p=0.000564; std_effect=1.346  | freq=0.619; models=5; sparse=True; sign_consistency=1.000 |
| Cingulate gyrus 25    | NOC2LP2      | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.801 | p=0.000872; std_effect=1.175  | freq=0.607; models=5; sparse=True; sign_consistency=1.000 |
| Cingulate gyrus 25    | GEMIN7       | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.8   | p=0.00103; std_effect=1.282   | freq=0.604; models=5; sparse=True; sign_consistency=1.000 |
| Cingulate gyrus 25    | C1orf174     | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.787 | p=0.00141; std_effect=1.257   | freq=0.559; models=5; sparse=True; sign_consistency=1.000 |
| Hippocampal Formation | TOR2A        | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.846 | p=0.00231; std_effect=0.968   | freq=0.694; models=5; sparse=True; sign_consistency=1.000 |
| Hippocampal Formation | IL4I1        | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.793 | p=0.00315; std_effect=0.936   | freq=0.569; models=5; sparse=True; sign_consistency=1.000 |
| Hippocampal Formation | ZNF663P      | up_in_MDD   | Moderate-confidence ML-prioritized Tier 2 |               0.787 | p=0.000182; std_effect=1.164  | freq=0.489; models=5; sparse=True; sign_consistency=1.000 |
| Hippocampal Formation | HLA-F-AS1    | up_in_MDD   | Moderate-confidence ML-prioritized Tier 2 |               0.782 | p=0.000117; std_effect=1.138  | freq=0.440; models=5; sparse=True; sign_consistency=1.000 |
| Hippocampal Formation | TRMT2A       | down_in_MDD | High-confidence ML-prioritized Tier 2     |               0.758 | p=0.000714; std_effect=-1.055 | freq=0.508; models=5; sparse=True; sign_consistency=0.894 |
| Insula                | TMEM161B-DT  | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.89  | p=0.000306; std_effect=1.036  | freq=0.731; models=5; sparse=True; sign_consistency=1.000 |
| Insula                | BMS1P2       | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.874 | p=0.00021; std_effect=1.037   | freq=0.722; models=5; sparse=True; sign_consistency=1.000 |
| Insula                | BMS1P2-AGAP9 | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.846 | p=0.000209; std_effect=1.034  | freq=0.634; models=5; sparse=True; sign_consistency=1.000 |
| Insula                | ARL14EP-DT   | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.802 | p=0.0019; std_effect=0.891    | freq=0.562; models=5; sparse=True; sign_consistency=1.000 |
| Insula                | LOC101929710 | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.8   | p=0.00548; std_effect=0.804   | freq=0.579; models=5; sparse=True; sign_consistency=1.000 |
| OFC BA11              | QPRT         | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.833 | p=6e-06; std_effect=1.019     | freq=0.547; models=5; sparse=True; sign_consistency=1.000 |
| OFC BA11              | LOC105373386 | down_in_MDD | High-confidence ML-prioritized Tier 2     |               0.787 | p=0.00183; std_effect=-0.823  | freq=0.550; models=5; sparse=True; sign_consistency=0.956 |
| OFC BA11              | LOC105378913 | down_in_MDD | High-confidence ML-prioritized Tier 2     |               0.781 | p=0.0011; std_effect=-0.857   | freq=0.596; models=5; sparse=True; sign_consistency=0.966 |
| OFC BA11              | RAB24        | up_in_MDD   | Moderate-confidence ML-prioritized Tier 2 |               0.765 | p=0.00432; std_effect=0.822   | freq=0.478; models=5; sparse=True; sign_consistency=1.000 |
| OFC BA11              | MYBPC1       | up_in_MDD   | Moderate-confidence ML-prioritized Tier 2 |               0.756 | p=0.000277; std_effect=0.950  | freq=0.384; models=5; sparse=True; sign_consistency=1.000 |
| Ventral Striatum      | MORN5        | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.871 | p=0.000291; std_effect=0.920  | freq=0.724; models=5; sparse=True; sign_consistency=1.000 |
| Ventral Striatum      | FSHB         | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.837 | p=0.00131; std_effect=0.906   | freq=0.721; models=5; sparse=True; sign_consistency=1.000 |
| Ventral Striatum      | MRPL2        | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.832 | p=0.00111; std_effect=0.894   | freq=0.659; models=5; sparse=True; sign_consistency=1.000 |
| Ventral Striatum      | MOSPD3       | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.825 | p=0.000574; std_effect=0.997  | freq=0.634; models=5; sparse=True; sign_consistency=1.000 |
| Ventral Striatum      | POF1B        | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.821 | p=0.00293; std_effect=0.879   | freq=0.631; models=5; sparse=True; sign_consistency=1.000 |
| dlPFC BA8/9           | PPT2         | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.831 | p=0.000941; std_effect=0.947  | freq=0.634; models=5; sparse=True; sign_consistency=1.000 |
| dlPFC BA8/9           | NOL6         | down_in_MDD | High-confidence ML-prioritized Tier 2     |               0.821 | p=0.00428; std_effect=-0.804  | freq=0.631; models=5; sparse=True; sign_consistency=0.930 |
| dlPFC BA8/9           | LINC00269    | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.809 | p=0.000842; std_effect=0.948  | freq=0.551; models=5; sparse=True; sign_consistency=1.000 |
| dlPFC BA8/9           | LOC107984015 | down_in_MDD | High-confidence ML-prioritized Tier 2     |               0.802 | p=0.000534; std_effect=-1.001 | freq=0.620; models=5; sparse=True; sign_consistency=0.941 |
| dlPFC BA8/9           | LINC02897    | up_in_MDD   | High-confidence ML-prioritized Tier 2     |               0.788 | p=0.00146; std_effect=0.895   | freq=0.551; models=5; sparse=True; sign_consistency=1.000 |


| region_or_subregion   |   n_Tier2_genes |   n_high_confidence |   n_moderate_confidence |   n_with_sparse_support |   n_with_v2_support |   median_ML_priority_score | top_10_gene_symbols                                                                                         |
|:----------------------|----------------:|--------------------:|------------------------:|------------------------:|--------------------:|---------------------------:|:------------------------------------------------------------------------------------------------------------|
| Cingulate gyrus 25    |             186 |                   8 |                      29 |                     123 |                 162 |                   0.557238 | LOC107985795;LOC105377644;NOC2LP2;GEMIN7;C1orf174;EXOC3-AS1;PCBD1;LOC105372507;MEOX2;CD247                  |
| Hippocampal Formation |             127 |                   3 |                      25 |                      93 |                  91 |                   0.535342 | TOR2A;IL4I1;ZNF663P;HLA-F-AS1;TRMT2A;ARL14EP-DT;COQ4;CKAP2L;ADRA2B;TMEM161B-DT                              |
| Insula                |              70 |                   6 |                      18 |                      65 |                  68 |                   0.611851 | TMEM161B-DT;BMS1P2;BMS1P2-AGAP9;ARL14EP-DT;LOC101929710;POU6F2-AS2;BTN2A3P;LOC101928268;LINC00242;LINC01691 |
| OFC BA11              |              70 |                   3 |                      28 |                      69 |                  70 |                   0.638499 | QPRT;LOC105373386;LOC105378913;RAB24;MYBPC1;POLR3E;LOC105377975;TRPV3;C9orf163;NOL6                         |
| Ventral Striatum      |             146 |                   8 |                      18 |                      97 |                 116 |                   0.531131 | MORN5;FSHB;MRPL2;MOSPD3;POF1B;LINC01361;BMS1P2;KRT5;BMS1P2-AGAP9;OR4F3                                      |
| dlPFC BA8/9           |              72 |                   7 |                      12 |                      58 |                  71 |                   0.598543 | PPT2;NOL6;LINC00269;LOC107984015;LINC02897;LINC01990;LOC105371009;EMC6;CARF;MSANTD2-AS1                     |

## 12. Recommended figures
Recommended Tier 2 figures are in `outputs/length_gt_100_clean_pipeline/images/tier2`.
Recommended enrichment figures are in `outputs/length_gt_100_clean_pipeline/images/enrichment`.
Recommended ML figures are in `outputs/length_gt_100_clean_pipeline/images/ml_v1` and `outputs/length_gt_100_clean_pipeline/images/ml_prioritized`, and come from the extended ML v1 prioritization branch.

Corrected count-plot caption: **Extended ML-prioritized Tier 2 candidate genes. Bars show high-confidence genes and additional moderate-confidence genes, excluding high-confidence genes.** Corrected plot path: `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended\plots\Extended_ML_prioritized_Tier2_gene_counts_corrected.png`.

## 13. Conclusions
The clean path preserves the exploratory Tier 2 biological signal while using Tier 4 ranked signatures only as ML support for stable prioritization.

## 14. Limitations
Tier 2 selection used the full dataset, so downstream ML performance is optimistic. There is no external validation cohort in this clean pipeline. Candidate genes are not validated biomarkers.

## 15. Next steps
Use extended ML-prioritized Tier 2 genes for enrichment follow-up, AHBA/module mapping, and lower-dimensional pathway or module-level ML. Do not start deep learning until stable lower-dimensional features or external validation data are available.

## 16. Output files
Primary clean outputs are under `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline`. The reproducibility manifest is `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\pipeline_manifest.json`.

## 17. Glossary
LMM: linear mixed model.
Tier 2: nominal exploratory candidate genes with p_value_diagnosis < 0.01 and abs standardized effect size > 0.8.
ML v1: exploratory classical ML using Tier 2 expression features.
