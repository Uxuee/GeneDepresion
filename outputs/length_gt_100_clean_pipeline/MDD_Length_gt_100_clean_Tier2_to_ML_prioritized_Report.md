# Length > 100 MDD transcriptomics: Tier 2 candidate genes, enrichment, and Tier2-only ML prioritization

## 1. Executive summary
This clean pipeline organizes the completed Length > 100 analysis into a reproducible Tier2-only path from filtered/normalized gene expression through LMM-derived Tier 2 candidate genes, enrichment, ML v1 using Tier 2 features only, and final Tier2FixedOnly ML-prioritized genes.

Tier 2 genes are exploratory candidate genes, not FDR-significant discoveries. ML was used as a Tier2-guided gene-prioritization tool. Classification performance is reported only descriptively and should not be interpreted as unbiased diagnostic prediction. The final prioritized genes should be interpreted as internally supported candidate genes, not validated biomarkers.

## 2. Clean pipeline overview
The clean path is: raw counts/expression inputs -> Length > 100 gene filtering -> GeTMM -> log2(GeTMM + 1) -> broad-region LMM -> standardized effect sizes -> evidence tiers -> frontal subregion fixed-effect sensitivity -> final Tier 2 candidate pools -> Tier 2 enrichment -> ML v1 using Tier2_fixed only -> Tier2FixedOnly ML-prioritized genes -> report-ready images.

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
Tier 4 was explored previously but is excluded from this clean Tier2-only pipeline and report story.

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

## 10. ML v1 using Tier 2 genes only
In ML, samples are rows and genes are features. ML v1 used Tier 2 genes as features. Training/testing was by cross-validation within the same dataset. Every sample can be held out as test in some fold, but the test data are not external. Because Tier 2 was defined using the full dataset, performance is exploratory.

## 11. Tier2-only ML-prioritized genes
The clean prioritization uses only `feature_set == Tier2_fixed` feature-importance outputs. Tier4 and strict-univariate feature sets are excluded.

| region_or_subregion   |   n_Tier2_genes |   n_high_confidence |   n_moderate_confidence |   median_ML_priority_score |
|:----------------------|----------------:|--------------------:|------------------------:|---------------------------:|
| Cingulate gyrus 25    |             186 |                 113 |                      45 |                   0.803641 |
| Hippocampal Formation |             127 |                  81 |                      15 |                   0.820703 |
| Insula                |              70 |                  56 |                       1 |                   0.86121  |
| OFC BA11              |              70 |                  62 |                       0 |                   0.882313 |
| Ventral Striatum      |             146 |                  83 |                      15 |                   0.825653 |
| dlPFC BA8/9           |              72 |                  47 |                       1 |                   0.855496 |

## 12. Recommended figures
Recommended Tier 2 figures are in `outputs/length_gt_100_clean_pipeline/images/tier2`.
Recommended enrichment figures are in `outputs/length_gt_100_clean_pipeline/images/enrichment`.
Recommended ML figures are in `outputs/length_gt_100_clean_pipeline/images/ml_v1` and `outputs/length_gt_100_clean_pipeline/images/ml_prioritized`.

## 13. Conclusions
The clean Tier2-only path preserves the exploratory biological signal while separating candidate-gene discovery, enrichment interpretation, and ML-based prioritization.

## 14. Limitations
Tier 2 selection used the full dataset, so downstream ML performance is optimistic. There is no external validation cohort in this clean pipeline. Candidate genes are not validated biomarkers.

## 15. Next steps
Use prioritized Tier2FixedOnly genes for enrichment follow-up, AHBA/module mapping, and lower-dimensional pathway or module-level ML. Do not start deep learning until stable lower-dimensional features or external validation data are available.

## 16. Output files
Primary clean outputs are under `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline`. The reproducibility manifest is `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\pipeline_manifest.json`.

## 17. Glossary
LMM: linear mixed model.
Tier 2: nominal exploratory candidate genes with p_value_diagnosis < 0.01 and abs standardized effect size > 0.8.
ML v1: exploratory classical ML using Tier 2 expression features.
