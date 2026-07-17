# Published AHBA Module-Score ML Report

Published AHBA modules were used as external predefined biological feature sets. Module-score ML was exploratory and was used to assess whether AHBA-informed module-level expression patterns carried MDD/Control signal, not to establish a validated diagnostic classifier.

## Methods

- Main feature set: all 32 published AHBA modules.
- Module scores: within each CV fold, each gene was z-scored using training samples only; test samples were transformed using training means and SDs; module score was the mean standardized expression of available genes in the module.
- Models: majority baseline, logistic regression L2/L1/elastic net, linear SVM, and exploratory random forest.
- CV repeats/seeds: 20; grouped stratified folds were used when BrainID repeats occurred within a region.

## Main all32-module performance

| region                | model                     | cv_type                        |   n_samples |   n_MDD |   n_Control |   balanced_accuracy_mean |   balanced_accuracy_std |   roc_auc_mean |   roc_auc_std | leakage_status                    |
|:----------------------|:--------------------------|:-------------------------------|------------:|--------:|------------:|-------------------------:|------------------------:|---------------:|--------------:|:----------------------------------|
| Frontal_Cortex        | linear_SVM                | StratifiedGroupKFold           |         122 |      63 |          59 |                 0.58668  |               0.106298  |       0.632764 |      0.128532 | predefined_external_all32_modules |
| Frontal_Cortex        | logistic_L2               | StratifiedGroupKFold           |         122 |      63 |          59 |                 0.56584  |               0.102079  |       0.589336 |      0.12644  | predefined_external_all32_modules |
| Frontal_Cortex        | elastic_net               | StratifiedGroupKFold           |         122 |      63 |          59 |                 0.551672 |               0.103437  |       0.579871 |      0.125265 | predefined_external_all32_modules |
| Frontal_Cortex        | logistic_L1               | StratifiedGroupKFold           |         122 |      63 |          59 |                 0.537226 |               0.0960907 |       0.558339 |      0.119148 | predefined_external_all32_modules |
| Frontal_Cortex        | random_forest_exploratory | StratifiedGroupKFold           |         122 |      63 |          59 |                 0.532475 |               0.115493  |       0.553464 |      0.13546  | predefined_external_all32_modules |
| Frontal_Cortex        | majority_baseline         | StratifiedGroupKFold           |         122 |      63 |          59 |                 0.495    |               0.05      |       0.5      |      0        | predefined_external_all32_modules |
| Hippocampal_Formation | logistic_L2               | StratifiedKFold_repeated_seeds |          43 |      24 |          19 |                 0.540083 |               0.164577  |       0.483833 |      0.175832 | predefined_external_all32_modules |
| Hippocampal_Formation | elastic_net               | StratifiedKFold_repeated_seeds |          43 |      24 |          19 |                 0.537667 |               0.162107  |       0.486917 |      0.193618 | predefined_external_all32_modules |
| Hippocampal_Formation | linear_SVM                | StratifiedKFold_repeated_seeds |          43 |      24 |          19 |                 0.5085   |               0.142631  |       0.46375  |      0.179785 | predefined_external_all32_modules |
| Hippocampal_Formation | logistic_L1               | StratifiedKFold_repeated_seeds |          43 |      24 |          19 |                 0.5085   |               0.151174  |       0.472875 |      0.198205 | predefined_external_all32_modules |
| Hippocampal_Formation | majority_baseline         | StratifiedKFold_repeated_seeds |          43 |      24 |          19 |                 0.5      |               0         |       0.5      |      0        | predefined_external_all32_modules |
| Hippocampal_Formation | random_forest_exploratory | StratifiedKFold_repeated_seeds |          43 |      24 |          19 |                 0.37675  |               0.145076  |       0.320917 |      0.160602 | predefined_external_all32_modules |
| Insula                | linear_SVM                | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.630917 |               0.138324  |       0.642483 |      0.176715 | predefined_external_all32_modules |
| Insula                | random_forest_exploratory | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.592583 |               0.1494    |       0.639858 |      0.180509 | predefined_external_all32_modules |
| Insula                | logistic_L2               | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.56625  |               0.141076  |       0.570667 |      0.186149 | predefined_external_all32_modules |
| Insula                | elastic_net               | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.53075  |               0.13097   |       0.528267 |      0.177291 | predefined_external_all32_modules |
| Insula                | majority_baseline         | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.5      |               0         |       0.5      |      0        | predefined_external_all32_modules |
| Insula                | logistic_L1               | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.494417 |               0.136762  |       0.49555  |      0.176531 | predefined_external_all32_modules |
| Ventral_Striatum      | linear_SVM                | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.595333 |               0.141175  |       0.630817 |      0.173329 | predefined_external_all32_modules |
| Ventral_Striatum      | logistic_L1               | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.560833 |               0.160791  |       0.563317 |      0.185709 | predefined_external_all32_modules |
| Ventral_Striatum      | logistic_L2               | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.55875  |               0.133167  |       0.58155  |      0.173889 | predefined_external_all32_modules |
| Ventral_Striatum      | elastic_net               | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.554333 |               0.136805  |       0.5675   |      0.182409 | predefined_external_all32_modules |
| Ventral_Striatum      | random_forest_exploratory | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.53225  |               0.140679  |       0.577575 |      0.157796 | predefined_external_all32_modules |
| Ventral_Striatum      | majority_baseline         | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.5      |               0         |       0.5      |      0        | predefined_external_all32_modules |

## Best model per region/feature set

| region                | feature_set                                | model                     |   balanced_accuracy_mean |   roc_auc_mean | leakage_status                                   |
|:----------------------|:-------------------------------------------|:--------------------------|-------------------------:|---------------:|:-------------------------------------------------|
| Frontal_Cortex        | all32_predefined_modules                   | linear_SVM                |                 0.58668  |       0.632764 | predefined_external_all32_modules                |
| Frontal_Cortex        | post_enrichment_FDR005_modules_sensitivity | majority_baseline         |                 0.495    |       0.5      | post_enrichment_selected_sensitivity_exploratory |
| Frontal_Cortex        | post_enrichment_FDR010_modules_sensitivity | random_forest_exploratory |                 0.523569 |       0.51604  | post_enrichment_selected_sensitivity_exploratory |
| Frontal_Cortex        | post_enrichment_M26_only_sensitivity       | linear_SVM                |                 0.50116  |       0.506133 | post_enrichment_selected_sensitivity_exploratory |
| Hippocampal_Formation | all32_predefined_modules                   | logistic_L2               |                 0.540083 |       0.483833 | predefined_external_all32_modules                |
| Hippocampal_Formation | post_enrichment_FDR005_modules_sensitivity | logistic_L2               |                 0.538833 |       0.560708 | post_enrichment_selected_sensitivity_exploratory |
| Hippocampal_Formation | post_enrichment_FDR010_modules_sensitivity | elastic_net               |                 0.54175  |       0.555667 | post_enrichment_selected_sensitivity_exploratory |
| Hippocampal_Formation | post_enrichment_M26_only_sensitivity       | linear_SVM                |                 0.580667 |       0.536833 | post_enrichment_selected_sensitivity_exploratory |
| Insula                | all32_predefined_modules                   | linear_SVM                |                 0.630917 |       0.642483 | predefined_external_all32_modules                |
| Insula                | post_enrichment_FDR005_modules_sensitivity | linear_SVM                |                 0.59475  |       0.62815  | post_enrichment_selected_sensitivity_exploratory |
| Insula                | post_enrichment_FDR010_modules_sensitivity | linear_SVM                |                 0.609667 |       0.626317 | post_enrichment_selected_sensitivity_exploratory |
| Insula                | post_enrichment_M26_only_sensitivity       | linear_SVM                |                 0.64575  |       0.7107   | post_enrichment_selected_sensitivity_exploratory |
| Ventral_Striatum      | all32_predefined_modules                   | linear_SVM                |                 0.595333 |       0.630817 | predefined_external_all32_modules                |
| Ventral_Striatum      | post_enrichment_FDR005_modules_sensitivity | random_forest_exploratory |                 0.648917 |       0.742592 | post_enrichment_selected_sensitivity_exploratory |
| Ventral_Striatum      | post_enrichment_FDR010_modules_sensitivity | logistic_L1               |                 0.60475  |       0.663367 | post_enrichment_selected_sensitivity_exploratory |
| Ventral_Striatum      | post_enrichment_M26_only_sensitivity       | random_forest_exploratory |                 0.588167 |       0.590425 | post_enrichment_selected_sensitivity_exploratory |

## Limitations

This is exploratory module-level ML. All 32-module models use external predefined AHBA modules and avoid feature-selection leakage, but they are still evaluated only within this dataset. Post-enrichment selected module subsets are sensitivity analyses and are more optimistic because the modules were selected after enrichment.