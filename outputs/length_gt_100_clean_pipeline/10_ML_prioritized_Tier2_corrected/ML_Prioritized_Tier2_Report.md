# ML-prioritized Tier 2 genes for MDD/control separation

## Purpose

This analysis identifies Tier 2 exploratory DE candidate genes repeatedly prioritized by ML models for separating MDD vs Control.

## Important Interpretation

This is gene prioritization, not diagnostic classifier validation. ML v1 is used as Tier2-guided feature prioritization / exploratory classification screening. ML v2 is used as stricter support when a Tier 2 gene was selected in leakage-controlled folds.

## Input Data

Files used:

- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\DE_candidate_gene_sets\Tier2_working_exploratory\Ventral_Striatum__Tier2_all.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\DE_candidate_gene_sets\Tier2_working_exploratory\Insula__Tier2_all.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\DE_candidate_gene_sets\Tier2_working_exploratory\Hippocampal_Formation__Tier2_all.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\DE_candidate_gene_sets\Tier2_working_exploratory\Cingulate_gyrus_25__Tier2_all.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\DE_candidate_gene_sets\Tier2_working_exploratory\dlPFC_BA8_9__Tier2_all.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\DE_candidate_gene_sets\Tier2_working_exploratory\OFC_BA11__Tier2_all.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_pilot_Tier2\feature_importance\fold_level_feature_importance.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_strict_v2\feature_stability\selected_genes_per_fold_all_models.csv`

Missing optional files:

- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_pilot_Tier2\feature_importance\feature_importance_summary.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_pilot_Tier2\ML_performance_summary.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_pilot_Tier2\best_model_by_region.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_strict_v2\feature_stability\selection_frequency_all_models.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_strict_v2\ML_strict_v2_performance_summary.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_strict_v2\best_strict_model_by_region.csv`

## Methods

The candidate pool was restricted to existing Tier 2 genes. Final v1 scoring used Tier2_fixed only; historical Tier4 v1 sets did not enter the final score. Strict v2 used nested train-only feature-selection/top-k evidence. Random forest was retained as non-directional selection/rank evidence, the majority-class baseline was excluded from gene-level evidence, and missing sign consistency could not satisfy the high-confidence sign gate. Predictive-performance metrics such as ROC-AUC did not enter the gene score. The component weights were 0.40/0.15/0.10/0.15/0.10/0.05/0.05.

## Results By Region/Subregion

| region_or_subregion   |   n_Tier2_genes |   n_high_confidence |   n_moderate_confidence |   n_with_sparse_support |   n_with_v2_support |   median_ML_priority_score | top_10_gene_symbols                                                                                      |
|:----------------------|----------------:|--------------------:|------------------------:|------------------------:|--------------------:|---------------------------:|:---------------------------------------------------------------------------------------------------------|
| Cingulate gyrus 25    |             186 |                   6 |                      30 |                     123 |                 162 |                   0.481605 | LOC107985795;LOC105377644;NOC2LP2;GEMIN7;C1orf174;PCBD1;EXOC3-AS1;LOC105372507;DLGAP4-AS1;MEOX2          |
| Hippocampal Formation |             127 |                   2 |                      24 |                      93 |                  91 |                   0.459393 | TOR2A;IL4I1;ZNF663P;HLA-F-AS1;TRMT2A;COQ4;CKAP2L;ARL14EP-DT;TMEM161B-DT;ADRA2B                           |
| Insula                |              70 |                   5 |                      15 |                      65 |                  68 |                   0.534709 | TMEM161B-DT;BMS1P2;BMS1P2-AGAP9;LOC101929710;ARL14EP-DT;POU6F2-AS2;BTN2A3P;SMIM10;LOC101928268;LINC01691 |
| OFC BA11              |              70 |                   3 |                      24 |                      69 |                  70 |                   0.561776 | QPRT;LOC105373386;LOC105378913;RAB24;POLR3E;TRPV3;LOC105377975;MYBPC1;C9orf163;NOL6                      |
| Ventral Striatum      |             146 |                   7 |                      16 |                      97 |                 116 |                   0.465759 | MORN5;FSHB;MRPL2;MOSPD3;POF1B;LINC01361;KRT5;BMS1P2;BMS1P2-AGAP9;LOC105371067                            |
| dlPFC BA8/9           |              72 |                   7 |                      10 |                      58 |                  71 |                   0.524953 | PPT2;NOL6;LOC107984015;LINC00269;LINC02897;LOC105371009;LINC01990;EMC6;CARF;PPT2-EGFL8                   |

## Top Genes By Region/Subregion

### Cingulate gyrus 25

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| LOC107985795  |            0.794806 | High-confidence ML-prioritized Tier 2     |                         0.6256 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC105377644  |            0.777733 | High-confidence ML-prioritized Tier 2     |                         0.5968 |                     5 | up_in_MDD       | up_in_MDD      |
| NOC2LP2       |            0.753504 | High-confidence ML-prioritized Tier 2     |                         0.5872 |                     5 | up_in_MDD       | up_in_MDD      |
| GEMIN7        |            0.752134 | High-confidence ML-prioritized Tier 2     |                         0.5776 |                     5 | up_in_MDD       | up_in_MDD      |
| C1orf174      |            0.738035 | High-confidence ML-prioritized Tier 2     |                         0.5392 |                     5 | up_in_MDD       | up_in_MDD      |
| PCBD1         |            0.725991 | Moderate-confidence ML-prioritized Tier 2 |                         0.4976 |                     5 | up_in_MDD       | up_in_MDD      |
| EXOC3-AS1     |            0.717291 | Moderate-confidence ML-prioritized Tier 2 |                         0.408  |                     5 | up_in_MDD       | up_in_MDD      |
| LOC105372507  |            0.715616 | High-confidence ML-prioritized Tier 2     |                         0.508  |                     5 | up_in_MDD       | up_in_MDD      |
| DLGAP4-AS1    |            0.702139 | Moderate-confidence ML-prioritized Tier 2 |                         0.4704 |                     5 | down_in_MDD     | down_in_MDD    |
| MEOX2         |            0.701122 | Moderate-confidence ML-prioritized Tier 2 |                         0.4496 |                     5 | up_in_MDD       | up_in_MDD      |

### Hippocampal Formation

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| TOR2A         |            0.817642 | High-confidence ML-prioritized Tier 2     |                         0.6944 |                     5 | up_in_MDD       | up_in_MDD      |
| IL4I1         |            0.743743 | High-confidence ML-prioritized Tier 2     |                         0.5488 |                     5 | up_in_MDD       | up_in_MDD      |
| ZNF663P       |            0.740996 | Moderate-confidence ML-prioritized Tier 2 |                         0.4864 |                     5 | up_in_MDD       | up_in_MDD      |
| HLA-F-AS1     |            0.724195 | Moderate-confidence ML-prioritized Tier 2 |                         0.4256 |                     5 | up_in_MDD       | up_in_MDD      |
| TRMT2A        |            0.704332 | Moderate-confidence ML-prioritized Tier 2 |                         0.472  |                     5 | down_in_MDD     | down_in_MDD    |
| COQ4          |            0.701104 | Moderate-confidence ML-prioritized Tier 2 |                         0.4576 |                     5 | down_in_MDD     | down_in_MDD    |
| CKAP2L        |            0.699789 | Moderate-confidence ML-prioritized Tier 2 |                         0.3984 |                     5 | down_in_MDD     | down_in_MDD    |
| ARL14EP-DT    |            0.699264 | Moderate-confidence ML-prioritized Tier 2 |                         0.464  |                     5 | up_in_MDD       | up_in_MDD      |
| TMEM161B-DT   |            0.695619 | Moderate-confidence ML-prioritized Tier 2 |                         0.4736 |                     5 | up_in_MDD       | up_in_MDD      |
| ADRA2B        |            0.692397 | Moderate-confidence ML-prioritized Tier 2 |                         0.3992 |                     5 | down_in_MDD     | down_in_MDD    |

### Insula

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| TMEM161B-DT   |            0.863462 | High-confidence ML-prioritized Tier 2     |                         0.7312 |                     5 | up_in_MDD       | up_in_MDD      |
| BMS1P2        |            0.838361 | High-confidence ML-prioritized Tier 2     |                         0.7024 |                     5 | up_in_MDD       | up_in_MDD      |
| BMS1P2-AGAP9  |            0.804612 | High-confidence ML-prioritized Tier 2     |                         0.6176 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC101929710  |            0.758471 | High-confidence ML-prioritized Tier 2     |                         0.5792 |                     5 | up_in_MDD       | up_in_MDD      |
| ARL14EP-DT    |            0.752132 | High-confidence ML-prioritized Tier 2     |                         0.5488 |                     5 | up_in_MDD       | up_in_MDD      |
| POU6F2-AS2    |            0.71675  | Moderate-confidence ML-prioritized Tier 2 |                         0.4016 |                     5 | down_in_MDD     | down_in_MDD    |
| BTN2A3P       |            0.70692  | Moderate-confidence ML-prioritized Tier 2 |                         0.4968 |                     5 | up_in_MDD       | up_in_MDD      |
| SMIM10        |            0.669643 | Moderate-confidence ML-prioritized Tier 2 |                         0.384  |                     5 | up_in_MDD       | up_in_MDD      |
| LOC101928268  |            0.668309 | Moderate-confidence ML-prioritized Tier 2 |                         0.344  |                     5 | up_in_MDD       | up_in_MDD      |
| LINC01691     |            0.666752 | Moderate-confidence ML-prioritized Tier 2 |                         0.4384 |                     5 | down_in_MDD     | down_in_MDD    |

### OFC BA11

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| QPRT          |            0.789042 | High-confidence ML-prioritized Tier 2     |                         0.5456 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC105373386  |            0.746031 | High-confidence ML-prioritized Tier 2     |                         0.5456 |                     5 | down_in_MDD     | down_in_MDD    |
| LOC105378913  |            0.736479 | High-confidence ML-prioritized Tier 2     |                         0.5776 |                     5 | down_in_MDD     | down_in_MDD    |
| RAB24         |            0.71449  | Moderate-confidence ML-prioritized Tier 2 |                         0.4768 |                     5 | up_in_MDD       | up_in_MDD      |
| POLR3E        |            0.705837 | Moderate-confidence ML-prioritized Tier 2 |                         0.4496 |                     5 | down_in_MDD     | down_in_MDD    |
| TRPV3         |            0.698844 | Moderate-confidence ML-prioritized Tier 2 |                         0.4768 |                     5 | down_in_MDD     | down_in_MDD    |
| LOC105377975  |            0.692251 | Moderate-confidence ML-prioritized Tier 2 |                         0.4224 |                     5 | up_in_MDD       | up_in_MDD      |
| MYBPC1        |            0.691249 | Moderate-confidence ML-prioritized Tier 2 |                         0.3664 |                     5 | up_in_MDD       | up_in_MDD      |
| C9orf163      |            0.689395 | Moderate-confidence ML-prioritized Tier 2 |                         0.4808 |                     5 | up_in_MDD       | up_in_MDD      |
| NOL6          |            0.685334 | Moderate-confidence ML-prioritized Tier 2 |                         0.4288 |                     5 | down_in_MDD     | down_in_MDD    |

### Ventral Striatum

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| MORN5         |            0.832806 | High-confidence ML-prioritized Tier 2     |                         0.6976 |                     5 | up_in_MDD       | up_in_MDD      |
| FSHB          |            0.806363 | High-confidence ML-prioritized Tier 2     |                         0.7104 |                     5 | up_in_MDD       | up_in_MDD      |
| MRPL2         |            0.794063 | High-confidence ML-prioritized Tier 2     |                         0.6496 |                     5 | up_in_MDD       | up_in_MDD      |
| MOSPD3        |            0.783282 | High-confidence ML-prioritized Tier 2     |                         0.6208 |                     5 | up_in_MDD       | up_in_MDD      |
| POF1B         |            0.779861 | High-confidence ML-prioritized Tier 2     |                         0.6176 |                     5 | up_in_MDD       | up_in_MDD      |
| LINC01361     |            0.75403  | High-confidence ML-prioritized Tier 2     |                         0.6112 |                     5 | down_in_MDD     | down_in_MDD    |
| KRT5          |            0.734177 | High-confidence ML-prioritized Tier 2     |                         0.5728 |                     5 | down_in_MDD     | down_in_MDD    |
| BMS1P2        |            0.731951 | Moderate-confidence ML-prioritized Tier 2 |                         0.4752 |                     5 | up_in_MDD       | up_in_MDD      |
| BMS1P2-AGAP9  |            0.712512 | Moderate-confidence ML-prioritized Tier 2 |                         0.4576 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC105371067  |            0.696534 | Moderate-confidence ML-prioritized Tier 2 |                         0.5072 |                     5 | up_in_MDD       | up_in_MDD      |

### dlPFC BA8/9

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| PPT2          |            0.795594 | High-confidence ML-prioritized Tier 2     |                         0.6304 |                     5 | up_in_MDD       | up_in_MDD      |
| NOL6          |            0.79136  | High-confidence ML-prioritized Tier 2     |                         0.6288 |                     5 | down_in_MDD     | down_in_MDD    |
| LOC107984015  |            0.768869 | High-confidence ML-prioritized Tier 2     |                         0.6064 |                     5 | down_in_MDD     | down_in_MDD    |
| LINC00269     |            0.765225 | High-confidence ML-prioritized Tier 2     |                         0.544  |                     5 | up_in_MDD       | up_in_MDD      |
| LINC02897     |            0.746004 | High-confidence ML-prioritized Tier 2     |                         0.5488 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC105371009  |            0.738269 | High-confidence ML-prioritized Tier 2     |                         0.5456 |                     5 | down_in_MDD     | down_in_MDD    |
| LINC01990     |            0.733693 | High-confidence ML-prioritized Tier 2     |                         0.5184 |                     5 | up_in_MDD       | up_in_MDD      |
| EMC6          |            0.715964 | Moderate-confidence ML-prioritized Tier 2 |                         0.4616 |                     5 | up_in_MDD       | up_in_MDD      |
| CARF          |            0.689796 | Moderate-confidence ML-prioritized Tier 2 |                         0.4864 |                     5 | up_in_MDD       | up_in_MDD      |
| PPT2-EGFL8    |            0.682434 | Moderate-confidence ML-prioritized Tier 2 |                         0.448  |                     5 | up_in_MDD       | up_in_MDD      |

## Interpretation

Regions with more high/moderate genes have more internally supported ML-prioritized Tier 2 candidates. These genes are candidates for follow-up, not validated biomarkers.

## Recommended Next Step

Use these ML-prioritized Tier 2 genes for enrichment revisit, overlap with previous Tier 2 enrichment terms, AHBA/WGCNA module mapping, and pathway/module-level ML.

## Warnings

Do not claim these are validated biomarkers. Do not claim the classifier is externally validated. Use language such as ML-prioritized candidate genes, exploratory gene prioritization, internally supported, and requires validation.
