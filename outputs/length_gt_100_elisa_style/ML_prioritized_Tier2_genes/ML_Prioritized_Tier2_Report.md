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
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_pilot_Tier2\feature_importance\feature_importance_summary.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_pilot_Tier2\ML_performance_summary.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_pilot_Tier2\best_model_by_region.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_strict_v2\feature_stability\selected_genes_per_fold_all_models.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_strict_v2\feature_stability\selection_frequency_all_models.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_strict_v2\ML_strict_v2_performance_summary.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_elisa_style\ML_strict_v2\best_strict_model_by_region.csv`

Missing optional files:

None.

## Methods

The candidate pool was restricted to existing Tier 2 genes. The script extracted ML feature-importance events, recomputed conservative selected flags, summarized selection frequency, model support, sparse-model support, random-forest support, rank stability, sign consistency, LMM/ML sign agreement, and an ML priority score scaled from 0 to 1.

## Results By Region/Subregion

| region_or_subregion   |   n_Tier2_genes |   n_high_confidence |   n_moderate_confidence |   n_with_sparse_support |   n_with_v2_support |   median_ML_priority_score | top_10_gene_symbols                                                                                         |
|:----------------------|----------------:|--------------------:|------------------------:|------------------------:|--------------------:|---------------------------:|:------------------------------------------------------------------------------------------------------------|
| Cingulate gyrus 25    |             186 |                   8 |                      29 |                     123 |                 162 |                   0.557238 | LOC107985795;LOC105377644;NOC2LP2;GEMIN7;C1orf174;EXOC3-AS1;PCBD1;LOC105372507;MEOX2;CD247                  |
| Hippocampal Formation |             127 |                   3 |                      25 |                      93 |                  91 |                   0.535342 | TOR2A;IL4I1;ZNF663P;HLA-F-AS1;TRMT2A;ARL14EP-DT;COQ4;CKAP2L;ADRA2B;TMEM161B-DT                              |
| Insula                |              70 |                   6 |                      18 |                      65 |                  68 |                   0.611851 | TMEM161B-DT;BMS1P2;BMS1P2-AGAP9;ARL14EP-DT;LOC101929710;POU6F2-AS2;BTN2A3P;LOC101928268;LINC00242;LINC01691 |
| OFC BA11              |              70 |                   3 |                      28 |                      69 |                  70 |                   0.638499 | QPRT;LOC105373386;LOC105378913;RAB24;MYBPC1;POLR3E;LOC105377975;TRPV3;C9orf163;NOL6                         |
| Ventral Striatum      |             146 |                   8 |                      18 |                      97 |                 116 |                   0.531131 | MORN5;FSHB;MRPL2;MOSPD3;POF1B;LINC01361;BMS1P2;KRT5;BMS1P2-AGAP9;OR4F3                                      |
| dlPFC BA8/9           |              72 |                   7 |                      12 |                      58 |                  71 |                   0.598543 | PPT2;NOL6;LINC00269;LOC107984015;LINC02897;LINC01990;LOC105371009;EMC6;CARF;MSANTD2-AS1                     |

## Top Genes By Region/Subregion

### Cingulate gyrus 25

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| LOC107985795  |            0.831105 | High-confidence ML-prioritized Tier 2     |                       0.632267 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC105377644  |            0.823215 | High-confidence ML-prioritized Tier 2     |                       0.619467 |                     5 | up_in_MDD       | up_in_MDD      |
| NOC2LP2       |            0.800813 | High-confidence ML-prioritized Tier 2     |                       0.6072   |                     5 | up_in_MDD       | up_in_MDD      |
| GEMIN7        |            0.799979 | High-confidence ML-prioritized Tier 2     |                       0.604267 |                     5 | up_in_MDD       | up_in_MDD      |
| C1orf174      |            0.786783 | High-confidence ML-prioritized Tier 2     |                       0.5592   |                     5 | up_in_MDD       | up_in_MDD      |
| EXOC3-AS1     |            0.783291 | Moderate-confidence ML-prioritized Tier 2 |                       0.44     |                     5 | up_in_MDD       | up_in_MDD      |
| PCBD1         |            0.776746 | High-confidence ML-prioritized Tier 2     |                       0.514933 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC105372507  |            0.766872 | High-confidence ML-prioritized Tier 2     |                       0.520667 |                     5 | up_in_MDD       | up_in_MDD      |
| MEOX2         |            0.753867 | Moderate-confidence ML-prioritized Tier 2 |                       0.456267 |                     5 | up_in_MDD       | up_in_MDD      |
| CD247         |            0.752885 | High-confidence ML-prioritized Tier 2     |                       0.503733 |                     5 | up_in_MDD       | up_in_MDD      |

### Hippocampal Formation

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| TOR2A         |            0.846446 | High-confidence ML-prioritized Tier 2     |                       0.6944   |                     5 | up_in_MDD       | up_in_MDD      |
| IL4I1         |            0.793366 | High-confidence ML-prioritized Tier 2     |                       0.5688   |                     5 | up_in_MDD       | up_in_MDD      |
| ZNF663P       |            0.786581 | Moderate-confidence ML-prioritized Tier 2 |                       0.489067 |                     5 | up_in_MDD       | up_in_MDD      |
| HLA-F-AS1     |            0.78188  | Moderate-confidence ML-prioritized Tier 2 |                       0.440267 |                     5 | up_in_MDD       | up_in_MDD      |
| TRMT2A        |            0.757771 | High-confidence ML-prioritized Tier 2     |                       0.508    |                     5 | down_in_MDD     | down_in_MDD    |
| ARL14EP-DT    |            0.756321 | Moderate-confidence ML-prioritized Tier 2 |                       0.48     |                     5 | up_in_MDD       | up_in_MDD      |
| COQ4          |            0.751067 | Moderate-confidence ML-prioritized Tier 2 |                       0.460267 |                     5 | down_in_MDD     | down_in_MDD    |
| CKAP2L        |            0.7486   | Moderate-confidence ML-prioritized Tier 2 |                       0.3984   |                     5 | down_in_MDD     | down_in_MDD    |
| ADRA2B        |            0.74782  | Moderate-confidence ML-prioritized Tier 2 |                       0.4012   |                     5 | down_in_MDD     | down_in_MDD    |
| TMEM161B-DT   |            0.747577 | Moderate-confidence ML-prioritized Tier 2 |                       0.478933 |                     5 | up_in_MDD       | up_in_MDD      |

### Insula

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| TMEM161B-DT   |            0.889593 | High-confidence ML-prioritized Tier 2     |                       0.7312   |                     5 | up_in_MDD       | up_in_MDD      |
| BMS1P2        |            0.873534 | High-confidence ML-prioritized Tier 2     |                       0.7224   |                     5 | up_in_MDD       | up_in_MDD      |
| BMS1P2-AGAP9  |            0.845983 | High-confidence ML-prioritized Tier 2     |                       0.6336   |                     5 | up_in_MDD       | up_in_MDD      |
| ARL14EP-DT    |            0.801768 | High-confidence ML-prioritized Tier 2     |                       0.562133 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC101929710  |            0.799627 | High-confidence ML-prioritized Tier 2     |                       0.5792   |                     5 | up_in_MDD       | up_in_MDD      |
| POU6F2-AS2    |            0.773825 | Moderate-confidence ML-prioritized Tier 2 |                       0.418933 |                     5 | down_in_MDD     | down_in_MDD    |
| BTN2A3P       |            0.760365 | High-confidence ML-prioritized Tier 2     |                       0.508133 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC101928268  |            0.739137 | Moderate-confidence ML-prioritized Tier 2 |                       0.366667 |                     5 | up_in_MDD       | up_in_MDD      |
| LINC00242     |            0.727954 | Moderate-confidence ML-prioritized Tier 2 |                       0.349067 |                     5 | up_in_MDD       | up_in_MDD      |
| LINC01691     |            0.727472 | Moderate-confidence ML-prioritized Tier 2 |                       0.453067 |                     5 | down_in_MDD     | down_in_MDD    |

### OFC BA11

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| QPRT          |            0.83251  | High-confidence ML-prioritized Tier 2     |                       0.546933 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC105373386  |            0.787298 | High-confidence ML-prioritized Tier 2     |                       0.5496   |                     5 | down_in_MDD     | down_in_MDD    |
| LOC105378913  |            0.781435 | High-confidence ML-prioritized Tier 2     |                       0.596267 |                     5 | down_in_MDD     | down_in_MDD    |
| RAB24         |            0.764628 | Moderate-confidence ML-prioritized Tier 2 |                       0.478133 |                     5 | up_in_MDD       | up_in_MDD      |
| MYBPC1        |            0.75602  | Moderate-confidence ML-prioritized Tier 2 |                       0.383733 |                     5 | up_in_MDD       | up_in_MDD      |
| POLR3E        |            0.755025 | Moderate-confidence ML-prioritized Tier 2 |                       0.4496   |                     5 | down_in_MDD     | down_in_MDD    |
| LOC105377975  |            0.751652 | Moderate-confidence ML-prioritized Tier 2 |                       0.4384   |                     5 | up_in_MDD       | up_in_MDD      |
| TRPV3         |            0.744135 | Moderate-confidence ML-prioritized Tier 2 |                       0.483467 |                     5 | down_in_MDD     | down_in_MDD    |
| C9orf163      |            0.742143 | Moderate-confidence ML-prioritized Tier 2 |                       0.488133 |                     5 | up_in_MDD       | up_in_MDD      |
| NOL6          |            0.729095 | Moderate-confidence ML-prioritized Tier 2 |                       0.4528   |                     5 | down_in_MDD     | down_in_MDD    |

### Ventral Striatum

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| MORN5         |            0.871387 | High-confidence ML-prioritized Tier 2     |                       0.724267 |                     5 | up_in_MDD       | up_in_MDD      |
| FSHB          |            0.836799 | High-confidence ML-prioritized Tier 2     |                       0.721067 |                     5 | up_in_MDD       | up_in_MDD      |
| MRPL2         |            0.832149 | High-confidence ML-prioritized Tier 2     |                       0.658933 |                     5 | up_in_MDD       | up_in_MDD      |
| MOSPD3        |            0.824715 | High-confidence ML-prioritized Tier 2     |                       0.634133 |                     5 | up_in_MDD       | up_in_MDD      |
| POF1B         |            0.820559 | High-confidence ML-prioritized Tier 2     |                       0.630933 |                     5 | up_in_MDD       | up_in_MDD      |
| LINC01361     |            0.790042 | High-confidence ML-prioritized Tier 2     |                       0.625867 |                     5 | down_in_MDD     | down_in_MDD    |
| BMS1P2        |            0.783982 | Moderate-confidence ML-prioritized Tier 2 |                       0.489867 |                     5 | up_in_MDD       | up_in_MDD      |
| KRT5          |            0.770843 | High-confidence ML-prioritized Tier 2     |                       0.583467 |                     5 | down_in_MDD     | down_in_MDD    |
| BMS1P2-AGAP9  |            0.761357 | Moderate-confidence ML-prioritized Tier 2 |                       0.458933 |                     5 | up_in_MDD       | up_in_MDD      |
| OR4F3         |            0.749262 | Moderate-confidence ML-prioritized Tier 2 |                       0.409333 |                     5 | up_in_MDD       | up_in_MDD      |

### dlPFC BA8/9

| gene_symbol   |   ML_priority_score | ML_priority_class                         |   selection_frequency_combined |   model_support_count | lmm_direction   | ml_direction   |
|:--------------|--------------------:|:------------------------------------------|-------------------------------:|----------------------:|:----------------|:---------------|
| PPT2          |            0.831354 | High-confidence ML-prioritized Tier 2     |                       0.6344   |                     5 | up_in_MDD       | up_in_MDD      |
| NOL6          |            0.820772 | High-confidence ML-prioritized Tier 2     |                       0.631467 |                     5 | down_in_MDD     | down_in_MDD    |
| LINC00269     |            0.809094 | High-confidence ML-prioritized Tier 2     |                       0.550667 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC107984015  |            0.802391 | High-confidence ML-prioritized Tier 2     |                       0.619733 |                     5 | down_in_MDD     | down_in_MDD    |
| LINC02897     |            0.787699 | High-confidence ML-prioritized Tier 2     |                       0.551467 |                     5 | up_in_MDD       | up_in_MDD      |
| LINC01990     |            0.779863 | High-confidence ML-prioritized Tier 2     |                       0.519733 |                     5 | up_in_MDD       | up_in_MDD      |
| LOC105371009  |            0.778374 | High-confidence ML-prioritized Tier 2     |                       0.5496   |                     5 | down_in_MDD     | down_in_MDD    |
| EMC6          |            0.775586 | Moderate-confidence ML-prioritized Tier 2 |                       0.490267 |                     5 | up_in_MDD       | up_in_MDD      |
| CARF          |            0.737026 | Moderate-confidence ML-prioritized Tier 2 |                       0.491733 |                     5 | up_in_MDD       | up_in_MDD      |
| MSANTD2-AS1   |            0.735452 | Moderate-confidence ML-prioritized Tier 2 |                       0.433067 |                     5 | up_in_MDD       | up_in_MDD      |

## Interpretation

Regions with more high/moderate genes have more internally supported ML-prioritized Tier 2 candidates. These genes are candidates for follow-up, not validated biomarkers.

## Recommended Next Step

Use these ML-prioritized Tier 2 genes for enrichment revisit, overlap with previous Tier 2 enrichment terms, AHBA/WGCNA module mapping, and pathway/module-level ML.

## Warnings

Do not claim these are validated biomarkers. Do not claim the classifier is externally validated. Use language such as ML-prioritized candidate genes, exploratory gene prioritization, internally supported, and requires validation.
