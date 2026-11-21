> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Experimental whole-transcriptome ML–RRA prioritisation

## Status and boundary

This analysis is an experimental prioritisation redesign. It does not replace or alter the canonical corrected 671-row ML priority result. The old final score combined hand-specified, correlated components as `0.40F + 0.15M + 0.10S + 0.15R + 0.10D + 0.05L + 0.05V`. This branch replaces only the experimental aggregation strategy with a literature-established rank-aggregation framework. It does not delete, relabel, benchmark against, or make current changes to the frozen legacy results. The six original analysis units are Ventral Striatum, Insula, Hippocampal Formation, CG25, dlPFC BA8/9, and OFC BA11; pooled broad Frontal Cortex is excluded.

## What was done

Each normalized whole-transcriptome matrix was analysed independently. Five repeats of stratified five-fold outer cross-validation produced 25 held-out evaluations per model family. Within every outer-training split, three-fold inner cross-validation jointly selected the feature-selection method (ANOVA F or absolute standardized mean difference), the number of retained genes (25, 50, 100, 250, or 500), and compact model hyperparameters. Median imputation, scaling for linear models, feature selection, and tuning were all fit only on training data.

The five model families were L1 logistic regression, Elastic Net logistic regression, L2 logistic regression, linear SVM, and Random Forest. Fold rankings used absolute linear coefficients or Random Forest feature importance. Signs were exported for linear models as annotation but did not enter aggregation.

## Exact settings

| Region | Samples | MDD | Control | Genes |
| --- | --- | --- | --- | --- |
| Ventral Striatum | 48 | 26 | 22 | 23917 |
| Insula | 48 | 26 | 22 | 23090 |
| Hippocampal Formation | 43 | 24 | 19 | 24371 |
| CG25 | 27 | 12 | 15 | 23941 |
| dlPFC BA8/9 | 48 | 26 | 22 | 23941 |
| OFC BA11 | 47 | 25 | 22 | 23941 |

- Outer CV: five stratified folds × five repeats; seeds 2026, 2027, 2028, 2029, and 2030.
- Inner CV: three stratified folds, balanced-accuracy tuning.
- Selectors: ANOVA F and absolute standardized mean difference, fit on inner-training samples only.
- Candidate k: 25, 50, 100, 250, and 500.
- L1/L2 logistic and linear SVM: `C ∈ {0.01, 0.1, 1}`.
- Elastic Net logistic: `C ∈ {0.01, 0.1, 1}`, `l1_ratio = 0.5`.
- Random Forest: 200 trees, `max_depth ∈ {3, None}`, `min_samples_leaf = 2`.
- Median imputation for all models; standard scaling for linear models; balanced class weights.
- Python 3.12.14, NumPy 2.3.5, pandas 3.0.1, scikit-learn 1.7.2; R 4.6.1 and RobustRankAggreg 1.2.1.

## Leakage controls

The outer test samples were never used to fit imputation, scaling, supervised selection, selector choice, k choice, or model parameters. Those operations were repeated inside each inner-training fold. After tuning, the selected pipeline was refit only on the full outer-training partition and evaluated once on the held-out outer fold. Tier2, Tier4, the legacy 671 pairs, and old confidence labels were not available to fitting or aggregation.

## Robust rank aggregation

Level 1 used `RobustRankAggreg::aggregateRanks(method = "RRA")` to aggregate the 25 outer-fold rankings separately for every region and model. Level 2 applied the same canonical R implementation to the five model-level rankings per region. Held-out predictive metrics did not weight either RRA level. RRA scores quantify unusually consistent high ranking relative to the specified gene universe; they are not MDD-association p-values.

## Output summary

- Outer model fits: 750
- Fold-ranking rows: 67,907
- Model-level RRA rows: 27,916
- Consensus RRA rows: 11,792
- Canonical outputs modified: no

No high/moderate confidence threshold was invented for the new consensus. The main result is the continuous within-region `consensus_rank` and `consensus_rra_score`, accompanied by transparent per-model ranks.

## Held-out prediction results

| region | model | n_outer_fits | balanced_accuracy_mean | balanced_accuracy_sd | f1_mean | roc_auc_mean | precision_mean | recall_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CG25 | elastic_net | 25 | 0.483 | 0.099 | 0.400 | 0.464 | 0.333 | 0.553 |
| CG25 | linear_svm | 25 | 0.530 | 0.144 | 0.487 | 0.604 | 0.466 | 0.580 |
| CG25 | logistic_l1 | 25 | 0.467 | 0.110 | 0.126 | 0.422 | 0.123 | 0.133 |
| CG25 | logistic_l2 | 25 | 0.517 | 0.160 | 0.470 | 0.536 | 0.449 | 0.540 |
| CG25 | random_forest | 25 | 0.577 | 0.185 | 0.431 | 0.547 | 0.487 | 0.433 |
| Hippocampal Formation | elastic_net | 25 | 0.442 | 0.104 | 0.499 | 0.465 | 0.444 | 0.604 |
| Hippocampal Formation | linear_svm | 25 | 0.469 | 0.177 | 0.513 | 0.515 | 0.529 | 0.532 |
| Hippocampal Formation | logistic_l1 | 25 | 0.480 | 0.121 | 0.248 | 0.482 | 0.295 | 0.236 |
| Hippocampal Formation | logistic_l2 | 25 | 0.501 | 0.166 | 0.527 | 0.495 | 0.581 | 0.542 |
| Hippocampal Formation | random_forest | 25 | 0.425 | 0.138 | 0.539 | 0.440 | 0.501 | 0.614 |
| Insula | elastic_net | 25 | 0.550 | 0.117 | 0.585 | 0.566 | 0.567 | 0.631 |
| Insula | linear_svm | 25 | 0.554 | 0.142 | 0.557 | 0.574 | 0.600 | 0.547 |
| Insula | logistic_l1 | 25 | 0.549 | 0.173 | 0.492 | 0.565 | 0.573 | 0.464 |
| Insula | logistic_l2 | 25 | 0.562 | 0.161 | 0.571 | 0.586 | 0.605 | 0.563 |
| Insula | random_forest | 25 | 0.590 | 0.123 | 0.636 | 0.562 | 0.629 | 0.673 |
| OFC BA11 | elastic_net | 25 | 0.521 | 0.169 | 0.546 | 0.501 | 0.555 | 0.576 |
| OFC BA11 | linear_svm | 25 | 0.522 | 0.175 | 0.546 | 0.532 | 0.564 | 0.560 |
| OFC BA11 | logistic_l1 | 25 | 0.486 | 0.139 | 0.429 | 0.494 | 0.444 | 0.432 |
| OFC BA11 | logistic_l2 | 25 | 0.532 | 0.148 | 0.551 | 0.586 | 0.574 | 0.560 |
| OFC BA11 | random_forest | 25 | 0.559 | 0.159 | 0.578 | 0.580 | 0.564 | 0.616 |
| Ventral Striatum | elastic_net | 25 | 0.509 | 0.121 | 0.514 | 0.520 | 0.477 | 0.584 |
| Ventral Striatum | linear_svm | 25 | 0.512 | 0.148 | 0.527 | 0.516 | 0.519 | 0.557 |
| Ventral Striatum | logistic_l1 | 25 | 0.545 | 0.157 | 0.415 | 0.556 | 0.441 | 0.413 |
| Ventral Striatum | logistic_l2 | 25 | 0.496 | 0.136 | 0.529 | 0.490 | 0.538 | 0.543 |
| Ventral Striatum | random_forest | 25 | 0.436 | 0.136 | 0.515 | 0.441 | 0.482 | 0.581 |
| dlPFC BA8/9 | elastic_net | 25 | 0.523 | 0.111 | 0.535 | 0.510 | 0.525 | 0.572 |
| dlPFC BA8/9 | linear_svm | 25 | 0.579 | 0.130 | 0.604 | 0.564 | 0.622 | 0.611 |
| dlPFC BA8/9 | logistic_l1 | 25 | 0.517 | 0.165 | 0.469 | 0.520 | 0.453 | 0.499 |
| dlPFC BA8/9 | logistic_l2 | 25 | 0.631 | 0.133 | 0.657 | 0.657 | 0.671 | 0.664 |
| dlPFC BA8/9 | random_forest | 25 | 0.610 | 0.127 | 0.653 | 0.629 | 0.652 | 0.691 |

These metrics describe MDD/control classification on held-out folds. Good F1, balanced accuracy, or ROC-AUC does not by itself establish a biologically valid gene ranking.

## Top 20 consensus genes per region

### Ventral Striatum

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 254956 | MORN5 | 6.389e-22 | 5 | 1 | 1 | 1 | 1 | 1 |
| 2 | 2488 | FSHB | 6.542e-19 | 5 | 4 | 4 | 3 | 3 | 3 |
| 3 | 51069 | MRPL2 | 6.542e-19 | 5 | 2 | 2 | 2 | 2 | 4 |
| 4 | 64598 | MOSPD3 | 1.074e-17 | 5 | 3 | 3 | 7 | 4 | 2 |
| 5 | 79983 | POF1B | 2.094e-17 | 5 | 7 | 5 | 5 | 8 | 5 |
| 6 | 101927498 | LINC01361 | 1.029e-16 | 5 | 9 | 6 | 6 | 6 | 11 |
| 7 | 3852 | KRT5 | 4.852e-16 | 5 | 5 | 7 | 4 | 9 | 15 |
| 8 | 642826 | BMS1P2 | 1.582e-15 | 5 | 6 | 10 | 9 | 10 | 19 |
| 9 | 26683 | OR4F3 | 2.044e-15 | 5 | 8 | 8 | 14 | 13 | 20 |
| 10 | 5743 | PTGS2 | 2.044e-15 | 5 | 19 | 20 | 19 | 11 | 17 |
| 11 | 100505894 | TMEM161B-DT | 5.087e-15 | 5 | 16 | 12 | 13 | 7 | 24 |
| 12 | 96626 | LIMS3 | 5.087e-15 | 5 | 24 | 19 | 24 | 24 | 8 |
| 13 | 54433 | GAR1 | 7.591e-15 | 5 | 11 | 21 | 15 | 18 | 26 |
| 14 | 594840 | ZFAT-AS1 | 9.168e-15 | 5 | 27 | 13 | 12 | 16 | 7 |
| 15 | 105370512 | LOC105370512 | 1.1e-14 | 5 | 12 | 14 | 22 | 23 | 28 |
| 16 | 9683 | N4BP1 | 2.144e-14 | 5 | 10 | 15 | 27 | 32 | 32 |
| 17 | 105371067 | LOC105371067 | 2.903e-14 | 5 | 15 | 28 | 18 | 34 | 12 |
| 18 | 80705 | TSGA10 | 3.863e-14 | 5 | 18 | 17 | 32 | 36 | 10 |
| 19 | 197021 | LCTL | 6.542e-14 | 5 | 21 | 16 | 11 | 12 | 40 |
| 20 | 105372764 | LOC105372764 | 7.402e-14 | 5 | 17 | 11 | 39 | 39 | 41 |
### Insula

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 642826 | BMS1P2 | 2.438e-20 | 5 | 2 | 2 | 2 | 2 | 1 |
| 2 | 100505894 | TMEM161B-DT | 7.801e-19 | 5 | 1 | 1 | 1 | 1 | 4 |
| 3 | 1803 | DPP4 | 5.924e-18 | 5 | 6 | 6 | 4 | 4 | 3 |
| 4 | 113939912 | BMS1P2-AGAP9 | 2.496e-17 | 5 | 8 | 3 | 3 | 3 | 2 |
| 5 | 105376607 | ARL14EP-DT | 1.896e-16 | 5 | 3 | 4 | 12 | 12 | 10 |
| 6 | 107985195 | LINC01691 | 7.988e-16 | 5 | 9 | 8 | 16 | 15 | 15 |
| 7 | 101928268 | LOC101928268 | 1.44e-15 | 5 | 17 | 12 | 14 | 18 | 16 |
| 8 | 100689074 | POU6F2-AS2 | 3.111e-15 | 5 | 21 | 5 | 8 | 6 | 8 |
| 9 | 101929710 | LOC101929710 | 3.111e-15 | 5 | 4 | 7 | 7 | 7 | 21 |
| 10 | 153768 | PRELID2 | 6.066e-15 | 5 | 7 | 10 | 5 | 5 | 24 |
| 11 | 80054 | CEBPA-DT | 7.44e-15 | 5 | 16 | 15 | 20 | 19 | 25 |
| 12 | 78994 | PRR14 | 1.563e-14 | 5 | 29 | 18 | 13 | 13 | 13 |
| 13 | 54718 | BTN2A3P | 5.283e-14 | 5 | 5 | 9 | 27 | 20 | 37 |
| 14 | 128153 | SPATA17 | 6.036e-14 | 5 | 32 | 16 | 22 | 37 | 38 |
| 15 | 388569 | ZNF324B | 6.036e-14 | 5 | 15 | 11 | 36 | 38 | 20 |
| 16 | 7010 | TEK | 1.12e-13 | 5 | 43 | 36 | 24 | 23 | 32 |
| 17 | 4188 | MDFI | 1.747e-13 | 5 | 38 | 14 | 29 | 22 | 47 |
| 18 | 107986284 | LOC107986284 | 2.152e-13 | 5 | 28 | 49 | 32 | 44 | 42 |
| 19 | 401288 | LINC00242 | 3.498e-13 | 5 | 24 | 23 | 18 | 14 | 54 |
| 20 | 387990 | TOMM20L | 5.924e-13 | 5 | 10 | 25 | 53 | 60 | 59 |
### Hippocampal Formation

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 27433 | TOR2A | 1.413e-19 | 5 | 1 | 1 | 1 | 1 | 3 |
| 2 | 65083 | NOL6 | 1.817e-18 | 5 | 5 | 4 | 4 | 5 | 1 |
| 3 | 27037 | TRMT2A | 4.522e-18 | 5 | 4 | 6 | 2 | 2 | 2 |
| 4 | 259307 | IL4I1 | 9.774e-18 | 5 | 3 | 3 | 3 | 4 | 7 |
| 5 | 100505894 | TMEM161B-DT | 5.816e-17 | 5 | 2 | 2 | 5 | 3 | 10 |
| 6 | 79918 | SETD6 | 5.816e-17 | 5 | 10 | 10 | 8 | 9 | 8 |
| 7 | 285830 | HLA-F-AS1 | 3.128e-16 | 5 | 9 | 5 | 10 | 11 | 14 |
| 8 | 105371295 | LINC02137 | 1.861e-15 | 5 | 19 | 20 | 18 | 12 | 16 |
| 9 | 729652 | TEDC2-AS1 | 2.997e-15 | 5 | 7 | 11 | 22 | 14 | 13 |
| 10 | 11341 | SCRG1 | 8.345e-15 | 5 | 21 | 13 | 16 | 24 | 27 |
| 11 | 154092 | LINC01010 | 8.345e-15 | 5 | 22 | 24 | 27 | 16 | 9 |
| 12 | 151 | ADRA2B | 1.193e-14 | 5 | 15 | 26 | 9 | 18 | 29 |
| 13 | 170487 | ACTL10 | 1.413e-14 | 5 | 11 | 30 | 19 | 21 | 19 |
| 14 | 348378 | SHISAL2A | 1.413e-14 | 5 | 14 | 9 | 15 | 15 | 30 |
| 15 | 105377202 | LOC105377202 | 2.276e-14 | 5 | 6 | 33 | 17 | 27 | 18 |
| 16 | 100130934 | ZNF663P | 5.247e-14 | 5 | 20 | 7 | 23 | 22 | 39 |
| 17 | 80758 | PRR7 | 5.955e-14 | 5 | 40 | 22 | 11 | 10 | 6 |
| 18 | 101928283 | LOC101928283 | 1.334e-13 | 5 | 46 | 39 | 25 | 47 | 36 |
| 19 | 51117 | COQ4 | 1.643e-13 | 5 | 8 | 15 | 13 | 19 | 49 |
| 20 | 101927931 | LINC01480 | 2.007e-13 | 5 | 43 | 51 | 40 | 25 | 21 |
### CG25

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 116349 | EXOC3-AS1 | 2.034e-20 | 5 | 1 | 2 | 1 | 1 | 1 |
| 2 | 105377644 | LOC105377644 | 6.51e-19 | 5 | 3 | 4 | 2 | 2 | 2 |
| 3 | 107985795 | LOC107985795 | 1.987e-18 | 5 | 2 | 3 | 4 | 5 | 5 |
| 4 | 79760 | GEMIN7 | 1.987e-18 | 5 | 5 | 5 | 3 | 3 | 4 |
| 5 | 401010 | NOC2LP2 | 4.943e-18 | 5 | 4 | 1 | 5 | 4 | 6 |
| 6 | 339448 | C1orf174 | 2.083e-17 | 5 | 8 | 6 | 6 | 6 | 3 |
| 7 | 64773 | PCED1A | 1.024e-16 | 5 | 10 | 8 | 9 | 9 | 11 |
| 8 | 100037417 | DDTL | 2.36e-16 | 5 | 6 | 13 | 10 | 13 | 8 |
| 9 | 115290 | FBXO17 | 1.201e-15 | 5 | 13 | 18 | 14 | 12 | 17 |
| 10 | 26071 | RTL8A | 2.034e-15 | 5 | 15 | 15 | 20 | 17 | 20 |
| 11 | 645 | BLVRB | 2.596e-15 | 5 | 21 | 7 | 12 | 7 | 16 |
| 12 | 2671 | GFER | 3.276e-15 | 5 | 22 | 10 | 7 | 8 | 7 |
| 13 | 10974 | ADIRF | 7.553e-15 | 5 | 12 | 19 | 26 | 26 | 26 |
| 14 | 581 | BAX | 2.488e-14 | 5 | 17 | 21 | 23 | 33 | 27 |
| 15 | 5092 | PCBD1 | 3.844e-14 | 5 | 32 | 31 | 28 | 31 | 36 |
| 16 | 113230 | MISP3 | 2.417e-13 | 5 | 52 | 16 | 11 | 16 | 21 |
| 17 | 374659 | HDDC3 | 2.417e-13 | 5 | 40 | 52 | 38 | 38 | 30 |
| 18 | 79094 | CHAC1 | 3.199e-13 | 5 | 14 | 39 | 16 | 20 | 55 |
| 19 | 79812 | MMRN2 | 8.583e-13 | 5 | 67 | 64 | 25 | 21 | 31 |
| 20 | 109616972 | SNORA68B | 1.23e-12 | 5 | 72 | 54 | 62 | 50 | 69 |
### dlPFC BA8/9

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 83460 | EMC6 | 6.51e-19 | 5 | 4 | 1 | 2 | 1 | 1 |
| 2 | 100996279 | LINC00269 | 4.943e-18 | 5 | 1 | 3 | 5 | 6 | 4 |
| 3 | 9374 | PPT2 | 4.943e-18 | 5 | 2 | 2 | 4 | 3 | 6 |
| 4 | 107984015 | LOC107984015 | 3.754e-17 | 5 | 3 | 4 | 9 | 7 | 7 |
| 5 | 388759 | LINC02897 | 6.357e-17 | 5 | 6 | 5 | 8 | 10 | 10 |
| 6 | 105373490 | LOC105373490 | 1.024e-16 | 5 | 7 | 8 | 11 | 8 | 3 |
| 7 | 414062 | CCL3L3 | 2.034e-15 | 5 | 18 | 13 | 10 | 12 | 20 |
| 8 | 100532746 | PPT2-EGFL8 | 3.276e-15 | 5 | 9 | 10 | 14 | 13 | 22 |
| 9 | 101928419 | LOC101928419 | 4.092e-15 | 5 | 23 | 6 | 13 | 14 | 14 |
| 10 | 107983974 | LOC107983974 | 1.304e-14 | 5 | 15 | 11 | 19 | 19 | 29 |
| 11 | 57529 | RTL9 | 1.545e-14 | 5 | 24 | 27 | 23 | 30 | 17 |
| 12 | 101927659 | SMAD1-AS2 | 2.488e-14 | 5 | 14 | 16 | 30 | 33 | 28 |
| 13 | 121274 | ZNF641 | 2.488e-14 | 5 | 11 | 12 | 22 | 23 | 33 |
| 14 | 147 | ADRA1B | 4.408e-14 | 5 | 37 | 7 | 29 | 31 | 26 |
| 15 | 105369340 | PPP1R14B-AS1 | 7.365e-14 | 5 | 41 | 25 | 27 | 41 | 34 |
| 16 | 105371009 | LOC105371009 | 9.345e-14 | 5 | 5 | 9 | 12 | 11 | 43 |
| 17 | 105379086 | LOC105379086 | 1.309e-13 | 5 | 46 | 26 | 6 | 9 | 8 |
| 18 | 57061 | HYMAI | 1.458e-13 | 5 | 47 | 14 | 3 | 4 | 2 |
| 19 | 112268140 | LOC112268140 | 1.987e-13 | 5 | 31 | 19 | 35 | 50 | 36 |
| 20 | 283738 | NTRK3-AS1 | 3.199e-13 | 5 | 55 | 32 | 20 | 21 | 21 |
### OFC BA11

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 65083 | NOL6 | 6.51e-19 | 5 | 1 | 1 | 4 | 1 | 1 |
| 2 | 5468 | PPARG | 4.943e-18 | 5 | 2 | 2 | 2 | 5 | 6 |
| 3 | 101929549 | GAU1 | 1.024e-16 | 5 | 8 | 6 | 11 | 7 | 3 |
| 4 | 4604 | MYBPC1 | 9.026e-16 | 5 | 9 | 4 | 16 | 17 | 9 |
| 5 | 105378913 | LOC105378913 | 1.201e-15 | 5 | 3 | 3 | 5 | 2 | 18 |
| 6 | 105373386 | LOC105373386 | 1.574e-15 | 5 | 14 | 18 | 12 | 4 | 19 |
| 7 | 2273 | FHL1 | 2.596e-15 | 5 | 4 | 21 | 9 | 6 | 15 |
| 8 | 51296 | SLC15A3 | 4.092e-15 | 5 | 11 | 19 | 23 | 8 | 2 |
| 9 | 100873982 | ABCC5-AS1 | 1.094e-14 | 5 | 13 | 11 | 28 | 12 | 24 |
| 10 | 23475 | QPRT | 1.094e-14 | 5 | 10 | 5 | 1 | 15 | 28 |
| 11 | 53917 | RAB24 | 1.304e-14 | 5 | 29 | 27 | 25 | 23 | 29 |
| 12 | 389023 | DPP10-AS1 | 1.545e-14 | 5 | 6 | 15 | 20 | 22 | 30 |
| 13 | 706 | TSPO | 1.545e-14 | 5 | 26 | 16 | 30 | 19 | 16 |
| 14 | 4814 | NINJ1 | 2.133e-14 | 5 | 7 | 12 | 8 | 11 | 32 |
| 15 | 105371659 | LOC105371659 | 3.844e-14 | 5 | 32 | 7 | 6 | 9 | 36 |
| 16 | 105377975 | LOC105377975 | 5.037e-14 | 5 | 22 | 25 | 13 | 21 | 38 |
| 17 | 6948 | TCN2 | 5.037e-14 | 5 | 38 | 30 | 38 | 36 | 11 |
| 18 | 55620 | STAP2 | 1.458e-13 | 5 | 24 | 35 | 47 | 29 | 45 |
| 19 | 105370472 | LOC105370472 | 1.62e-13 | 5 | 23 | 48 | 29 | 44 | 31 |
| 20 | 101927366 | LOC101927366 | 1.987e-13 | 5 | 15 | 50 | 17 | 33 | 13 |

## Interpretation cautions

This is internal resampling-based prioritisation, not external validation and not biomarker validation. Small regional sample sizes limit performance estimates and ranking stability. Genes missing from sparse fold/model rankings are handled by the RRA partial-list procedure rather than being assigned fabricated coefficients. Model performance is reported descriptively and is never used to up-weight a model's gene evidence.

The same cohort is reused across correlated outer resamples and model families, so the rankings are not independent. There is no external validation, empirical label-permutation calibration, or current benchmark against the handmade legacy score, Stability Selection, or simpler baselines.

## Next step

The next phase is a separate benchmark against the old handmade score, Stability Selection, simpler baseline rankings, and label-permutation/null analyses. None of those comparisons was performed in this implementation phase.
