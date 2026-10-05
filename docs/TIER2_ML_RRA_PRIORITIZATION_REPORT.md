# Primary Tier2 → ML → RRA gene-prioritisation report

## A. Scientific objective

The gene-wise diagnosis models first define a prespecified exploratory Tier2 candidate universe. The primary ML question is then: among those candidates, which genes receive the most stable multivariate support across five model families and repeated resampling? The result is described as **ML-prioritized exploratory Tier2 candidates**, not validated biomarkers or diagnostic genes.

The canonical upstream source contains 671 region–gene pairs and was generated before historical ML scoring. Only region, gene ID, and gene symbol were imported. Historical ML scores, confidence labels, selection frequencies, Tier4, WGCNA, AHBA, and all-gene RRA ranks did not enter fitting or RRA.

## B. Difference from the frozen all-gene sensitivity analysis

The frozen all-gene analysis asks whether ML can discover useful genes from the full transcriptome using strict train-only feature selection. This primary analysis asks which genes are most consistently prioritized **within the already screened Tier2 universe**. Accordingly, it has no ANOVA/SMD selector choice and no nested top-k grid. The all-gene analysis remains frozen as a separate sensitivity analysis.

## C. Important limitation

Tier2 membership was defined from the complete cohort using Diagnosis. The candidate universe is therefore not independent of Diagnosis, and outer-fold classifier metrics may be optimistic. They are descriptive/internal sanity checks, not unbiased diagnostic-validation estimates. The primary endpoint is within-study gene prioritisation. There is no external validation, and no gene is a validated biomarker.

## D. Exact ML procedure

| Region | Samples | MDD | Control | Tier2 genes |
| --- | --- | --- | --- | --- |
| Ventral Striatum | 48 | 26 | 22 | 146 |
| Insula | 48 | 26 | 22 | 70 |
| Hippocampal Formation | 43 | 24 | 19 | 127 |
| CG25 | 27 | 12 | 15 | 186 |
| dlPFC BA8/9 | 48 | 26 | 22 | 72 |
| OFC BA11 | 47 | 25 | 22 | 70 |

- Outer resampling: five stratified folds × five repeats; seeds 2026–2030.
- Inner tuning: three stratified folds using balanced accuracy.
- Preprocessing: training-only median imputation; training-only standardization for linear models.
- Candidate universe: all region-specific Tier2 genes, fixed throughout each regional analysis.
- L1/L2 logistic and linear SVM: `C ∈ {0.01, 0.1, 1}`.
- Elastic Net logistic: `C ∈ {0.01, 0.1, 1}`, `l1_ratio = 0.5`.
- Random Forest: 200 trees, `max_depth ∈ {3, None}`, `min_samples_leaf = 2`.
- L1/Elastic Net rank nonzero absolute coefficients; L2/SVM rank all Tier2 genes; Random Forest ranks nonzero feature importance and is nondirectional.
- Empty sparse-model lists are retained as empty evidence.
- A most-frequent dummy classifier is reported only as a descriptive reference and contributes no gene evidence.
- Level 1 RRA aggregates 25 outer-fit lists per region/model. Level 2 aggregates the five model-level lists per region.
- R 4.6.1; `RobustRankAggreg` 1.2.1; `aggregateRanks(method = "RRA")`.
- Performance and coefficient sign never weight RRA.

## E. Candidate accounting

| region | total_tier2_candidates | candidates_with_any_fold_ranking_evidence | candidates_absent_from_every_fold_ranking | candidates_in_consensus_rra |
| --- | --- | --- | --- | --- |
| Ventral Striatum | 146 | 146 | 0 | 146 |
| Insula | 70 | 70 | 0 | 70 |
| Hippocampal Formation | 127 | 127 | 0 | 127 |
| CG25 | 186 | 186 | 0 | 186 |
| dlPFC BA8/9 | 72 | 72 | 0 | 72 |
| OFC BA11 | 70 | 70 | 0 | 70 |

Tier2 membership and ML ranking support are separate concepts. A Tier2 gene remains a candidate even if a sparse model never assigns it nonzero evidence.

## Descriptive internal classifier performance

| region | model | n_outer_fits | contributes_gene_ranking | balanced_accuracy_mean | balanced_accuracy_sd | f1_mean | roc_auc_mean | precision_mean | recall_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CG25 | dummy_most_frequent | 25 | False | 0.500 | 0.000 | 0.000 | 0.500 | 0.000 | 0.000 |
| CG25 | elastic_net | 25 | True | 0.803 | 0.136 | 0.739 | 0.920 | 0.907 | 0.660 |
| CG25 | linear_svm | 25 | True | 0.937 | 0.106 | 0.917 | 1.000 | 1.000 | 0.873 |
| CG25 | logistic_l1 | 25 | True | 0.697 | 0.161 | 0.620 | 0.869 | 0.643 | 0.660 |
| CG25 | logistic_l2 | 25 | True | 0.903 | 0.107 | 0.880 | 0.993 | 0.950 | 0.860 |
| CG25 | random_forest | 25 | True | 0.793 | 0.150 | 0.687 | 0.947 | 0.827 | 0.640 |
| Hippocampal Formation | dummy_most_frequent | 25 | False | 0.500 | 0.000 | 0.716 | 0.500 | 0.558 | 1.000 |
| Hippocampal Formation | elastic_net | 25 | True | 0.855 | 0.122 | 0.858 | 0.935 | 0.895 | 0.840 |
| Hippocampal Formation | linear_svm | 25 | True | 0.912 | 0.106 | 0.915 | 0.974 | 0.954 | 0.890 |
| Hippocampal Formation | logistic_l1 | 25 | True | 0.824 | 0.107 | 0.816 | 0.894 | 0.900 | 0.778 |
| Hippocampal Formation | logistic_l2 | 25 | True | 0.892 | 0.117 | 0.888 | 0.965 | 0.952 | 0.850 |
| Hippocampal Formation | random_forest | 25 | True | 0.761 | 0.176 | 0.784 | 0.878 | 0.777 | 0.812 |
| Insula | dummy_most_frequent | 25 | False | 0.500 | 0.000 | 0.702 | 0.500 | 0.542 | 1.000 |
| Insula | elastic_net | 25 | True | 0.938 | 0.089 | 0.940 | 0.991 | 0.951 | 0.939 |
| Insula | linear_svm | 25 | True | 0.958 | 0.062 | 0.961 | 0.992 | 0.965 | 0.963 |
| Insula | logistic_l1 | 25 | True | 0.854 | 0.081 | 0.845 | 0.949 | 0.925 | 0.800 |
| Insula | logistic_l2 | 25 | True | 0.954 | 0.062 | 0.957 | 0.977 | 0.965 | 0.955 |
| Insula | random_forest | 25 | True | 0.854 | 0.126 | 0.862 | 0.937 | 0.875 | 0.871 |
| OFC BA11 | dummy_most_frequent | 25 | False | 0.500 | 0.000 | 0.695 | 0.500 | 0.533 | 1.000 |
| OFC BA11 | elastic_net | 25 | True | 0.852 | 0.135 | 0.862 | 0.940 | 0.864 | 0.880 |
| OFC BA11 | linear_svm | 25 | True | 0.926 | 0.077 | 0.934 | 0.982 | 0.926 | 0.952 |
| OFC BA11 | logistic_l1 | 25 | True | 0.738 | 0.122 | 0.737 | 0.862 | 0.800 | 0.736 |
| OFC BA11 | logistic_l2 | 25 | True | 0.931 | 0.078 | 0.938 | 0.976 | 0.933 | 0.952 |
| OFC BA11 | random_forest | 25 | True | 0.877 | 0.115 | 0.883 | 0.925 | 0.895 | 0.880 |
| Ventral Striatum | dummy_most_frequent | 25 | False | 0.500 | 0.000 | 0.702 | 0.500 | 0.542 | 1.000 |
| Ventral Striatum | elastic_net | 25 | True | 0.941 | 0.100 | 0.948 | 0.954 | 0.937 | 0.963 |
| Ventral Striatum | linear_svm | 25 | True | 0.934 | 0.094 | 0.937 | 0.950 | 0.955 | 0.925 |
| Ventral Striatum | logistic_l1 | 25 | True | 0.915 | 0.108 | 0.925 | 0.943 | 0.926 | 0.929 |
| Ventral Striatum | logistic_l2 | 25 | True | 0.923 | 0.086 | 0.923 | 0.937 | 0.963 | 0.895 |
| Ventral Striatum | random_forest | 25 | True | 0.746 | 0.165 | 0.774 | 0.874 | 0.771 | 0.803 |
| dlPFC BA8/9 | dummy_most_frequent | 25 | False | 0.500 | 0.000 | 0.702 | 0.500 | 0.542 | 1.000 |
| dlPFC BA8/9 | elastic_net | 25 | True | 0.899 | 0.099 | 0.911 | 0.943 | 0.895 | 0.941 |
| dlPFC BA8/9 | linear_svm | 25 | True | 0.927 | 0.075 | 0.932 | 0.964 | 0.941 | 0.932 |
| dlPFC BA8/9 | logistic_l1 | 25 | True | 0.828 | 0.136 | 0.855 | 0.892 | 0.835 | 0.893 |
| dlPFC BA8/9 | logistic_l2 | 25 | True | 0.889 | 0.112 | 0.895 | 0.951 | 0.907 | 0.895 |
| dlPFC BA8/9 | random_forest | 25 | True | 0.810 | 0.116 | 0.829 | 0.882 | 0.820 | 0.849 |

Because Tier2 was selected from the full cohort using Diagnosis, these values are optimistic/internal and must not be interpreted as unbiased diagnostic validation.

## Top 20 Tier2-RRA genes by region

### Ventral Striatum

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 254956 | MORN5 | 2.412e-09 | 5 | 1 | 1 | 1 | 2 | 1 |
| 2 | 79983 | POF1B | 2.355e-07 | 5 | 3 | 3 | 5 | 4 | 4 |
| 3 | 2488 | FSHB | 5.861e-07 | 5 | 6 | 5 | 4 | 6 | 5 |
| 4 | 51069 | MRPL2 | 8.707e-07 | 5 | 2 | 2 | 2 | 1 | 19 |
| 5 | 64598 | MOSPD3 | 2.47e-06 | 5 | 4 | 4 | 3 | 3 | 8 |
| 6 | 3852 | KRT5 | 7.537e-06 | 5 | 9 | 8 | 9 | 9 | 10 |
| 7 | 107985657 | LOC107985657 | 1.214e-05 | 5 | 7 | 6 | 7 | 5 | 11 |
| 8 | 594840 | ZFAT-AS1 | 2.798e-05 | 5 | 13 | 11 | 6 | 7 | 3 |
| 9 | 101927498 | LINC01361 | 7.903e-05 | 5 | 8 | 10 | 8 | 8 | 16 |
| 10 | 642826 | BMS1P2 | 0.000107 | 5 | 5 | 7 | 13 | 17 | 7 |
| 11 | 101928161 | NCKAP5-AS2 | 0.0003884 | 5 | 18 | 15 | 12 | 11 | 22 |
| 12 | 96626 | LIMS3 | 0.0004851 | 5 | 16 | 16 | 17 | 18 | 23 |
| 13 | 26683 | OR4F3 | 0.001081 | 5 | 10 | 13 | 19 | 27 | 21 |
| 14 | 105372764 | LOC105372764 | 0.001297 | 5 | 28 | 19 | 18 | 21 | 26 |
| 15 | 113939912 | BMS1P2-AGAP9 | 0.00146 | 5 | 37 | 9 | 10 | 12 | 13 |
| 16 | 105371067 | LOC105371067 | 0.001546 | 5 | 15 | 14 | 11 | 10 | 29 |
| 17 | 125 | ADH1B | 0.005972 | 5 | 38 | 21 | 15 | 15 | 31 |
| 18 | 151 | ADRA2B | 0.008732 | 5 | 41 | 28 | 39 | 40 | 36 |
| 19 | 6103 | RPGR | 0.01134 | 5 | 59 | 22 | 14 | 20 | 9 |
| 20 | 441509 | GLRA4 | 0.01346 | 5 | 11 | 12 | 23 | 14 | 114 |
### Insula

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 100505894 | TMEM161B-DT | 7.229e-07 | 5 | 1 | 1 | 1 | 3 | 2 |
| 2 | 101929710 | LOC101929710 | 2.313e-05 | 5 | 2 | 2 | 3 | 1 | 6 |
| 3 | 642826 | BMS1P2 | 5e-05 | 5 | 6 | 5 | 2 | 7 | 1 |
| 4 | 105376607 | ARL14EP-DT | 9.748e-05 | 5 | 3 | 3 | 5 | 5 | 8 |
| 5 | 54718 | BTN2A3P | 0.001105 | 5 | 8 | 10 | 10 | 13 | 12 |
| 6 | 105379042 | LOC105379042 | 0.002259 | 5 | 15 | 7 | 9 | 4 | 14 |
| 7 | 26063 | DECR2 | 0.0023 | 5 | 4 | 4 | 7 | 2 | 52 |
| 8 | 107985195 | LINC01691 | 0.003875 | 5 | 5 | 8 | 6 | 8 | 36 |
| 9 | 101929741 | METTL14-DT | 0.004224 | 5 | 10 | 14 | 11 | 17 | 13 |
| 10 | 113939912 | BMS1P2-AGAP9 | 0.004224 | 5 | 17 | 12 | 4 | 14 | 3 |
| 11 | 100689074 | POU6F2-AS2 | 0.01215 | 5 | 12 | 21 | 12 | 19 | 16 |
| 12 | 80054 | CEBPA-DT | 0.02369 | 5 | 18 | 16 | 24 | 20 | 18 |
| 13 | 101928268 | LOC101928268 | 0.03535 | 5 | 14 | 25 | 23 | 26 | 20 |
| 14 | 80036 | TRPM3 | 0.03535 | 5 | 13 | 9 | 19 | 10 | 26 |
| 15 | 116143 | DNAAF10 | 0.04269 | 5 | 16 | 6 | 8 | 6 | 27 |
| 16 | 171177 | RHOV | 0.04368 | 5 | 11 | 13 | 15 | 11 | 48 |
| 17 | 388569 | ZNF324B | 0.04368 | 5 | 9 | 15 | 13 | 12 | 29 |
| 18 | 401288 | LINC00242 | 0.05576 | 5 | 7 | 11 | 16 | 9 | 55 |
| 19 | 10090 | UST | 0.06102 | 5 | 20 | 29 | 28 | 29 | 22 |
| 20 | 102724078 | LOC102724078 | 0.06102 | 5 | 23 | 26 | 29 | 21 | 28 |
### Hippocampal Formation

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 27433 | TOR2A | 9.55e-08 | 5 | 1 | 1 | 1 | 1 | 7 |
| 2 | 100505894 | TMEM161B-DT | 1.177e-06 | 5 | 3 | 6 | 3 | 4 | 3 |
| 3 | 105376607 | ARL14EP-DT | 8.936e-06 | 5 | 4 | 5 | 6 | 6 | 9 |
| 4 | 105370868 | DIS3L-AS1 | 5.619e-05 | 5 | 11 | 13 | 13 | 10 | 8 |
| 5 | 259307 | IL4I1 | 5.817e-05 | 5 | 2 | 2 | 5 | 3 | 22 |
| 6 | 79918 | SETD6 | 8.139e-05 | 5 | 12 | 9 | 7 | 14 | 5 |
| 7 | 151 | ADRA2B | 0.0001198 | 5 | 6 | 4 | 4 | 5 | 34 |
| 8 | 101929077 | LINC02334 | 0.0001587 | 5 | 10 | 16 | 11 | 15 | 6 |
| 9 | 100130934 | ZNF663P | 0.000286 | 5 | 15 | 12 | 18 | 16 | 17 |
| 10 | 27037 | TRMT2A | 0.000286 | 5 | 8 | 10 | 17 | 18 | 4 |
| 11 | 150468 | CKAP2L | 0.0003747 | 5 | 19 | 17 | 14 | 9 | 10 |
| 12 | 51117 | COQ4 | 0.0005948 | 5 | 9 | 3 | 2 | 2 | 29 |
| 13 | 126208 | ZNF787 | 0.0006181 | 5 | 7 | 15 | 19 | 17 | 21 |
| 14 | 107987272 | LOC107987272 | 0.001478 | 5 | 18 | 7 | 9 | 13 | 25 |
| 15 | 79713 | IGFLR1 | 0.00252 | 5 | 13 | 11 | 8 | 7 | 117 |
| 16 | 285830 | HLA-F-AS1 | 0.005663 | 5 | 5 | 8 | 16 | 11 | 33 |
| 17 | 3430 | IFI35 | 0.005923 | 5 | 14 | 21 | 33 | 27 | 15 |
| 18 | 80758 | PRR7 | 0.006876 | 5 | 16 | 30 | 23 | 34 | 24 |
| 19 | 154092 | LINC01010 | 0.007949 | 5 | 20 | 24 | 35 | 35 | 18 |
| 20 | 54496 | PRMT7 | 0.007949 | 5 | 21 | 14 | 20 | 21 | 35 |
### CG25

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 107985795 | LOC107985795 | 3.775e-07 | 5 | 1 | 2 | 2 | 7 | 6 |
| 2 | 79760 | GEMIN7 | 3.617e-06 | 5 | 5 | 6 | 7 | 11 | 10 |
| 3 | 105372507 | LOC105372507 | 5.255e-06 | 5 | 3 | 3 | 4 | 3 | 52 |
| 4 | 401010 | NOC2LP2 | 5.255e-06 | 5 | 4 | 4 | 3 | 2 | 14 |
| 5 | 919 | CD247 | 8.339e-06 | 5 | 9 | 10 | 13 | 6 | 2 |
| 6 | 105377644 | LOC105377644 | 1.208e-05 | 5 | 2 | 1 | 5 | 14 | 8 |
| 7 | 100861412 | FSBP | 4.244e-05 | 5 | 12 | 12 | 8 | 16 | 18 |
| 8 | 26071 | RTL8A | 4.864e-05 | 5 | 7 | 5 | 1 | 5 | 128 |
| 9 | 23783 | ANKRD62P1-PARP4P3 | 0.0001999 | 5 | 8 | 7 | 10 | 1 | 71 |
| 10 | 5092 | PCBD1 | 0.0004108 | 5 | 10 | 8 | 11 | 12 | 41 |
| 11 | 101926987 | DLGAP4-AS1 | 0.0004607 | 5 | 15 | 16 | 29 | 27 | 5 |
| 12 | 64773 | PCED1A | 0.0004607 | 5 | 23 | 18 | 24 | 29 | 11 |
| 13 | 4223 | MEOX2 | 0.0007541 | 5 | 14 | 13 | 9 | 9 | 44 |
| 14 | 27098 | CLUL1 | 0.001617 | 5 | 17 | 14 | 12 | 4 | 129 |
| 15 | 101927947 | LOC101927947 | 0.00178 | 5 | 32 | 38 | 37 | 33 | 34 |
| 16 | 101927655 | ZASP | 0.002026 | 5 | 31 | 39 | 36 | 32 | 19 |
| 17 | 115290 | FBXO17 | 0.002026 | 5 | 11 | 24 | 20 | 22 | 39 |
| 18 | 339448 | C1orf174 | 0.002026 | 5 | 6 | 15 | 19 | 39 | 7 |
| 19 | 101927051 | PDXP-DT | 0.002602 | 5 | 20 | 17 | 26 | 41 | 12 |
| 20 | 116349 | EXOC3-AS1 | 0.004144 | 5 | 18 | 23 | 35 | 45 | 4 |
### dlPFC BA8/9

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 100996279 | LINC00269 | 8.075e-06 | 5 | 5 | 4 | 3 | 5 | 2 |
| 2 | 65083 | NOL6 | 8.075e-06 | 5 | 2 | 1 | 1 | 2 | 5 |
| 3 | 9374 | PPT2 | 1.455e-05 | 5 | 1 | 2 | 2 | 1 | 7 |
| 4 | 101929607 | LINC01990 | 2.009e-05 | 5 | 3 | 3 | 6 | 4 | 6 |
| 5 | 107984015 | LOC107984015 | 8.468e-05 | 5 | 4 | 6 | 5 | 8 | 3 |
| 6 | 388759 | LINC02897 | 0.0004162 | 5 | 8 | 5 | 4 | 11 | 8 |
| 7 | 83460 | EMC6 | 0.00139 | 5 | 10 | 10 | 13 | 14 | 4 |
| 8 | 100287015 | MCPH1-DT | 0.003472 | 5 | 6 | 7 | 8 | 7 | 28 |
| 9 | 100532746 | PPT2-EGFL8 | 0.003669 | 5 | 17 | 15 | 10 | 15 | 14 |
| 10 | 105371009 | LOC105371009 | 0.01672 | 5 | 7 | 11 | 12 | 9 | 40 |
| 11 | 202781 | PAXIP1-DT | 0.01672 | 5 | 12 | 12 | 11 | 3 | 59 |
| 12 | 100507283 | MSANTD2-AS1 | 0.02273 | 5 | 11 | 8 | 9 | 13 | 29 |
| 13 | 107984438 | LOC107984438 | 0.03708 | 5 | 14 | 18 | 27 | 20 | 22 |
| 14 | 3178 | HNRNPA1 | 0.03708 | 5 | 13 | 21 | 20 | 27 | 13 |
| 15 | 79800 | CARF | 0.03925 | 5 | 15 | 9 | 15 | 6 | 70 |
| 16 | 105369340 | PPP1R14B-AS1 | 0.05013 | 5 | 9 | 14 | 16 | 12 | 38 |
| 17 | 147 | ADRA1B | 0.1174 | 5 | 34 | 23 | 34 | 28 | 20 |
| 18 | 101928942 | PRKG2-AS1 | 0.1792 | 5 | 21 | 30 | 35 | 24 | 37 |
| 19 | 105379086 | LOC105379086 | 0.1792 | 5 | 26 | 37 | 26 | 30 | 27 |
| 20 | 101927003 | PKIA-AS1 | 0.2048 | 5 | 16 | 25 | 33 | 38 | 11 |
### OFC BA11

| consensus_rank | gene_id | gene_symbol | consensus_rra_score | model_appearance_count | logistic_l1_rank | elastic_net_rank | logistic_l2_rank | linear_svm_rank | random_forest_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 105373386 | LOC105373386 | 8.145e-05 | 5 | 3 | 2 | 2 | 1 | 39 |
| 2 | 105378913 | LOC105378913 | 8.145e-05 | 5 | 1 | 3 | 1 | 3 | 8 |
| 3 | 53917 | RAB24 | 0.0006136 | 5 | 4 | 5 | 4 | 4 | 18 |
| 4 | 162514 | TRPV3 | 0.001105 | 5 | 12 | 8 | 13 | 8 | 7 |
| 5 | 158055 | C9orf163 | 0.001257 | 5 | 6 | 6 | 3 | 2 | 17 |
| 6 | 23475 | QPRT | 0.001257 | 5 | 2 | 1 | 6 | 5 | 16 |
| 7 | 79571 | GCC1 | 0.0023 | 5 | 5 | 7 | 5 | 7 | 68 |
| 8 | 105377975 | LOC105377975 | 0.003119 | 5 | 8 | 14 | 14 | 16 | 11 |
| 9 | 55620 | STAP2 | 0.009222 | 5 | 9 | 10 | 9 | 6 | 63 |
| 10 | 92283 | ZNF461 | 0.01863 | 5 | 10 | 4 | 12 | 10 | 69 |
| 11 | 101929549 | GAU1 | 0.01915 | 5 | 11 | 17 | 15 | 23 | 2 |
| 12 | 388574 | RPL23AP87 | 0.05576 | 5 | 16 | 13 | 11 | 12 | 58 |
| 13 | 100505894 | TMEM161B-DT | 0.06102 | 5 | 18 | 9 | 26 | 11 | 29 |
| 14 | 100873982 | ABCC5-AS1 | 0.06102 | 5 | 22 | 29 | 17 | 21 | 10 |
| 15 | 100506380 | RNF32-DT | 0.07229 | 5 | 19 | 19 | 23 | 15 | 30 |
| 16 | 55718 | POLR3E | 0.08682 | 5 | 7 | 11 | 18 | 14 | 64 |
| 17 | 100507283 | MSANTD2-AS1 | 0.09982 | 5 | 21 | 28 | 22 | 18 | 32 |
| 18 | 105371729 | LOC105371729 | 0.1062 | 5 | 13 | 12 | 19 | 9 | 38 |
| 19 | 4814 | NINJ1 | 0.1062 | 5 | 15 | 16 | 8 | 19 | 60 |
| 20 | 100129482 | ZNF37BP | 0.1164 | 5 | 17 | 18 | 33 | 20 | 28 |

RRA scores are aggregation/ranking scores, not MDD-association p-values.

## F. Model agreement and recurrence

Model-level ranks and outer-fit recurrence are provided as transparent stability diagnostics. Agreement reflects repeated analyses of the same cohort and is not independent replication.

## Tier2-RRA versus frozen all-gene RRA

| region | n_tier2_candidates | n_present_in_allgene_rra | proportion_present_in_allgene_rra | spearman_rho_common_genes | spearman_p_value_descriptive | median_allgene_rank_percentile | top10_overlap_count | top10_overlap_proportion_of_tier2_top | top10_jaccard | top25_overlap_count | top25_overlap_proportion_of_tier2_top | top25_jaccard | top50_overlap_count | top50_overlap_proportion_of_tier2_top | top50_jaccard |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CG25 | 186 | 169 | 0.909 | 0.685 | 1.028017417528295e-24 | 0.116 | 5 | 0.5 | 0.3333333333333333 | 14 | 0.56 | 0.3888888888888889 | 23 | 0.46 | 0.2987012987012987 |
| Hippocampal Formation | 127 | 119 | 0.937 | 0.575 | 7.683000029332312e-12 | 0.081 | 5 | 0.5 | 0.3333333333333333 | 13 | 0.52 | 0.3513513513513513 | 27 | 0.54 | 0.3698630136986301 |
| Insula | 70 | 70 | 1.000 | 0.522 | 3.5350564113604225e-06 | 0.033 | 6 | 0.6 | 0.4285714285714285 | 13 | 0.52 | 0.3513513513513513 | 26 | 0.52 | 0.3513513513513513 |
| OFC BA11 | 70 | 70 | 1.000 | 0.392 | 0.0007830988069511 | 0.016 | 3 | 0.3 | 0.1764705882352941 | 13 | 0.52 | 0.3513513513513513 | 33 | 0.66 | 0.4925373134328358 |
| Ventral Striatum | 146 | 134 | 0.918 | 0.562 | 1.6921185754155214e-12 | 0.125 | 8 | 0.8 | 0.6666666666666666 | 13 | 0.52 | 0.3513513513513513 | 19 | 0.38 | 0.2345679012345679 |
| dlPFC BA8/9 | 72 | 72 | 1.000 | 0.632 | 2.6072647088455127e-09 | 0.034 | 6 | 0.6 | 0.4285714285714285 | 16 | 0.64 | 0.4705882352941176 | 30 | 0.6 | 0.4285714285714285 |

This is a **within-study convergence/sensitivity comparison**. Neither ranking influenced the other, but both use the same cohort, so agreement is not independent validation. Direct top-10/top-25/top-50 overlaps and rank correlations are descriptive.

## G. Interpretation and next steps

The primary output is a continuous region-specific ranking of ML-prioritized exploratory Tier2 candidates. No High/Moderate/Low thresholds were created. Label permutations, Stability Selection, handmade-score benchmarking, LASSO-only baselines, limma/voom, DESeq2, and external validation remain later phases and were not run here.
