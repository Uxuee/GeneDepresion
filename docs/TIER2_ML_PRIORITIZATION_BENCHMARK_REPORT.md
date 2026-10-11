# Tier2 ML prioritisation method benchmark

## 1. Scientific question

This benchmark asks whether the frozen five-model Tier2–RRA ranking adds useful prioritisation stability beyond two simpler sparse alternatives: LASSO alone and classical half-sampling Stability Selection. The endpoint is gene prioritisation—not classifier accuracy.

Tier2 was defined using Diagnosis in the complete cohort. Therefore all results are exploratory, within-study comparisons. They are not external validation, and the genes are not validated biomarkers.

## 2. Frozen primary Tier2–RRA method

The primary result was read without modification from frozen commit `df3596ed208bb7472e166e5fc0330b69453a9e31`. It combines L1, Elastic Net, L2, linear SVM, and Random Forest rankings across 25 repeated outer splits using two-level Robust Rank Aggregation. No benchmark output was allowed to change this ranking.

| Region | Tier2 candidates |
| --- | --- |
| CG25 | 186 |
| Hippocampal Formation | 127 |
| Insula | 70 |
| OFC BA11 | 70 |
| Ventral Striatum | 146 |
| dlPFC BA8/9 | 72 |

## 3. LASSO-alone benchmark

The exact 150 L1-logistic fits already present in the frozen primary run were reused: 25 fits per region, training-only median imputation and scaling, and training-only three-fold tuning over `C = 0.01, 0.1, 1`. The primary quantity is selection frequency across 25 fits. Tie-breakers are mean selected rank, mean absolute coefficient, and stable gene ID. No weighted score was created.

## 4. Stability Selection benchmark

No compatible Stability Selection package was installed in the frozen environment. The prospective design therefore used a transparent classical procedure: 100 deterministic stratified half-samples per region, L1 logistic regression, training-only preprocessing, and three-fold tuning over the same C grid within every half-sample. The continuous selection probability is primary. A predeclared 0.60 threshold is descriptive only.

| region | genes_probability_ge_0_60 |
| --- | --- |
| CG25 | 0 |
| Hippocampal Formation | 1 |
| Insula | 2 |
| OFC BA11 | 0 |
| Ventral Striatum | 5 |
| dlPFC BA8/9 | 3 |

The low threshold counts in some regions are a result, not a reason to change the threshold. In particular, half-sampling the small CG25 cohort is demanding and produces sparse, variable selections.

## 5. Fair comparison design

All methods received exactly the same 671 Tier2 region–gene pairs and the same regional samples. No Tier4, historical score, all-gene rank, WGCNA, or AHBA information entered fitting. Never-selected genes were assigned a shared worst-rank block for correlation; they were not silently dropped. Top-set comparisons used top 10, top 25, and `ceil(20% × regional Tier2 universe)`.

## 6. Rank agreement

| method_left | method_right | spearman_min | spearman_median | spearman_max | kendall_min | kendall_median | kendall_max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| LASSO | Stability Selection | 0.674 | 0.821 | 0.897 | 0.552 | 0.683 | 0.748 |
| Tier2-RRA | LASSO | 0.654 | 0.745 | 0.834 | 0.489 | 0.577 | 0.686 |
| Tier2-RRA | Stability Selection | 0.602 | 0.638 | 0.688 | 0.439 | 0.472 | 0.521 |

RRA–LASSO Spearman agreement ranged from 0.654 to 0.834. RRA–Stability Selection agreement ranged from 0.602 to 0.688. These are moderate-to-strong within-study convergence signals, not independent replication.

## 7. Top-set overlap

| method_left | method_right | top_definition | jaccard_min | jaccard_median | jaccard_max |
| --- | --- | --- | --- | --- | --- |
| LASSO | Stability Selection | top10 | 0.429 | 0.667 | 1.000 |
| LASSO | Stability Selection | top20percent | 0.520 | 0.634 | 0.793 |
| LASSO | Stability Selection | top25 | 0.429 | 0.640 | 0.786 |
| Tier2-RRA | LASSO | top10 | 0.538 | 0.818 | 0.818 |
| Tier2-RRA | LASSO | top20percent | 0.500 | 0.699 | 0.875 |
| Tier2-RRA | LASSO | top25 | 0.562 | 0.667 | 0.852 |
| Tier2-RRA | Stability Selection | top10 | 0.538 | 0.667 | 0.818 |
| Tier2-RRA | Stability Selection | top20percent | 0.579 | 0.657 | 0.750 |
| Tier2-RRA | Stability Selection | top25 | 0.471 | 0.562 | 0.724 |

The frozen RRA top 10 shared 7–9 genes with LASSO and 7–9 genes with Stability Selection across regions. Thus the strongest candidates are not artifacts of one aggregation method, although ordering and lower-ranked candidates differ.

## 8. Resampling stability

| method | stability_metric | top_definition | regional_mean_min | regional_mean_median | regional_mean_max |
| --- | --- | --- | --- | --- | --- |
| LASSO | jaccard | top10 | 0.184 | 0.296 | 0.503 |
| LASSO | jaccard | top20percent | 0.222 | 0.369 | 0.418 |
| LASSO | jaccard | top25 | 0.222 | 0.386 | 0.422 |
| LASSO | kendall_rank | all_candidates | 0.260 | 0.427 | 0.543 |
| LASSO | spearman_rank | all_candidates | 0.299 | 0.467 | 0.567 |
| Stability Selection | jaccard | top10 | 0.048 | 0.143 | 0.195 |
| Stability Selection | jaccard | top20percent | 0.048 | 0.166 | 0.209 |
| Stability Selection | jaccard | top25 | 0.048 | 0.168 | 0.209 |
| Stability Selection | kendall_rank | all_candidates | 0.065 | 0.168 | 0.264 |
| Stability Selection | spearman_rank | all_candidates | 0.069 | 0.179 | 0.275 |
| Tier2-RRA | jaccard | top10 | 0.212 | 0.321 | 0.476 |
| Tier2-RRA | jaccard | top20percent | 0.299 | 0.411 | 0.469 |
| Tier2-RRA | jaccard | top25 | 0.346 | 0.456 | 0.499 |
| Tier2-RRA | kendall_rank | all_candidates | 0.431 | 0.576 | 0.673 |
| Tier2-RRA | spearman_rank | all_candidates | 0.486 | 0.631 | 0.716 |

Split-level five-model RRA had higher mean all-candidate rank stability than LASSO alone in every region and substantially higher stability than half-sampling Stability Selection. Stability Selection was especially variable for CG25 and OFC BA11. This does not invalidate its continuous final ranking; it shows that individual half-sample selections are sensitive to the small cohorts.

## 9. Coverage and sparsity

| region | method | tier2_universe_size | genes_with_any_support | proportion_with_any_support | support_mean | support_median | support_minimum | support_maximum |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Ventral Striatum | Tier2-RRA | 146.000 | 146.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Insula | Tier2-RRA | 70.000 | 70.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Hippocampal Formation | Tier2-RRA | 127.000 | 127.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| CG25 | Tier2-RRA | 186.000 | 186.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| dlPFC BA8/9 | Tier2-RRA | 72.000 | 72.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| OFC BA11 | Tier2-RRA | 70.000 | 70.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| CG25 | LASSO | 186.000 | 71.000 | 0.382 | 0.073 | 0.000 | 0.000 | 0.800 |
| dlPFC BA8/9 | LASSO | 72.000 | 49.000 | 0.681 | 0.241 | 0.080 | 0.000 | 1.000 |
| Hippocampal Formation | LASSO | 127.000 | 58.000 | 0.457 | 0.140 | 0.000 | 0.000 | 1.000 |
| Insula | LASSO | 70.000 | 54.000 | 0.771 | 0.267 | 0.140 | 0.000 | 0.960 |
| OFC BA11 | LASSO | 70.000 | 64.000 | 0.914 | 0.298 | 0.220 | 0.000 | 0.960 |
| Ventral Striatum | LASSO | 146.000 | 70.000 | 0.479 | 0.121 | 0.000 | 0.000 | 1.000 |
| Ventral Striatum | Stability Selection | 146.000 | 113.000 | 0.774 | 0.089 | 0.040 | 0.000 | 0.890 |
| Insula | Stability Selection | 70.000 | 68.000 | 0.971 | 0.178 | 0.130 | 0.000 | 0.750 |
| Hippocampal Formation | Stability Selection | 127.000 | 104.000 | 0.819 | 0.088 | 0.040 | 0.000 | 0.710 |
| CG25 | Stability Selection | 186.000 | 136.000 | 0.731 | 0.037 | 0.020 | 0.000 | 0.330 |
| dlPFC BA8/9 | Stability Selection | 72.000 | 68.000 | 0.944 | 0.166 | 0.120 | 0.000 | 0.680 |
| OFC BA11 | Stability Selection | 70.000 | 70.000 | 1.000 | 0.181 | 0.160 | 0.010 | 0.480 |

RRA ranks the full Tier2 universe because dense L2/SVM models contribute evidence. LASSO is deliberately sparse. Stability Selection gives more genes at least occasional support, but generally low selection probabilities. Coverage is descriptive and is not scored as better or worse.

## 10. Descriptive predictive metrics

| region | n_outer_fits | runtime_seconds | balanced_accuracy_mean | f1_mean | roc_auc_mean | precision_mean | recall_mean | balanced_accuracy_sd | f1_sd | roc_auc_sd | precision_sd | recall_sd | source | performance_interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CG25 | 25.000 | 1.181 | 0.697 | 0.620 | 0.869 | 0.643 | 0.660 | 0.161 | 0.268 | 0.149 | 0.308 | 0.328 | reused_exact_frozen_primary_logistic_l1_events | optimistic within-study diagnostics, not unbiased out-of-sample estimates |
| Hippocampal Formation | 25.000 | 1.249 | 0.824 | 0.816 | 0.894 | 0.900 | 0.778 | 0.107 | 0.139 | 0.100 | 0.107 | 0.187 | reused_exact_frozen_primary_logistic_l1_events | optimistic within-study diagnostics, not unbiased out-of-sample estimates |
| Insula | 25.000 | 1.229 | 0.854 | 0.845 | 0.949 | 0.925 | 0.800 | 0.081 | 0.101 | 0.082 | 0.096 | 0.160 | reused_exact_frozen_primary_logistic_l1_events | optimistic within-study diagnostics, not unbiased out-of-sample estimates |
| OFC BA11 | 25.000 | 1.073 | 0.738 | 0.737 | 0.862 | 0.800 | 0.736 | 0.122 | 0.138 | 0.146 | 0.148 | 0.229 | reused_exact_frozen_primary_logistic_l1_events | optimistic within-study diagnostics, not unbiased out-of-sample estimates |
| Ventral Striatum | 25.000 | 1.184 | 0.915 | 0.925 | 0.943 | 0.926 | 0.929 | 0.108 | 0.095 | 0.097 | 0.112 | 0.096 | reused_exact_frozen_primary_logistic_l1_events | optimistic within-study diagnostics, not unbiased out-of-sample estimates |
| dlPFC BA8/9 | 25.000 | 1.033 | 0.828 | 0.855 | 0.892 | 0.835 | 0.893 | 0.136 | 0.116 | 0.123 | 0.150 | 0.136 | reused_exact_frozen_primary_logistic_l1_events | optimistic within-study diagnostics, not unbiased out-of-sample estimates |

These are optimistic within-study diagnostics, not unbiased out-of-sample estimates. They do not decide which prioritisation method is preferable.

## 11. All-gene sensitivity comparison

| region | method | n_tier2_candidates | n_present_in_allgene_rra | spearman_rho_common | top10_intersection | top10_jaccard | top10_overlap_fraction_method | top10_overlap_fraction_allgene | top25_intersection | top25_jaccard | top25_overlap_fraction_method | top25_overlap_fraction_allgene | top20percent_intersection | top20percent_jaccard | top20percent_overlap_fraction_method | top20percent_overlap_fraction_allgene |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Ventral Striatum | Tier2-RRA | 146.000 | 134.000 | 0.562 | 8.000 | 0.667 | 0.800 | 0.800 | 13.000 | 0.351 | 0.520 | 0.520 | 14.000 | 0.304 | 0.467 | 0.467 |
| Insula | Tier2-RRA | 70.000 | 70.000 | 0.522 | 6.000 | 0.429 | 0.600 | 0.600 | 13.000 | 0.351 | 0.520 | 0.520 | 10.000 | 0.556 | 0.714 | 0.714 |
| Hippocampal Formation | Tier2-RRA | 127.000 | 119.000 | 0.575 | 5.000 | 0.333 | 0.500 | 0.500 | 13.000 | 0.351 | 0.520 | 0.520 | 14.000 | 0.368 | 0.538 | 0.538 |
| CG25 | Tier2-RRA | 186.000 | 169.000 | 0.685 | 5.000 | 0.333 | 0.500 | 0.500 | 14.000 | 0.389 | 0.560 | 0.560 | 19.000 | 0.333 | 0.500 | 0.500 |
| dlPFC BA8/9 | Tier2-RRA | 72.000 | 72.000 | 0.632 | 6.000 | 0.429 | 0.600 | 0.600 | 16.000 | 0.471 | 0.640 | 0.640 | 6.000 | 0.250 | 0.400 | 0.400 |
| OFC BA11 | Tier2-RRA | 70.000 | 70.000 | 0.392 | 3.000 | 0.176 | 0.300 | 0.300 | 13.000 | 0.351 | 0.520 | 0.520 | 6.000 | 0.273 | 0.429 | 0.429 |
| CG25 | LASSO | 186.000 | 169.000 | 0.713 | 6.000 | 0.429 | 0.600 | 0.600 | 15.000 | 0.429 | 0.600 | 0.600 | 22.000 | 0.407 | 0.579 | 0.579 |
| dlPFC BA8/9 | LASSO | 72.000 | 72.000 | 0.630 | 5.000 | 0.333 | 0.500 | 0.500 | 14.000 | 0.389 | 0.560 | 0.560 | 6.000 | 0.250 | 0.400 | 0.400 |
| Hippocampal Formation | LASSO | 127.000 | 119.000 | 0.609 | 5.000 | 0.333 | 0.500 | 0.500 | 13.000 | 0.351 | 0.520 | 0.520 | 14.000 | 0.368 | 0.538 | 0.538 |
| Insula | LASSO | 70.000 | 70.000 | 0.669 | 6.000 | 0.429 | 0.600 | 0.600 | 14.000 | 0.389 | 0.560 | 0.560 | 7.000 | 0.333 | 0.500 | 0.500 |
| OFC BA11 | LASSO | 70.000 | 70.000 | 0.662 | 3.000 | 0.176 | 0.300 | 0.300 | 15.000 | 0.429 | 0.600 | 0.600 | 5.000 | 0.217 | 0.357 | 0.357 |
| Ventral Striatum | LASSO | 146.000 | 134.000 | 0.682 | 9.000 | 0.818 | 0.900 | 0.900 | 14.000 | 0.389 | 0.560 | 0.560 | 17.000 | 0.395 | 0.567 | 0.567 |
| Ventral Striatum | Stability Selection | 146.000 | 134.000 | 0.727 | 9.000 | 0.818 | 0.900 | 0.900 | 14.000 | 0.389 | 0.560 | 0.560 | 17.000 | 0.395 | 0.567 | 0.567 |
| Insula | Stability Selection | 70.000 | 70.000 | 0.718 | 5.000 | 0.333 | 0.500 | 0.500 | 14.000 | 0.389 | 0.560 | 0.560 | 9.000 | 0.474 | 0.643 | 0.643 |
| Hippocampal Formation | Stability Selection | 127.000 | 119.000 | 0.709 | 5.000 | 0.333 | 0.500 | 0.500 | 14.000 | 0.389 | 0.560 | 0.560 | 14.000 | 0.368 | 0.538 | 0.538 |
| CG25 | Stability Selection | 186.000 | 169.000 | 0.624 | 7.000 | 0.538 | 0.700 | 0.700 | 12.000 | 0.316 | 0.480 | 0.480 | 17.000 | 0.288 | 0.447 | 0.447 |
| dlPFC BA8/9 | Stability Selection | 72.000 | 72.000 | 0.740 | 4.000 | 0.250 | 0.400 | 0.400 | 16.000 | 0.471 | 0.640 | 0.640 | 9.000 | 0.429 | 0.600 | 0.600 |
| OFC BA11 | Stability Selection | 70.000 | 70.000 | 0.656 | 4.000 | 0.250 | 0.400 | 0.400 | 16.000 | 0.471 | 0.640 | 0.640 | 7.000 | 0.333 | 0.500 | 0.500 |

This comparison uses the frozen all-gene RRA only after benchmark rankings were fixed. It is within-study sensitivity/convergence, not independent validation.

## 12. Top genes by method

| region | method | rank | gene_id | gene_symbol | primary_evidence |
| --- | --- | --- | --- | --- | --- |
| Ventral Striatum | Tier2-RRA | 1 | 254956 | MORN5 | 2.412e-09 |
| Ventral Striatum | Tier2-RRA | 2 | 79983 | POF1B | 2.355e-07 |
| Ventral Striatum | Tier2-RRA | 3 | 2488 | FSHB | 5.861e-07 |
| Ventral Striatum | Tier2-RRA | 4 | 51069 | MRPL2 | 8.707e-07 |
| Ventral Striatum | Tier2-RRA | 5 | 64598 | MOSPD3 | 2.47e-06 |
| Ventral Striatum | LASSO | 1 | 254956 | MORN5 | 1 |
| Ventral Striatum | LASSO | 2 | 51069 | MRPL2 | 1 |
| Ventral Striatum | LASSO | 3 | 64598 | MOSPD3 | 1 |
| Ventral Striatum | LASSO | 4 | 2488 | FSHB | 1 |
| Ventral Striatum | LASSO | 5 | 642826 | BMS1P2 | 1 |
| Ventral Striatum | Stability Selection | 1 | 254956 | MORN5 | 0.89 |
| Ventral Striatum | Stability Selection | 2 | 51069 | MRPL2 | 0.85 |
| Ventral Striatum | Stability Selection | 3 | 2488 | FSHB | 0.66 |
| Ventral Striatum | Stability Selection | 4 | 79983 | POF1B | 0.65 |
| Ventral Striatum | Stability Selection | 5 | 64598 | MOSPD3 | 0.64 |
| Insula | Tier2-RRA | 1 | 100505894 | TMEM161B-DT | 7.229e-07 |
| Insula | Tier2-RRA | 2 | 101929710 | LOC101929710 | 2.313e-05 |
| Insula | Tier2-RRA | 3 | 642826 | BMS1P2 | 5e-05 |
| Insula | Tier2-RRA | 4 | 105376607 | ARL14EP-DT | 9.748e-05 |
| Insula | Tier2-RRA | 5 | 54718 | BTN2A3P | 0.001105 |
| Insula | LASSO | 1 | 100505894 | TMEM161B-DT | 0.96 |
| Insula | LASSO | 2 | 101929710 | LOC101929710 | 0.96 |
| Insula | LASSO | 3 | 105376607 | ARL14EP-DT | 0.96 |
| Insula | LASSO | 4 | 26063 | DECR2 | 0.88 |
| Insula | LASSO | 5 | 107985195 | LINC01691 | 0.84 |
| Insula | Stability Selection | 1 | 100505894 | TMEM161B-DT | 0.75 |
| Insula | Stability Selection | 2 | 101929710 | LOC101929710 | 0.64 |
| Insula | Stability Selection | 3 | 105376607 | ARL14EP-DT | 0.55 |
| Insula | Stability Selection | 4 | 26063 | DECR2 | 0.55 |
| Insula | Stability Selection | 5 | 107985195 | LINC01691 | 0.53 |
| Hippocampal Formation | Tier2-RRA | 1 | 27433 | TOR2A | 9.55e-08 |
| Hippocampal Formation | Tier2-RRA | 2 | 100505894 | TMEM161B-DT | 1.177e-06 |
| Hippocampal Formation | Tier2-RRA | 3 | 105376607 | ARL14EP-DT | 8.936e-06 |
| Hippocampal Formation | Tier2-RRA | 4 | 105370868 | DIS3L-AS1 | 5.619e-05 |
| Hippocampal Formation | Tier2-RRA | 5 | 259307 | IL4I1 | 5.817e-05 |
| Hippocampal Formation | LASSO | 1 | 27433 | TOR2A | 1 |
| Hippocampal Formation | LASSO | 2 | 259307 | IL4I1 | 0.92 |
| Hippocampal Formation | LASSO | 3 | 100505894 | TMEM161B-DT | 0.8 |
| Hippocampal Formation | LASSO | 4 | 285830 | HLA-F-AS1 | 0.8 |
| Hippocampal Formation | LASSO | 5 | 105376607 | ARL14EP-DT | 0.76 |
| Hippocampal Formation | Stability Selection | 1 | 27433 | TOR2A | 0.71 |
| Hippocampal Formation | Stability Selection | 2 | 259307 | IL4I1 | 0.5 |
| Hippocampal Formation | Stability Selection | 3 | 100505894 | TMEM161B-DT | 0.47 |
| Hippocampal Formation | Stability Selection | 4 | 27037 | TRMT2A | 0.43 |
| Hippocampal Formation | Stability Selection | 5 | 151 | ADRA2B | 0.41 |
| CG25 | Tier2-RRA | 1 | 107985795 | LOC107985795 | 3.775e-07 |
| CG25 | Tier2-RRA | 2 | 79760 | GEMIN7 | 3.617e-06 |
| CG25 | Tier2-RRA | 3 | 105372507 | LOC105372507 | 5.255e-06 |
| CG25 | Tier2-RRA | 4 | 401010 | NOC2LP2 | 5.255e-06 |
| CG25 | Tier2-RRA | 5 | 919 | CD247 | 8.339e-06 |
| CG25 | LASSO | 1 | 107985795 | LOC107985795 | 0.8 |
| CG25 | LASSO | 2 | 105377644 | LOC105377644 | 0.72 |
| CG25 | LASSO | 3 | 105372507 | LOC105372507 | 0.72 |
| CG25 | LASSO | 4 | 401010 | NOC2LP2 | 0.68 |
| CG25 | LASSO | 5 | 79760 | GEMIN7 | 0.6 |
| CG25 | Stability Selection | 1 | 107985795 | LOC107985795 | 0.33 |
| CG25 | Stability Selection | 2 | 79760 | GEMIN7 | 0.25 |
| CG25 | Stability Selection | 3 | 339448 | C1orf174 | 0.23 |
| CG25 | Stability Selection | 4 | 105377644 | LOC105377644 | 0.19 |
| CG25 | Stability Selection | 5 | 100861412 | FSBP | 0.18 |
| dlPFC BA8/9 | Tier2-RRA | 1 | 100996279 | LINC00269 | 8.075e-06 |
| dlPFC BA8/9 | Tier2-RRA | 2 | 65083 | NOL6 | 8.075e-06 |
| dlPFC BA8/9 | Tier2-RRA | 3 | 9374 | PPT2 | 1.455e-05 |
| dlPFC BA8/9 | Tier2-RRA | 4 | 101929607 | LINC01990 | 2.009e-05 |
| dlPFC BA8/9 | Tier2-RRA | 5 | 107984015 | LOC107984015 | 8.468e-05 |
| dlPFC BA8/9 | LASSO | 1 | 9374 | PPT2 | 1 |
| dlPFC BA8/9 | LASSO | 2 | 65083 | NOL6 | 0.96 |
| dlPFC BA8/9 | LASSO | 3 | 107984015 | LOC107984015 | 0.88 |
| dlPFC BA8/9 | LASSO | 4 | 100287015 | MCPH1-DT | 0.84 |
| dlPFC BA8/9 | LASSO | 5 | 105371009 | LOC105371009 | 0.84 |
| dlPFC BA8/9 | Stability Selection | 1 | 65083 | NOL6 | 0.68 |
| dlPFC BA8/9 | Stability Selection | 2 | 9374 | PPT2 | 0.66 |
| dlPFC BA8/9 | Stability Selection | 3 | 100996279 | LINC00269 | 0.64 |
| dlPFC BA8/9 | Stability Selection | 4 | 105371009 | LOC105371009 | 0.51 |
| dlPFC BA8/9 | Stability Selection | 5 | 107984015 | LOC107984015 | 0.49 |
| OFC BA11 | Tier2-RRA | 1 | 105373386 | LOC105373386 | 8.145e-05 |
| OFC BA11 | Tier2-RRA | 2 | 105378913 | LOC105378913 | 8.145e-05 |
| OFC BA11 | Tier2-RRA | 3 | 53917 | RAB24 | 0.0006136 |
| OFC BA11 | Tier2-RRA | 4 | 162514 | TRPV3 | 0.001105 |
| OFC BA11 | Tier2-RRA | 5 | 158055 | C9orf163 | 0.001257 |
| OFC BA11 | LASSO | 1 | 105378913 | LOC105378913 | 0.96 |
| OFC BA11 | LASSO | 2 | 23475 | QPRT | 0.92 |
| OFC BA11 | LASSO | 3 | 105373386 | LOC105373386 | 0.8 |
| OFC BA11 | LASSO | 4 | 53917 | RAB24 | 0.76 |
| OFC BA11 | LASSO | 5 | 158055 | C9orf163 | 0.76 |
| OFC BA11 | Stability Selection | 1 | 105378913 | LOC105378913 | 0.48 |
| OFC BA11 | Stability Selection | 2 | 92283 | ZNF461 | 0.4 |
| OFC BA11 | Stability Selection | 3 | 158055 | C9orf163 | 0.39 |
| OFC BA11 | Stability Selection | 4 | 23475 | QPRT | 0.37 |
| OFC BA11 | Stability Selection | 5 | 55620 | STAP2 | 0.37 |

## 13. Interpretation

The evidence best matches **Scenario E, with an important Scenario B nuance**. LASSO alone recovers much of the RRA shortlist and therefore provides a strong, simpler benchmark. However, split-level RRA rankings are consistently more stable than LASSO rankings, indicating that the ensemble contributes useful consensus beyond one sparse model. Stability Selection converges substantially with the final RRA top genes but is not more stable under the prespecified half-sampling design; it therefore supports rather than challenges retaining RRA.

This is not proof of biological truth. All methods reuse the same small cohort and the same Diagnosis-informed Tier2 screen.

## 14. Limitations

- Tier2 was defined in the complete cohort, so downstream performance metrics are optimistic.
- Resampling observations are not independent biological replications.
- CG25 has only 27 samples; half-samples are especially small.
- Stability Selection probabilities depend on the prespecified C grid and half-sampling design.
- Dense-model RRA and sparse methods answer slightly different questions about coverage.
- No external cohort was tested.

## 15. Recommendation for the next phase

Retain Tier2–RRA as the primary continuous consensus ranking, present LASSO as the transparent simpler benchmark, and present Stability Selection as a supportive sensitivity analysis. Emphasize that LASSO recovers most leading candidates, while RRA provides measurably greater resampling stability. The next prospective phase should be Tier2 conditional label-permutation/null testing. That phase was not run here. Limma/voom, DESeq2, and external validation were also not run.

## Runtime and complexity

| method | final_fits | inner_evaluation_fits | runtime_seconds | models | tuning_burden | interpretability |
| --- | --- | --- | --- | --- | --- | --- |
| Tier2-RRA | 750.000 | 6300.000 | 489.813 | 5.000 | five model families; 14 parameter settings total per split | consensus ranking; model-level ranks retained |
| LASSO | 150.000 | 1350.000 | 6.949 | 1.000 | three C values in 3-fold CV per outer split | selection frequency with transparent tie-breakers; exact frozen events reused |
| Stability Selection | 600.000 | 5400.000 | 20.051 | 1.000 | three C values in 3-fold CV per half-sample | continuous selection probability |
