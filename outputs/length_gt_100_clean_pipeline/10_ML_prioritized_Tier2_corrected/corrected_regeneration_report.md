# Corrected ML-prioritisation regeneration report

## Verification

- v1 SHA-256: 406ba687a9f3f4244de0d296550e1b2fcfb043f504c175c56328b447133cde1e.
- v2 SHA-256: 7b33e7d0a6eb2ff2efaffdc223aced7f4722f18236a2902cc9102257b2215ef8.
- Corrected CSV reconstruction maximum difference: 9.71e-17.
- Frozen legacy reconstruction maximum difference: 9.71e-17.
- Frozen legacy class matches: 671/671.
- The frozen CSV used legacy first-three weights 0.30/0.20/0.15; the corrected CSV uses the mandated 0.40/0.15/0.10.

## Counts

- Old v2 denominator: 180.
- Corrected v2 denominator values: [150].
- Old high/moderate/low: 35/130/506.
- Corrected high/moderate/low: 30/119/522.

## Confidence transitions

| transition           |   n_gene_region_rows |
|:---------------------|---------------------:|
| high_to_high         |                   30 |
| high_to_moderate     |                    5 |
| low_to_low           |                  506 |
| moderate_to_low      |                   16 |
| moderate_to_moderate |                  114 |

## Isolated correction counts

| scenario               |   high |   moderate |   low |
|:-----------------------|-------:|-----------:|------:|
| legacy_reconstructed   |     35 |        130 |   506 |
| verified_weights_only  |     34 |        131 |   506 |
| RF_sign_fix_only       |     35 |        130 |   506 |
| baseline_removal_only  |     31 |        118 |   522 |
| missing_sign_gate_only |     35 |        130 |   506 |
| fully_corrected        |     30 |        119 |   522 |

## Isolated correction effects

| scenario               |   n_scores_changed |   maximum_absolute_score_change |   n_confidence_labels_changed |
|:-----------------------|-------------------:|--------------------------------:|------------------------------:|
| verified_weights_only  |                644 |                       0.09712   |                             1 |
| RF_sign_fix_only       |                253 |                       0.05      |                             0 |
| baseline_removal_only  |                644 |                       0.0260653 |                            20 |
| missing_sign_gate_only |                  0 |                       0         |                             0 |
| fully_corrected        |                644 |                       0.107462  |                            21 |

## Real-gene arithmetic

| analysis_level    | region_or_subregion   |   gene_id | gene_symbol   | old_class   | corrected_class   |   v1_selected_events |   v1_denominator |   selection_frequency_v1 |   v2_selected_events |   v2_denominator |   selection_frequency_v2 |   combined_selection_frequency |   model_support_scaled |   sparse_model_support_binary |   importance_rank_score |   corrected_sign_consistency |   sign_consistency_score |   lmm_evidence_scaled |   v2_support_score |   reconstructed_score |   csv_score |   absolute_arithmetic_difference | confidence_label                      |
|:------------------|:----------------------|----------:|:--------------|:------------|:------------------|---------------------:|-----------------:|-------------------------:|---------------------:|-----------------:|-------------------------:|-------------------------------:|-----------------------:|------------------------------:|------------------------:|-----------------------------:|-------------------------:|----------------------:|-------------------:|----------------------:|------------:|---------------------------------:|:--------------------------------------|
| frontal_subregion | Cingulate gyrus 25    | 107985795 | LOC107985795  | high        | high              |                  108 |              125 |                    0.864 |                   70 |              150 |                 0.466667 |                         0.6256 |                      1 |                             1 |                0.799251 |                            1 |                        1 |              0.49357  |                  1 |              0.794806 |    0.794806 |                                0 | High-confidence ML-prioritized Tier 2 |
| four_group        | Insula                |    153768 | PRELID2       | moderate    | low               |                   34 |              125 |                    0.272 |                   30 |              150 |                 0.2      |                         0.2288 |                      1 |                             1 |                0.34827  |                            1 |                        1 |              0.376898 |                  1 |              0.562605 |    0.562605 |                                0 | Exploratory ML-prioritized Tier 2     |

These remain exploratory prioritisation outputs, not validated biomarkers.
