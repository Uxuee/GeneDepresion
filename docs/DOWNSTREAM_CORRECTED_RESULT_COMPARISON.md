# Downstream Corrected Result Comparison

## Scope

This report isolates changes caused by replacing the legacy 35-high/130-moderate ranking with the canonical 30-high/119-moderate ranking. No expression normalisation, diagnosis model, tier, ML training/scoring, AHBA module, WGCNA module, eigengene, connectivity, or module-level ML model was rebuilt.

## Membership

| version   | region                |   high |   moderate |   high_plus_moderate |
|:----------|:----------------------|-------:|-----------:|---------------------:|
| legacy    | Cingulate gyrus 25    |      8 |         29 |                   37 |
| corrected | Cingulate gyrus 25    |      6 |         30 |                   36 |
| legacy    | Hippocampal Formation |      3 |         25 |                   28 |
| corrected | Hippocampal Formation |      2 |         24 |                   26 |
| legacy    | Insula                |      6 |         18 |                   24 |
| corrected | Insula                |      5 |         15 |                   20 |
| legacy    | OFC BA11              |      3 |         28 |                   31 |
| corrected | OFC BA11              |      3 |         24 |                   27 |
| legacy    | Ventral Striatum      |      8 |         18 |                   26 |
| corrected | Ventral Striatum      |      7 |         16 |                   23 |
| legacy    | dlPFC BA8/9           |      7 |         12 |                   19 |
| corrected | dlPFC BA8/9           |      7 |         10 |                   17 |

## AHBA

Across the directly comparable ML-dependent tests, 7 rows changed nominal/FDR significance category. The corrected analysis has no published-AHBA module enrichment surviving BH-FDR < 0.10. The former Ventral Striatum up-in-MDD M26 result moved from FDR<0.10 to nominal p<0.05 after TSGA10 left the overlap. These changes are annotation/context changes, not evidence for or against a validated disease mechanism.

| candidate_set      | scope                                         | published_module   |   overlap_genes_old |   overlap_genes_corrected |   p_value_hypergeometric_old |   p_value_hypergeometric_corrected |   fdr_BH_old |   fdr_BH_corrected | old_significance   | corrected_significance   | genes_lost   |   genes_gained |
|:-------------------|:----------------------------------------------|:-------------------|--------------------:|--------------------------:|-----------------------------:|-----------------------------------:|-------------:|-------------------:|:-------------------|:-------------------------|:-------------|---------------:|
| high_confidence    | all_regions                                   | M19                |                   3 |                         3 |                   0.0548849  |                          0.0247726 |    1         |           0.500357 | none               | nominal_p<0.05           | nan          |            nan |
| high_confidence    | all_regions                                   | M24                |                   3 |                         3 |                   0.0682259  |                          0.0312723 |    1         |           0.500357 | none               | nominal_p<0.05           | nan          |            nan |
| high_plus_moderate | region:Insula|direction:down_in_MDD           | M18                |                   1 |                         0 |                   0.0317989  |                          1         |    1         |           1        | nominal_p<0.05     | none                     | ZKSCAN2      |            nan |
| high_plus_moderate | region:Insula|direction:down_in_MDD           | M6                 |                   1 |                         1 |                   0.0741442  |                          0.0377844 |    1         |           1        | none               | nominal_p<0.05           | nan          |            nan |
| high_plus_moderate | region:OFC BA11|direction:up_in_MDD           | M31                |                   2 |                         2 |                   0.0547404  |                          0.0358404 |    1         |           0.89289  | none               | nominal_p<0.05           | nan          |            nan |
| high_plus_moderate | region:Ventral Striatum|direction:down_in_MDD | M18                |                   1 |                         0 |                   0.0317989  |                          1         |    0.508783  |           1        | nominal_p<0.05     | none                     | ZKSCAN2      |            nan |
| high_plus_moderate | region:Ventral Striatum|direction:up_in_MDD   | M26                |                   3 |                         2 |                   0.00294979 |                          0.028013  |    0.0943933 |           0.882283 | FDR<0.10           | nominal_p<0.05           | TSGA10       |            nan |

## WGCNA

Underlying module assignments were preserved byte-for-byte. 10 broad/frontal enrichment rows changed nominal/FDR category. Broad Ventral Striatum MEdarkred high+moderate enrichment moved from FDR<0.10 to nominal after PRELID2 and TMEM269-DT (gene IDs 153768 and 728975) left its overlap. In the original frontal subregions, no enrichment row changed category; corrected OFC BA11 MEtan retains FDR<0.10 enrichment for the global high set and global high+moderate set.

| analysis_unit   | region                | module           | gene_set                                |   overlap_count_old |   overlap_count_corrected |   p_value_old |   p_value_corrected |   fdr_BH_old |   fdr_BH_corrected | significance_category_old   | significance_category_corrected   | genes_lost    |   genes_gained |
|:----------------|:----------------------|:-----------------|:----------------------------------------|--------------------:|--------------------------:|--------------:|--------------------:|-------------:|-------------------:|:----------------------------|:----------------------------------|:--------------|---------------:|
| broad_regions   | Frontal Cortex        | MEcyan           | high_plus_moderate_ML_genes             |                   7 |                         6 |   0.03886     |          0.0661765  |    0.336786  |           0.430147 | nominal_p<0.05              | none                              | 79812         |            nan |
| broad_regions   | Frontal Cortex        | MEsalmon         | high_confidence_ML_genes                |                   2 |                         1 |   0.039055    |          0.229016   |    1         |           1        | nominal_p<0.05              | none                              | 5092          |            nan |
| broad_regions   | Hippocampal Formation | MEcoral1         | high_plus_moderate_ML_genes             |                   2 |                         1 |   0.026855    |          0.205681   |    0.617664  |           1        | nominal_p<0.05              | none                              | 105371659     |            nan |
| broad_regions   | Insula                | MEblue           | same_region_high_plus_moderate_ML_genes |                   8 |                         6 |   0.0178176   |          0.088267   |    0.516709  |           1        | nominal_p<0.05              | none                              | 153768;3909   |            nan |
| broad_regions   | Insula                | MEbrown          | high_plus_moderate_ML_genes             |                  17 |                        16 |   0.0547781   |          0.0412037  |    1         |           1        | none                        | nominal_p<0.05                    | 79812         |            nan |
| broad_regions   | Ventral Striatum      | MEdarkolivegreen | high_plus_moderate_ML_genes             |                   3 |                         2 |   0.0497691   |          0.166787   |    0.993685  |           1        | nominal_p<0.05              | none                              | 101927003     |            nan |
| broad_regions   | Ventral Striatum      | MEdarkred        | high_plus_moderate_ML_genes             |                   8 |                         6 |   0.000127594 |          0.00217014 |    0.0063797 |           0.108507 | FDR<0.10                    | nominal_p<0.05                    | 153768;728975 |            nan |
| broad_regions   | Ventral Striatum      | MEdarkred        | same_region_high_plus_moderate_ML_genes |                   2 |                         1 |   0.0184074   |          0.167278   |    0.920368  |           1        | nominal_p<0.05              | none                              | 728975        |            nan |
| broad_regions   | Ventral Striatum      | MEhoneydew1      | high_confidence_ML_genes                |                   1 |                         1 |   0.0501293   |          0.0412002  |    1         |           1        | none                        | nominal_p<0.05                    | nan           |            nan |
| broad_regions   | Ventral Striatum      | MEsienna3        | high_plus_moderate_ML_genes             |                   2 |                         2 |   0.0596211   |          0.0483675  |    0.993685  |           1        | none                        | nominal_p<0.05                    | nan           |            nan |

## Stable core and networks

The WGCNA structural core was not recalculated because its selection uses module membership, expression presence, and kME/hub rank—not ML priority. Versioned tables preserve every structural-core row and refresh only `is_high_confidence_ML` and `is_high_plus_moderate_ML`. Base network membership and edges are preserved; corrected node-annotation tables identify changed ML labels. Visual-only network rerenders were generated in the versioned Phase 2 folder. The main-report frontal rerender contains only the original frontal subregions (dlPFC BA8/9, OFC BA11, and CG25); pooled broad Frontal Cortex is excluded from that figure.

## Other downstream analyses

- Tier2 enrichment: unchanged input; not rerun.
- Tier4 ranked GSEA: unchanged input; not rerun.
- AHBA module-score ML: fixed published modules, independent of the shortlist; preserved.
- MDD module/pathway-level ML2: fixed WGCNA module scores; preserved. Only ML membership annotations changed.
- Preliminary manuscript package and dirty report copies: not overwritten. They remain queued for a separate publication refresh after human review of the changes above.

## Interpretation safeguard

Crossing a nominal or FDR threshold after changing the candidate list is reported as a sensitivity outcome, not a new biological discovery. The corrected results remain exploratory and require external validation.
