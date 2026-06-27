# High/Moderate ML-Priority Count Definition Audit

## Question

Are moderate-confidence ML-prioritized Tier 2 counts inclusive or exclusive of high-confidence genes?

## Answer

The current final CSV files show that **moderate-confidence is exclusive of high-confidence**. High-confidence gene-region pairs are not also present in the moderate-confidence file.

- Total ranked Tier 2 gene-region pairs: 671
- High-confidence gene-region pairs: 35
- Moderate-confidence gene-region pairs, excluding high: 130
- High + moderate total: 165
- Exploratory/lower-priority Tier 2 rows: 506
- High/moderate overlap rows: 0

## Interpretation

Report the counts as: **35 high-confidence and 130 additional moderate-confidence ML-prioritized Tier 2 gene-region pairs, for 165 high/moderate candidates total**.

Use **high-confidence genes**, not **highly expressed genes**. High-confidence refers to the ML-priority class, not expression level.

## Existing Figure Check

The existing plot `MLprioritized_counts_by_region.png` was generated from `n_high_confidence` and `n_moderate_confidence`, but it draws both bars from zero rather than stacking high + additional moderate. Because moderate-confidence is exclusive, the existing plot should not be interpreted as a stacked high/moderate total.

Existing plotted counts match the final CSV tables numerically, but the older visualization is visually ambiguous for totals. A corrected stacked plot was created:

`D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended\plots\Extended_ML_prioritized_Tier2_gene_counts_corrected.png`

Recommended caption: **Extended ML-prioritized Tier 2 candidate genes. Bars show high-confidence genes and additional moderate-confidence genes, excluding high-confidence genes.**

## Per-Region Clean Count Definitions

| region_or_subregion   |   n_total_ranked_Tier2 |   n_high_confidence |   n_moderate_excluding_high |   n_high_plus_moderate_total |   n_exploratory_or_lower_priority |
|:----------------------|-----------------------:|--------------------:|----------------------------:|-----------------------------:|----------------------------------:|
| Cingulate gyrus 25    |                    186 |                   8 |                          29 |                           37 |                               149 |
| Hippocampal Formation |                    127 |                   3 |                          25 |                           28 |                                99 |
| Insula                |                     70 |                   6 |                          18 |                           24 |                                46 |
| OFC BA11              |                     70 |                   3 |                          28 |                           31 |                                39 |
| Ventral Striatum      |                    146 |                   8 |                          18 |                           26 |                               120 |
| dlPFC BA8/9           |                     72 |                   7 |                          12 |                           19 |                                53 |

## Files Created

- `high_moderate_count_audit_by_region.csv`
- `high_moderate_gene_membership_audit.csv`
- `high_moderate_count_definition_report.md`
- `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended\final_clean_high_moderate_summary_by_region.csv`
- `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended\plots\Extended_ML_prioritized_Tier2_gene_counts_corrected.png`
