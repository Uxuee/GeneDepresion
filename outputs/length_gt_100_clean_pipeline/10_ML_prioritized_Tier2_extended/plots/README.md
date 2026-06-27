# Step 10 Plot Notes

## ML Priority Model Support Heatmap Audit

`MLpriority_model_support_heatmap.png` is a renamed copy of `heatmap_top_stable_ML_prioritized_genes.png`. It does **not** show all high-confidence or all high + moderate-confidence ML-prioritized Tier 2 genes. It shows the top 10 gene-region pairs per region/subregion by `ML_priority_score`, pivoted into a gene-symbol by region heatmap.

Recommended caption/title: **Top stable ML-prioritized Tier 2 genes by region** or **Top 10 stable ML-prioritized Tier 2 gene-region pairs per region**. Use "high-confidence genes" only for the ML-priority class; do not use "highly expressed genes" for this plot.

Audit report: `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended\heatmap_count_audit\heatmap_count_audit_report.md`
Corrected plot: `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended\heatmap_count_audit\MLpriority_model_support_heatmap_corrected.png`

## High/Moderate Count Definition

Moderate-confidence counts are **exclusive** of high-confidence counts in the current final CSV files. Report counts as 35 high-confidence and 130 additional moderate-confidence gene-region pairs, for 165 high/moderate candidates total.

`MLprioritized_counts_by_region.png` uses the correct source counts but is visually ambiguous because high and moderate bars are drawn from zero. Prefer `Extended_ML_prioritized_Tier2_gene_counts_corrected.png` for reports.

Recommended caption: **Extended ML-prioritized Tier 2 candidate genes. Bars show high-confidence genes and additional moderate-confidence genes, excluding high-confidence genes.**

Do not use "highly expressed genes" for this figure; use "high-confidence genes."
