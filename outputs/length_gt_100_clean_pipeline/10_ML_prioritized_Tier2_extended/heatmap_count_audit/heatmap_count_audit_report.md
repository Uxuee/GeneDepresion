# Heatmap Count Audit: ML Priority Model Support Heatmap

Audited plot: `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended\plots\MLpriority_model_support_heatmap.png`

## Source Trace

`MLpriority_model_support_heatmap.png` is a byte-identical copy of `heatmap_top_stable_ML_prioritized_genes.png`. The plot was generated in `code/extract_stable_ML_prioritized_Tier2_genes.py` from the in-memory final ranked Tier 2 table. The clean pipeline then copied it to `MLpriority_model_support_heatmap.png` in `code/length_gt_100_clean_pipeline/pipeline_lib.py`.

The plotting code uses:

```python
heat_source = result.sort_values("ML_priority_score", ascending=False).groupby("region_or_subregion").head(10)
heat = heat_source.pivot_table(index="gene_symbol", columns="region_or_subregion", values="ML_priority_score", aggfunc="max", fill_value=0)
```

## What The Heatmap Actually Shows

The heatmap shows the **top 10 stable ML-prioritized Tier 2 gene-region pairs per region/subregion**, ranked by `ML_priority_score`. It is not a plot of all high-confidence genes and it is not a plot of all high + moderate-confidence genes. Duplicate gene symbols appearing in more than one region are collapsed into one heatmap row by the pivot table, with one column per region/subregion.

## Count Audit

| Quantity | Count |
|---|---:|
| Final total ranked Tier 2 gene-region pairs | 671 |
| Final high-confidence count | 35 |
| Final moderate-confidence count | 130 |
| Final high + moderate count | 165 |
| Heatmap gene-region pairs before pivot | 60 |
| Heatmap unique gene symbols / matrix rows | 55 |
| Heatmap high-confidence gene-region pairs | 34 |
| Heatmap moderate-confidence gene-region pairs | 26 |
| Heatmap exploratory/other gene-region pairs | 0 |

Expected extended-strategy final counts match current clean files: **True**.

The heatmap matches its actual coded subset rule: **True**. It does not match the full high/moderate-confidence counts because it was never designed to show all 165 high + moderate-confidence gene-region pairs.

## Agreement With Current Final Counts

The plot agrees with the current ranking table in the sense that its rows come from the current top 10 per region/subregion by `ML_priority_score`. It should not be interpreted as showing all 35 high-confidence genes or all 130 moderate-confidence genes. It is a selected top-stable visualization.

## Recommendation For Main Report

The heatmap can be used in the main report only if the title/caption clearly states the subset. Recommended wording:

**Top stable ML-prioritized Tier 2 genes by region**

or

**Top 10 stable ML-prioritized Tier 2 gene-region pairs per region**

Avoid wording that implies all high-confidence or all high + moderate-confidence genes are shown. Use **high-confidence genes**, not **highly expressed genes**, because confidence here is the ML-priority class and not expression level.

## Corrected Plot

A corrected copy with an explicit title and caption was created here:

`D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended\heatmap_count_audit\MLpriority_model_support_heatmap_corrected.png`

## Files Created

- `heatmap_gene_list.csv`: reconstructed heatmap gene-region pairs and inclusion reason.
- `heatmap_vs_final_counts_summary.csv`: final counts versus heatmap subset counts.
- `heatmap_count_audit_report.md`: this report.
- `MLpriority_model_support_heatmap_corrected.png`: corrected plot with explicit subset title/caption.
