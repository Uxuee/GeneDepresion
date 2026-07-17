# MDD WGCNA Cluster Visualization Plots

These plots were generated from the saved Step 12 MDD WGCNA outputs. The full WGCNA pipeline was not rerun.

## What Was Created

- Module eigengene correlation heatmaps: show which modules have similar sample-level eigengene patterns.
- Sample PCA plots: show whether samples separate in module-eigengene space by Diagnosis.
- Top-module expression heatmaps: show row-standardized expression of top hub genes across samples.
- Top-module hub networks: show small correlation networks among hub genes, with Tier2/ML overlap genes highlighted when present.
- Module summary bubbleplots: combine module size, Diagnosis evidence, Tier2/ML enrichment, and highlighted AHBA overlap.

## Generated vs Skipped Plots

The plot audit found 32 PNG plots, 32 matching PDF plots, two CSV summary files, and one README. PCA/module-clustering plots were generated successfully. t-SNE was attempted but skipped because the required sklearn functionality failed in this environment, so no t-SNE plots are included.

Generated:

- Module eigengene correlation heatmaps were generated for each processed region.
- PCA/module-eigengene sample plots were generated for each processed region.
- Top-module expression heatmaps were generated for prioritized non-grey modules.
- Module hub-network plots were generated for selected prioritized modules.
- Module summary bubbleplots were generated for each processed region.

Skipped:

- t-SNE plots were attempted but skipped because sklearn failed in this environment.
- No t-SNE image files are included in this folder or listed as generated in the plot audit.

## Files

- Script: `code/create_mdd_wgcna_cluster_visualizations.py`
- Plot inventory: `cluster_visualization_plot_inventory.csv`
- Best modules table: `best_modules_summary_for_plotting.csv`

## Important Limitations

- MEgrey is excluded from the main biological plots because it is the unassigned gene bin.
- UMAP was not available in the current Python environment. t-SNE was attempted as a fallback, but scikit-learn failed with a threadpoolctl error, so only PCA sample projections are included.
- A true cleaned or zoomed WGCNA gene dendrogram cannot be recreated from the saved CSVs alone because the original gene dendrogram/TOM object was not saved. The existing Step6 dendrogram PNGs remain available in the parent `plots/` folder.
- No module eigengene survived BH-FDR < 0.10 for Diagnosis. Raw Diagnosis p < 0.05 should be treated as exploratory.

## Most Presentation-Friendly Plots

1. `module_summary_bubbleplot_Ventral_Striatum.png` for the clearest compact evidence summary.
2. `top_module_expression_heatmap_Ventral_Striatum_MElavenderblush3.png` for the strongest multi-evidence module.
3. `module_network_Ventral_Striatum_MElavenderblush3.png` for a simple hub-gene module view.
4. `sample_PCA_module_eigengenes_*` plots to show that module eigengenes do not produce a clean corrected Diagnosis separation.
5. `module_eigengene_correlation_heatmap_*` plots to explain relationships among discovered modules.

## Region Notes

### Frontal Cortex

Top plotted modules: MEdarkgreen (score 5), MEblack (score 5), MEsalmon (score 4), MEdarkmagenta (score 3), MEblue (score 3).
Diagnosis signal: 3 module(s) had raw p < 0.05, but none survived Diagnosis BH-FDR < 0.10.
Tier2 enrichment: 1 module(s) had Tier2 FDR < 0.05.
High+moderate ML enrichment: 0 module(s) had FDR < 0.10.
AHBA overlap: 23 module(s) had best AHBA overlap FDR < 0.05; highlighted M26/M28/M2/M6/M7 evidence is shown in the bubbleplots and triage table.

### Hippocampal Formation

Top plotted modules: MEgreen (score 4), MEthistle3 (score 3), MEdarkslateblue (score 3), MEmagenta (score 3), MEdarkgreen (score 3).
Diagnosis signal: 0 module(s) had raw p < 0.05, but none survived Diagnosis BH-FDR < 0.10.
Tier2 enrichment: 1 module(s) had Tier2 FDR < 0.05.
High+moderate ML enrichment: 0 module(s) had FDR < 0.10.
AHBA overlap: 52 module(s) had best AHBA overlap FDR < 0.05; highlighted M26/M28/M2/M6/M7 evidence is shown in the bubbleplots and triage table.

### Insula

Top plotted modules: MEorange (score 5), MEdarkorange (score 3), MEbrown (score 3), MElightyellow (score 3), MEmagenta (score 3).
Diagnosis signal: 2 module(s) had raw p < 0.05, but none survived Diagnosis BH-FDR < 0.10.
Tier2 enrichment: 0 module(s) had Tier2 FDR < 0.05.
High+moderate ML enrichment: 0 module(s) had FDR < 0.10.
AHBA overlap: 22 module(s) had best AHBA overlap FDR < 0.05; highlighted M26/M28/M2/M6/M7 evidence is shown in the bubbleplots and triage table.

### Ventral Striatum

Top plotted modules: MElavenderblush3 (score 8), MElightcyan1 (score 6), MEdarkred (score 6), MElightgreen (score 5), MEskyblue2 (score 3).
Diagnosis signal: 3 module(s) had raw p < 0.05, but none survived Diagnosis BH-FDR < 0.10.
Tier2 enrichment: 2 module(s) had Tier2 FDR < 0.05.
High+moderate ML enrichment: 1 module(s) had FDR < 0.10.
AHBA overlap: 33 module(s) had best AHBA overlap FDR < 0.05; highlighted M26/M28/M2/M6/M7 evidence is shown in the bubbleplots and triage table.

## Interpretation

The plots support prioritized exploratory follow-up rather than a clean MDD/control module-level separation. Ventral Striatum modules, especially MElavenderblush3, are the most biologically interesting because they combine raw Diagnosis evidence, Tier2 enrichment, highlighted AHBA overlap, and interpretable hub genes.
