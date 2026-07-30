# MDD modules main-report figure placement — original groups only

Use only these analysis units in the main report: Ventral Striatum, Insula, Hippocampal Formation, dlPFC BA8/9, OFC BA11, and CG25 (Cingulate gyrus 25). Do not use broad pooled Frontal Cortex panels or rows.

## Recommended figure order

1. **Module construction workflow.** Use the existing deterministic WGCNA workflow schematic. The workflow is methodological and is not a pooled-region result.
2. **Module colors and module counts.** Use region/subregion-specific dendrogram and module-size panels for the six retained units. Do not include the pooled Frontal Cortex dendrogram or counts.
3. **Candidate-gene enrichment inside MDD-derived modules.** Use Tier2/ML enrichment panels or tables filtered to the six retained units.
4. **Main MDD module results by original region/subregion.** Present results for Ventral Striatum, Insula, Hippocampal Formation, dlPFC BA8/9, OFC BA11, and CG25 only.
   
   Use the final 10-module prioritized table as the main result. Do **not** place `frontal_subregion_module_summary_bubble_plot.png`, its PDF, or the later legend-only variant. Replace that former figure slot with the report section **“Module prioritization strategy”**, which explains the integrated evidence workflow from unsupervised WGCNA through eigengene testing, enrichment, annotation, module-level ML, and final exploratory selection.
5. **Biological interpretation of prioritized MDD-derived modules.** Place the corrected figures in this order:
   - `biological_category_heatmap_prioritized_modules_original_groups_only.png`
   - `top_GO_terms_dotplot_by_module_original_groups_only.png`
   - `hub_gene_annotation_plot_top_modules_original_groups_only.png`
   - `cell_type_enrichment_dotplot_original_groups_only.png`
   
   Use the matching PDFs for manuscript assembly. The cell-type figure contains only modules with valid marker-enrichment rows; absence of a panel is not evidence of no cell-type biology.
6. **Network/cluster visualizations of prioritized modules.** Use saved module-network panels for retained original units only.
7. **Module-level ML summary.** Use only a figure/table whose rows have been explicitly filtered to the six retained units; do not place an existing pooled Frontal Cortex panel.
8. **External validation plan.** Present the validation workflow/plan without treating pooled Frontal Cortex as a main result.

## Caption wording

“Prioritized MDD-derived modules from the original region/subregion analyses. Broad pooled Frontal Cortex was excluded; frontal results are shown separately for dlPFC BA8/9, OFC BA11, and CG25.”
