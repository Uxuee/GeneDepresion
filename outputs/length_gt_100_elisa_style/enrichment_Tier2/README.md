# Tier 2 Enrichment Analysis

This folder contains over-representation enrichment analysis for the Length > 100 Tier 2 exploratory DE candidate gene sets.

Analyses include GO Biological Process, GO Molecular Function, GO Cellular Component, KEGG, Reactome when available, and MSigDB collections when available.

Tables are in `tables/`.
Plots are in `plots/`.
Unmapped input and background genes are in `unmapped_genes/`.

Backgrounds were defined as all genes tested in the same region/subregion analysis, using the packaged background universe files from `DE_candidate_gene_sets/`.

To rerun:

```powershell
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" code\run_Tier2_enrichment_length_gt_100.R
```

Tier 2 genes are exploratory DE candidate genes and should not be described as FDR-significant DE genes.
