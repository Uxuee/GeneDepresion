# MDD WGCNA Module Report

Generated: 2026-07-16 19:15:37 KST

Run scope: Limited validation run. Active STEP6 controls: STEP6_MAX_REGIONS=1; STEP6_GENE_FILTER_TARGETS=5000; STEP6_PARAM_GRID_LIMIT=1; STEP6_VARIANTS=unadjusted_only

## What was done

This analysis built WGCNA coexpression modules directly from the normalized MDD/control GeTMM-log2 expression matrices.

The modules were not built from DEG-only genes, Tier 2-only genes, or ML-prioritized-only genes. Genes were filtered by expression quality and variability/MAD, then existing Tier 2 and ML-prioritized gene sets were tested for enrichment after module construction.

Signed WGCNA networks were used as the main analysis. The script tried robust bicor correlation first and used Pearson only if bicor failed.

## Sample-size and usability overview

- Regions processed: Frontal Cortex
- Parameter-grid runs recorded: 1
- Usable parameter-grid runs: 1
- Regions with at least one usable module solution: Frontal Cortex

A solution was marked usable when it had at least 3 non-grey modules, the largest module was below 70% of genes, grey genes were below 50%, and the region had at least 15 samples.

## Diagnosis-associated modules

No module eigengene reached BH-FDR < 0.10 for Diagnosis in the saved best solutions.

## Tier 2 and ML-prioritized enrichment

- Frontal Cortex MEred / all_Tier2_genes: overlap=13, p=0.0031, FDR=0.017
- Frontal Cortex MEgrey / all_Tier2_genes: overlap=21, p=0.00876, FDR=0.0321
- Frontal Cortex MEgrey / high_confidence_ML_genes: overlap=4, p=0.00318, FDR=0.0349
- Frontal Cortex MEgrey / high_plus_moderate_ML_genes: overlap=10, p=7.13e-05, FDR=0.000785
- Frontal Cortex MEyellow / all_Tier2_genes: overlap=20, p=0.00087, FDR=0.00957

## Published AHBA module overlap

The AHBA overlap step tests each MDD-derived module against published AHBA modules and highlights M26, M28, M2, M6, and M7.

- Frontal Cortex MEturquoise vs M26: overlap=25, p=1.56e-06, FDR=3.93e-05
- Frontal Cortex MEturquoise vs M28: overlap=18, p=0.000846, FDR=0.011
- Frontal Cortex MEbrown vs M26: overlap=30, p=6.41e-07, FDR=1.88e-05
- Frontal Cortex MEblue vs M6: overlap=106, p=1.48e-18, FDR=6.51e-17
- Frontal Cortex MEblack vs M7: overlap=10, p=0.000959, FDR=0.0121
- Frontal Cortex MEmagenta vs M2: overlap=5, p=0.00216, FDR=0.0245
- Frontal Cortex MEmagenta vs M6: overlap=12, p=0.000105, FDR=0.00194

## Limitations

- This is exploratory and should be interpreted as module discovery, not validation.
- Some regions have modest sample sizes, so module stability may be limited.
- Residualized analyses are sensitivity analyses and may use Age/Sex if those covariates are available.
- No external validation cohort has been used yet.

## Main output files

- `mdd_wgcna_parameter_sensitivity.csv`
- `best_solution_gene_to_module_assignments.csv`
- `best_solution_module_size_table.csv`
- `best_solution_module_eigengenes.csv`
- `module_trait_association_results.csv`
- `module_enrichment_for_Tier2_and_ML_genes.csv`
- `mdd_modules_vs_published_AHBA_modules_enrichment.csv`
- `module_hub_genes_top50.csv`
- `mdd_wgcna_manifest.json`
