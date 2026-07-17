# Notebook: Clear Length > 100 MDD Workflow

This notebook-style document walks through the current clear workflow in:

`D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline`

It is written as a readable analysis notebook, not as a rerun log.

## 1. Goal

The goal is to keep one clear MDD/control transcriptomics workflow:

- one normalization path: GeTMM followed by `log2(GeTMM + 1)`
- one main candidate-gene path: Tier 2 exploratory evidence
- one main ML interpretation: ML-prioritized Tier 2 genes
- clear follow-up branches for AHBA contextualization and MDD WGCNA modules

## 2. Input Folder

The clean workflow root is:

```text
D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline
```

The manifest is:

```text
D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\pipeline_manifest.json
```

The manifest records the scripts, reused input files, final report path, ML strategy, and final output counts.

## 3. Normalized Expression

Expression matrices are here:

```text
outputs\length_gt_100_clean_pipeline\02_normalized_expression
```

They are expected to be genes by samples:

```text
rows = genes
columns = samples
values = log2(GeTMM + 1)
```

These matrices are the correct inputs for module-level work such as WGCNA. WGCNA should use many expressed/informative genes, not only Tier 2 genes.

## 4. Statistical Modeling

Broad-region statistical results are here:

```text
outputs\length_gt_100_clean_pipeline\03_broad_region_LMM
```

Standardized effect sizes are here:

```text
outputs\length_gt_100_clean_pipeline\04_standardized_effect_size
```

The standardized diagnosis effect is:

```text
beta_diagnosis / SD(expression of gene)
```

This helps compare effects across genes.

## 5. Evidence Tiers

Evidence tiers are here:

```text
outputs\length_gt_100_clean_pipeline\05_evidence_tiers
```

Definitions:

- Tier 1: Storey q < 0.10 and abs standardized effect size > 0.8
- Tier 2: p_value_diagnosis < 0.01 and abs standardized effect size > 0.8
- Tier 3: p_value_diagnosis < 0.01 and abs beta_diagnosis > 0.2
- Tier 4: ranked support layer, not a final separate biological gene list

Broad-region summary:

| Brain region | Genes tested | Tier 1 | Tier 2 | Tier 3 |
| --- | ---: | ---: | ---: | ---: |
| Frontal_Cortex | 23,941 | 0 | 1 | 99 |
| Hippocampal_Formation | 24,371 | 0 | 127 | 70 |
| Insula | 23,090 | 0 | 70 | 92 |
| Ventral_Striatum | 23,917 | 4 | 146 | 167 |

## 6. Frontal Subregion Sensitivity

Frontal subregion outputs are here:

```text
outputs\length_gt_100_clean_pipeline\06_frontal_subregion_sensitivity
```

This step is used because broad Frontal Cortex combines OFC, dlPFC, and Cg25. The frontal subregion sensitivity model helps separate these anatomical subregions.

## 7. Final Tier 2 Candidate Gene Sets

Final Tier 2 candidate gene sets are here:

```text
outputs\length_gt_100_clean_pipeline\07_Tier2_candidate_gene_sets
```

Important table:

```text
Tier2_candidate_genes_all_regions.csv
```

Final candidate counts:

| Region or subregion | Tier 2 genes | Up in MDD | Down in MDD |
| --- | ---: | ---: | ---: |
| Cingulate gyrus 25 | 186 | 137 | 49 |
| Hippocampal Formation | 127 | 80 | 47 |
| Insula | 70 | 37 | 33 |
| OFC BA11 | 70 | 38 | 32 |
| Ventral Striatum | 146 | 76 | 70 |
| dlPFC BA8/9 | 72 | 24 | 48 |

## 8. Tier 2 Enrichment

Tier 2 enrichment outputs are here:

```text
outputs\length_gt_100_clean_pipeline\08_Tier2_enrichment
```

This step asks whether the exploratory Tier 2 genes are overrepresented in known pathways or gene sets. It is biological interpretation, not independent validation.

## 9. ML v1 Prioritization

Main ML outputs are here:

```text
outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended
```

The main final table is:

```text
final_ranked_ML_prioritized_Tier2_genes.csv
```

The final ranking contains:

- 671 ranked Tier 2 genes
- 35 high-confidence ML-prioritized Tier 2 genes
- 130 moderate-confidence ML-prioritized Tier 2 genes

Interpretation:

ML is used to prioritize genes already supported by the transcriptomics evidence. It is not external validation and should not be presented as a validated diagnostic classifier.

## 10. Report-Ready Images

Figures are collected here:

```text
outputs\length_gt_100_clean_pipeline\images
```

Useful subfolders:

- `images\tier2`
- `images\enrichment`
- `images\ml_v1`
- `images\ml_prioritized`

These are the figures to use in a report, presentation, or GitHub README.

## 11. AHBA Contextualization

Published AHBA module work is here:

```text
outputs\length_gt_100_clean_pipeline\11_AHBA_modules_enrichment
```

This branch asks how MDD candidate genes relate to published healthy adult human brain coexpression modules.

Do not overwrite these outputs when running MDD-derived WGCNA modules.

## 12. MDD WGCNA Modules

MDD-derived WGCNA work is here:

```text
outputs\length_gt_100_clean_pipeline\12_MDD_WGCNA_modules
```

Script:

```text
code\Step6_MDD_WGCNA_modules.R
```

Important method rule:

```text
Build WGCNA modules from normalized expression data first.
Then test Tier 2 / ML-prioritized gene enrichment inside modules.
```

The current saved WGCNA output is a limited validation run, not the final full grid.

## 13. How to Explain the Workflow in One Paragraph

This project analyzes MDD/control RNA-seq data using a clean Length > 100 gene-filtered workflow. Expression is normalized with GeTMM and transformed as `log2(GeTMM + 1)`. Statistical models estimate MDD diagnosis effects by brain region and frontal subregion. Genes are assigned exploratory evidence tiers, with Tier 2 used as the main candidate-gene pool. Enrichment analyses provide biological context, and ML v1 is used only to prioritize Tier 2 genes by internal feature-support stability. AHBA and MDD WGCNA analyses are follow-up module-level context, not replacements for the main candidate-gene workflow.

## 14. What Not to Overclaim

Do not say:

- these are validated biomarkers
- ML proves diagnostic prediction
- Tier 2 genes are FDR-significant discoveries
- WGCNA modules should be built only from Tier 2 genes

Safer wording:

- exploratory candidate genes
- internally prioritized Tier 2 genes
- module-level context
- follow-up hypotheses for validation
