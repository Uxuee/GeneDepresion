# Step 3 Length-Cutoff Sensitivity Report

This report summarizes the full Step 3 sensitivity run across the gene-length
cutoffs created by:

```text
code/Step1_length_cutoff_sensitivity_GeTMM.R
```

The full Step 3 script was:

```text
code/Step3_run_all_DGE_methods_on_length_cutoffs.R
```

Local output folder:

```text
outputs/step3_length_sensitivity
```

Generated output files are not committed to GitHub.

## Main Transcriptomics Methods

These are the primary methods to interpret:

- `dream` for Frontal Cortex
- `limma_voom` for Insula, Ventral Striatum, and Hippocampal Formation

| Length cutoff | Brain area | Method | Genes tested | Raw p < 0.05 | Raw p < 0.01 | FDR < 0.05 | FDR < 0.10 | Minimum FDR | Top gene |
|---:|---|---|---:|---:|---:|---:|---:|---:|---|
| 20 | Frontal_Cortex | dream | 24,286 | 1,490 | 282 | 0 | 0 | 0.12651821 | 5147 |
| 50 | Frontal_Cortex | dream | 24,286 | 1,497 | 288 | 0 | 1 | 0.07904222 | 5147 |
| 100 | Frontal_Cortex | dream | 23,941 | 1,458 | 277 | 0 | 0 | 0.12709347 | 5147 |
| 200 | Frontal_Cortex | dream | 23,739 | 1,447 | 275 | 0 | 0 | 0.12248101 | 5147 |
| 500 | Frontal_Cortex | dream | 23,308 | 1,399 | 265 | 0 | 1 | 0.08908718 | 5147 |
| 20 | Insula | limma_voom | 23,375 | 941 | 193 | 0 | 0 | 0.83519450 | 10090 |
| 50 | Insula | limma_voom | 23,374 | 941 | 193 | 0 | 0 | 0.83511097 | 10090 |
| 100 | Insula | limma_voom | 23,090 | 943 | 194 | 0 | 0 | 0.84302246 | 10090 |
| 200 | Insula | limma_voom | 22,924 | 939 | 197 | 0 | 0 | 0.73597461 | 10090 |
| 500 | Insula | limma_voom | 22,541 | 931 | 188 | 0 | 0 | 0.74124317 | 10090 |
| 20 | Ventral_Striatum | limma_voom | 24,275 | 2,472 | 559 | 0 | 2 | 0.09301074 | 107987147 |
| 50 | Ventral_Striatum | limma_voom | 24,273 | 2,473 | 559 | 0 | 2 | 0.09304780 | 107987147 |
| 100 | Ventral_Striatum | limma_voom | 23,917 | 2,423 | 533 | 0 | 2 | 0.08157989 | 107987147 |
| 200 | Ventral_Striatum | limma_voom | 23,722 | 2,352 | 516 | 0 | 3 | 0.08422257 | 107987147 |
| 500 | Ventral_Striatum | limma_voom | 23,268 | 2,375 | 537 | 0 | 2 | 0.08684684 | 107987147 |
| 20 | Hippocampal_Formation | limma_voom | 24,700 | 800 | 117 | 0 | 0 | 0.99995020 | 65083 |
| 50 | Hippocampal_Formation | limma_voom | 24,698 | 800 | 117 | 0 | 0 | 0.99996464 | 65083 |
| 100 | Hippocampal_Formation | limma_voom | 24,371 | 795 | 117 | 0 | 0 | 0.99992971 | 65083 |
| 200 | Hippocampal_Formation | limma_voom | 24,182 | 786 | 112 | 0 | 0 | 0.99914675 | 65083 |
| 500 | Hippocampal_Formation | limma_voom | 23,726 | 774 | 111 | 0 | 0 | 0.99960647 | 65083 |

## Top-Gene Stability

This compares each cutoff's top 100 genes by raw p-value with the same method's
baseline `Length > 20` top 100 genes.

| Brain area | Method | Length cutoff | Top 100 overlap vs Length > 20 |
|---|---|---:|---:|
| Frontal_Cortex | dream | 50 | 96% |
| Frontal_Cortex | dream | 100 | 99% |
| Frontal_Cortex | dream | 200 | 97% |
| Frontal_Cortex | dream | 500 | 93% |
| Insula | limma_voom | 50 | 100% |
| Insula | limma_voom | 100 | 98% |
| Insula | limma_voom | 200 | 96% |
| Insula | limma_voom | 500 | 97% |
| Ventral_Striatum | limma_voom | 50 | 100% |
| Ventral_Striatum | limma_voom | 100 | 97% |
| Ventral_Striatum | limma_voom | 200 | 91% |
| Ventral_Striatum | limma_voom | 500 | 97% |
| Hippocampal_Formation | limma_voom | 50 | 100% |
| Hippocampal_Formation | limma_voom | 100 | 100% |
| Hippocampal_Formation | limma_voom | 200 | 93% |
| Hippocampal_Formation | limma_voom | 500 | 96% |

The top genes are very stable across cutoffs.

## Comparator Models

The normalized-expression comparator models were run as diagnostics.

Important caution:

- Frontal fixed-effect LM produced many FDR-significant genes.
- That model ignores donor random effects, so it should not be treated as the
  main Frontal Cortex result.
- Frontal LMM gave a small number of exploratory FDR hits, but it also produced
  many singular fits.

Frontal LMM singular-fit counts:

| Length cutoff | Singular fits |
|---:|---:|
| 20 | 1,133 |
| 50 | 1,132 |
| 100 | 1,017 |
| 200 | 923 |
| 500 | 885 |

Effect-size correlations between transcriptomics methods and comparator models
were high:

- Insula limma-voom vs LM: about 0.90 to 0.93
- Ventral Striatum limma-voom vs LM: about 0.90 to 0.93
- Hippocampal Formation limma-voom vs LM: about 0.88 to 0.91
- Frontal Cortex dream vs LMM: about 0.90 to 0.94

This means the direction and relative size of diagnosis effects are broadly
consistent across model families, even though the FDR calls differ.

## Interpretation

The stricter length filters do reduce the number of genes tested, but they do
not create a clear, strong improvement in FDR for the main transcriptomics
models.

For the main models:

- Frontal Cortex dream keeps the same top gene, `5147`, across all cutoffs.
- Ventral Striatum limma-voom keeps the same top gene, `107987147`, across all
  cutoffs.
- Insula keeps the same top gene, `10090`, across all cutoffs.
- Hippocampal Formation keeps the same top gene, `65083`, across all cutoffs.
- Ventral Striatum remains the only region with consistent exploratory
  FDR < 0.10 findings.

There is no evidence here that `Length > 500` gives a biologically cleaner or
more convincing result than the baseline. It removes genes, but the top genes
remain similar and the FDR pattern does not substantially improve.

The most conservative interpretation is:

- The baseline `Length > 20` does not appear to be causing an obvious collapse
  of the result.
- `Length > 100` or `Length > 200` could be reasonable sensitivity cutoffs to
  mention, because they remove more short genes while preserving strong top-gene
  stability.
- Do not choose a cutoff only because it gives one or two more exploratory FDR
  hits.

## Key Output CSV Files

The main local summary files are:

```text
outputs/step3_length_sensitivity/Step3_length_cutoff_all_methods_summary.csv
outputs/step3_length_sensitivity/Step3_top_gene_stability_vs_Length_gt_20_all_methods.csv
outputs/step3_length_sensitivity/Step3_cross_method_top100_overlap.csv
outputs/step3_length_sensitivity/Step3_cross_method_effect_size_correlation.csv
```
