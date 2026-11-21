> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Upstream statistical-method benchmark: frozen design

## Scientific question

This benchmark asks whether the **already frozen** Tier2 candidates are broadly supported when the same RNA-seq samples and eligible genes are analysed with limma-voom and DESeq2. It is a within-study method-concordance analysis, not independent replication and not a procedure for creating a replacement Tier2.

## Audited reference analysis

The frozen analysis uses six original regions only: Ventral Striatum, Insula, Hippocampal Formation, CG25, dlPFC BA8/9, and OFC BA11. Pooled Frontal Cortex is excluded from the benchmark conclusions.

For Ventral Striatum, Insula, and Hippocampal Formation, the frozen model is ordinary least squares on `log2(GeTMM + 1)`:

`expression_value ~ Diagnosis + Age + Sex + RIN`

The original code sets `Diagnosis` levels to `CTRL`, then `MDD`; therefore `DiagnosisMDD` and a positive diagnosis beta mean higher expression in MDD. Samples missing any model covariate are removed before all gene fits. BH adjustment is performed separately across every tested gene in each region.

For CG25, dlPFC BA8/9, and OFC BA11, the frozen frontal-subregion sensitivity model uses the same fixed-effect formula after subsetting the broad frontal matrix by original anatomical region. Each donor appears once within a subregion, so a donor random intercept is not estimable. The diagnosis coefficient has the same MDD-minus-Control interpretation, and BH adjustment is regional.

The standardized diagnosis effect is the frozen diagnosis beta divided by the gene-wise sample standard deviation of `log2(GeTMM + 1)` expression in the applicable region/subregion. It is not a residual-standardized effect.

## Audited normalization and filtering path

The canonical path is:

1. raw integer NCBI count matrix;
2. duplicated GeneIDs summed;
3. NCBI GRCh38.p13 gene-length annotation matched;
4. gene length greater than 100 bp retained;
5. within-analysis-region expression filter applied;
6. reads per kilobase calculated;
7. TMM normalization factors calculated;
8. GeTMM produced;
9. `log2(GeTMM + 1)` used by the frozen LM.

The benchmark methods will instead receive the matching **raw integer counts** after the same frozen length/expression eligibility filter. GeTMM, RPK, and transformed values will not enter limma-voom or DESeq2.

## Frozen reference sources

- Broad-region all-gene LM and standardized effects: `D:/genes/GeneDepresion/outputs/length_gt_100_clean_pipeline/04_standardized_effect_size/LMM_all_genes_with_standardized_effect_size_4groups.csv`.
- Frontal-subregion all-gene LM and standardized effects: immutable Git blob `f700cb323a47b24825b11bc1bba33cc07e6430bb` from historical commit `ac37ec9` in `D:/genes/GeneDepresion`.
- Frozen Tier2: `outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv`.
- Raw counts: `D:/genes/genedata/GSE102556_raw_counts_GRCh38.p13_NCBI.tsv` and the exact frozen region-filtered count matrices under `D:/genes/GeneDepresion/outputs/step1_length_sensitivity/Length_gt_100/03_filtered_raw_counts/keep_all_samples/`.
- Metadata: the corresponding frozen regional metadata plus `metadata_clean_matched_all_samples_with_brain_area.csv` under the same Step1 output root.
- Annotation: `D:/genes/genedata/Human.GRCh38.p13.annot.tsv`.

## Validation targets before fitting

The matched comparison must contain exactly the same complete-case samples and eligible genes as the frozen LM:

| Region | Expected samples | Expected eligible genes | Expected frozen Tier2 pairs |
|---|---:|---:|---:|
| Ventral Striatum | 48 | 23,917 | 146 |
| Insula | 48 | 23,090 | 70 |
| Hippocampal Formation | 43 | 24,371 | 127 |
| CG25 | 27 | 23,941 | 186 |
| dlPFC BA8/9 | 48 | 23,941 | 72 |
| OFC BA11 | 47 | 23,941 | 70 |

If sample identities, covariates, raw-count integrality, or gene identities fail to match, fitting must stop.

## Frozen benchmark specifications

### limma-voom

For each region: raw counts → `edgeR::DGEList` → `edgeR::calcNormFactors(method="TMM")` → `limma::voom` → `limma::lmFit` → `limma::eBayes`. The design is `~ Diagnosis + Age + Sex + RIN`, with `CTRL` as reference. No quality weights or post-result tuning will be used.

### DESeq2

For each region: raw integer counts → `DESeqDataSetFromMatrix` with design `~ Age + Sex + RIN + Diagnosis` → standard `DESeq` fitting → explicit `Diagnosis: MDD versus CTRL` contrast. Native `padj` will be retained, and a separate BH adjustment of raw p-values across the full shared eligible regional universe will be stored.

The documented denominator for the explicit comparison is the same eligible tested gene universe used by the frozen LM and limma-voom; DESeq2's method-native independent-filtering `padj` remains a separate field.

## Prospective comparison definitions

- **Primary genome-wide rank:** two-sided raw-p rank, ascending raw p with GeneID as deterministic tie-breaker.
- **Signed evidence sensitivity statistic:** `sign(effect) * -log10(raw p)`.
- Direction agreement is based only on coefficient sign.
- Pairwise Spearman and Kendall correlations use primary regional ranks.
- Top-set overlap is predeclared at top 10, 25, 50, 1%, and 5%.
- Tier2 enrichment is predeclared at benchmark top 1%, 5%, and 10%, with BH correction across region × method × threshold tests.
- No RRA, consensus tier, weighted score, or replacement Tier2 will be created.

## Interpretation boundaries

All methods use the same cohort, and Tier2 was selected with the frozen LM. Agreement is therefore within-study statistical-method concordance. Effect magnitudes are not directly interchangeable across LM, voom, and DESeq2 scales. Cross-method agreement cannot turn a non-FDR result into a significant disease association, provide external validation, or establish a biomarker.
