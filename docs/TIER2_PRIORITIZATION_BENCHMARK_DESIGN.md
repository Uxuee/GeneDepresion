> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Prospective design: Tier2 gene-prioritisation benchmark

## Frozen inputs and comparison boundary

This benchmark starts from frozen commit `df3596ed208bb7472e166e5fc0330b69453a9e31`. The primary Tier2–ML–RRA result in `18_Tier2_ML_RRA_prioritization` is read-only. The benchmark will not retrain, rewrite, or tune that ranking.

The sole candidate source is `outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv` (SHA-256 `1672ddf95a4f4c9e44ae42e7accc208c311e2e1da8533e1808dcf91bbf7b67a1`). Only region, stable gene ID, and gene symbol are used. It contains 671 unique region–gene pairs: Ventral Striatum 146, Insula 70, Hippocampal Formation 127, CG25 186, dlPFC BA8/9 72, and OFC BA11 70. These counts are validation expectations, not a replacement for the source table.

The six normalized expression matrices, sample identities, and Diagnosis labels are exactly those recorded in the frozen primary manifest: 48, 48, 43, 27, 48, and 47 samples respectively. No pooled Frontal Cortex, Tier4, historical ML score/label, WGCNA, AHBA, or all-gene rank enters benchmark fitting.

## Frozen primary method

The primary analysis used five stratified outer folds repeated with seeds 2026–2030 (25 splits per region/model), training-only median imputation, training-only scaling for linear models, three-fold inner tuning by balanced accuracy, and five models: L1 logistic, Elastic Net logistic, L2 logistic, linear SVM, and Random Forest. Two-level RRA used R 4.6.1 and `RobustRankAggreg` 1.2.1: 25 fold lists were aggregated within each region/model, then the five model lists were aggregated per region. Predictive performance and coefficient direction did not weight RRA.

## Benchmark A: LASSO alone

The frozen primary run already contains the exact requested L1-logistic fits on the same 25 outer splits, with the same training-only preprocessing, inner tuning, model grid (`C = 0.01, 0.1, 1`), tolerance, and samples. Reusing those immutable event-level records is mathematically identical to rerunning LASSO and avoids needless duplicate stochastic computation.

For each region, genes are ordered by:

1. descending selection frequency across exactly 25 outer fits;
2. ascending mean rank among fits where selected;
3. descending mean absolute coefficient among fits where selected;
4. stable gene ID.

The output also records selection count, median selected rank, coefficient sign consistency, and descriptive outer-fold metrics. Never-selected genes receive frequency zero and a shared worst-rank block for rank correlations; stable gene ID is used only to make exported row order deterministic, not to imply evidence differences.

## Benchmark B: classical subsampling Stability Selection

### Package audit

The project runtime contains scikit-learn 1.7.2 but no installed Python `stability-selection` package. The isolated R library contains only `RobustRankAggreg` 1.2.1 and no `stabs` installation. Adding an unpinned package would weaken the frozen environment. Therefore this phase uses a transparent implementation of the classical repeated half-sampling procedure with sparse logistic regression, rather than substituting a different algorithm.

### Prespecified variant

- 100 deterministic, stratified half-samples per region.
- Seeds 31001–31100; the region is included in deterministic seed derivation.
- Each half-sample contains `floor(n_class / 2)` observations from each Diagnosis class.
- Within every half-sample, median imputation and standardization are learned only from that half-sample.
- L1 logistic regularisation is selected by three-fold stratified CV within the half-sample, optimizing balanced accuracy over the same `C = 0.01, 0.1, 1` grid.
- A gene is selected when `abs(coef) > 1e-12`, matching the validated Tier2 pipeline.
- The primary ranking is descending selection probability (`selected subsamples / 100`), then descending mean absolute coefficient when selected, then stable gene ID.
- Never-selected genes receive probability zero and a shared worst-rank block for correlations.
- A selection-probability threshold of 0.60 is reported only as a predeclared descriptive reference. The continuous probability ranking is primary; the threshold will not be tuned from biological results.

This design requires 600 final sparse fits plus 5,400 inner-CV evaluation fits. It is deterministic and retains all zero-selection events explicitly.

## Fair comparison and prespecified interpretation

All three methods use the identical 671-pair Tier2 universe. Pairwise comparisons use top 10, top 25, and top 20% (`ceil(0.20 × regional universe)`). Top-set metrics report intersection, Jaccard, and overlap relative to each list. Spearman and Kendall correlations use all regional candidates; never-selected sparse-method genes share the worst rank.

Resampling stability uses the 25 frozen LASSO outer lists, 100 Stability Selection subsamples, and 25 split-level five-model RRA consensus lists reconstructed only for stability diagnostics with the same RRA package. Sparse lists are never padded for Jaccard calculations. Rank correlations use the shared-worst-rank rule.

The report will apply the five interpretation scenarios supplied before analysis. It will not create a composite winner score. Classifier performance is descriptive only because Tier2 was defined using Diagnosis in the complete cohort. Agreement is within-study convergence, not external validation or proof of biological truth.
