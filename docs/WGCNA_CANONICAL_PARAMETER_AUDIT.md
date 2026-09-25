# WGCNA canonical parameter audit

## Scope and authority

This audit recovers the module-construction method from the executable canonical scripts, not from report prose:

- `code/Step6_MDD_WGCNA_modules.R`
- `code/Step6_Frontal_Subregion_WGCNA_modules.R`

The retained parameter tables and manifests were used to identify the solution actually selected for each region. The existing assignments remain canonical. Phase 3B writes only to `outputs/length_gt_100_clean_pipeline/17_WGCNA_reconstruction_validation/` and is a noncanonical validation candidate.

## Exact construction algorithm

1. Read the normalized matrix with genes in rows and samples in columns, transpose to samples × genes, coerce expression to numeric, remove missing/blank gene identifiers, and retain the first row when a gene identifier is duplicated.
2. Match expression samples to metadata by sample ID. Regions with fewer than 15 matched samples are not eligible.
3. Create an unadjusted `log2(GeTMM + 1)` branch. The historical workflow also created a residualized sensitivity branch when usable technical covariates were available. That branch used `model.matrix` and `lm.fit` with the available non-diagnosis fields among `BrainBank`, `SeqBatch`, `RIN`, `rin`, `Age`, `age`, `Sex`, and `gender`. Diagnosis was explicitly not regressed out.
4. Remove genes containing any missing expression value. Compute base-R `mad(..., na.rm=TRUE)` and `var(..., na.rm=TRUE)`, retain genes with finite MAD and variance and variance greater than zero, then apply `WGCNA::goodSamplesGenes`. Rank genes by decreasing MAD. Broad-region candidates were the top 5,000, 10,000, and 15,000 genes; an all-usable branch was added only when no more than 20,000 genes remained. Frontal-subregion candidates used the top 5,000 genes.
5. Test powers 1–20 with `WGCNA::pickSoftThreshold`, `networkType="signed"`, and bicor using `maxPOutliers=0.05` and `use="pairwise.complete.obs"`. Select the first power whose signed scale-free fit is at least 0.80. If no power reaches 0.80, select the power with the largest signed R², resolving a tie in favor of the smaller power. If bicor fails, retry with Pearson correlation and pairwise-complete observations.
6. Build a signed adjacency with the selected power and the selected correlation method. Calculate a signed TOM using `WGCNA::TOMsimilarity(..., TOMType="signed")`; define dissimilarity as `1 - TOM`; cluster genes using average-linkage hierarchical clustering.
7. Evaluate all 18 combinations of `deepSplit={2,3,4}`, `minModuleSize={20,30,50}`, and `mergeCutHeight={0.15,0.25}`. Dynamic cutting used `dynamicTreeCut::cutreeDynamic`, the TOM dissimilarity, `pamRespectsDendro=FALSE`, and the grid-specific minimum cluster size. Close modules were merged with `WGCNA::mergeCloseModules` at the grid-specific cut height.
8. Mark a candidate usable when it has at least three non-grey modules, its largest non-grey module contains less than 70% of genes, its grey fraction is below 50%, and at least 15 samples are present.
9. Rank candidates by: usable before unusable; larger retained gene universe; more non-grey modules; then smaller largest-module fraction. The first row is the saved best solution. No diagnosis statistic, candidate-gene enrichment, ML result, or functional annotation enters this selection.
10. Name modules from WGCNA color labels, with `MEgrey` reserved for unassigned genes. Module color names are arbitrary labels and cannot be compared across runs without gene-overlap alignment.
11. Calculate module eigengenes with `WGCNA::moduleEigengenes` and order them with `WGCNA::orderMEs`.
12. Calculate kME in the saved workflow as the Pearson correlation between each module gene's expression profile and the corresponding eigengene, using pairwise-complete observations. Rank hubs by decreasing `abs(kME)`, then by decreasing absolute gene–Diagnosis correlation; retain the first 50 per non-grey module.

## Selected historical solutions

| Region | Samples | Selected branch | MAD genes | Power | Power rule | Correlation | deepSplit | minModuleSize | mergeCutHeight | Non-grey modules | Grey fraction |
|---|---:|---|---:|---:|---|---|---:|---:|---:|---:|---:|
| Ventral Striatum | 48 | residualized technical covariates | 15,000 | 11 | first R² ≥ 0.80 | bicor | 4 | 20 | 0.15 | 50 | 0 |
| Insula | 48 | unadjusted | 15,000 | 14 | first R² ≥ 0.80 | bicor | 4 | 20 | 0.15 | 29 | 0 |
| Hippocampal Formation | 43 | unadjusted | 15,000 | 10 | first R² ≥ 0.80 | bicor | 4 | 20 | 0.15 | 69 | 0 |
| CG25 / Cingulate gyrus 25 | 27 | unadjusted | 5,000 | 17 | maximum R² because none reached 0.80 | bicor | 4 | 20 | 0.15 | 27 | 0.0024 |
| dlPFC BA8/9 | 48 | residualized technical covariates | 5,000 | 15 | first R² ≥ 0.80 | bicor | 3 | 20 | 0.15 | 21 | 0.0058 |
| OFC BA11 | 47 | residualized technical covariates | 5,000 | 9 | first R² ≥ 0.80 | bicor | 4 | 20 | 0.15 | 18 | 0 |

The pooled broad Frontal Cortex network is excluded from the six-unit main validation. It remains a historically separate supplementary analysis and is not interchangeable with dlPFC BA8/9, OFC BA11, or CG25.

## Diagnosis-model implementation found in code

Module construction is unsupervised: Diagnosis is not used in MAD filtering, correlation, adjacency, TOM construction, tree cutting, merging, or solution selection.

After construction, the saved scripts tested eigengenes with `lm`. The executable code built the formula from all usable fields among `Age`, `age`, `Sex`, `gender`, `BrainBank`, `SeqBatch`, `RIN`, `rin`, `PMI`, and `pmi`, after `Diagnosis_numeric`. Thus, where both normalized and lowercase aliases existed, both names appeared in the formula even though they represented duplicate information. If sample size was below `5 + 2 × number_of_covariate_fields`, the script fell back to `ME ~ Diagnosis_numeric`. This is more literal than the simplified report wording `ME ~ Diagnosis + Age + Sex + RIN + PMI`; Phase 3B reproduces the executable implementation for validation and records every fitted formula. BH correction was applied within region across non-grey module Diagnosis p-values.

## Historical versus Phase-3B environment

The retained historical manifest records R 4.4.2, WGCNA 1.74, dynamicTreeCut 1.63-1, limma 3.62.2, readr 2.1.5, dplyr 1.1.4, tidyr 1.3.1, and stringr 1.5.1. Other loaded packages are listed in the historical manifest. A complete historical lockfile was not retained.

Phase 3B uses a newly controlled R 4.6.1 library. Its actual versions are recorded in `environment/phase3b/package_versions.csv`, `environment/phase3b/sessionInfo.txt`, and the real `environment/phase3b/renv.lock`. This environment is a reconstruction environment; it is not represented as the unknown historical environment.

## Equivalence decision fixed before enrichment

- `EXACT`: ARI = 1 and 100% gene-assignment agreement after label-invariant mapping.
- `EFFECTIVELY_EQUIVALENT`: ARI ≥ 0.90, median best-module Jaccard ≥ 0.90, mapped gene-assignment agreement ≥ 0.90, and median absolute aligned eigengene correlation ≥ 0.95.
- `MATERIAL_DIFFERENCE`: a completed reconstruction that does not satisfy those thresholds.
- `FAILED_TO_RECONSTRUCT`: the independent run cannot be completed from frozen inputs without guessing.

These thresholds are encoded in `module_comparison/predeclared_equivalence_thresholds.json` before functional or ML enrichment results are calculated.
