> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Tier2 ML conditional-null pilot report

## 1. Scientific question

Given the already-fixed Tier2 candidate universe, does the original MDD/Control labeling produce more coherent and reproducible downstream ML gene prioritisation than labels shuffled within each region?

This 20-permutation pilot completed successfully for all six analysis units. It reran the full five-model Tier2-RRA pipeline and the exact LASSO-alone benchmark under each permuted label vector. No observed fit, preprocessing object, hyperparameter, coefficient, importance, ranking, or RRA result was reused.

## 2. Why the null is conditional

**This is not an end-to-end global-null test.** Tier2 was selected using the original Diagnosis labels in the full cohort, so the fixed 671-pair universe is already enriched for observed-label association. The experiment asks whether downstream ML remains structured with respect to those original labels after the candidate set is fixed. It cannot answer whether the full expression → gene-wise model → Tier2 → ML pipeline would produce the result under a global null.

This design is therefore expected to favor the observed labels. Any separation is conditional evidence for label-dependent prioritisation, not independent validation, a biomarker claim, or a resolved gene-level significance test.

## 3. Frozen Tier2 candidate universe

The exact 671 region–gene pairs were preserved: Ventral Striatum 146, Insula 70, Hippocampal Formation 127, CG25 186, dlPFC BA8/9 72, and OFC BA11 70. Expression values, sample IDs, gene identities, and class counts were unchanged. Only labels were shuffled within region using deterministic seeds 41001–41020.

## 4. Permutation procedure

For every seed the pipeline regenerated five-fold × five-repeat stratified outer CV, training-only median imputation/scaling, three-fold inner tuning on balanced accuracy, all five models, 750 feature-model fits, 150 DummyClassifier prediction references, fold gene rankings, two-level Robust Rank Aggregation, split-level RRA, and the LASSO selection-frequency ranking. The DummyClassifier never contributed gene evidence.

All 20 permutations completed; none failed. A completed seed was rerun with `--resume`, verified by checksums, and skipped without recomputation. Single-seed ML wall times ranged from 515.1 to 957.6 seconds (median 800.2; summed per-seed wall time 15856.5 seconds). Parallel execution affects wall time only, not method or seeds.

Pilot empirical probabilities use `(1 + null as-or-more-extreme) / 21`; the minimum possible value is 0.0476.

## 5. Tier2-RRA null results

RRA evidence was regionally heterogeneous. Ventral Striatum showed the clearest separation from shuffled labels across rank correlation and all top-set overlaps. Hippocampal Formation, dlPFC BA8/9, and CG25 showed separation for some broader top-set endpoints but not uniformly across rank and top-10 measures. Insula was suggestive but did not reach the pilot floor. OFC BA11 observed RRA resampling stability was similar to or worse than the shuffled-label reference.

| Region | Mean rank Spearman: observed vs null | Pilot p | Top-10 Jaccard: observed vs null | Pilot p | Top-20% Jaccard: observed vs null | Pilot p |
|---|---:|---:|---:|---:|---:|---:|
| Ventral Striatum | 0.716 vs 0.624 | 0.0476 | 0.476 vs 0.260 | 0.0476 | 0.469 vs 0.344 | 0.0476 |
| Insula | 0.617 vs 0.553 | 0.1905 | 0.358 vs 0.297 | 0.2381 | 0.405 vs 0.322 | 0.1429 |
| Hippocampal Formation | 0.696 vs 0.620 | 0.0952 | 0.285 vs 0.265 | 0.3810 | 0.451 vs 0.338 | 0.0476 |
| CG25 | 0.644 vs 0.640 | 0.4286 | 0.212 vs 0.261 | 0.8095 | 0.413 vs 0.345 | 0.0476 |
| dlPFC BA8/9 | 0.601 vs 0.554 | 0.1905 | 0.357 vs 0.289 | 0.0952 | 0.410 vs 0.326 | 0.0476 |
| OFC BA11 | 0.486 vs 0.587 | 1.0000 | 0.246 vs 0.300 | 0.7619 | 0.299 vs 0.328 | 0.5714 |

Consensus-score concentration was also mixed. Several broader top-set median RRA-score summaries were stronger than the null, but the minimum RRA score alone was usually not unusual. RRA scores are aggregation scores, not MDD association p-values. Counts supported by at least three models were uninformative because the dense models supplied broad evidence to every candidate.

## 6. LASSO null results

Observed LASSO selection-frequency concentration exceeded the null in every region, especially for the leading 10–25 genes. Resample ranking stability was clearer in Ventral Striatum and dlPFC, moderate in Hippocampal Formation and Insula, and not better than null in CG25 or OFC BA11.

| Region | Max frequency: observed vs null | Pilot p | Mean rank Spearman: observed vs null | Pilot p | Top-10 Jaccard: observed vs null | Pilot p |
|---|---:|---:|---:|---:|---:|---:|
| Ventral Striatum | 1.000 vs 0.506 | 0.0476 | 0.567 vs 0.355 | 0.0476 | 0.503 vs 0.219 | 0.0476 |
| Insula | 0.960 vs 0.446 | 0.0476 | 0.479 vs 0.341 | 0.1905 | 0.314 vs 0.255 | 0.2381 |
| Hippocampal Formation | 1.000 vs 0.488 | 0.0476 | 0.455 vs 0.321 | 0.1429 | 0.278 vs 0.215 | 0.2857 |
| CG25 | 0.800 vs 0.462 | 0.1905 | 0.310 vs 0.327 | 0.4762 | 0.184 vs 0.220 | 0.6190 |
| dlPFC BA8/9 | 1.000 vs 0.450 | 0.0476 | 0.492 vs 0.375 | 0.0476 | 0.335 vs 0.268 | 0.0952 |
| OFC BA11 | 0.960 vs 0.384 | 0.0476 | 0.299 vs 0.338 | 0.5714 | 0.195 vs 0.253 | 0.8095 |

Thus the simple sparse model supports label-dependent concentration, but it does not establish uniformly reproducible rankings in all regions.

## 7. Ranking stability

The primary endpoint does not support one project-wide yes/no statement. Ventral Striatum provides clear conditional stability evidence. Hippocampal Formation and dlPFC provide partial evidence, depending on top-set breadth. CG25 and Insula are mixed. OFC BA11 fails the predeclared ranking-stability criterion despite strong predictive separation.

This regional heterogeneity is scientifically important and should not be hidden by pooling endpoints or constructing a weighted score.

## 8. Observed top-gene behavior

Leading observed genes were usually ranked better than their own null median, but several reappeared near the top under multiple shuffled labels. Gene-specific probabilities are exploratory only; with 20 permutations they have coarse resolution and no BH/FDR interpretation.

| Region | Observed rank-1 gene | Null top-10 recurrence | Null median rank (range) | Exploratory rank-tail p |
|---|---|---:|---:|---:|
| Ventral Striatum | MORN5 | 4/20 | 32.5 (2–89) | 0.0476 |
| Insula | TMEM161B-DT | 10/20 | 10.5 (1–33) | 0.0952 |
| Hippocampal Formation | TOR2A | 0/20 | 93.5 (12–99) | 0.0476 |
| CG25 | LOC107985795 | 2/20 | 38.5 (2–90) | 0.0476 |
| dlPFC BA8/9 | LINC00269 | 11/20 | 10.0 (2–17) | 0.0476 |
| OFC BA11 | LOC105373386 | 8/20 | 17.5 (1–38) | 0.1429 |

Full top-25 rank and RRA-score distributions are preserved in `RRA_observed_topgene_null_behavior.csv`; the corresponding LASSO distributions are in `LASSO_observed_topgene_null_behavior.csv`.

## 9. Secondary prediction diagnostics

Observed best-of-five balanced accuracy exceeded every one of the 20 corresponding null maxima in all six regions. This is a secondary diagnostic and **must not be used as evidence of an unbiased diagnostic classifier**, because Tier2 was selected from the same cohort using real Diagnosis.

| Region | Observed best-of-five BA | Null mean | Null 95% interval | Pilot p |
|---|---:|---:|---:|---:|
| Ventral Striatum | 0.941 | 0.557 | 0.471–0.710 | 0.0476 |
| Insula | 0.958 | 0.546 | 0.481–0.663 | 0.0476 |
| Hippocampal Formation | 0.912 | 0.553 | 0.468–0.703 | 0.0476 |
| CG25 | 0.937 | 0.581 | 0.447–0.707 | 0.0476 |
| dlPFC BA8/9 | 0.927 | 0.535 | 0.477–0.654 | 0.0476 |
| OFC BA11 | 0.931 | 0.520 | 0.465–0.652 | 0.0476 |

## 10. Interpretation using the predeclared scenarios

No single scenario describes every region:

- Ventral Striatum is closest to **Scenario A**: RRA and LASSO stability/concentration both exceed the conditional null.
- Hippocampal Formation and dlPFC BA8/9 provide partial Scenario A support, with stronger separation for broader top sets and LASSO concentration than for every RRA endpoint.
- CG25 and Insula are mixed: leading-gene concentration exists, but ranking-stability endpoints are inconsistent.
- OFC BA11 is closest to **Scenario E**: prediction strongly exceeds null, while RRA and LASSO resampling stability do not. Prediction alone cannot justify OFC gene prioritisation.

Overall, the pilot provides conditional evidence for label-dependent prioritisation in some regions, not uniform validation of the entire Tier2-RRA ranking.

## 11. Limitations

1. Tier2 was selected with the observed labels, so the null is conditional and favoring of the observed result is expected.
2. Twenty permutations give a minimum empirical probability of 0.0476 and cannot resolve smaller tails.
3. Multiple separate endpoints and regions were inspected; the pilot probabilities are descriptive and are not multiplicity-adjusted claims.
4. Small regional sample sizes and correlated expression features can destabilize gene-level rankings.
5. RRA scores are rank-aggregation evidence, not disease-association p-values.
6. Prediction metrics are optimistic within-study diagnostics and not external validation.

## 12. Recommendation

**Stop at the completed 20-permutation pilot for now and proceed to the upstream statistical-method benchmark while retaining this analysis as conditional descriptive evidence.** Extending to 99 would improve tail resolution, but it would not repair the central candidate-selection bias and would require substantial additional computation. A 99-permutation extension should be reserved for a later explicit decision if precise conditional tails become necessary.

The primary report should retain cautious, region-specific language. It should not claim uniformly stable ML prioritisation, and OFC BA11 gene-level ML conclusions should be softened substantially. No validated MDD genes, biomarkers, or independent validation are established here.

## Resume command

The completed pilot can be verified/resumed without recomputing finished seeds:

```bat
.runtime\python\python.exe code\tier2_conditional_null\04_run_pilot.py --start 1 --end 20 --workers 3 --resume
```
