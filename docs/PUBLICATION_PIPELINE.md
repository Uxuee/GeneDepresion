# Publication Pipeline

ML-prioritised genes are internally supported exploratory candidates, not validated biomarkers.

1. **Raw RNA-seq and metadata:** GSE102556 expression/count inputs and sample annotations supply diagnosis, anatomy, and covariates.
2. **Filtering and normalisation:** the main clean analysis uses the Length >100 definition, GeTMM normalisation, and log2(GeTMM + 1) expression.
3. **QC/covariates:** sample, expression, metadata, and region-specific checks precede diagnosis modelling.
4. **Gene-wise models:** saved models estimate covariate-adjusted Diagnosis effects and p-values. This is the statistical-association layer.
5. **Standardised effects:** Diagnosis coefficients are scaled for comparable magnitude while preserving direction.
6. **Evidence tiers:** Tier 1/2/3 are rule-based groups. Tier 2 is the exploratory candidate population. Tier 4 is an all-tested-gene ranking used in historical sensitivity/performance runs; it is not a final list and does not enter the corrected score.
7. **Enrichment:** Tier 2 over-representation analyses supply biological context, not individual-gene validation.
8. **ML v1:** the final score uses `Tier2_fixed` only: five models × five repeats × five folds = 125 events per region.
9. **Strict ML v2:** two train-only selectors × five top-k values × five models × three folds = 150 non-baseline events per region.
10. **Corrected prioritisation:** deterministic post-processing combines the events using `ML_PRIORITY_SCORE_FINAL.md`. Predictive metrics do not enter the score. Result: 671 pairs = 30 high + 119 moderate + 522 low.
11. **Biological follow-up:** enrichment, AHBA, WGCNA, hub genes, and module-level ML are exploratory context and validation targets. Legacy-list-dependent outputs require regeneration.

| Layer | Question | Permitted interpretation |
|---|---|---|
| Gene-wise model | Association with Diagnosis after adjustment? | Within-dataset statistical association |
| Tier 2 | Meets exploratory effect/p-value rules? | Exploratory candidate |
| ML score | Repeated internal ML support? | Internally prioritised candidate |
| Predictive performance | Internal CV performance? | Exploratory performance only |
| Enrichment/AHBA/WGCNA | Biological overlap/context? | Annotation, not validation |
| External validation | Reproduces independently? | Not complete |

Exact tier cutoffs, covariate formulas, and source tables remain defined by their saved scripts and outputs; this cleanup does not revise them.
