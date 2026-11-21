> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Upstream statistical-method benchmark report

## 1. Why Tier2 needed benchmarking

Tier2 was defined using a gene-wise LM on `log2(GeTMM + 1)`. This analysis asks whether its directions and genome-wide positions remain similar when the same samples and eligible genes are analysed using limma-voom and DESeq2. This is **within-study statistical-method concordance**, not independent replication.

## 2. Frozen current LM and Tier2 definition

**SOURCE-DERIVED FACT:** The frozen model was `expression ~ Diagnosis + Age + Sex + RIN`, with Control as reference and positive diagnosis effects meaning higher expression in MDD. Tier2 remains exactly 671 region–gene pairs selected by raw diagnosis p < 0.01 and absolute standardized diagnosis effect > 0.8. No membership was changed.

The LM used `log2(GeTMM + 1)`. Its standardized effect equals the diagnosis beta divided by the gene-wise expression SD. BH correction was regional across the full eligible LM universe.

## 3. Matched samples and genes

**SOURCE-DERIVED FACT:** Matching was exact. No benchmark sample or gene was added or removed relative to its regional frozen LM universe.

| Region                | Matched samples | Matched eligible genes |
| --------------------- | --------------- | ---------------------- |
| Ventral Striatum      | 48              | 23917                  |
| Insula                | 48              | 23090                  |
| Hippocampal Formation | 43              | 24371                  |
| CG25                  | 27              | 23941                  |
| dlPFC BA8/9           | 48              | 23941                  |
| OFC BA11              | 47              | 23941                  |

The exact filtered raw-count source is `D:/genes/GeneDepresion/outputs/step1_length_sensitivity/Length_gt_100/03_filtered_raw_counts/keep_all_samples/`. A same-sized but non-identical later Step1 export was rejected before fitting.

## 4. limma-voom method

Raw integer counts were passed through `edgeR::DGEList`, TMM normalization, standard `limma::voom`, `lmFit`, and `eBayes`, using `~ Diagnosis + Age + Sex + RIN`. No quality weights or result-driven tuning were used. The diagnosis coefficient is MDD minus Control.

## 5. DESeq2 method

Raw integer counts were fitted with standard DESeq2 using `~ Age + Sex + RIN + Diagnosis` and the explicit MDD-versus-Control contrast. Native DESeq2 `padj` and a separate BH correction across the same eligible regional universe are both retained.

| component             | value                                         |
| --------------------- | --------------------------------------------- |
| R                     | R version 4.4.2 (2024-10-31 ucrt)             |
| limma                 | 3.62.2                                        |
| edgeR                 | 4.4.2                                         |
| DESeq2                | 1.46.0                                        |
| limma_normalization   | TMM via edgeR::calcNormFactors; standard voom |
| DESeq2_size_factors   | DESeq2 median-of-ratios estimateSizeFactors   |
| DESeq2_dispersion_fit | standard DESeq() defaults                     |
| contrast              | Diagnosis MDD - CTRL                          |

## 6. Genome-wide direction agreement

**SOURCE-DERIVED FACT:** Genome-wide direction agreement between the LM and benchmark methods was 92.9%–96.9%. All 671 frozen Tier2 pairs had the same direction in all three methods.

| region                | Frozen LM vs DESeq2 | Frozen LM vs limma-voom |
| --------------------- | ------------------- | ----------------------- |
| Ventral Striatum      | 95.3%               | 96.0%                   |
| Insula                | 95.0%               | 95.4%                   |
| Hippocampal Formation | 93.4%               | 93.0%                   |
| CG25                  | 96.0%               | 96.9%                   |
| dlPFC BA8/9           | 93.2%               | 93.5%                   |
| OFC BA11              | 96.0%               | 95.5%                   |

**INTERPRETATION:** Direction is highly robust to statistical-model choice in this cohort.

## 7. Genome-wide rank agreement

The primary rank is two-sided raw p ascending, with GeneID as deterministic tie-breaker. Exact coefficient magnitudes were not compared as if they were interchangeable.

| region                | Frozen LM vs DESeq2 | Frozen LM vs limma-voom |
| --------------------- | ------------------- | ----------------------- |
| Ventral Striatum      | 0.955               | 0.965                   |
| Insula                | 0.946               | 0.948                   |
| Hippocampal Formation | 0.904               | 0.901                   |
| CG25                  | 0.954               | 0.966                   |
| dlPFC BA8/9           | 0.912               | 0.927                   |
| OFC BA11              | 0.965               | 0.952                   |

**SOURCE-DERIVED FACT:** LM-versus-limma and LM-versus-DESeq2 Spearman correlations were 0.901–0.966 across regions. Kendall correlations are retained in `method_rank_correlations.csv`.

**INTERPRETATION:** Genome-wide ordering is strongly concordant, although the top few positions are not identical.

## 8. Top-set overlap

Top 10, 25, 50, 1%, and 5% pairwise intersections and Jaccard indices are retained in `method_topset_overlap.csv`. They show partial rather than identical leading-gene lists, which is expected from different count and variance models.

## 9. Direct audit of Tier2

The complete 671-row table `tier2_statistical_method_audit.csv` preserves LM, limma, and DESeq2 effects, p-values, adjusted values, ranks, percentiles, and direction-agreement flags separately. No consensus score or new tier was constructed.

## 10. Tier2 enrichment in benchmark rankings

| Region                | Method     | Tier2 in top 1% | Expected | Enrichment ratio | Tier2 median percentile |
| --------------------- | ---------- | --------------- | -------- | ---------------- | ----------------------- |
| Ventral Striatum      | limma-voom | 108.0           | 1.47     | 73.7             | 0.49                    |
| Ventral Striatum      | DESeq2     | 99.0            | 1.47     | 67.6             | 0.60                    |
| Insula                | limma-voom | 69.0            | 0.70     | 98.5             | 0.22                    |
| Insula                | DESeq2     | 67.0            | 0.70     | 95.7             | 0.28                    |
| Hippocampal Formation | limma-voom | 116.0           | 1.27     | 91.2             | 0.32                    |
| Hippocampal Formation | DESeq2     | 95.0            | 1.27     | 74.7             | 0.57                    |
| CG25                  | limma-voom | 145.0           | 1.86     | 77.8             | 0.43                    |
| CG25                  | DESeq2     | 136.0           | 1.86     | 72.9             | 0.53                    |
| dlPFC BA8/9           | limma-voom | 66.0            | 0.72     | 91.4             | 0.26                    |
| dlPFC BA8/9           | DESeq2     | 67.0            | 0.72     | 92.8             | 0.25                    |
| OFC BA11              | limma-voom | 62.0            | 0.70     | 88.4             | 0.18                    |
| OFC BA11              | DESeq2     | 64.0            | 0.70     | 91.2             | 0.21                    |

**SOURCE-DERIVED FACT:** Regional Tier2 median positions were within the strongest 0.18%–0.60% of both benchmark rankings. Top-1% enrichment ratios ranged from 67.6-fold to 98.5-fold, and every predeclared top-set enrichment remained below BH 0.05.

**INTERPRETATION:** Tier2 is strongly concentrated near the top of both benchmarks. Because Tier2 and all benchmark methods use the same cohort, this is expected to some degree and must not be called external validation.

## 11. FDR results

| Region                | Method                  | Raw p<0.01 | Adjusted<0.05 | Adjusted<0.10 | Minimum adjusted |
| --------------------- | ----------------------- | ---------- | ------------- | ------------- | ---------------- |
| Ventral Striatum      | Frozen LM               | 591        | 0             | 0             | 0.124            |
| Ventral Striatum      | limma-voom              | 533        | 0             | 2             | 0.082            |
| Ventral Striatum      | DESeq2_same_universe_BH | 825        | 8             | 30            | 0.021            |
| Ventral Striatum      | DESeq2_native_padj      | 825        | 11            | 38            | 0.008            |
| Insula                | Frozen LM               | 197        | 0             | 0             | 0.934            |
| Insula                | limma-voom              | 194        | 0             | 0             | 0.843            |
| Insula                | DESeq2_same_universe_BH | 262        | 0             | 0             | 0.158            |
| Insula                | DESeq2_native_padj      | 262        | 0             | 0             | 0.158            |
| Hippocampal Formation | Frozen LM               | 142        | 0             | 0             | 0.997            |
| Hippocampal Formation | limma-voom              | 117        | 0             | 0             | 1.000            |
| Hippocampal Formation | DESeq2_same_universe_BH | 211        | 0             | 0             | 0.272            |
| Hippocampal Formation | DESeq2_native_padj      | 211        | 0             | 0             | 0.272            |
| CG25                  | Frozen LM               | 203        | 0             | 0             | 0.593            |
| CG25                  | limma-voom              | 174        | 0             | 0             | 0.999            |
| CG25                  | DESeq2_same_universe_BH | 293        | 0             | 1             | 0.065            |
| CG25                  | DESeq2_native_padj      | 293        | 0             | 1             | 0.065            |
| dlPFC BA8/9           | Frozen LM               | 265        | 0             | 0             | 0.838            |
| dlPFC BA8/9           | limma-voom              | 325        | 0             | 0             | 0.496            |
| dlPFC BA8/9           | DESeq2_same_universe_BH | 318        | 0             | 0             | 0.218            |
| dlPFC BA8/9           | DESeq2_native_padj      | 318        | 0             | 0             | 0.218            |
| OFC BA11              | Frozen LM               | 298        | 0             | 0             | 0.144            |
| OFC BA11              | limma-voom              | 281        | 1             | 1             | 0.041            |
| OFC BA11              | DESeq2_same_universe_BH | 341        | 1             | 1             | 0.005            |
| OFC BA11              | DESeq2_native_padj      | 341        | 1             | 1             | 0.005            |

**SOURCE-DERIVED FACT:** The frozen LM produced no BH<0.05 gene. limma-voom produced one in OFC BA11 (`QPRT`). The explicit same-universe DESeq2 BH produced eight in Ventral Striatum and one in OFC BA11; native DESeq2 `padj` produced eleven and one, respectively. No other region produced BH/padj<0.05 genes.

The genes reaching BH<0.05 in limma or the explicit DESeq2 comparison were:

| Region           | Gene         | Frozen Tier2 | LM BH | limma BH | DESeq2 same-universe BH | DESeq2 native padj | Three-way direction |
| ---------------- | ------------ | ------------ | ----- | -------- | ----------------------- | ------------------ | ------------------- |
| Ventral Striatum | YTHDF2       | False        | 0.133 | 0.332    | 0.043                   | 0.024              | True                |
| Ventral Striatum | TSGA10       | True         | 0.124 | 0.082    | 0.021                   | 0.008              | True                |
| Ventral Striatum | PTPN4        | True         | 0.178 | 0.332    | 0.048                   | 0.029              | True                |
| Ventral Striatum | KIAA1191     | False        | 0.124 | 0.426    | 0.043                   | 0.024              | True                |
| Ventral Striatum | NPTX2        | True         | 0.208 | 0.332    | 0.048                   | NA                 | True                |
| Ventral Striatum | BMS1P1       | True         | 0.124 | 0.141    | 0.022                   | 0.008              | True                |
| Ventral Striatum | LOC107987147 | True         | 0.124 | 0.082    | 0.040                   | NA                 | True                |
| Ventral Striatum | LRRC74B      | False        | 0.376 | 0.428    | 0.040                   | NA                 | True                |
| OFC BA11         | QPRT         | True         | 0.144 | 0.041    | 0.005                   | 0.005              | True                |

All have the same effect direction across the three methods. These are method-specific same-cohort findings. They do not replace Tier2, and they should be independently checked before any manuscript claim.

## 12. Region-specific interpretation

- **Ventral Striatum:** very high directional/rank concordance; DESeq2 alone produces eight same-universe BH<0.05 genes. This is supportive but method-specific.
- **Insula:** high concordance and strong Tier2 concentration; no method produces an FDR-supported gene.
- **Hippocampal Formation:** the lowest, but still strong, LM rank concordance (about 0.90); no FDR-supported gene.
- **CG25:** very high LM/limma and LM/DESeq2 concordance despite the smaller sample size; no BH<0.05 gene.
- **dlPFC BA8/9:** high direction/rank concordance and strong Tier2 ranking concentration; no BH<0.05 gene.
- **OFC BA11:** high concordance; `QPRT` is rank 1 in all three methods and reaches BH<0.05 in limma and DESeq2, but not the frozen LM.

## 13. Current Tier2-RRA top-gene audit

| Region                | Tier2-RRA #1 gene | LM effect | LM p     | LM rank | limma effect | limma p  | limma rank | DESeq2 effect | DESeq2 p | DESeq2 rank | Three-way direction |
| --------------------- | ----------------- | --------- | -------- | ------- | ------------ | -------- | ---------- | ------------- | -------- | ----------- | ------------------- |
| Ventral Striatum      | MORN5             | 0.500     | 2.91e-04 | 26      | 0.577        | 4.97e-04 | 25         | 0.630         | 5.05e-05 | 15          | yes                 |
| Insula                | TMEM161B-DT       | 0.344     | 3.06e-04 | 7       | 0.356        | 1.43e-04 | 3          | 0.366         | 3.69e-05 | 5           | yes                 |
| Hippocampal Formation | TOR2A             | 0.160     | 0.002    | 26      | 0.260        | 0.005    | 43         | 0.267         | 0.005    | 109         | yes                 |
| CG25                  | LOC107985795      | 0.280     | 9.63e-04 | 23      | 0.558        | 0.002    | 24         | 0.570         | 0.013    | 349         | yes                 |
| dlPFC BA8/9           | LINC00269         | 0.150     | 8.42e-04 | 13      | 0.652        | 0.002    | 41         | 0.676         | 8.34e-04 | 36          | yes                 |
| OFC BA11              | LOC105373386      | -0.036    | 0.002    | 45      | -0.471       | 0.002    | 49         | -0.471        | 0.010    | 343         | yes                 |

**SOURCE-DERIVED FACT:** All six current Tier2-RRA rank-one genes retain the same direction across methods. Their benchmark ranks range from 3 to 49 under limma and 5 to 349 under DESeq2; CG25 `LOC107985795` and OFC `LOC105373386` show the largest DESeq2 rank weakening.

## 14. Limitations

1. All methods use the same cohort; this is not replication.
2. Tier2 was selected using the frozen LM, so enrichment in correlated benchmark rankings is not independent.
3. Regional sample sizes are small, particularly CG25.
4. LM, limma-voom, and DESeq2 use different normalization and distributional assumptions.
5. Effect magnitudes are on related but non-identical scales.
6. Method agreement does not establish a biological mechanism or biomarker.
7. No external cohort was included.
8. Cross-method agreement cannot convert an absent LM FDR result into LM significance.

## 15. Recommendation

**INTERPRETATION:** The result is closest to prospective **Scenario A**, with the important same-cohort limitation: effect directions, genome-wide rankings, and Tier2 concentration are strongly concordant across methods. Tier2 is robust enough to retain unchanged as an exploratory candidate screen.

Do not create a consensus tier or modify Tier2. Report the limma/DESeq2 FDR findings as method-specific sensitivity results, not replacements for the frozen analysis. Continue to describe Tier2 and downstream ML as exploratory prioritisation requiring external validation.
