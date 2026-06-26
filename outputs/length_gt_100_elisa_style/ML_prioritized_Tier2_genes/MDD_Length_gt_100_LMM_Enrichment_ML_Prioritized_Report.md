# Length > 100 MDD Transcriptomics: LMM, Enrichment, and ML-prioritized Tier 2 Genes

## 1. Executive summary

We analyzed Length > 100 genes using GeTMM-normalized log2 expression. The main gene-level model was a linear mixed model (LMM), run gene by gene, to estimate MDD vs Control expression effects.

Strict corrected single-gene discoveries were limited, so the analysis moved to an evidence-tiered exploratory interpretation. Tier 2 exploratory candidate genes were defined using nominal association evidence and standardized effect size:

- `p_value_diagnosis < 0.01`
- `abs(standardized_effect_size_diagnosis) > 0.8`

Tier 2 genes were then used for enrichment and ML-guided prioritization. ML v1 is interpreted as Tier2-guided ML gene prioritization, not validated classifier performance. Strict ML v2 is interpreted as a leakage-controlled prediction check.

The final stable ML-prioritized extraction ranked 671 Tier 2 gene-region pairs and identified 35 high-confidence and 130 moderate-confidence ML-prioritized Tier 2 genes.

**The final output should be interpreted as ML-prioritized candidate genes, not validated MDD biomarkers.**

## 2. Analysis workflow

```text
raw counts
-> Length > 100 filter
-> GeTMM normalization
-> log2(GeTMM + 1)
-> LMM
-> standardized effect size
-> evidence tiers
-> Tier 2 enrichment
-> ML v1 Tier2-guided prioritization
-> strict ML v2 prediction check
-> stable ML-prioritized Tier 2 genes
-> enrichment revisit / AHBA / module-level ML
```

## 3. LMM stage

The LMM stage was run gene by gene. For each gene, the model estimated the MDD vs Control expression effect while accounting for available covariates and repeated-subject structure where possible.

Main quantities:

- `beta_diagnosis`: the estimated MDD effect size. Positive values mean higher expression in MDD; negative values mean lower expression in MDD.
- `p_value_diagnosis`: nominal statistical evidence for the diagnosis association.
- `FDR` / `q-value`: corrected evidence after accounting for many genes being tested.

Strict corrected discoveries were very limited, so the analysis moved to evidence-tiered exploratory interpretation rather than relying only on FDR-significant single genes.

## 4. Standardized effect size

The standardized effect size was computed as:

```text
standardized_effect_size_diagnosis = beta_diagnosis / SD(expression of the gene)
```

This matters because raw beta values are not equally comparable across genes. A raw beta can look small for one gene but still be large relative to that gene's natural variability. The standardized effect size asks whether the MDD effect is large compared with the variability of that gene.

The sign is preserved:

- positive standardized effect size = higher in MDD
- negative standardized effect size = lower in MDD

## 5. Evidence tiers

The analysis uses four evidence tiers.

| Tier | Definition | Interpretation |
|---|---|---|
| Tier 1 | `Storey q < 0.10` and `abs(standardized_effect_size_diagnosis) > 0.8` | Strict corrected/high-confidence anchors |
| Tier 2 | `p_value_diagnosis < 0.01` and `abs(standardized_effect_size_diagnosis) > 0.8` | Main exploratory DE candidate genes |
| Tier 3 | `p_value_diagnosis < 0.01` and `abs(beta_diagnosis) > 0.2` | Old raw-beta comparison baseline |
| Tier 4 | Full ranked gene signatures | Used for ranked analyses such as GSEA and feature-ranking comparisons |

Main Tier 2 counts:

| Region/subregion | Tier 2 genes |
|---|---:|
| Ventral Striatum | 146 |
| Insula | 70 |
| Hippocampal Formation | 127 |
| Cingulate gyrus 25 | 186 |
| dlPFC BA8/9 | 72 |
| OFC BA11 | 70 |
| broad Frontal Cortex | 1, documentation only |

Broad Frontal Cortex was not used as the primary frontal set because it had only one Tier 2 gene after standardization.

## 6. Frontal subregion sensitivity

Broad Frontal Cortex was split into three frontal subregions:

- Cingulate gyrus 25
- dlPFC BA8/9
- OFC BA11

Cingulate gyrus 25 had the strongest frontal Tier 2 signal. This suggested that the broad Frontal Cortex grouping was diluted by subregion heterogeneity.

Caveat: the subregion model was not identical to the grouped model because the random effect structure was limited after splitting. The frontal subregion analysis should therefore be interpreted as sensitivity analysis rather than a perfect replacement for the grouped model.

## 7. P-value correction sensitivity

Multiple p-value correction strategies were tested. Standard across-gene correction remained conservative, and only limited corrected discoveries survived.

Larger counts appeared under sensitivity correction scopes, but those were not used as primary claims because they depend on the correction scope. Therefore Tier 2 was retained as the working exploratory candidate set.

**Tier 2 is not an FDR-significant gene list. It is an exploratory candidate set combining nominal evidence and standardized effect size.**

## 8. Enrichment stage

Enrichment asked what pathways or biological processes were overrepresented among Tier 2 genes. The enrichment run used the correct tested-gene background for each region or subregion. Broad Frontal Cortex was not interpreted because it had only one Tier 2 gene.

Enrichment is exploratory and supports biological interpretation; it is not proof that a pathway is definitively altered in MDD. The enrichment report also recommended ranked GSEA later using Tier 4 signatures.

Enrichment should now be revisited using:

- high-confidence ML-prioritized Tier 2 genes
- moderate-confidence ML-prioritized Tier 2 genes
- overlap with previous Tier 2 enrichment terms

## 9. ML v1: Tier2-guided ML prioritization

ML v1 used fixed Tier 2 and fixed Tier 4 features selected or ranked from the full dataset.

**ML v1 is interpreted as Tier2-guided ML gene prioritization, not validated classifier performance.**

This makes ML v1 useful for identifying which candidate genes help separate MDD vs Control in the available dataset. It is not an unbiased prediction test because the feature sets were selected or ranked before cross-validation using the full dataset.

ML v1 produced very high classification metrics. Those metrics are not the main result. The useful output is feature importance and model support for Tier 2 genes.

## 10. ML v2: strict leakage-controlled prediction check

ML v2 selected features inside cross-validation. The held-out test folds were not used for:

- feature selection
- scaling
- imputation
- hyperparameter tuning
- model fitting

ML v2 is therefore a stricter internal prediction check. Performance dropped compared with v1. Cingulate gyrus 25 retained the strongest moderate signal, while other regions were weak or near chance.

This means raw-gene diagnostic prediction is unstable. It does not mean Tier 2 genes are useless for biological prioritization.

**ML v2 acts as a validation check on raw-gene predictability, while ML v1 and the stable feature extraction support gene prioritization.**

## 11. Stable ML-prioritized Tier 2 gene extraction

The final ML-prioritized extraction used existing Tier 2 genes as the candidate pool and summarized repeated ML support across the available ML outputs.

Candidate pool: 671 Tier 2 gene-region pairs.

Main result:

- 35 high-confidence ML-prioritized Tier 2 genes
- 130 moderate-confidence ML-prioritized Tier 2 genes
- 671 total Tier 2 genes ranked

High-confidence counts:

| Region/subregion      |   High-confidence genes |
|:----------------------|------------------------:|
| Cingulate gyrus 25    |                       8 |
| Ventral Striatum      |                       8 |
| dlPFC BA8/9           |                       7 |
| Insula                |                       6 |
| Hippocampal Formation |                       3 |
| OFC BA11              |                       3 |

Moderate-confidence counts:

| Region/subregion      |   Moderate-confidence genes |
|:----------------------|----------------------------:|
| Cingulate gyrus 25    |                          29 |
| Ventral Striatum      |                          18 |
| dlPFC BA8/9           |                          12 |
| Insula                |                          18 |
| Hippocampal Formation |                          25 |
| OFC BA11              |                          28 |

Full summary:

| Region/subregion      |   Tier 2 genes |   High-confidence |   Moderate-confidence |   Sparse-model support |   Strict v2 support |   Median ML priority score |
|:----------------------|---------------:|------------------:|----------------------:|-----------------------:|--------------------:|---------------------------:|
| Cingulate gyrus 25    |            186 |                 8 |                    29 |                    123 |                 162 |                      0.557 |
| Ventral Striatum      |            146 |                 8 |                    18 |                     97 |                 116 |                      0.531 |
| dlPFC BA8/9           |             72 |                 7 |                    12 |                     58 |                  71 |                      0.599 |
| Insula                |             70 |                 6 |                    18 |                     65 |                  68 |                      0.612 |
| Hippocampal Formation |            127 |                 3 |                    25 |                     93 |                  91 |                      0.535 |
| OFC BA11              |             70 |                 3 |                    28 |                     69 |                  70 |                      0.638 |

Top 5 genes per region/subregion by `ML_priority_score` are shown below. The full CSV contains all 671 ranked Tier 2 gene-region pairs.

### Cingulate gyrus 25

| Gene         |   Beta |   p value |   Std effect | LMM direction   |   ML selection freq |   Model support | Sparse support   |   Sign consistency |   LMM/ML sign agree |   ML priority | Priority class                        |
|:-------------|-------:|----------:|-------------:|:----------------|--------------------:|----------------:|:-----------------|-------------------:|--------------------:|--------------:|:--------------------------------------|
| LOC107985795 |  0.28  |  0.000963 |         1.28 | up_in_MDD       |               0.632 |               5 | True             |                  1 |                   1 |         0.831 | High-confidence ML-prioritized Tier 2 |
| LOC105377644 |  0.261 |  0.000564 |         1.35 | up_in_MDD       |               0.619 |               5 | True             |                  1 |                   1 |         0.823 | High-confidence ML-prioritized Tier 2 |
| NOC2LP2      |  0.105 |  0.000872 |         1.18 | up_in_MDD       |               0.607 |               5 | True             |                  1 |                   1 |         0.801 | High-confidence ML-prioritized Tier 2 |
| GEMIN7       |  0.242 |  0.00103  |         1.28 | up_in_MDD       |               0.604 |               5 | True             |                  1 |                   1 |         0.8   | High-confidence ML-prioritized Tier 2 |
| C1orf174     |  0.156 |  0.00141  |         1.26 | up_in_MDD       |               0.559 |               5 | True             |                  1 |                   1 |         0.787 | High-confidence ML-prioritized Tier 2 |

### Ventral Striatum

| Gene   |   Beta |   p value |   Std effect | LMM direction   |   ML selection freq |   Model support | Sparse support   |   Sign consistency |   LMM/ML sign agree |   ML priority | Priority class                        |
|:-------|-------:|----------:|-------------:|:----------------|--------------------:|----------------:|:-----------------|-------------------:|--------------------:|--------------:|:--------------------------------------|
| MORN5  |  0.5   |  0.000291 |        0.92  | up_in_MDD       |               0.724 |               5 | True             |                  1 |                   1 |         0.871 | High-confidence ML-prioritized Tier 2 |
| FSHB   |  0.368 |  0.00131  |        0.906 | up_in_MDD       |               0.721 |               5 | True             |                  1 |                   1 |         0.837 | High-confidence ML-prioritized Tier 2 |
| MRPL2  |  0.148 |  0.00111  |        0.894 | up_in_MDD       |               0.659 |               5 | True             |                  1 |                   1 |         0.832 | High-confidence ML-prioritized Tier 2 |
| MOSPD3 |  0.213 |  0.000574 |        0.997 | up_in_MDD       |               0.634 |               5 | True             |                  1 |                   1 |         0.825 | High-confidence ML-prioritized Tier 2 |
| POF1B  |  0.124 |  0.00293  |        0.879 | up_in_MDD       |               0.631 |               5 | True             |                  1 |                   1 |         0.821 | High-confidence ML-prioritized Tier 2 |

### dlPFC BA8/9

| Gene         |   Beta |   p value |   Std effect | LMM direction   |   ML selection freq |   Model support | Sparse support   |   Sign consistency |   LMM/ML sign agree |   ML priority | Priority class                        |
|:-------------|-------:|----------:|-------------:|:----------------|--------------------:|----------------:|:-----------------|-------------------:|--------------------:|--------------:|:--------------------------------------|
| PPT2         |  0.107 |  0.000941 |        0.947 | up_in_MDD       |               0.634 |               5 | True             |              1     |                   1 |         0.831 | High-confidence ML-prioritized Tier 2 |
| NOL6         | -0.12  |  0.00428  |       -0.804 | down_in_MDD     |               0.631 |               5 | True             |              0.93  |                   1 |         0.821 | High-confidence ML-prioritized Tier 2 |
| LINC00269    |  0.15  |  0.000842 |        0.948 | up_in_MDD       |               0.551 |               5 | True             |              1     |                   1 |         0.809 | High-confidence ML-prioritized Tier 2 |
| LOC107984015 | -0.345 |  0.000534 |       -1     | down_in_MDD     |               0.62  |               5 | True             |              0.941 |                   1 |         0.802 | High-confidence ML-prioritized Tier 2 |
| LINC02897    |  0.118 |  0.00146  |        0.895 | up_in_MDD       |               0.551 |               5 | True             |              1     |                   1 |         0.788 | High-confidence ML-prioritized Tier 2 |

### Insula

| Gene         |   Beta |   p value |   Std effect | LMM direction   |   ML selection freq |   Model support | Sparse support   |   Sign consistency |   LMM/ML sign agree |   ML priority | Priority class                        |
|:-------------|-------:|----------:|-------------:|:----------------|--------------------:|----------------:|:-----------------|-------------------:|--------------------:|--------------:|:--------------------------------------|
| TMEM161B-DT  |  0.344 |  0.000306 |        1.04  | up_in_MDD       |               0.731 |               5 | True             |                  1 |                   1 |         0.89  | High-confidence ML-prioritized Tier 2 |
| BMS1P2       |  0.311 |  0.00021  |        1.04  | up_in_MDD       |               0.722 |               5 | True             |                  1 |                   1 |         0.874 | High-confidence ML-prioritized Tier 2 |
| BMS1P2-AGAP9 |  0.344 |  0.000209 |        1.03  | up_in_MDD       |               0.634 |               5 | True             |                  1 |                   1 |         0.846 | High-confidence ML-prioritized Tier 2 |
| ARL14EP-DT   |  0.564 |  0.0019   |        0.891 | up_in_MDD       |               0.562 |               5 | True             |                  1 |                   1 |         0.802 | High-confidence ML-prioritized Tier 2 |
| LOC101929710 |  0.494 |  0.00548  |        0.804 | up_in_MDD       |               0.579 |               5 | True             |                  1 |                   1 |         0.8   | High-confidence ML-prioritized Tier 2 |

### Hippocampal Formation

| Gene      |    Beta |   p value |   Std effect | LMM direction   |   ML selection freq |   Model support | Sparse support   |   Sign consistency |   LMM/ML sign agree |   ML priority | Priority class                            |
|:----------|--------:|----------:|-------------:|:----------------|--------------------:|----------------:|:-----------------|-------------------:|--------------------:|--------------:|:------------------------------------------|
| TOR2A     |  0.16   |  0.00231  |        0.968 | up_in_MDD       |               0.694 |               5 | True             |              1     |                   1 |         0.846 | High-confidence ML-prioritized Tier 2     |
| IL4I1     |  0.0987 |  0.00315  |        0.936 | up_in_MDD       |               0.569 |               5 | True             |              1     |                   1 |         0.793 | High-confidence ML-prioritized Tier 2     |
| ZNF663P   |  0.343  |  0.000182 |        1.16  | up_in_MDD       |               0.489 |               5 | True             |              1     |                   1 |         0.787 | Moderate-confidence ML-prioritized Tier 2 |
| HLA-F-AS1 |  0.252  |  0.000117 |        1.14  | up_in_MDD       |               0.44  |               5 | True             |              1     |                   1 |         0.782 | Moderate-confidence ML-prioritized Tier 2 |
| TRMT2A    | -0.139  |  0.000714 |       -1.06  | down_in_MDD     |               0.508 |               5 | True             |              0.894 |                   1 |         0.758 | High-confidence ML-prioritized Tier 2     |

### OFC BA11

| Gene         |    Beta |   p value |   Std effect | LMM direction   |   ML selection freq |   Model support | Sparse support   |   Sign consistency |   LMM/ML sign agree |   ML priority | Priority class                            |
|:-------------|--------:|----------:|-------------:|:----------------|--------------------:|----------------:|:-----------------|-------------------:|--------------------:|--------------:|:------------------------------------------|
| QPRT         |  0.281  |  6e-06    |        1.02  | up_in_MDD       |               0.547 |               5 | True             |              1     |                   1 |         0.833 | High-confidence ML-prioritized Tier 2     |
| LOC105373386 | -0.0365 |  0.00183  |       -0.823 | down_in_MDD     |               0.55  |               5 | True             |              0.956 |                   1 |         0.787 | High-confidence ML-prioritized Tier 2     |
| LOC105378913 | -0.185  |  0.0011   |       -0.857 | down_in_MDD     |               0.596 |               5 | True             |              0.966 |                   1 |         0.781 | High-confidence ML-prioritized Tier 2     |
| RAB24        |  0.117  |  0.00432  |        0.822 | up_in_MDD       |               0.478 |               5 | True             |              1     |                   1 |         0.765 | Moderate-confidence ML-prioritized Tier 2 |
| MYBPC1       |  0.266  |  0.000277 |        0.95  | up_in_MDD       |               0.384 |               5 | True             |              1     |                   1 |         0.756 | Moderate-confidence ML-prioritized Tier 2 |

## 12. How ML priority was computed

The ML priority score combined several forms of support:

- selection frequency across folds, repeats, models, and feature sets
- support across model families
- sparse-model support from LASSO / logistic L1 and Elastic Net
- coefficient or feature-importance magnitude
- rank stability
- sign consistency for signed linear models
- agreement between ML direction and LMM direction
- LMM evidence score
- strict v2 support when available

High-confidence genes are Tier 2 genes with stronger ML priority score, repeated selection, support from more than one model, sparse-model support, and consistent direction when signed coefficients are available.

Moderate-confidence genes have weaker but still meaningful ML support.

## 13. Interpretation of the new result

The project now has a refined candidate hierarchy:

1. Tier 1 genes: strict corrected anchors, mainly Ventral Striatum.
2. Tier 2 genes: exploratory MDD-associated candidate genes.
3. ML-prioritized Tier 2 genes: Tier 2 genes most repeatedly useful for MDD/control separation.
4. High-confidence ML-prioritized genes: strongest candidates for follow-up.
5. Moderate-confidence ML-prioritized genes: secondary candidates.

This is stronger than using all Tier 2 genes equally because it combines LMM evidence and ML support.

## 14. Current conclusions

The analysis does not support a validated diagnostic classifier. The analysis does support an exploratory MDD candidate-gene prioritization workflow.

Ventral Striatum and Cingulate gyrus 25 remain important, but for different reasons:

- Ventral Striatum has strict Tier 1 anchors and 8 high-confidence ML-prioritized genes.
- Cingulate gyrus 25 has the strongest frontal Tier 2 signal and 8 high-confidence ML-prioritized genes.

dlPFC BA8/9 and Insula also contribute meaningful ML-prioritized candidates. Hippocampal Formation and OFC BA11 have smaller high-confidence sets but should not be discarded.

## 15. Recommended next steps

1. Enrichment revisit: compare high/moderate ML-prioritized genes with previous Tier 2 enrichment terms.
2. Ranked GSEA: use Tier 4 signed ranked signatures.
3. AHBA/WGCNA module mapping: test whether ML-prioritized Tier 2 genes map to AHBA coexpression modules.
4. Pathway/module-level ML: use pathway scores or module eigengenes instead of raw genes.
5. External validation search: look for another MDD postmortem RNA-seq dataset or compatible validation resource.

## 16. Warnings / limitations

- Tier 2 genes are exploratory, not FDR-significant.
- ML-prioritized genes are candidate genes, not validated biomarkers.
- ML v1 performance is optimistic for prediction because features were selected or ranked from the full dataset.
- ML v2 is stricter but still internal validation only.
- Sample sizes are small.
- Postmortem transcriptomics has technical and biological noise.
- No external validation has been performed yet.

## 17. File and output summary

Key output folders:

- `outputs/length_gt_100_elisa_style/evidence_tiers_4groups/`
- `outputs/length_gt_100_elisa_style/frontal_subregion_sensitivity/`
- `outputs/length_gt_100_elisa_style/pvalue_correction_sensitivity/`
- `outputs/length_gt_100_elisa_style/DE_candidate_gene_sets/`
- `outputs/length_gt_100_elisa_style/enrichment_Tier2/`
- `outputs/length_gt_100_elisa_style/ML_pilot_Tier2/`
- `outputs/length_gt_100_elisa_style/ML_strict_v2/`
- `outputs/length_gt_100_elisa_style/ML_prioritized_Tier2_genes/`

Key new files:

- `outputs/length_gt_100_elisa_style/ML_prioritized_Tier2_genes/ML_prioritized_Tier2_genes_all_ranked.csv`
- `outputs/length_gt_100_elisa_style/ML_prioritized_Tier2_genes/ML_prioritized_Tier2_genes_high_confidence.csv`
- `outputs/length_gt_100_elisa_style/ML_prioritized_Tier2_genes/ML_prioritized_Tier2_genes_moderate_confidence.csv`
- `outputs/length_gt_100_elisa_style/ML_prioritized_Tier2_genes/ML_prioritized_Tier2_gene_summary_by_region.csv`
- `outputs/length_gt_100_elisa_style/ML_prioritized_Tier2_genes/model_support_summary.csv`
- `outputs/length_gt_100_elisa_style/ML_prioritized_Tier2_genes/gene_lists/`
- `outputs/length_gt_100_elisa_style/ML_prioritized_Tier2_genes/ML_Prioritized_Tier2_Report.md`

This integrated report is saved in two places:

- `docs/MDD_Length_gt_100_LMM_Enrichment_ML_Prioritized_Report.md`
- `outputs/length_gt_100_elisa_style/ML_prioritized_Tier2_genes/MDD_Length_gt_100_LMM_Enrichment_ML_Prioritized_Report.md`

## 18. Short glossary

| Term | Simple meaning |
|---|---|
| LMM | Linear mixed model; a regression model that can include fixed effects and repeated-subject/random-effect structure. |
| `beta_diagnosis` | Estimated MDD vs Control expression effect. Positive means higher in MDD; negative means lower in MDD. |
| `p_value_diagnosis` | Nominal evidence that diagnosis is associated with expression for one gene. |
| FDR/q-value | Multiple-testing-corrected p-value used when many genes are tested. |
| Standardized effect size | MDD effect divided by gene expression variability, making genes more comparable. |
| Tier 1 | Strict corrected high-confidence anchors. |
| Tier 2 | Exploratory candidate genes selected by nominal p-value and standardized effect size. |
| Tier 3 | Older raw-beta exploratory comparison set. |
| Tier 4 | Full ranked signatures for ranked analyses and feature-ranking comparisons. |
| Enrichment | A test asking whether a gene set is overrepresented in known pathways or biological categories. |
| ML v1 | Tier2-guided exploratory ML used for gene prioritization, not unbiased prediction. |
| ML v2 | Strict leakage-controlled internal prediction check. |
| Feature-selection leakage | When test data influence which features are selected, making prediction performance optimistic. |
| ML-prioritized Tier 2 genes | Tier 2 genes repeatedly supported by ML feature importance or selection. |
| AHBA | Allen Human Brain Atlas; a resource for spatial brain gene-expression mapping. |
| WGCNA/module eigengene | A coexpression module method; module eigengenes summarize groups of coexpressed genes. |
