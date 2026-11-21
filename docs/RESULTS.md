# Verified existing results

No numerical output was recalculated. Sources are the frozen benchmark/synthesis manifests and the existing full-Tier2/Labonte summary CSV, inventoried with SHA-256.

- Frozen Tier2: 671 region-gene pairs; primary Tier2-RRA preserves all 671.
- Corrected weighted benchmark: 30 high, 119 moderate, 522 low; labels stay with that method.
- Frozen regional LM: zero BH <0.05 genes. Matched benchmarks: one limma-voom and nine DESeq2 same-universe BH findings. These are method-specific records, not ten independently replicated genes and not additions to Tier2.
- All 671 Tier2 pairs have concordant effect direction across LM, limma-voom and DESeq2. Genome-wide raw-p Spearman correlations with LM span 0.901–0.966. Concordance is within this cohort.
- Conditional-null pilot: strongest regional support in Ventral Striatum; partial in Hippocampal Formation and dlPFC; mixed in Insula and CG25; OFC predictive performance did not establish stable gene ranking. With 20 permutations the minimum plus-one empirical quantity is 1/21.

## Completed direct Labonte comparison

| Region | Tier2 pairs | Labonte-union matches | No observed union match | Match % |
|---|---:|---:|---:|---:|
| CG25 | 186 | 92 | 94 | 49.5 |
| OFC BA11 | 70 | 35 | 35 | 50.0 |
| dlPFC BA8/9 | 72 | 25 | 47 | 34.7 |
| Insula | 70 | 38 | 32 | 54.3 |
| Ventral Striatum | 146 | 69 | 77 | 47.3 |
| Hippocampal Formation | 127 | 41 | 86 | 32.3 |

300 of 671 pairs (44.7%) appear in the mapped sex-union nominal DEG lists; 371 have no observed match. Regional proportions vary from 32.3% to 54.3%. These are gene-region pairs, not 671 unique genes. Original versus mapped list counts and exact symbols are retained in `external_validation/rra_validation/tier2_branch_primary/simple_matching/tier2_vs_labonte_direct_summary.csv` and its report.

This disagreement is descriptive. Different preprocessing, normalization, covariates, sex-specific versus combined analyses, and statistical models mean that a stricter threshold does not imply a subset of Labonte's lists. Absence from a significance-filtered published list does not establish an untested gene's nonsignificance. This is same-dataset comparison, not independent validation.

## Completed external-reference work

Downloaded source tables and standardized data exist under `external_validation/reference_studies/`. Corrected weighted comparisons, whole-transcriptome RRA sensitivity comparisons and primary Tier2-RRA comparisons occupy separate directories. Top 10% and next 20% boundaries and all existing intersections remain unchanged. MDD-specific reference counts and the direct all-Tier2 comparison are completed local products, not planned analyses. Statistical comparison artifacts are preserved as previously computed, not rerun here.

Goes BH-FDR, Ramaker nominal p and Lanz selected-table/MDD-specific caveats must be read from their existing reference reports; coverage matches must not be presented as significant-gene replication. Independent-cohort validation of prioritized biomarkers remains unfulfilled.
