# Insula: final evidence synthesis

## Statistical-method concordance

All 70 frozen Tier2 pairs retained the same direction across the LM, limma-voom, and DESeq2. The regional LM–limma and LM–DESeq2 genome-wide rank correlations are documented in the frozen statistical benchmark. No new statistical consensus score was calculated.

## Tier2-RRA, LASSO, and Stability Selection

The primary order remains the frozen Tier2-RRA order. LASSO and Stability Selection are displayed as separate same-cohort benchmark evidence and do not alter that order.

## All-gene sensitivity

All-gene RRA ranks are shown only as within-study sensitivity evidence. Missing all-gene matches remain `NA`; they are not interpreted as negative findings.

## Conditional-null result

Mixed stability evidence; suggestive but inconsistent across endpoints. This is the frozen 20-permutation pilot interpretation, not a new gene-level significance threshold.

## Leading candidates

| RRA rank | Gene         | LM rank | limma rank | DESeq2 rank | LASSO rank | LASSO frequency | Stability rank | Stability probability | All-gene RRA rank |
| -------- | ------------ | ------- | ---------- | ----------- | ---------- | --------------- | -------------- | --------------------- | ----------------- |
| 1        | TMEM161B-DT  | 7       | 3          | 5           | 1          | 0.96            | 1              | 0.75                  | 2.0               |
| 2        | LOC101929710 | 97      | 59         | 350         | 2          | 0.96            | 2              | 0.64                  | 9.0               |
| 3        | BMS1P2       | 5       | 6          | 4           | 6          | 0.80            | 13             | 0.29                  | 1.0               |
| 4        | ARL14EP-DT   | 28      | 57         | 13          | 3          | 0.96            | 3              | 0.55                  | 5.0               |
| 5        | BTN2A3P      | 54      | 27         | 127         | 8          | 0.72            | 7              | 0.41                  | 13.0              |
| 6        | LOC105379042 | 103     | 91         | 133         | 14         | 0.60            | 8              | 0.41                  | 64.0              |
| 7        | DECR2        | 42      | 78         | 208         | 4          | 0.88            | 4              | 0.55                  | 133.0             |
| 8        | LINC01691    | 38      | 84         | 113         | 5          | 0.84            | 5              | 0.53                  | 6.0               |
| 9        | METTL14-DT   | 87      | 158        | 156         | 9          | 0.72            | 12             | 0.32                  | 34.0              |
| 10       | BMS1P2-AGAP9 | 4       | 5          | 3           | 20         | 0.32            | 9              | 0.40                  | 4.0               |

The regional Tier2-RRA leader is **TMEM161B-DT**. Its LM/limma/DESeq2 ranks are 7/3/5; its LASSO and Stability Selection ranks are 1/1. A number-one ML rank does not establish biological causality or diagnostic validity.

## Caveats

- All statistical and ML evidence comes from the same cohort.
- Tier2 was selected using the frozen LM.
- RRA scores are rank-aggregation quantities, not MDD-association p-values.
- Conditional-null resolution is limited to 20 permutations.
- WGCNA is contextual annotation only; AHBA is not used as gene-level replication.
- No result here is a validated biomarker or external replication.
