# Ventral Striatum: final evidence synthesis

## Statistical-method concordance

All 146 frozen Tier2 pairs retained the same direction across the LM, limma-voom, and DESeq2. The regional LM–limma and LM–DESeq2 genome-wide rank correlations are documented in the frozen statistical benchmark. No new statistical consensus score was calculated.

## Tier2-RRA, LASSO, and Stability Selection

The primary order remains the frozen Tier2-RRA order. LASSO and Stability Selection are displayed as separate same-cohort benchmark evidence and do not alter that order.

## All-gene sensitivity

All-gene RRA ranks are shown only as within-study sensitivity evidence. Missing all-gene matches remain `NA`; they are not interpreted as negative findings.

## Conditional-null result

Clearest conditional RRA/LASSO stability support. This is the frozen 20-permutation pilot interpretation, not a new gene-level significance threshold.

## Leading candidates

| RRA rank | Gene         | LM rank | limma rank | DESeq2 rank | LASSO rank | LASSO frequency | Stability rank | Stability probability | All-gene RRA rank |
| -------- | ------------ | ------- | ---------- | ----------- | ---------- | --------------- | -------------- | --------------------- | ----------------- |
| 1        | MORN5        | 26      | 25         | 15          | 1          | 1.00            | 1              | 0.89                  | 1.0               |
| 2        | POF1B        | 193     | 43         | 115         | 6          | 0.96            | 4              | 0.65                  | 5.0               |
| 3        | FSHB         | 94      | 12         | 64          | 4          | 1.00            | 3              | 0.66                  | 2.0               |
| 4        | MRPL2        | 78      | 189        | 256         | 2          | 1.00            | 2              | 0.85                  | 3.0               |
| 5        | MOSPD3       | 42      | 106        | 244         | 3          | 1.00            | 5              | 0.64                  | 4.0               |
| 6        | KRT5         | 225     | 230        | 381         | 8          | 0.80            | 9              | 0.31                  | 7.0               |
| 7        | LOC107985657 | 424     | 301        | 1366        | 7          | 0.80            | 7              | 0.37                  | 28.0              |
| 8        | ZFAT-AS1     | 435     | 555        | 747         | 13         | 0.40            | 12             | 0.25                  | 14.0              |
| 9        | LINC01361    | 134     | 284        | 793         | 9          | 0.76            | 6              | 0.37                  | 6.0               |
| 10       | BMS1P2       | 18      | 6          | 17          | 5          | 1.00            | 8              | 0.36                  | 8.0               |

The regional Tier2-RRA leader is **MORN5**. Its LM/limma/DESeq2 ranks are 26/25/15; its LASSO and Stability Selection ranks are 1/1. A number-one ML rank does not establish biological causality or diagnostic validity.

## Caveats

- All statistical and ML evidence comes from the same cohort.
- Tier2 was selected using the frozen LM.
- RRA scores are rank-aggregation quantities, not MDD-association p-values.
- Conditional-null resolution is limited to 20 permutations.
- WGCNA is contextual annotation only; AHBA is not used as gene-level replication.
- No result here is a validated biomarker or external replication.
