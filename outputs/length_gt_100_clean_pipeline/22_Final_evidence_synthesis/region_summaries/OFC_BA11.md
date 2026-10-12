# OFC BA11: final evidence synthesis

## Statistical-method concordance

All 70 frozen Tier2 pairs retained the same direction across the LM, limma-voom, and DESeq2. The regional LM–limma and LM–DESeq2 genome-wide rank correlations are documented in the frozen statistical benchmark. No new statistical consensus score was calculated.

## Tier2-RRA, LASSO, and Stability Selection

The primary order remains the frozen Tier2-RRA order. LASSO and Stability Selection are displayed as separate same-cohort benchmark evidence and do not alter that order.

## All-gene sensitivity

All-gene RRA ranks are shown only as within-study sensitivity evidence. Missing all-gene matches remain `NA`; they are not interpreted as negative findings.

## Conditional-null result

Scenario E-like: prediction exceeded null, but gene-ranking stability did not. This is the frozen 20-permutation pilot interpretation, not a new gene-level significance threshold.

## Leading candidates

| RRA rank | Gene         | LM rank | limma rank | DESeq2 rank | LASSO rank | LASSO frequency | Stability rank | Stability probability | All-gene RRA rank |
| -------- | ------------ | ------- | ---------- | ----------- | ---------- | --------------- | -------------- | --------------------- | ----------------- |
| 1        | LOC105373386 | 45      | 49         | 343         | 3          | 0.80            | 7              | 0.35                  | 6.0               |
| 2        | LOC105378913 | 24      | 42         | 118         | 1          | 0.96            | 1              | 0.48                  | 5.0               |
| 3        | RAB24        | 110     | 310        | 220         | 4          | 0.76            | 15             | 0.29                  | 11.0              |
| 4        | TRPV3        | 46      | 47         | 56          | 10         | 0.68            | 9              | 0.33                  | 33.0              |
| 5        | C9orf163     | 62      | 56         | 221         | 5          | 0.76            | 3              | 0.39                  | 26.0              |
| 6        | QPRT         | 1       | 1          | 1           | 2          | 0.92            | 4              | 0.37                  | 10.0              |
| 7        | GCC1         | 92      | 292        | 255         | 8          | 0.68            | 6              | 0.37                  | 27.0              |
| 8        | LOC105377975 | 57      | 80         | 107         | 6          | 0.72            | 10             | 0.32                  | 16.0              |
| 9        | STAP2        | 103     | 507        | 237         | 11         | 0.64            | 5              | 0.37                  | 18.0              |
| 10       | ZNF461       | 17      | 55         | 60          | 9          | 0.68            | 2              | 0.40                  | 44.0              |

The regional Tier2-RRA leader is **LOC105373386**. Its LM/limma/DESeq2 ranks are 45/49/343; its LASSO and Stability Selection ranks are 3/7. A number-one ML rank does not establish biological causality or diagnostic validity.

## Caveats

- All statistical and ML evidence comes from the same cohort.
- Tier2 was selected using the frozen LM.
- RRA scores are rank-aggregation quantities, not MDD-association p-values.
- Conditional-null resolution is limited to 20 permutations.
- WGCNA is contextual annotation only; AHBA is not used as gene-level replication.
- No result here is a validated biomarker or external replication.
