# Hippocampal Formation: final evidence synthesis

## Statistical-method concordance

All 127 frozen Tier2 pairs retained the same direction across the LM, limma-voom, and DESeq2. The regional LM–limma and LM–DESeq2 genome-wide rank correlations are documented in the frozen statistical benchmark. No new statistical consensus score was calculated.

## Tier2-RRA, LASSO, and Stability Selection

The primary order remains the frozen Tier2-RRA order. LASSO and Stability Selection are displayed as separate same-cohort benchmark evidence and do not alter that order.

## All-gene sensitivity

All-gene RRA ranks are shown only as within-study sensitivity evidence. Missing all-gene matches remain `NA`; they are not interpreted as negative findings.

## Conditional-null result

Partial support, strongest for broader top-set endpoints. This is the frozen 20-permutation pilot interpretation, not a new gene-level significance threshold.

## Leading candidates

| RRA rank | Gene        | LM rank | limma rank | DESeq2 rank | LASSO rank | LASSO frequency | Stability rank | Stability probability | All-gene RRA rank |
| -------- | ----------- | ------- | ---------- | ----------- | ---------- | --------------- | -------------- | --------------------- | ----------------- |
| 1        | TOR2A       | 26      | 43         | 109         | 1          | 1.00            | 1              | 0.71                  | 1.0               |
| 2        | TMEM161B-DT | 36      | 14         | 33          | 3          | 0.80            | 3              | 0.47                  | 5.0               |
| 3        | ARL14EP-DT  | 60      | 35         | 30          | 5          | 0.76            | 8              | 0.30                  | 70.0              |
| 4        | DIS3L-AS1   | 46      | 47         | 255         | 11         | 0.60            | 9              | 0.28                  | 34.0              |
| 5        | IL4I1       | 40      | 340        | 2581        | 2          | 0.92            | 2              | 0.50                  | 4.0               |
| 6        | SETD6       | 16      | 100        | 250         | 12         | 0.60            | 12             | 0.26                  | 6.0               |
| 7        | ADRA2B      | 38      | 370        | 59          | 9          | 0.68            | 5              | 0.41                  | 12.0              |
| 8        | LINC02334   | 123     | 174        | 237         | 10         | 0.60            | 11             | 0.26                  | 32.0              |
| 9        | ZNF663P     | 3       | 3          | 12          | 14         | 0.48            | 14             | 0.25                  | 16.0              |
| 10       | TRMT2A      | 6       | 28         | 69          | 8          | 0.72            | 4              | 0.43                  | 3.0               |

The regional Tier2-RRA leader is **TOR2A**. Its LM/limma/DESeq2 ranks are 26/43/109; its LASSO and Stability Selection ranks are 1/1. A number-one ML rank does not establish biological causality or diagnostic validity.

## Caveats

- All statistical and ML evidence comes from the same cohort.
- Tier2 was selected using the frozen LM.
- RRA scores are rank-aggregation quantities, not MDD-association p-values.
- Conditional-null resolution is limited to 20 permutations.
- WGCNA is contextual annotation only; AHBA is not used as gene-level replication.
- No result here is a validated biomarker or external replication.
