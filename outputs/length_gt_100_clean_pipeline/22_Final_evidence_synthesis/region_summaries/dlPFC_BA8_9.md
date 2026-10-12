# dlPFC BA8/9: final evidence synthesis

## Statistical-method concordance

All 72 frozen Tier2 pairs retained the same direction across the LM, limma-voom, and DESeq2. The regional LM–limma and LM–DESeq2 genome-wide rank correlations are documented in the frozen statistical benchmark. No new statistical consensus score was calculated.

## Tier2-RRA, LASSO, and Stability Selection

The primary order remains the frozen Tier2-RRA order. LASSO and Stability Selection are displayed as separate same-cohort benchmark evidence and do not alter that order.

## All-gene sensitivity

All-gene RRA ranks are shown only as within-study sensitivity evidence. Missing all-gene matches remain `NA`; they are not interpreted as negative findings.

## Conditional-null result

Partial support, with stronger broader-top-set and LASSO evidence. This is the frozen 20-permutation pilot interpretation, not a new gene-level significance threshold.

## Leading candidates

| RRA rank | Gene         | LM rank | limma rank | DESeq2 rank | LASSO rank | LASSO frequency | Stability rank | Stability probability | All-gene RRA rank |
| -------- | ------------ | ------- | ---------- | ----------- | ---------- | --------------- | -------------- | --------------------- | ----------------- |
| 1        | LINC00269    | 13      | 41         | 36          | 7          | 0.80            | 3              | 0.64                  | 2.0               |
| 2        | NOL6         | 108     | 424        | 475         | 2          | 0.96            | 1              | 0.68                  | 21.0              |
| 3        | PPT2         | 16      | 79         | 103         | 1          | 1.00            | 2              | 0.66                  | 3.0               |
| 4        | LINC01990    | 162     | 176        | 172         | 6          | 0.80            | 7              | 0.45                  | 27.0              |
| 5        | LOC107984015 | 7       | 15         | 43          | 3          | 0.88            | 5              | 0.49                  | 4.0               |
| 6        | LINC02897    | 26      | 39         | 150         | 8          | 0.76            | 6              | 0.47                  | 5.0               |
| 7        | EMC6         | 2       | 17         | 26          | 10         | 0.68            | 11             | 0.38                  | 1.0               |
| 8        | MCPH1-DT     | 27      | 81         | 86          | 4          | 0.84            | 10             | 0.39                  | 36.0              |
| 9        | PPT2-EGFL8   | 48      | 190        | 197         | 16         | 0.40            | 13             | 0.33                  | 8.0               |
| 10       | LOC105371009 | 38      | 118        | 152         | 5          | 0.84            | 4              | 0.51                  | 16.0              |

The regional Tier2-RRA leader is **LINC00269**. Its LM/limma/DESeq2 ranks are 13/41/36; its LASSO and Stability Selection ranks are 7/3. A number-one ML rank does not establish biological causality or diagnostic validity.

## Caveats

- All statistical and ML evidence comes from the same cohort.
- Tier2 was selected using the frozen LM.
- RRA scores are rank-aggregation quantities, not MDD-association p-values.
- Conditional-null resolution is limited to 20 permutations.
- WGCNA is contextual annotation only; AHBA is not used as gene-level replication.
- No result here is a validated biomarker or external replication.
