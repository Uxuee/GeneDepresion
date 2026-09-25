# Publication audit numerical checks

Machine-readable evidence: `docs/PUBLICATION_AUDIT_VERIFICATION.json`. Verification code: `audit/publication_audit_verification.py`. Canonical baseline files were read from commit `ca25ec35af359a4d409abf0f728833d2ff54761c` when the working-tree copy was deleted or modified.

## Summary

| Check | Reconstructed | Canonical | Absolute difference | Result |
|---|---:|---:|---:|---|
| GeTMM log value | 0.6842735183893409 | 0.6842735183893410 | 1.11e-16 | PASS |
| Positive diagnosis beta | 0.3437465610993611 | 0.3437465610993620 | 9.44e-16 | PASS |
| Negative diagnosis beta | -0.3161932734289576 | -0.3161932734289570 | 6.11e-16 | PASS |
| Insula BH vector | full vector | full vector | max 1.11e-16 | PASS |
| Standardized effect | 0.7010544075559330 | 0.7010544075559327 | 2.22e-16 | PASS |
| Tier1/2/3 memberships | 95,319 rows | 95,319 rows | 0 mismatches | PASS |
| Tier4 score example | 5.380829883783952 | 5.380829883783952 | 0 | PASS formula only |
| Tier2 GO p | 0.000240982637150624 | 0.000240982637150600 | 2.36e-17 | PASS |
| ML all scores | 671 rows | 671 rows | max 9.71e-17 | PASS |
| AHBA M26 p | 0.02801303788284664 | 0.02801303788284660 | 4.16e-17 | PASS |
| WGCNA MEdarkred p | 0.00217013628974135 | 0.00217013628974130 | 5.25e-17 | PASS |

## 1. GeTMM/log normalization

Sources: raw GSE102556 count table, NCBI annotation length, canonical VS normalization-factor table, and canonical normalized expression matrix. Implementation: `Step1_length_cutoff_sensitivity_GeTMM.R`.

Example: Entrez `100287102`, sample `GSM2740799`.

```text
raw count                       = 5
gene length                     = 1652 bp = 1.652 kb
RPK                             = 5 / 1.652 = 3.026634382566586
sum RPK in sample               = 4778316.68550817
TMM factor                      = 1.04369398492787
effective library              = 4778316.68550817 * 1.04369398492787
                                = 4987100.38274535
GeTMM                           = RPK / effective library * 1e6
                                = 0.6068926129977864
log2(GeTMM + 1)                = 0.6842735183893409
canonical                       = 0.6842735183893410
absolute difference             = 1.1102230246251565e-16
```

Result: PASS.

## 2. Diagnosis model coefficients

Sources: Insula normalized matrix, saved metadata, and canonical gene-effect table. Formula: `expression ~ Diagnosis + Age + Sex + RIN`, with Control and Female as references; n=48, residual df=43.

| Gene | Direction | Reconstructed beta | Canonical beta | Reconstructed SE | Canonical SE | t | Reconstructed p | Canonical p |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| BMS1P2-AGAP9 / 113939912 | MDD-up | 0.343746561099361 | 0.343746561099362 | 0.0848350442977467 | 0.0848350442977468 | 4.0519406095 | 0.000208838853615768 | 0.000208838853615700 |
| UST / 10090 | MDD-down | -0.316193273428958 | -0.316193273428957 | 0.0738470496685215 | 0.0738470496685216 | -4.2817319696 | 0.000102019507615156 | 0.000102019507615100 |

Result: PASS.

## 3. BH correction

The complete 23,090-value Insula diagnosis-p vector was adjusted with the standard Benjamini–Hochberg step-up procedure and compared row-for-row with the canonical output. Maximum absolute difference: `1.11e-16`; minimum adjusted p: `0.9336862754678356`. Family: all filtered/tested genes in Insula. Result: PASS.

## 4. Standardized diagnosis effect and Tier check

Source row: pooled Frontal Cortex Entrez `5147`.

```text
beta                            = 0.0968409325231598
within-unit expression SD       = 0.1381361153705227
beta / SD                       = 0.7010544075559330
canonical standardized effect   = 0.7010544075559327
absolute difference             = 2.22e-16
```

Its diagnosis p satisfies p<0.01, but `abs(effect)=0.7011` does not exceed 0.8; therefore it is not Tier2. Full broad-table recomputation found zero membership mismatches for Tier1, Tier2, or Tier3 across 95,319 rows. Result: PASS.

## 5. Tier4 score

For the same row:

```text
signed_score = sign(0.7010544) * -log10(p) = 5.380829883783952
stored score                                      = 5.380829883783952
```

Result: PASS for the formula. A canonical Tier4 ranking/GSEA package was not present, so ranking uniqueness, ties, zero-p handling, GSEA permutations, and GSEA FDR cannot be numerically certified.

## 6. Tier2 enrichment

Source: `08_Tier2_enrichment/tables/four_group__Hippocampal_Formation__down_in_MDD__GO_MF__full.csv`, GO:0017017.

```text
N (GO-mapped tested background) = 15,343
K (term genes in background)    = 12
n (mapped candidate genes)      = 30
k (overlap)                     = 2
p = P[X >= 2], X~Hypergeom(N,K,n)
reconstructed p                 = 0.00024098263715062357
canonical p                     = 0.00024098263715060000
absolute difference             = 2.36e-17
BH over 122 saved GO:MF rows    = 0.0265974725010149
canonical BH                    = 0.0265974725010151
```

Result: PASS. The denominator 30 comes from the parsed GeneRatio; each database uses the subset of the same-unit tested universe that maps into that database.

## 7. v1 and v2 selection frequencies

The v1 raw file checksum is `406ba687a9f3f4244de0d296550e1b2fcfb043f504c175c56328b447133cde1e`; v2 is `7b33e7d0a6eb2ff2efaffdc223aced7f4722f18236a2902cc9102257b2215ef8`.

For LOC107985795/CG25:

```text
v1 = 108 / 125 = 0.864
v2 =  70 / 150 = 0.4666666666666667
combined = 0.40*v1 + 0.60*v2 = 0.6256
```

PRELID2/Insula: v1 `34/125=0.272`; v2 `30/150=0.20`; combined `0.2288`, below the moderate gate of 0.25.

## 8. Complete ML-priority scores

Formula:

```text
0.40 combined frequency + 0.15 model support + 0.10 sparse support
+ 0.15 importance rank + 0.10 sign consistency + 0.05 LMM evidence
+ 0.05 v2 support
```

| Gene/unit | v1 events | v2 events | Combined | Model | Sparse | Rank | Sign | LMM | v2 support | Reconstructed score | Class |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| LOC107985795 / CG25 | 108 | 70 | 0.6256 | 1 | 1 | 0.799251273 | 1 | 0.493570483 | 1 | 0.7948062150587105 | high |
| PRELID2 / Insula | 34 | 30 | 0.2288 | 1 | 1 | 0.348269771 | 1 | 0.376898396 | 1 | 0.5626053853815516 | exploratory/low (frequency gate fails) |
| TMEM161B-DT / Insula | 116 | 90 | 0.7312 | 1 | 1 | 0.860214143 | 1 | 0.839000656 | 1 | 0.8634621542622630 | high |
| LOC105372507 / CG25 | 95 | 51 | 0.5080 | 1 | 1 | 0.685500841 | 1 | 0.191820320 | 1 | 0.7156161420787578 | high; near high boundary |
| GLRA4 / VS | 85 | 0 | 0.2720 | 1 | 1 | 0.563379036 | 1 | 0.179875500 | 0 | 0.5523006303639477 | moderate |

Every reconstructed score equals its CSV value to machine precision; maximum difference over all 671 rows is `9.71e-17`. Counts: 30 high, 119 moderate, 522 low. Result: PASS.

## 9. Published-AHBA enrichment

Corrected VS/up-in-MDD/M26 high+moderate set:

```text
N background                       = 12,042
K M26 genes in background          = 554
n candidate genes in background    = 6
k overlap                          = 2 (MORN5; RPGR)
hypergeometric p                   = 0.02801303788284664
canonical p                        = 0.02801303788284660
BH                                 = 0.8822833419054857
```

Result: PASS; nominal only. A null example is the corrected VS/down-in-MDD/M18 row after ZKSCAN2 loss: k=0, p=1, BH=1.

## 10. WGCNA enrichment

Corrected VS MEdarkred/global high+moderate:

```text
N network genes                    = 15,000
K module genes                     = 227
n corrected set genes in universe  = 88
k overlap                          = 6
one-sided Fisher p                 = 0.00217013628974135
canonical p                        = 0.00217013628974130
BH                                 = 0.1085068144870676
```

Result: PASS; nominal only. OFC BA11 MEtan/global high has overlap 2, p=0.001780636064698, BH=0.0320514491645654 and remains FDR<0.10.

## 11. Correction transitions

Legacy counts 35/130/506 changed to 30/119/522. Retained comparison reconstructs five high→moderate, sixteen moderate→low, zero upward transitions, and 650 unchanged rows. M26 lost TSGA10. MEdarkred lost PRELID2 and TMEM269-DT.
