# Final cross-method evidence synthesis

## 1. Purpose

This report assembles the frozen statistical, ML, sensitivity, benchmark, and conditional-null evidence for exactly 671 frozen Tier2 region–gene pairs. It does not rerun a model, change Tier2, or create a new score or confidence category.

## 2. Frozen analysis stages

The synthesis covers the original LM/Tier2 screen, all-gene ML-RRA sensitivity, primary Tier2 five-model ML-RRA, LASSO and Stability Selection benchmarks, the 20-permutation conditional-null pilot, LM/limma-voom/DESeq2 concordance, and limited WGCNA module context.

## 3. Why no additional scoring was introduced

The evidence streams answer different questions and use different scales. Combining them numerically would impose new, post-hoc weights. The frozen Tier2-RRA rank therefore remains the sole candidate order, while every other field remains visible separately.

## 4. Statistical robustness of Tier2

All 671 Tier2 pairs retained the same effect direction across the LM, limma-voom, and DESeq2. Genome-wide raw-p rank correlations were strong in all regions. Tier2 is therefore robust to statistical-model choice as an exploratory screen, but this is within-study concordance rather than replication. No frozen LM gene survived BH FDR.

## 5. ML prioritisation benchmark results

Tier2-RRA is an exploratory multivariate prioritisation. LASSO and Stability Selection often converge with its leading candidates, but neither supplies independent validation. The all-gene RRA comparison is a sensitivity analysis using the same cohort.

## 6. Conditional-null findings

- **Ventral Striatum:** clearest conditional RRA/LASSO stability support.
- **Hippocampal Formation:** partial support.
- **dlPFC BA8/9:** partial support.
- **Insula:** mixed evidence.
- **CG25:** mixed evidence.
- **OFC BA11:** Scenario E-like; predictive performance exceeded null but ranking stability did not.

The pilot used 20 permutations, so its smallest possible plus-one empirical quantity is 1/21. It supports region-specific caution, not validated gene calls.

## 7. All-gene sensitivity results

All-gene ranks are retained as separate sensitivity evidence. Strong convergence can support a narrative example, but a weak or missing all-gene rank does not automatically invalidate a Tier2 candidate.

## 8. Region-by-region evidence

Detailed summaries are in `22_Final_evidence_synthesis/region_summaries/`. They preserve the frozen regional interpretation and list the top ten candidates without reranking.

## 9. Leading candidate evidence profiles

| Region                | Leader       | LM rank | limma rank | DESeq2 rank | LASSO rank | LASSO frequency | Stability rank | Stability probability | All-gene RRA rank | Null top-10 / 20 |
| --------------------- | ------------ | ------- | ---------- | ----------- | ---------- | --------------- | -------------- | --------------------- | ----------------- | ---------------- |
| Ventral Striatum      | MORN5        | 26      | 25         | 15          | 1          | 1.0             | 1              | 0.89                  | 1.0               | 4.0              |
| Insula                | TMEM161B-DT  | 7       | 3          | 5           | 1          | 0.96            | 1              | 0.75                  | 2.0               | 10.0             |
| Hippocampal Formation | TOR2A        | 26      | 43         | 109         | 1          | 1.0             | 1              | 0.71                  | 1.0               | 0.0              |
| CG25                  | LOC107985795 | 23      | 24         | 349         | 1          | 0.8             | 1              | 0.33                  | 3.0               | 2.0              |
| dlPFC BA8/9           | LINC00269    | 13      | 41         | 36          | 7          | 0.8             | 3              | 0.64                  | 2.0               | 11.0             |
| OFC BA11              | LOC105373386 | 45      | 49         | 343         | 3          | 0.8             | 7              | 0.35                  | 6.0               | 8.0              |

The leaders are useful audit examples, not proven biological winners. MORN5 has the most favorable regional null context. TOR2A and LINC00269 retain partial regional support. TMEM161B-DT and LOC107985795 require mixed-evidence caveats. LOC105373386 requires the strongest caution because OFC ranking stability did not exceed the conditional null and its DESeq2 rank is weaker.

## 10. Method-specific FDR sensitivity findings

The separate file `method_specific_FDR_findings.csv` contains 10 method-specific records: one limma-voom finding and nine same-universe DESeq2 findings. These include OFC `QPRT` and the eight Ventral Striatum DESeq2 findings. They are not promoted into Tier2 or a new shortlist.

## 11. Biological-context annotations

Unambiguous Phase-3B WGCNA module assignments and module diagnosis p/BH values are appended as context only. Missing module assignments remain `NA`. No module survived diagnosis BH FDR. AHBA is not represented as a gene-level replication field because no unambiguous gene–region linkage suitable for that claim was available; it remains context only.

## 12. Main limitations

1. All evidence streams use the same cohort.
2. Tier2 was defined from the original LM, so downstream concordance is not independent.
3. Regional sample sizes are small, especially CG25.
4. Statistical and ML methods answer different questions and use different scales.
5. Conditional-null resolution is limited to 20 permutations.
6. WGCNA and AHBA are contextual, not independent MDD replication.
7. No external cohort is included.

## 13. What can be claimed now

- Tier2 is robust to statistical-model choice as an exploratory candidate screen.
- Tier2-RRA supplies an exploratory multivariate prioritisation of frozen Tier2 candidates.
- RRA broadly converges with LASSO and Stability Selection for several leading candidates.
- Ventral Striatum shows the clearest conditional evidence for label-dependent prioritisation.
- Several candidates show within-study statistical and ML convergence.

## 14. What cannot be claimed now

The project has not established a validated biomarker, diagnostic gene, causal MDD gene, clinically predictive model, proven mechanism, independent replication, or external validation. Method-specific FDR findings cannot be described as canonical LM significance.

## 15. Recommended next scientific step

The next major scientific step is preregistered external validation in an independent cohort, preserving the frozen Tier2 identities and primary Tier2-RRA ordering before looking at external outcomes.
