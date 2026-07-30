# Module prioritization strategy

1. WGCNA modules were built from normalized expression data after unsupervised MAD-based gene filtering.
2. Module eigengenes were tested with `ME ~ Diagnosis + Age + Sex + RIN + PMI`.
3. No module survived Diagnosis BH-FDR correction; raw p-values below 0.05 are nominal only.
4. Modules were evaluated for Tier2 candidate-gene enrichment.
5. Modules were evaluated for enrichment of ML-prioritized Tier2 genes.
6. GO/MSigDB enrichment and hub genes ranked by `abs(kME)` supported biological interpretation.
7. Cell-type marker enrichment was used only as limited supportive annotation.
8. Module-level ML tested whether module-score features carried exploratory MDD/control signal.
9. Ten modules were shortlisted by integrating these evidence streams rather than using any one score or threshold.
10. The modules are exploratory prioritization and external-validation targets, not validated biomarkers.

## Main-report instruction

Keep the final 10-module prioritized table as the main result. Do not include the frontal-subregion module-summary bubble plot. Its red highlighting was an implementation-specific follow-up flag and did not encode WGCNA module color, effect direction, Diagnosis FDR significance, or final module selection.

