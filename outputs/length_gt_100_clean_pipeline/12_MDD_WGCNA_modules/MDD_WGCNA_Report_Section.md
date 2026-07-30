# MDD WGCNA module analysis

## Interpretation summary

The main analysis units are Ventral Striatum, Insula, Hippocampal Formation, dlPFC BA8/9, OFC BA11, and CG25. Modules are presented as exploratory prioritization and external-validation targets. No module eigengene association with Diagnosis survived BH-FDR correction; therefore, none of the modules should be described as a validated disease-associated module or biomarker.

The frontal-subregion module-summary bubble plot is intentionally excluded from the main report. Its red highlighting reflected an implementation-specific follow-up score, not WGCNA module color, Diagnosis effect direction, Diagnosis FDR significance, or the final integrated module shortlist.

## Module prioritization strategy

The final shortlist was obtained through the following sequence:

1. **Unsupervised module construction.** WGCNA modules were built from normalized expression data after genes were ranked by median absolute deviation (MAD). MAD filtering used expression variability only and did not use Diagnosis labels.
2. **Module-eigengene testing.** Each module was summarized by its module eigengene and tested using:

   `ME ~ Diagnosis + Age + Sex + RIN + PMI`

3. **Corrected Diagnosis inference.** Diagnosis p-values were corrected within each analysis unit using BH-FDR. No module survived Diagnosis BH-FDR correction. Raw p-values below 0.05 are described only as nominal.
4. **Tier2 enrichment.** Modules were evaluated for over-representation of Tier2 candidate genes. This evidence was added after unsupervised module construction.
5. **ML-prioritized Tier2 enrichment.** Modules were separately evaluated for enrichment of ML-prioritized Tier2 genes, including high-confidence and high-plus-moderate-confidence sets.
6. **Biological interpretation.** GO/MSigDB enrichment and hub genes ranked by absolute module membership, `abs(kME)`, were used to interpret module biology. Broad theme labels are keyword-based annotation summaries rather than formal mechanistic results.
7. **Limited cell-type support.** Cell-type marker enrichment was used only as a small supportive annotation layer. A missing marker result does not demonstrate absence of a cell type.
8. **Module-level ML.** Module-score features were tested for exploratory MDD/control classification signal. This analysis asks whether coordinated module-level expression carries predictive information; it does not validate a diagnostic biomarker.
9. **Integrated shortlist.** The final table contains 10 modules selected by integrating nominal module–Diagnosis context, Tier2 enrichment, ML-prioritized enrichment, biological annotation, hub-gene context, and module-level ML evidence where available. No single evidence stream alone defined the shortlist.
10. **Interpretation boundary.** All shortlisted modules are exploratory prioritization and external-validation targets. They are not validated biomarkers, causal mechanisms, or confirmed disease modules.

## Final prioritized module shortlist

This table is the main module-level result.

| Region/subregion | Module | Main prioritization evidence | Interpretation |
|---|---|---|---|
| Ventral Striatum | MElavenderblush3 | Same-region Tier2 enrichment, AHBA M6 overlap, nominal Diagnosis association | High-priority exploratory module |
| Ventral Striatum | MEdarkred | High-plus-moderate ML-prioritized enrichment | Exploratory enrichment-prioritized module; neuronal label should be treated cautiously |
| Ventral Striatum | MElightgreen | Nominal Diagnosis association and neuronal/synaptic annotation | Exploratory neuronal/synaptic module |
| Ventral Striatum | MElightcyan1 | Same-region Tier2 enrichment, AHBA overlap, excitatory-neuron marker support | Exploratory Tier2/annotation-supported module |
| Insula | MEorange | AHBA M2 overlap, nominal Diagnosis association, neuronal/axon annotation | Exploratory module; nominal status is covariate-formula-sensitive |
| Hippocampal Formation | MEgreen | Very strong same-region Tier2 enrichment | Exploratory enrichment-prioritized module with mixed glial/cytoskeletal/ECM annotation |
| dlPFC BA8/9 | MEblue | Very strong same-subregion Tier2 enrichment; substantial preservation in the 15,000-gene sensitivity run | Exploratory Tier2-prioritized module; “angiogenesis-related” is safer than “vascular/endothelial” |
| dlPFC BA8/9 | MEorange | Strongest frontal nominal Diagnosis signal in the 5,000-gene analysis, with AHBA/Tier2 context | Exploratory module; compact module identity was gene-universe-sensitive |
| CG25 | MEwhite | Strong same-subregion Tier2 enrichment | Exploratory Tier2-prioritized RNA-processing/spliceosomal module |
| OFC BA11 | MEtan | High-confidence ML-prioritized enrichment | Exploratory ML-enrichment module; prioritization was gene-universe-sensitive |

## Evidence interpretation

### Module eigengene–Diagnosis results

No prioritized module survived Diagnosis BH-FDR correction. Nominal raw-p signals can help rank exploratory targets but are not corrected evidence of disease association.

### Tier2 and ML-prioritized enrichment

Tier2 and ML-prioritized gene sets were tested after module construction. Their enrichment helps identify modules containing more candidate genes than expected under the relevant background, but it does not establish module causality or Diagnosis association.

### GO/MSigDB, hub genes, and broad themes

GO/MSigDB results and top `abs(kME)` hub genes provide biological context. Broad biological-theme labels were assigned by ordered keyword rules over top annotation terms, available marker results, and hub-gene names. They are interpretive summaries and should not replace the underlying enrichment tables.

In particular:

- dlPFC MEblue should be described as having **angiogenesis-related annotation context**, not as a confirmed vascular/endothelial module.
- Ventral Striatum MEdarkred has **mixed evidence with one neuronal-related term** and should not be called clearly neuronal.
- CG25 MEwhite is more specifically described as **RNA-processing/spliceosomal**.
- Hippocampal MEgreen is more specifically described as **mixed glial/cytoskeletal/ECM-related**.

### Cell-type annotation

The cell-marker analysis used a limited marker list and hypergeometric enrichment. Only valid marker-enrichment rows should be displayed. Cell-type results are supporting annotations, not cell-composition estimates or definitive module identities.

### Module-level ML

Module-level ML evaluates whether module-derived features carry reproducible exploratory classification signal. It should be reported separately from eigengene–Diagnosis inference and should not be used to relabel a module as Diagnosis-significant.

## Sensitivity and limitations

- No module survived Diagnosis BH-FDR correction.
- Covariate sensitivity preserved Diagnosis-effect direction for the prioritized modules, but some nominal rankings were formula-sensitive.
- The 15,000-MAD-gene frontal sensitivity analysis broadly supported dlPFC Tier2 context but showed that exact module identity and OFC prioritization depend on the gene universe.
- MAD filtering is unsupervised and does not use Diagnosis labels.
- CG25 has a smaller sample size and should be interpreted especially cautiously.
- Adaptive parameter selection and multiple downstream annotation layers increase the need for external validation.

## Main-report presentation rule

Keep the 10-module table above as the primary result. Use the corrected biological-annotation and module-network figures for supporting interpretation. Do not include `frontal_subregion_module_summary_bubble_plot.png` or its PDF/legend variants in the main report.

