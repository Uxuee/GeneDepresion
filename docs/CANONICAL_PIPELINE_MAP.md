# Canonical pipeline map

This is the human-readable companion to `provenance/canonical_pipeline.json`.
The JSON registry is authoritative for stage status, inputs, outputs, hashes,
dependencies, seeds, environment references, and replacements.

## First-page stage index

| Stage | Status | Input → script → canonical output | Milestone |
|---|---|---|---|
| 01 source data | CANONICAL EXTERNAL INPUT | GEO/AHBA resources → acquisition methods → machine-local checked files | external |
| 02 filtering + GeTMM | CANONICAL | stage 01 → `Step1_length_cutoff_sensitivity_GeTMM.R` → Length >100 normalized data | `ca25ec3` |
| 03 QC/covariates | CANONICAL | stage 02 → `Step2_Covariate_Exploration_FrontalCortex.R` → QC/metadata outputs | `ca25ec3` |
| 04 diagnosis models | CANONICAL | stages 02–03 → registered gene-wise scripts → frozen diagnosis tables | `ca25ec3` |
| 05 standardized effects | CANONICAL | stage 04 → clean step 04 wrapper → standardized-effect table | `ca25ec3` |
| 06 tiers | CANONICAL | stage 05 → clean step 05 wrapper → Tier2 candidate table | `ca25ec3` |
| 07 Tier2 enrichment | CANONICAL | stage 06 → clean step 08 wrapper → Tier2 enrichment folder | `ca25ec3` |
| 08 Tier4 ranking/GSEA | CANONICAL ranking; GSEA EXCLUDED | stage 05 → recovered Git generator → exact ranked Git artifacts; GSEA not recovered | `7888e98` |
| 09 ML v1 | CANONICAL | stage 06 → frozen snapshot generator → fold-level event file | `fa3b53f` |
| 10 ML v2 | CANONICAL | stage 06 → frozen snapshot generator → selected-gene event file | `fa3b53f` |
| 11 corrected priority | CANONICAL | stages 06/09/10 → corrected scorer → corrected 671-row ranking | `fa3b53f` |
| 12 AHBA enrichment | CANONICAL | stages 01/11 → corrected downstream generator → `AHBA_corrected` | `ca25ec3` |
| 13 WGCNA | CANONICAL | stages 02/03 → Step 6 WGCNA → frozen module assignments | `ca25ec3` |
| 14 corrected WGCNA ML enrichment | CANONICAL | stages 11/13 → corrected downstream generator → corrected WGCNA ML annotations | `ca25ec3` |
| 15 structural core/network | CANONICAL | stages 11/13 → corrected downstream generator → corrected core/network annotations | `ca25ec3` |
| 15b WGCNA biological annotation | UNDER RECOVERY | stage 13 → Step 7 annotation → historical low-coverage annotation | none |
| 16 AHBA module-score ML | CANONICAL all-32 only | stages 01/02 → AHBA score-ML script → module-score ML folder | `ca25ec3` |
| 17 module/pathway ML2 | EXCLUDED | retained outputs; exact generator missing | none |
| 18 figures/tables | UNDER RECOVERY | verified canonical stages → controlled assembly not yet frozen | none |
| 19 manuscript/report | UNDER RECOVERY | stage 18 → publication refresh not yet performed | none |

The registry contains full paths and SHA-256 values; abbreviated commit labels
above are for readability only.

## Dependency diagram

```text
01 external data
  └─> 02 normalization ─> 03 QC ─> 04 diagnosis models ─> 05 effects ─> 06 Tier2
                                                │                         ├─> 07 enrichment
                                                └─> 08 Tier4 ranking      ├─> 09 ML v1 ─┐
                                                                          └─> 10 ML v2 ─┴─> 11 corrected priority
01 AHBA ────────────────────────────────────────────────────────────────────────────────┬─> 12 corrected AHBA
02 + 03 ─> 13 WGCNA ────────────────────────────────────────────────┬───────────────────┼─> 14 corrected WGCNA ML
                                                                  ├─> 15 structural core/network
                                                                  └─> 15b annotation (pending)
01 + 02 ─> 16 AHBA all-32 module-score ML
17 module/pathway ML2 (excluded)     18 figures/tables (pending) ─> 19 report (pending)
```

## Publication path

1. External GSE102556, annotation, and AHBA resources.
2. Length >100 filtering, GeTMM normalization, and expression QC.
3. Gene-wise diagnosis models and standardized effects.
4. Exploratory Tier2 candidates and pathway enrichment.
5. ML v1 and strict nested v2 event records.
6. Corrected ML-priority score: 671 gene-region pairs = 30 high, 119 moderate, 522 low. The canonical high+moderate set contains 149 pairs. Denominators are 125 (v1) and 150 (v2).
7. Corrected-list-dependent downstream outputs under `16_corrected_ML_downstream`.
8. MDD-derived WGCNA using original report units: Ventral Striatum, Insula, Hippocampal Formation, dlPFC BA8/9, OFC BA11, and CG25. Broad pooled Frontal Cortex is supplementary only.
9. Published AHBA modules provide healthy-brain context. The all-32 module-score feature set is active; FDR005, FDR010, and M26-only are historical post-selection sensitivities.
10. Publication assembly after all pending provenance and annotation limitations are carried into captions and text.

## Status decisions

- Tier4 ranked signatures: exactly recoverable from Git objects; the recovery manifest records every object and SHA-256.
- Tier4 GSEA: excluded from publication because no exact generator, settings, or result artifacts were recovered.
- WGCNA biological annotation: pending controlled recovery. Historical identifier mapping covered 1,796/4,267 rows (42.1%), so it remains limited supportive interpretation.
- MDD module-level ML2 (historical Step 13): excluded from publication because its exact generator was not recovered.
- Independent external validation: not complete.

## Scientific interpretation guardrails

There are no validated biomarkers. Tier2, ML-prioritized genes, WGCNA modules,
AHBA overlaps, and module-level prediction are exploratory. No MDD WGCNA module
survived Diagnosis BH-FDR correction. AHBA is context, not MDD validation.

Method terminology: Ventral Striatum, Insula, Hippocampal Formation, CG25,
dlPFC BA8/9, and OFC BA11 use ordinary gene-wise linear models in the relevant
region/subregion branch. The pooled Frontal Cortex branch uses a mixed-effects
model with `(1 | BrainID)`. Historical filenames containing “LMM” are retained
as provenance and do not change this distinction.

Use `python code/run_canonical_pipeline.py --status` to inspect the registry and
`python code/run_canonical_pipeline.py --validate` to validate it. Historical
steps 9–12 in the old clean runner require the explicit `--allow-legacy` flag.
