> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Reproducibility and Publication-Readiness Checklist

| Item | Status | Evidence/action |
|---|---|---|
| Normalisation | Complete | Saved Length >100 GeTMM/log2 matrices and scripts |
| Statistical modelling | Complete | Saved gene-wise outputs; preserve exact formulas |
| QC/covariates | Needs documentation | Consolidate QC exclusions, missingness, rationale |
| Tier construction | Complete | Scripts/tables exist; Tier 2 remains exploratory |
| Tier 2 enrichment | Complete scientifically; needs repository cleanup | Resolve deleted tracked outputs before release |
| ML v1 event provenance | Complete | Exact file restored; checksum verified |
| Strict ML v2 event provenance | Complete | Exact file restored; checksum verified |
| Corrected ML prioritisation | Complete | 671 rows; 30/119/522; tests/reconstructions pass |
| Ranked GSEA | Needs documentation | Identify canonical outputs and distinguish from ORA |
| AHBA/WGCNA mapping | Needs rerun | Preserve WGCNA modules; refresh legacy-list-dependent overlays |
| Pathway/module-level ML | Needs rerun | Review and refresh affected inputs/exports |
| External validation | Future validation | Independent validation not complete |
| Final figures | Needs rerun | Replace legacy ML-derived panels |
| Manuscript methods | Needs documentation | Use final score specification and limitations |
| Reproducible environment | Needs documentation | Freeze Python/R versions and packages |
| Data availability/licensing | Needs documentation | GEO, redistribution, derived-data policy |
| Git release state | Needs documentation | Resolve unrelated dirty/deleted files; record approved commit |

## Corrected-score checks

- [x] Five scoring regression tests pass.
- [x] Ranking has 671 rows and classes 30/119/522.
- [x] Non-missing v1 denominators are 125.
- [x] Non-missing v2 denominators are 150.
- [x] Majority baseline is removed before v2 scoring.
- [x] Random forest is excluded from signed evidence.
- [x] Missing sign cannot pass the high gate.
- [x] Representative arithmetic matches the CSV within floating-point precision.
