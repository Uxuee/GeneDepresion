# GeneDepresion

## Canonical publication workflow

This repository preserves the corrected and provenance-controlled MDD
transcriptomic workflow. Start with these files:

1. [Canonical pipeline map](docs/CANONICAL_PIPELINE_MAP.md)
2. [Machine-readable pipeline registry](provenance/canonical_pipeline.json)
3. [Local path configuration example](config/project_paths.example.yaml)
4. [External-data acquisition and placement](data/README.md)
5. [Environment and test setup](environment/README.md)
6. [Registry-driven runner](code/run_canonical_pipeline.py)

Basic inspection commands:

```text
python code/run_canonical_pipeline.py --status
python code/run_canonical_pipeline.py --validate
python code/run_canonical_pipeline.py --dry-run
python code/run_canonical_pipeline.py --check-external-data
```

Copy `config/project_paths.example.yaml` to
`config/project_paths.local.yaml`, edit the local paths, and keep that local
file untracked. Large or restricted source data are not redistributed through
Git; [data/README.md](data/README.md) explains how to acquire and verify them.

## Canonical result identity

The corrected ML-prioritisation result contains **671 gene-region pairs: 30
high, 119 moderate, and 522 low confidence** (149 high+moderate). The event
denominators are 125 for v1 and 150 for v2. Its canonical location is:

```text
outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/
```

The similarly named `10_ML_prioritized_Tier2_extended` directory is historical
provenance only and must not be used as the current result.

## Scientific guardrails

- The six original region/subregion gene-wise analyses—Ventral Striatum,
  Insula, Hippocampal Formation, CG25, dlPFC BA8/9, and OFC BA11—use
  fixed-effect linear models. The supplementary pooled Frontal Cortex branch
  uses a mixed-effects model to account for repeated donors.
- No gene-wise Diagnosis result survives BH-FDR correction.
- No WGCNA module-eigengene Diagnosis result survives BH-FDR correction.
- AHBA is healthy-brain context only. Corrected M26 evidence is nominal, not a
  validated disease result.
- Ventral Striatum MEdarkred corrected high+moderate ML enrichment is nominal
  (`q=0.1085`).
- OFC BA11 MEtan corrected ML-overlap enrichment is exploratory and
  BH-FDR-supported for the global high set; it remains separate from its
  eigengene Diagnosis result and functional annotation.
- Functional module annotations describe the genes grouped in a module; they
  do not validate an MDD mechanism.
- Tier4 ranked signatures were recovered from Git artifacts. Tier4 GSEA is
  excluded because its exact generator/settings were not recovered.
- Historical Step13 MDD module-level ML2 is excluded because its exact
  generator was not recovered.
- Independent external-cohort validation is incomplete. There are no validated
  biomarkers.

## Provenance milestones

| Milestone | Commit |
|---|---|
| Corrected ML prioritisation | `fa3b53f71a9373cd9856df127991abce89bf52df` |
| Corrected ML-dependent downstream analyses | `ca25ec35af359a4d409abf0f728833d2ff54761c` |
| Canonical provenance/control layer | `a483b4f6267befbb98ef49f23a27ffe6fa48ed69` |
| End-to-end publication audit evidence | `df46778` |
| Exact WGCNA reconstruction and annotation validation | `dd59b1dc2b64d681ec9e7b8f0092b41b005471ef` |
| WGCNA functional-enrichment interpretation audit | `0de77ca` |

Only artifacts identified by the canonical registry or a stage manifest should
be treated as current. Historical reports and outputs are retained for
provenance and are indexed in
[docs/CANONICAL_RESULTS_INDEX.md](docs/CANONICAL_RESULTS_INDEX.md) and
[docs/STALE_PUBLICATION_REFERENCE_INVENTORY.md](docs/STALE_PUBLICATION_REFERENCE_INVENTORY.md).

## Testing

The canonical validator and Python suite are intentionally lightweight:

```text
python audit/validate_canonical_pipeline.py
python -m unittest discover -s tests -p "test*.py" -v
```

Phase‑3B R validation uses the reconstructed environment documented in
[environment/README.md](environment/README.md).
