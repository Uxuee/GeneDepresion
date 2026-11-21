# Reproducibility and safe inspection

Consolidated root: `D:/genes/GeneDepresion`. The original `D:/genes/GeneDepresion_RRA` checkout is retained unchanged. The GitHub destination is `https://github.com/Uxuee/GeneDepresion.git`; `Uxuee/DGEAnalysis` is excluded.

The approved consolidation creates no new numerical ranking. Frozen output hashes are checked against the pre-integration inventory and original RRA checkout. The historical registry is preserved byte-for-byte in `provenance/consolidation/canonical_pipeline.before.json`. Local-only checkpoint/permutation/matched-input products are copied only when absent, checksum-verified and never automatically staged.

## Inspection commands

```text
python audit/verify_consolidation.py
python code/run_canonical_pipeline.py --status
python code/run_canonical_pipeline.py --dry-run
python audit/validate_canonical_pipeline.py --allow-missing-local-artifacts
```

The last command is the historical guardrail validator in non-writing mode; the consolidation verifier additionally checks primary RRA identity. A passing integrity check is not a raw-data rerun or scientific validation. Avoid `code/validate_final_evidence_synthesis.py` as an unrestricted command: its main routine regenerates derived synthesis; lightweight callers must disable determinism regeneration.

## Dependencies and availability

Use `config/project_paths.example.yaml`, `data/README.md`, `environment/README.md`, and `environment/ml_rra/README.md`. Runtime availability is reported separately from file availability. The current registry does not register execution commands. Frozen normalized matrices, matched metadata, filtered integer counts and original annotation must match the input manifests before any separately approved rerun. Step 2 still auto-selects its first matching metadata file. Older R/Python paths and `D:/genes/GeneDepresion_RRA` manifest locations are historical provenance, not instructions to overwrite results.

The frozen regional LM is reconciled in `STEP3_PROVENANCE_RECONCILIATION.md`; exact generator reconstruction is not certified. Output directories and historical filenames are retained to protect dependencies. Numerical stage IDs in the original registry remain compatibility identifiers; new stages use descriptive names.

## Publication scope

All 227 original external-validation files are preserved and backed up; none is bulk-staged. The inventory separates source data, scripts, results, reports and intermediates. Code/report files are candidates for public release after content review, not blanket redistributable declarations. Original supplements, mappings, large tables, sample-level information and raw/normalized matrices remain local until license/privacy review. Newly written aggregate documentation contains no sample-level data.

Public fresh clones do not include every local artifact. The output manifest marks availability and hashes; source hashes cannot supply missing datasets. Do not publish credentials, runtime packages, cached libraries or the local backup.

## Local verification environment

The existing `D:/genes/GeneDepresion_RRA/.runtime/python/python.exe` supplies numpy, pandas, scipy, scikit-learn and matplotlib for the lightweight checks. The default `D:/miniconda/python.exe` lacks several of these dependencies. No packages were installed. Seven missing frozen ML input tables (163,586,183 bytes) were copied from the retained RRA checkout and hash-verified; all eight inspected ML input paths now resolve. See `provenance/consolidation/local_input_copies.json`.

The two pre-consolidation checkouts had 360 shared-file byte differences, all attributable solely to line endings. Original upstream bytes were preserved rather than replaced with the other checkout version. Exact imported RRA bytes and checkout rules preserve frozen stage hashes.
