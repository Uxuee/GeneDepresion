> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Phase 3A provenance recovery and project stabilization

## Outcome

Phase 3A created a central, machine-readable canonical registry and placed
explicit publication boundaries around unresolved historical branches. No WGCNA,
ML training, enrichment, differential-expression model, or scientific result was
rerun or altered.

## Recovered provenance

- Tier4 ranked signatures: exact generator and seven exact Git blobs recovered; hashes and row counts frozen.
- Corrected ML ranking: canonical SHA-256 frozen as `6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1`; counts remain 671 total, 30 high, 119 moderate, 522 low.
- Raw ML events: v1 and v2 SHA-256 values remain frozen; denominators remain 125 and 150.
- External local inputs: observed hashes recorded without committing large files.
- AHBA module-score ML: all-32 retained as the active feature set; FDR005/FDR010/M26-only labeled historical post-selection sensitivity.

## Unresolved or excluded

- No exact Tier4 GSEA implementation/settings/results were recovered: excluded from publication.
- No exact MDD module-level ML2 generator was recovered: its historical outputs are excluded from publication.
- WGCNA biological annotation had only 42.1% mapping in the historical run and lacks a fully pinned executable R environment: pending controlled recovery, not rerun here.
- External validation is incomplete.

## Control layer created

- `provenance/canonical_pipeline.json`: authoritative stage registry.
- `code/run_canonical_pipeline.py`: status/validation/registry execution gate.
- `audit/validate_canonical_pipeline.py`: schema, Git-object, checksum, and exclusion checks.
- `code/project_paths.py` plus a versioned example config: no new hard-coded machine path is required.
- `code/provenance_manifest.py`: reusable hash/manifest helper.
- guardrail tests preventing the legacy ML ranking and missing-generator Step 13 from becoming canonical.

## Worktree safety

All unrelated dirty files present at the Phase 3A baseline were left untouched.
The deleted tracked outputs were not restored, and the preliminary manuscript
PDF/package was not staged or regenerated. No commit or push was performed.

## Tier4/GSEA

Exact Tier4 ranked artifacts and the exact ranking method were found in Git;
they were recovered by identity, not regenerated. Ranking status is canonical.
No GSEA generator, outputs, package choice, gene-set version, permutation type,
permutation count, seed, size constraints, or correction scope was recoverable.
GSEA is excluded from publication pending method resolution.

## Step 13 ML2

The module-level ML2 generator was not recovered, so there is no source
commit/blob/path or valid reproduction comparison to report. The retained
outputs remain historical/unverified and are excluded from publication.

## AHBA sensitivity sets

The all-32 published-module model remains canonical because it is independent of
the corrected shortlist. FDR005, FDR010, and M26-only are historical
post-selection sensitivities and are removed from the active publication path.
No new ML run was needed to establish that corrected shortlist enrichment has no
BH-FDR <0.10 module.

## WGCNA biological annotation

Historical mapping coverage was 1,796/4,267 rows (42.1%); 2,471 rows were
unmapped. `org.Hs.eg.db` was unavailable, and the historical script's saved
mapping table used the `msigdbr` mapping. R is not currently available in the
audited shell, so a controlled versioned reannotation was not performed.
Annotations remain limited support; modules, eigengenes, enrichments, and
structural-core results were not changed.

## Data provenance

Observed GSE102556 counts, GRCh38 annotation, and Hawrylycz Table 6 checksums are
now recorded. Large inputs still require manual or helper-assisted acquisition;
AHBA donor retrieval details, licenses, and several database snapshots remain
UNKNOWN / NOT RECORDED where evidence was absent.

## Environment

Recovered exact versions are documented in `environment/README.md`. The existing
`environment.yml` is incomplete and unpinned. No fake lockfile was created. R
was unavailable; Python verification used a temporary local runtime that is not
a versioned project dependency.

## Paths

New canonical code uses central config/environment resolution. Historical hard-
coded `D:/` paths remain intentionally untouched as provenance. A repository-
wide migration of historical scripts was outside scope. The local config is
gitignored.

## Canonical pipeline

The authoritative registry is `provenance/canonical_pipeline.json`; the entry
point is `code/run_canonical_pipeline.py`; validation is performed by
`audit/validate_canonical_pipeline.py`. Tier4 GSEA and stage 17 are excluded;
WGCNA annotation and publication assembly remain pending recovery.

## Tests

Commands executed:

```text
C:/Users/palom/AppData/Local/Temp/codex-python-ml-audit/python.exe -m unittest tests.test_ml_priority_scoring tests.test_canonical_ml_priority_result tests.test_corrected_downstream_consistency -v
C:/Users/palom/AppData/Local/Temp/codex-python-ml-audit/python.exe -m unittest tests.test_phase3a_guardrails_unittest -v
C:/Users/palom/AppData/Local/Temp/codex-python-ml-audit/python.exe audit/validate_canonical_pipeline.py
```

The existing 16 tests passed. The seven new guardrail tests passed after the
registry stage-id correction. The validator passed with no errors or warnings.

## Publication blockers remaining

- Tier4 GSEA lacks recoverable method and result provenance.
- MDD module-level ML2 lacks its exact generator.
- WGCNA biological annotation needs a pinned, versioned reannotation if it is to
  carry more than limited contextual weight.
- Current-facing reports/manuscript figures still require a later controlled
  refresh from canonical outputs.
- Independent external validation has not been completed.
