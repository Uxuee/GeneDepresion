> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# GitHub release-readiness assessment

Assessment updated: 2026-09-26

Branch: `codex/clean-work-path-review-images`

Validated release-lineage base: `bcdfa135a1772caeb23db6ae2531243b0a76de27`

## Recommendation

**READY_TO_PUSH**

The canonical scientific state, provenance records, validation tools, and fresh-clone workflow are internally consistent. Ref-specific inspection established that the two oversized generated CSVs were never part of this release branch or any fetched remote ref. They existed only on the separate unpublished local `main` lineage, which was rewritten narrowly to omit them. Neither the release lineage nor any scientific milestone required rewriting.

No push has been performed. Private backup refs containing the old local-only history must remain local.

## Canonical scientific state

- Corrected ML prioritisation: 671 gene-region pairs: 30 high, 119 moderate, and 522 low confidence.
- High + moderate set: 149 gene-region pairs.
- v1 evidence denominator: 125 events per region.
- v2 evidence denominator: 150 non-baseline events per region.
- Canonical corrected-ranking SHA-256: `6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1`.
- Phase 2 corrected ML-dependent downstream analyses are rooted at `outputs/length_gt_100_clean_pipeline/16_corrected_ML_downstream/`.
- Phase 3B independently reconstructed all six report WGCNA networks with adjusted Rand index 1.0; established canonical assignments remain canonical.
- Phase 3C audits biological interpretation and annotation without replacing module assignments.

The historical `10_ML_prioritized_Tier2_extended` directory is provenance only. Tier4 GSEA and Step 13 module-level ML2 remain explicitly excluded from publication because required provenance was not recoverable. Broad pooled Frontal Cortex is not a main report analysis unit.

## Large generated artifacts

The local Git database contained two generated full-detail sensitivity tables of 442,234,871 and 248,786,630 bytes. Their exact audit is in `docs/LARGE_FILE_HISTORY_AUDIT.md`.

- Release HEAD before rewrite: neither blob reachable.
- Release upstream: neither blob reachable.
- `origin/main`: neither blob reachable.
- Old local `main`: both blobs reachable.
- Rewritten local `main`: neither blob reachable.

Only local `main` was rewritten. Important release milestone hashes are unchanged; the correspondence is recorded in `provenance/git_history_rewrite_map.csv`.

The release branch has no blob above 100 MB. Its largest reachable blob is 7,876,750 bytes, and its 1,149 unique reachable blobs total 157,470,283 bytes. The private local object store remains larger because the required recovery refs deliberately preserve the old history.

## Validation

- Python test discovery: **54/54 passed**.
- Phase-3B R reconstruction checks: **38/38 passed** in the configured source environment.
- Canonical validator: passed with zero errors and zero warnings.
- Tier4 Git-artifact verification: passed for seven retained artifacts.
- Canonical ranking checksum: verified.
- Phase-2, Phase-3B, and Phase-3C manifests: verified by the test and validation suites.
- Canonical counts and v1/v2 denominators: verified.
- `git diff --check`: passed.

## Fresh-clone validation

A fresh short-path clone from committed release history was tested without untracked workspace files:

- all 54 Python tests passed;
- canonical validator passed with zero errors and warnings;
- five expected informational notices identified canonical inputs deliberately kept local or external;
- canonical status and dry-run checks passed;
- Tier4 retained-artifact verification passed;
- external-data preflight produced the documented actionable missing-input status;
- no reachable blob exceeded 100 MB; and
- the clone remained clean after read-only verification.

`Rscript` is not installed system-wide and the local Phase-3B R library is intentionally not committed, so the 38 R checks were run in the configured source environment rather than the clean clone. The committed lockfile, versions, session information, scripts, manifests, and checksums preserve that reconstruction evidence.

## External data and environment setup

- Copy `config/project_paths.example.yaml` to a local ignored configuration and supply machine-specific locations.
- Run `python code/run_canonical_pipeline.py --check-external-data --config <local-config>` before stages requiring external data.
- Acquisition instructions and destinations are documented in `data/README.md` and `provenance/external_data_manifest.csv`.
- Environment restoration is documented in `environment/README.md` and `environment/phase3b/`.
- Raw counts, annotation data, AHBA donor files, and Hawrylycz Table 6 are not redistributed in Git.

## Remaining limitations

- This audit validates the canonical saved state and reconstruction evidence; it is not a full rerun from external raw data.
- Some external inputs require user-supplied public/licensed downloads and local configuration.
- Historical reports may contain superseded 35-high/130-moderate ML results and are clearly separated from canonical navigation.
- Private backup refs and `refs/original/refs/heads/main` contain pre-rewrite local-only history and must not be pushed.
- The preserved pre-release stash contains unrelated drafts, manuscript edits, and historical deletions and must not be applied to the release branch.

The branch is ready for an explicitly targeted push after user approval.
