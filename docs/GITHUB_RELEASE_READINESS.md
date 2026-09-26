# GitHub release-readiness assessment

Assessment date: 2026-09-26

Branch: `codex/clean-work-path-review-images`

Assessed commit: `9483973ffa11119663a77fe6d804f0e1061d3c70`

## Recommendation

**NOT_READY_TO_PUSH**

The canonical scientific state, provenance records, validation tools, and fresh-clone workflow are internally consistent. However, this branch contains two historical Git blobs larger than GitHub's 100 MB per-file limit. They were introduced after the current remote branch point, so a normal first push is expected to be rejected. Resolving this requires an explicit choice between removing those blobs from branch history or migrating them to Git LFS. No history rewrite was performed during this assessment.

## Canonical scientific state

- Corrected ML prioritisation: 671 gene-region pairs: 30 high, 119 moderate, and 522 low confidence.
- High + moderate set: 149 gene-region pairs.
- v1 evidence denominator: 125 events per region.
- v2 evidence denominator: 150 non-baseline events per region.
- Canonical corrected-ranking SHA-256: `6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1`.
- Phase 2 corrected ML-dependent downstream analyses are rooted at `outputs/length_gt_100_clean_pipeline/16_corrected_ML_downstream/`.
- Phase 3B independently reconstructed all six report WGCNA networks with adjusted Rand index 1.0; the established canonical assignments remain canonical.
- Phase 3C audits biological interpretation and annotation. It does not replace module assignments.

The historical `10_ML_prioritized_Tier2_extended` directory is retained only as provenance and is not a current input. Tier4 GSEA and Step 13 module-level ML2 remain explicitly excluded from publication because their required provenance was not recoverable. Broad pooled Frontal Cortex is not a main report analysis unit.

## Milestones included

| Commit | Purpose |
|---|---|
| `fa3b53f` | Freeze corrected ML prioritisation and publication provenance |
| `ca25ec3` | Freeze corrected ML-dependent downstream analyses |
| `a483b4f` | Establish canonical pipeline provenance and publication boundaries |
| `df46778` | Preserve end-to-end publication audit evidence |
| `dd59b1d` | Validate canonical WGCNA reconstruction and module annotation |
| `0de77ca` | Audit WGCNA functional enrichment interpretation |
| `01a20ae` | Prepare canonical repository navigation and data preflight |
| `ebefabb` | Harden fresh-clone validation on Windows |
| `95e4513` | Preserve canonical artifact line endings in fresh clones |
| `9483973` | Keep verification commands read only by default |

## Fresh-clone validation

A clone from committed Git objects was tested at `9483973` in a separate short Windows path with `core.longpaths=true`. It did not inherit untracked files from the working repository.

Results:

- Python test discovery: **54/54 passed**.
- Canonical validator in fresh-clone mode: **passed with zero errors and zero warnings**.
- The validator emitted five informational notices for canonical inputs deliberately kept local or external: normalized expression, diagnosis-model output, effect-size output, the v1 event file, and the v2 event file.
- Canonical status command: passed.
- Canonical dry run: passed and refused excluded, pending, validation-only, and audit-only stages.
- Tier4 Git-artifact verification: passed for seven retained artifacts.
- External-data preflight: exited with its documented actionable missing-data status and identified every missing destination/configuration item.
- Example path configuration loading: passed.
- Clone working tree after read-only checks: clean.

The first clone attempt exposed Windows path-length limitations in two deeply nested WGCNA figure names. The documented mitigation is a short clone root plus `core.longpaths=true`. A later clone check exposed line-ending-sensitive historical checksums; path-specific Git attributes now preserve the bytes expected by the canonical manifests.

`Rscript` is not installed system-wide and the local Phase-3B R library is intentionally not committed, so the 38 R reconstruction checks could not be re-executed inside the clean clone. The same checks were run in the configured source environment and passed **38/38**. The committed `renv.lock`, package-version record, session information, scripts, manifests, and checksums preserve the R reconstruction evidence.

## External data and environment setup

- Copy `config/project_paths.example.yaml` to a local, ignored path configuration and supply machine-specific locations.
- Run `python code/run_canonical_pipeline.py --check-external-data --config <local-config>` before attempting stages that require external data.
- Acquisition instructions and expected destinations are documented in `data/README.md` and `provenance/external_data_manifest.csv`.
- Python/R setup and Phase-3B environment restoration are documented in `environment/README.md` and `environment/phase3b/`.
- Raw count data, annotation data, AHBA donor files, and the Hawrylycz Table 6 resource are not redistributed in Git.

## Push-blocking Git objects

Both files were introduced by commit `c895ea4` and are absent from the current remote branch tip:

| Size | Blob | Historical path |
|---:|---|---|
| 442,234,871 bytes | `b0116e8471f0d549c8569a695e0e18c32d4a9167` | `outputs/length_gt_100_elisa_style/pvalue_correction_sensitivity/four_group_pvalue_correction_sensitivity_full.csv` |
| 248,786,630 bytes | `1a1302695ea66f91aaf7bc203dcf317609e795c7` | `outputs/length_gt_100_elisa_style/pvalue_correction_sensitivity/frontal_subregion_pvalue_correction_sensitivity_full.csv` |

No file at the assessed HEAD exceeds 100 MB, but GitHub evaluates all objects introduced by a push. Deleting the paths only from the current tree would therefore not fix the push.

## Required decision before pushing

Choose and approve one of these history-changing approaches:

1. Remove the two generated full-detail CSV blobs from branch history and retain compact summaries/manifests; or
2. Migrate the files to Git LFS and confirm that LFS is appropriate for the repository and remote.

After that operation, repeat the fresh-clone checks, compare canonical checksums, inspect the rewritten milestone mapping, and only then push.

## Remaining limitations

- A full end-to-end scientific rerun from raw external data was not performed; the audit validates the canonical saved state and reconstruction evidence.
- Some external inputs require user-supplied licensed/public downloads and local configuration.
- Historical reports and legacy outputs remain provenance only and may contain superseded 35-high/130-moderate ML results; canonical navigation warns against using them.
- Unrelated local report drafts, a manuscript PDF edit, and historical output deletions are not part of this release assessment and must remain excluded from release commits.
