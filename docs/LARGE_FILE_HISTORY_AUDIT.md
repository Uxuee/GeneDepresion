# Large-file Git history audit

Audit date: 2026-09-26

## Decision

The two files are **GENERATED_LARGE_ARTIFACT_DO_NOT_TRACK** outputs. They are not source data. They are produced by `code/pvalue_correction_sensitivity_length_gt_100.R` from saved four-region and frontal-subregion standardized-effect/LMM tables. The same commit retains the generating script, parameters, compact summaries, candidate-gene subsets, overlap table, plots, and narrative report.

The initial release-readiness assessment found the objects in the local Git object database and incorrectly inferred that they were reachable from the release branch. A ref-specific audit after fetching origin established that they were reachable only from the separate, unpublished local `main` lineage. They were never reachable from:

- `codex/clean-work-path-review-images`;
- `origin/codex/clean-work-path-review-images`; or
- `origin/main`.

Therefore the validated release milestones did not require rewriting. Only the unpublished local `main` ref was rewritten.

## Objects

| Historical path | Blob | Size | Data rows | SHA-256 |
|---|---|---:|---:|---|
| `outputs/length_gt_100_elisa_style/pvalue_correction_sensitivity/four_group_pvalue_correction_sensitivity_full.csv` | `b0116e8471f0d549c8569a695e0e18c32d4a9167` | 442,234,871 bytes | 2,235,339 | `55da513de476830a150d782d3d867ecf12949d0dd1998454f683f582a7f112a5` |
| `outputs/length_gt_100_elisa_style/pvalue_correction_sensitivity/frontal_subregion_pvalue_correction_sensitivity_full.csv` | `1a1302695ea66f91aaf7bc203dcf317609e795c7` | 248,786,630 bytes | 1,132,166 | `33328a1fa5d0cb6820e1796a66ebd8c3b5b339f7b6da80a3e1252841b3fd2718` |

Both objects were introduced together in old commit `c895ea4e27fafce09899ae58e90092e0831953f4` (`Add p-value correction sensitivity audit`). No later old-main commit modified or deleted either path. The rewritten counterpart is `1ed14c54af600ddeeed3ee73a9581f707a2648ab`, with only those two generated paths omitted.

## Regeneration and retained evidence

The generating workflow:

- does not rerun the LMM;
- reads `standardized_effect_size_4groups/LMM_all_genes_with_standardized_effect_size_4groups.csv`;
- reads `frontal_subregion_sensitivity/frontal_subregion_LMM_all_genes_with_standardized_effect_size.csv`;
- applies the documented correction methods and scopes;
- writes the two full-detail expansion tables plus compact summaries.

Retained compact outputs include:

- `four_group_correction_sensitivity_summary.csv`;
- `frontal_subregion_correction_sensitivity_summary.csv`;
- `correction_vs_current_Tier2_overlap.csv`;
- candidate-gene subset CSVs;
- diagnostic plots; and
- `Pvalue_Correction_Sensitivity_Report.md`.

The original blob IDs, byte sizes, row counts, and SHA-256 checksums above identify the removed generated artifacts exactly. Private local backup refs preserve pre-rewrite recovery without placing the artifacts on the release branch.

## Reachability evidence

Before rewrite:

- current release HEAD: not reachable;
- tracked release upstream: not reachable;
- `origin/main`: not reachable;
- old local `main`: reachable;
- occurrence relative to remotes: local and unpublished only.

After rewrite:

- current release HEAD: no blob above 100 MB;
- rewritten local `main`: no blob above 100 MB;
- release upstream: no blob above 100 MB.

The largest blob reachable from the final release branch is 7,876,750 bytes. Its 1,149 unique reachable blobs total 157,470,283 bytes. The private backup refs intentionally keep the old objects recoverable in this local object database; they must not be pushed.
