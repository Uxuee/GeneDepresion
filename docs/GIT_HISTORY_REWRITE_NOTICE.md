# Git history rewrite notice

On 2026-09-26, the unpublished local `main` ref was rewritten to remove two reproducible generated sensitivity-analysis tables that exceeded GitHub's 100 MB per-file limit:

- `outputs/length_gt_100_elisa_style/pvalue_correction_sensitivity/four_group_pvalue_correction_sensitivity_full.csv` — 442,234,871 bytes;
- `outputs/length_gt_100_elisa_style/pvalue_correction_sensitivity/frontal_subregion_pvalue_correction_sensitivity_full.csv` — 248,786,630 bytes.

The files were generated expansions, not source data. Their R generator, inputs-by-reference, parameters, compact summaries, gene subsets, plots, and report were retained. Exact blob identifiers, checksums, and row counts are recorded in `docs/LARGE_FILE_HISTORY_AUDIT.md` and `provenance/large_file_history_audit.csv`.

## Scope and method

After fetching origin, ref-specific reachability checks showed that the blobs existed only on the separate local `main` lineage. They were not reachable from the release branch, its upstream, or `origin/main`. The release branch therefore required no history rewrite.

The local `main` ref alone was rewritten using Git's path-specific `filter-branch --index-filter`, removing only the two exact paths. `git filter-repo` was not installed. No wildcard CSV removal was used, no remote ref was changed, and nothing was pushed.

Private local backup refs were created before rewriting:

- `refs/archive/pre_history_rewrite_2026_09_26` for the release HEAD;
- `refs/archive/pre_history_rewrite_main_2026_09_26` for the old local-main tip.

These recovery refs are local only and must not be pushed.

## Commit correspondence

The sensitivity-audit introduction changed from `c895ea4e27fafce09899ae58e90092e0831953f4` to `1ed14c54af600ddeeed3ee73a9581f707a2648ab`. The local-main tip changed from `dabae4f268388fe79a07a9a8dfd9591570a9e385` to `59f4fdb33ebe6ece42e8c37b316723d3ad565f73`.

All corrected publication milestones from `fa3b53f` through `bcdfa13` kept their original hashes because that lineage never contained either artifact. The complete requested identity map is in `provenance/git_history_rewrite_map.csv`.

No scientific result, canonical output, module assignment, prioritisation score, confidence category, or interpretation changed.
