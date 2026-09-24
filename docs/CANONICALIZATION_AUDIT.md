# Canonicalisation Audit

## Repository state at audit start

- Branch: `codex/clean-work-path-review-images`
- HEAD: `f4f60c2d12cae1e9ee600fc3ef0e6c74bfb450ef`
- Worktree: already dirty, with pre-existing modified, deleted, and untracked files. This cleanup preserved unrelated changes and did not restore/delete them.
- Corrected scorer: `code/extract_stable_ML_prioritized_Tier2_genes.py`
- Corrected result: `outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/`
- Legacy result: `outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended/`
- Relevant tests: `tests/test_ml_priority_scoring.py` plus the new canonical integrity test.

## Stale or contradictory material found

1. The legacy directory contains old 35/130/506 numerical outputs but its README and method file had been partly updated with corrected rules. Both now carry prominent legacy notices; historical content remains intact.
2. `docs/ML_PRIORITY_SCORE_CODE_AUDIT.md` correctly describes the pre-fix implementation, including denominator 180 and the old RF/baseline/missing-sign behaviour, but can be misread as current without its later changelog. Treat it as a historical audit.
3. `docs/ML_PRIORITY_SCORING_CORRECTNESS_FIX_CHANGELOG.md` contains both pre-regeneration observations and the later correction outcome. It is an audit trail, not the canonical method.
4. The AHBA report explicitly uses 35 high and 130 moderate legacy pairs.
5. Several integrated, supervisor, WGCNA, and workflow documents refer to `10_ML_prioritized_Tier2_extended` as the active input.
6. Several scripts still load the legacy final ranking directly. These are enumerated in `DOWNSTREAM_CORRECTION_IMPACT.md` and were not executed or changed.
7. The clean pipeline README/library retains historical “extended ML” wording and Step 10 paths. It must be updated when the downstream workflow is deliberately switched and regenerated.
8. The preliminary manuscript package and collected images inherit legacy-list-derived outputs and are not publication-current.

## Canonical structure established

- `PUBLICATION_PIPELINE.md`: actual pipeline and evidence-layer distinctions.
- `DATA_AND_RESULT_PROVENANCE.md`: files, hashes, seeds, legacy status, and Git state.
- `ML_PRIORITY_SCORE_FINAL.md`: sole authoritative corrected-score method.
- `CANONICAL_RESULTS_INDEX.md`: current versus historical results by stage.
- `DOWNSTREAM_CORRECTION_IMPACT.md`: controlled regeneration queue.
- `REPRODUCIBILITY_CHECKLIST.md`: current readiness status.
- `PUBLICATION_CHANGELOG.md`: publication cleanup history.
- `canonical_manifest.json`: machine-readable corrected-result identity.

## Deliberately not changed

- No raw processing, normalisation, QC, statistical model, tier, enrichment, ML training, score, threshold, AHBA analysis, WGCNA construction, module-level ML, or figure generation was rerun.
- Downstream input paths were not redirected yet because their saved outputs still reflect the legacy list and require a controlled regeneration.
- Historical audits and old numerical outputs were not deleted.
- No commit or push was performed.
