# Research roadmap

## Completed locally

Conservative RRA history integration; preservation of weighted ML and source datasets; external-reference retrieval/harmonization and existing intersection work; descriptive-stage documentation; registry migration; integrity and path checks. Review `provenance/consolidation/verification.json` for measured verification results.

## Before publication

Review the local branch diff and scientific designation. Confirm external-validation release candidates individually for licensing, private/sample-level fields and size. Recheck GitHub remote heads and auth immediately before an explicitly authorized branch push. Do not merge main, force-push or create another repository automatically.

## Future scientific/technical work — separate authorization required

- Recover and independently verify the exact frozen LM generator and input lineage.
- Make environment/path configuration reproducible across machines without changing frozen results.
- Plan independent-cohort validation, including full-pipeline selection control.
- Resolve mapping and selected-list/FDR provenance limitations in reference comparisons.
- Keep WGCNA/AHBA context and method-specific FDR findings distinct from primary rankings.
- Recover missing Tier4 GSEA/module-ML provenance before considering publication of excluded branches.

No model rerun, new statistical testing or ranking regeneration is part of consolidation.
