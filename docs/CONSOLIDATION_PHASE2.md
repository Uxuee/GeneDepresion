# Local consolidation record and publication procedure

Branch: `codex/mdd-consolidation-phase2`. Protection commit `d4bf578`; RRA integration commit `316ac54`. All original source branches/checkouts are retained. No files are moved out of their existing paths; integration adds missing tracked files and copies missing local-only stage products after hashing. New registry/documentation establishes current navigation.

Backup: `D:/genes/consolidation_backups/20261121_phase2/`. External-validation inventory: 227 files, 335,282,984 bytes, including 142 visible untracked files. Historical comparison: 85 historical files; 47 identical same-path files, two differences, 36 historical-only files. Both historical and current files remain intact. Full inventories and hashes are under `provenance/consolidation/`.

No old validation artifact is treated as current merely because its filename matches. External validation remains untracked; scripts/reports are review candidates, while original source data, standardized tables and intermediates remain local pending redistribution review. Current aggregate documentation is suitable for review without publishing those datasets.

## Exact proposed publication procedure â€” not executed

From `D:/genes/GeneDepresion`, after explicit publication approval:

```powershell
git remote get-url origin
gh auth status
git ls-remote --symref origin HEAD 'refs/heads/*'
git fetch --no-tags origin
git status --short
git log --oneline origin/main..codex/mdd-consolidation-phase2
git diff --stat origin/main...codex/mdd-consolidation-phase2
python audit/verify_consolidation.py
```

Verify origin is exactly `https://github.com/Uxuee/GeneDepresion.git`. If remote main has advanced, inspect its changes and integrate them into the consolidation branch only after checking compatibility; repeat verification. Review every staged/released file, particularly external-validation candidates. Stop on unexpected remote, auth failure, missing dependencies or checksum mismatch.

Only after the branch is approved for publication:

```powershell
git push --set-upstream origin codex/mdd-consolidation-phase2
gh pr create --repo Uxuee/GeneDepresion --base main --head codex/mdd-consolidation-phase2 --title 'Consolidate MDD workflow and designate exploratory Tier2-RRA' --body-file docs/CONSOLIDATION_PR_BODY.md
```

Do not force-push. PR creation and main merge are separate review steps; neither is performed by this phase. No operation targets `Uxuee/DGEAnalysis`.

## Local preservation details

Published-main ancestry was incorporated on the consolidation branch in `18038b8`; main itself was not checked out or merged into. Missing local stage artifacts: 717 files, 520,028,954 bytes. Missing frozen ML input tables: seven files, 163,586,183 bytes. These remain local and ignored. All 360 differing shared pre-consolidation file copies were line-ending-only; upstream originals were retained. The checksum-protected `docs/DOWNSTREAM_CORRECTED_RESULT_COMPARISON.md` remains byte-identical and its old canonical wording is scoped by the new workflow index.
