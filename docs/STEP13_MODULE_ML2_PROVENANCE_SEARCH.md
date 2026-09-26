# Step 13 module-level ML2 provenance search

No exact model generator was recovered from the working tree, ignored or
untracked files, reachable Git history and branches, stashes, or unreachable
Git objects. Filename, code-text, and schema searches were performed.

The closest candidates were rejected for specific reasons:

- `code/Step8_Define_WGCNA_Module_ML2_Stable_Core_Genes.py` consumes model output to define a stable core; it does not train the missing model.
- `code/ahba_modules/06_run_published_ahba_module_score_ml.py` is a separate AHBA published-module analysis.
- Existing Step 13 manifests/reports describe results but do not reconstruct the executable implementation.

| Candidate | Last reachable commit/date | Relevance | Matches retained model manifest? | Exact generator? |
|---|---|---|---|---|
| `code/Step8_Define_WGCNA_Module_ML2_Stable_Core_Genes.py` | `ca25ec35…`, 2026-09-25 | Defines stable-core genes from prior model output | downstream only | no |
| `code/ahba_modules/06_run_published_ahba_module_score_ml.py` | `257c5ecb…`, 2026-07-17 | Uses ML on published AHBA module scores | different feature construction/analysis | no |
| `module_level_ml2/stable_core_gene_manifest.json` | `770b4a0f…`, 2026-07-17 | Describes retained downstream artifact | descriptive, not executable | no |
| other Step13 reports/tables | reachable output history | supplies metrics/schema clues | insufficient to reconstruct implementation | no |
| stashes | none present | no candidate | no | no |
| unreachable commits/blobs from `git fsck` | all safely content-searched | no exact implementation match | no | no |

The retained Step 13 tables are therefore historical/unverified and excluded
from the canonical publication path. They may be restored only after the exact
generator, inputs, grid, seed, and environment are recovered or after a new,
clearly versioned analysis is specified.
