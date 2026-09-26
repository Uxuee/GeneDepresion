# Publication reproducibility gaps

Baseline assessed: `ca25ec35af359a4d409abf0f728833d2ff54761c`. This is a dry-run assessment; no expensive pipeline stage was rerun.

## Fresh-clone result

A new researcher cannot currently clone the repository and regenerate the publication results. The **first failure** is Step 1: scripts require `D:/genes/genedata/GSE102556_raw_counts_GRCh38.p13_NCBI.tsv` and `Human.GRCh38.p13.annot.tsv`; these files are not tracked and no checksum-verified acquisition command is provided. Metadata alone may be fetched from GEO, but counts and annotation are still required.

## Irrecoverable or publication-blocking provenance gaps

| Gap | Consequence | Required resolution |
|---|---|---|
| Canonical Tier4 ranking/GSEA artifacts not found | cannot verify rankings, libraries, permutations, FDR, or claimed GSEA findings | recover exact artifacts or regenerate only this stage from frozen upstream effects with a versioned manifest |
| Step 13 module-level ML2 generator not found | cannot inspect or reproduce feature construction, nested operations, model fits, or permutations | recover exact script and dependencies; reproduce checksums, or remove Step 13 from submission |
| Active publication materials use legacy results | clean user can follow wrong canonical path | rebuild report package from corrected manifests and add automated stale-reference check |

## Missing dependencies and environment definition

`environment.yml` is unpinned and incomplete for the current pipeline. It does not fully specify later R/Python requirements such as WGCNA, dynamicTreeCut, lme4/lmerTest, qvalue, clusterProfiler, ReactomePA, msigdbr, enrichplot, pytest, and exact versions of annotation databases. No `renv.lock`, conda-lock file, container digest, or complete `sessionInfo()` accompanies every canonical stage.

The Phase-2 manifest records the temporary Python 3.12/pandas 3.0.1/numpy 2.3.5/scipy 1.18.1/matplotlib 3.11.2 environment, but that does not cover upstream R stages. The machine’s normal PowerShell runtime currently fails to load `hostfxr.dll` with `0x80070005`; tests were run with an ACL-safe temporary Python and `unittest`. R was unavailable in this audit environment.

Resolution: provide separately locked Python and R environments (or a container), include test tooling, record `sessionInfo()`/`pip freeze`, and verify a minimal fresh environment.

## Hard-coded/local paths

Examples include:

- `D:/genes/genedata` for raw RNA-seq inputs;
- `D:/Rlibs` for local R packages;
- absolute output paths embedded in manifests/reports;
- scripts that infer the project directory from a particular Windows layout.

Resolution: central configuration/environment variables plus relative paths in portable manifests. Preserve original absolute paths only as historical provenance fields.

## Missing raw-data acquisition records

| Resource | Present? | Gap |
|---|---|---|
| GSE102556 metadata | cached/committed derivative; GEO fallback in code | exact fetched series-matrix/version/date not fully fixed |
| GSE102556 raw counts | local only | acquisition URL/command, checksum, redistribution status |
| NCBI GRCh38.p13 annotation | local only | source URL, date, checksum, release metadata |
| AHBA six-donor microarray inputs | local processing outputs retained | raw acquisition/download manifest and terms not complete |
| Hawrylycz 2015 Supplementary Table 6 | local source and load report | original download URL/date/checksum and redistribution permission |
| MSigDB | accessed through msigdbr | release version and license/registered-access documentation |
| GO/KEGG/Reactome/org.Hs.eg.db | package-backed | exact releases/versions and download dates |

## Canonical-versus-historical ambiguity

- The historical numbered runner still maps Step 10 to `10_ML_prioritized_Tier2_extended`; its README warns that the directory is legacy, but scripts and multiple current-facing documents still point there.
- Raw v1/v2 event files were recoverable from Git objects and checksum-verified, but their ordinary working-tree provenance is fragile.
- Essential ML v1/v2 generators exist in a committed output `code_snapshot_2026-06-27_tier4_score_export`, not the normal `code/` path. A clean clone contains them only because the snapshot is tracked.
- Phase-2 outputs have a manifest and checksums, but exact old/new network edge preservation was not checked.

Resolution: a single machine-readable canonical index should name every current input/output and reject legacy paths; historical paths should be excluded from active runners.

## Database and identifier drift

Tier2 enrichment correctly uses same-analysis tested-gene backgrounds, but the mapped N/K/n values depend on annotation-package/database versions. WGCNA biological annotation ran without `org.Hs.eg.db` and mapped 1,796/4,267 rows through msigdbr only; 2,471 were unmapped. AHBA intersections operate on normalized symbols, which can change under alias database updates.

Resolution: record resource versions and mapping tables, checksums, loss counts, and explicit duplicate/alias rules alongside each result. Do not silently remap canonical historical results.

## Tests missing for publication-critical invariants

Existing 16 tests genuinely exercise scorer rules and canonical/downstream consistency, but coverage is absent for:

- raw count/length duplicate handling and a GeTMM fixture;
- LM/LMM formula, coding, convergence/failure counts, and coefficient extraction;
- complete Tier1–4 membership and GSEA ranking/output integrity;
- enrichment universe and correction-family invariants;
- WGCNA selected-parameter and module/eigengene/kME checksums;
- Step 13 generator behavior and grouped-fold leakage;
- AHBA post-enrichment sensitivity-set freshness;
- figure/table source-row versus plotted-value checks;
- repository-wide prohibition of active legacy paths/counts.

The audit script provides lightweight independent spot checks, but those checks are not wired into the normal CI suite.

## Easy documentation fixes

- Add data-acquisition commands, citations, licenses/terms, expected filenames, sizes, and SHA-256 values.
- Correct LM versus LMM terminology.
- State exact multiple-testing families and enrichment universes.
- Label pooled Frontal Cortex supplementary and distinct from CG25/dlPFC/OFC.
- State that AHBA is context, module ML is exploratory, and external validation is incomplete.
- Add the corrected 30/119/522 values and score formula everywhere current-facing.

## Reruns that are and are not required

- No rerun is needed for verified GeTMM, gene-wise effects, Tier1–3, corrected ML scoring, corrected published-AHBA enrichment, corrected WGCNA ML enrichment, or structural-core flags.
- Tier4/GSEA requires recovery or a narrowly scoped rerun only if its claims remain.
- AHBA post-enrichment-selected module-score sensitivity subsets require removal or a narrowly scoped corrected rerun; the all32 analysis need not be rerun for shortlist correction.
- Step 13 should not be rerun until its exact generator is recovered and audited.
- Publication figures/tables/reports that consume legacy outputs must be refreshed after the scientific content is frozen.
