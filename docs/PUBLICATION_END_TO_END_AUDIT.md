# Full end-to-end publication audit

Audit baseline: `ca25ec35af359a4d409abf0f728833d2ff54761c`
Upstream corrected-ML milestone: `fa3b53f71a9373cd9856df127991abce89bf52df`
Audit date: 2026-09-25
Scope: adversarial read-only audit; no scientific pipeline output was changed.

## Executive conclusion

The corrected ML ranking and the corrected shortlist-dependent Phase-2 tables are numerically reproducible from retained event records and pass all 16 relevant regression tests. Core upstream calculations checked here—GeTMM, two diagnosis coefficients, BH correction, standardized effects, Tier membership, one Tier2 enrichment, the final ML score, one AHBA enrichment, and WGCNA enrichment—agree with canonical values to floating-point precision.

The project is **not yet submission-ready**. The principal blockers are publication/provenance rather than evidence of a newly discovered arithmetic error:

1. active report and README files still present the obsolete 35-high/130-moderate/165-total result, old score weights, legacy paths, and pre-correction AHBA/WGCNA claims;
2. no retained canonical Tier4 ranking/GSEA result files were found, so Tier4/GSEA numerical claims cannot be independently verified or regenerated from the documented current workflow;
3. the generator for Step 13 module-level ML2 was not found, preventing implementation-level verification or clean-clone regeneration;
4. AHBA module-score ML post-enrichment sensitivity subsets were selected from legacy enrichment results and are stale after correction, although the main all-32-module analysis is independent of the shortlist;
5. raw GSE102556 counts/annotation and AHBA microarray downloads are local-only, paths are hard-coded, and the environment is incomplete/unpinned.

No gene-wise diagnosis model, WGCNA module-eigengene diagnosis test, or corrected published-AHBA ML-shortlist test supports a validated biomarker claim. All ML and module findings remain internal, exploratory prioritization.

## 1. Git and baseline verification

- Branch: `codex/clean-work-path-review-images`.
- HEAD: `ca25ec35af359a4d409abf0f728833d2ff54761c` (exactly the requested baseline).
- `ca25ec35...` is a descendant of `fa3b53f...`.
- Canonical ranking SHA-256: `6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1` (verified).
- The worktree was already dirty. Modified documents/scripts, deleted historical images/tables, an edited preliminary PDF, and untracked supervisor documents/caches were not used as canonical numerical evidence. Where a deleted tracked input was needed, the audit read it directly from the baseline commit with `git show`.
- Audit additions are limited to the four requested reports, `audit/publication_audit_verification.py`, and its JSON evidence file. Nothing was staged, committed, restored, or pushed.

## 2. Implemented dependency DAG

| Stage | Generator / canonical source | Canonical output | Determinism / seed | Corrected-shortlist dependency | Status |
|---|---|---|---|---|---|
| Raw inputs | GSE102556 counts/metadata; NCBI GRCh38.p13 annotation | local `D:/genes/genedata` files | deterministic input | no | local-only |
| Filtering + GeTMM | `code/Step1_length_cutoff_sensitivity_GeTMM.R` | Step 1 normalized matrices copied into clean pipeline | deterministic | no | canonical calculation; external inputs absent from Git |
| QC/covariates | Step 2 scripts and saved audit tables | QC/PCA/covariate outputs | deterministic | no | retained |
| Gene-wise models | `code/Step3_run_all_DGE_methods_on_length_cutoffs.R`; frontal sensitivity script | gene-wise effect tables | deterministic | no | retained |
| Standardized effects | `code/make_standardized_effect_size_4groups_length_gt_100.py` | standardized effect tables | deterministic | no | retained |
| Tiers | `code/make_evidence_tiers_4groups_length_gt_100.py`, packaging helpers | Tier1–4 fields; Tier2 sets | deterministic | no | Tier1–3 retained; Tier4 rank file absent |
| Tier2 enrichment | committed snapshot `run_Tier2_enrichment_length_gt_100.R`; clean Step 8 copies results | `08_Tier2_enrichment` | deterministic conditional on databases | no | retained; database versions incompletely pinned |
| Tier4 GSEA | historical scripts/references only | expected ranked-GSEA files | seed/version not recoverable | no | canonical outputs not found |
| ML v1 | committed snapshot `run_ML_pilot_Tier2_length_gt_100.py` | v1 raw fold importance | seeds 2026–2030 | Tier2 only in canonical scorer | raw events restored and checksum-verified |
| ML v2 | committed snapshot `run_ML_strict_v2_length_gt_100.py` | v2 selected-gene events | outer seed 31001; model seeds 11/13 | starts from expression; joined to Tier2 for score | raw events restored and checksum-verified |
| Corrected ML score | `extract_stable_ML_prioritized_Tier2_genes.py` via wrapper | `10_ML_prioritized_Tier2_corrected` | deterministic | defines shortlist | canonical 671/30/119/522 |
| Published AHBA enrichment | fixed Hawrylycz module membership; Phase-2 regenerator | corrected AHBA tables | deterministic | yes, annotations only | corrected canonical |
| WGCNA construction | Step 6 broad/frontal R scripts | assignments, eigengenes, kME | deterministic for given data/settings | no | preserved |
| WGCNA ML enrichment | Phase-2 regenerator | corrected broad/frontal enrichment | deterministic | yes | corrected canonical |
| Structural core | Step 8 stable-core script | structural core | deterministic | no for membership; yes for flags | membership preserved, flags refreshed |
| WGCNA network figures | network renderers | Phase-2 figures/tables | deterministic reconstruction | ML labels only | edges reconstructed, not proven byte-preserved |
| AHBA module-score ML | `code/ahba_modules/06_run_published_ahba_module_score_ml.py` | Step 11 ML | seeded CV | all32: no; sensitivity subsets: yes via old enrichment | main preserved; sensitivity subsets stale |
| Module-level ML2 | generator not found; retained outputs/manifest only | Step 13 | manifest says grouped 15-fold CV | direct inputs fixed WGCNA features | not clean-clone reproducible |
| Reports/figures | multiple copy/report scripts | docs and manuscript package | mixed | mixed | several active files stale |

Hidden dependencies include absolute `D:/` locations, the local `D:/Rlibs`, database contents supplied by installed R packages, Git-object restoration of raw ML events, and output-directory code snapshots that are essential to provenance.

## 3. Sample and donor accounting

The saved metadata contain no duplicate sample IDs and no missing Age, Sex, RIN, PMI, or pH in the seven audited units.

| Analysis unit | Samples | MDD | Control | Donors | Repeated donors within unit | Max samples/donor | Consequence |
|---|---:|---:|---:|---:|---:|---:|---|
| Ventral Striatum | 48 | 26 | 22 | 48 | 0 | 1 | sample CV equals donor CV |
| Insula | 48 | 26 | 22 | 48 | 0 | 1 | sample CV equals donor CV |
| Hippocampal Formation | 43 | 24 | 19 | 43 | 0 | 1 | sample CV equals donor CV |
| Cingulate gyrus 25 / CG25 | 27 | 12 | 15 | 27 | 0 | 1 | small n; sample CV equals donor CV |
| dlPFC BA8/9 | 48 | 26 | 22 | 48 | 0 | 1 | sample CV equals donor CV |
| OFC BA11 | 47 | 25 | 22 | 47 | 0 | 1 | sample CV equals donor CV |
| pooled Frontal Cortex | 122 | 63 | 59 | 48 | 47 | 3 | repeated-donor model/grouped CV required |

The canonical v1/v2 gene-level ML analyzes the six original units, not pooled Frontal Cortex; therefore no within-unit donor crosses train/test folds. The pooled frontal gene-wise analysis uses `(1|BrainID)`, and Step 13 records grouped CV. Cross-tissue use of the same donors remains relevant when results from different regions are compared as if independent studies: they are not independent cohorts.

## 4. Normalization and filtering

The implemented path is: collapse duplicate count GeneIDs by summation; match NCBI GRCh38.p13 annotation; require finite length greater than the canonical 100-bp cutoff; within each region retain genes whose 80th-percentile raw count is greater than 10; calculate RPK as count/(length kb); calculate edgeR TMM factors on RPK; divide RPK by `colSums(RPK) * TMM_factor` and multiply by one million; store `log2(GeTMM + 1)`.

This matches the documented implemented GeTMM variant. Filtering occurs before normalization and region-wise multiple testing, so the tested universe changes with the filter. A retained value was reconstructed with an absolute error of `1.11e-16` (full arithmetic in the numerical appendix).

Risk: the gene-length source is a local NCBI annotation file without a retained download checksum/version manifest beyond its filename; annotation duplicate handling is not documented as explicitly as count-row duplicate handling.

## 5. QC and covariates

| Unit(s) | Gene-wise model covariates | WGCNA eigengene model | Availability / concern |
|---|---|---|---|
| VS, Insula, HF, CG25, dlPFC, OFC | Diagnosis + Age + Sex + RIN | Diagnosis + Age + Sex + RIN + PMI | complete in retained samples |
| pooled Frontal | above + OriginalRegion + `(1|BrainID)` | Diagnosis + Age + Sex + RIN + PMI | repeated donors handled in gene model |

The saved WGCNA covariate audit found low VIF for the current model (maximum 1.13–1.19), residual df 21–116, strong RIN/PC relevance, broad Sex/Age relevance, and weaker but defensible PMI relevance. RIN differed nominally by diagnosis in VS and HF; pooled-frontal pH differed by diagnosis. pH was not included and alcohol was 42–48% missing/inconsistently coded. No evidence was found that the primary covariate formula was selected to optimize diagnosis p-values, but the rationale needs prospective wording rather than outcome-based justification.

## 6. Gene-wise models and multiple testing

- True mixed-effects analysis: pooled Frontal Cortex, `expression ~ Diagnosis + Age + Sex + RIN + OriginalRegion + (1|BrainID)`, ML fit.
- Fixed-effect linear models: VS, Insula, HF, CG25, dlPFC, OFC, `expression ~ Diagnosis + Age + Sex + RIN`; within each unit there is one sample per donor.
- Diagnosis reference is Control; positive beta means higher expression in MDD and negative beta means lower expression in MDD.
- Calling every saved table “LMM” is a documentation-method mismatch: six original units are ordinary gene-wise linear models.
- Two Insula coefficients, SEs, statistics, and p-values were reconstructed to approximately `1e-15` or better.
- BH is applied within region/subregion across tested genes after filtering. Insula’s 23,090-value BH vector reproduced with maximum absolute error `1.11e-16`; minimum BH was 0.933686. Storey q-values are used for Tier1. No canonical gene-wise or WGCNA module diagnosis claim should be called FDR-significant unless its own stated family passes the threshold.

## 7. Effects and tier definitions

- Standardized effect = diagnosis beta / SD of that gene’s normalized expression among samples in the same analysis unit.
- Tier1: Storey q < 0.10 and absolute standardized effect > 0.8.
- Tier2: raw diagnosis p < 0.01 and absolute standardized effect > 0.8.
- Tier3: raw diagnosis p < 0.01 and absolute beta > 0.2.
- Tier4 score: `sign(standardized effect) * -log10(p)`; absolute magnitude determines extremeness where top-k is used.

All 95,319 broad-region rows recomputed with zero Tier1/2/3 membership mismatches. Broad counts were Frontal 0/1/99, HF 0/127/70, Insula 0/70/92, and VS 4/146/167. Frontal subregion Tier2 contributes 186 CG25, 72 dlPFC, and 70 OFC rows, yielding the canonical 671 gene-region pairs. The Tier4 formula was numerically checked for one row, but no canonical Tier4 ranking/GSEA files were found.

## 8. Tier2 enrichment and Tier4 GSEA

Tier2 enrichment maps Entrez IDs with `org.Hs.eg.db`, uses the tested genes from the same region/subregion as `universe`, splits all/up/down sets, runs clusterProfiler GO/KEGG and optional ReactomePA/MSigDB collections, and applies BH within each saved database/gene-set result family. The checked GO:MF result used N=15,343, K=12, n=30, k=2; raw p and BH both reproduced.

The universe is therefore not an inappropriate all-human universe. Database-specific mapped universes differ (appropriately) because only genes represented by each resource enter a test. Package/database versions and download dates were not captured sufficiently for exact future regeneration.

Tier4/GSEA is a publication blocker if its claims remain in the manuscript: the score formula is recoverable and independent of the ML correction, but retained ranked inputs, GSEA results, software/library versions, permutation settings, and complete correction provenance were not found. “Unchanged because ML did not change Tier4” is logically plausible, not a numerical verification of a missing result.

## 9. ML v1, v2, and final priority score

### v1

- Canonical scoring filters raw events to `Tier2_fixed`; historical Tier4 top-k v1 events do not contribute.
- Five models × five folds × five repeat seeds = 125 events/region.
- Models: L1, elastic-net, L2 logistic regression, linear SVM, and exploratory RF.
- Imputation and scaling are inside training-fold pipelines where applicable.
- L1/elastic selection uses nonzero absolute coefficients; L2/SVM/RF use event-wise top 20% importance.
- The Tier2 feature pool itself was defined on the full MDD dataset, so v1 predictive performance is optimistic and not a clean out-of-sample validation. It remains suitable only as within-study prioritization.

### strict v2

- Two selectors (ANOVA F and absolute standardized mean difference) × five top-k values (25/50/100/250/500) × five feature models × three outer folds = 150 gene-evidence events/region.
- Selection is fit on outer-training data; inner tuning, imputation, and scaling are training-only.
- Majority baseline is excluded from gene evidence and denominator.
- With one sample/donor in each analyzed unit, no donor crosses outer folds.

### final score

The score reconstructed all 671 rows with maximum absolute difference `9.71e-17`. Counts are exactly 30 high, 119 moderate, and 522 exploratory/low. RF is nondirectional; missing sign evidence cannot pass the high gate; the 0.5 sign imputation is used only in the weighted score. ROC-AUC, accuracy, and balanced accuracy do not enter the gene score.

Selection frequency, model support, sparse support, rank, sign consistency, and v2 support are correlated summaries of shared CV events. They are not six independent replications and must not be worded that way.

## 10. AHBA audit

Published modules are fixed M1–M32 assignments from Hawrylycz et al. 2015 Supplementary Table 6 (15,627 unique gene symbols). They are distinct from the project-generated AHBA WGCNA attempt. The ML shortlist does not define the published modules; it only defines the candidate set tested against them.

The corrected Ventral Striatum/up-in-MDD/M26 result reproduces N=12,042, K=554, n=6, k=2 (`MORN5;RPGR`), p=0.0280130, BH=0.882283. No corrected published-AHBA ML-shortlist row has BH <0.10. The lost overlap is TSGA10. This is nominal context, not external validation.

The main AHBA module-score ML using all 32 externally defined modules does not depend on the 30/119 shortlist and may remain as exploratory within-dataset analysis. Its `post_enrichment_FDR005`, `FDR010`, and M26-only sensitivity feature sets were selected from legacy enrichment and are stale after correction; corrected AHBA enrichment has no FDR<0.10 module subset.

## 11. WGCNA construction, integrity, and enrichment

Construction used normalized log2 GeTMM expression and unsupervised MAD ranking. Diagnosis was not used to construct modules. Broad units used the largest available tested filter (typically 15,000 MAD genes); original frontal subregions used 5,000. Networks are signed; bicor with maxPOutliers 0.05 is attempted first, with Pearson fallback; powers 1–20 seek the first signed R2 >=0.80, otherwise the best R2. Dynamic tree cut grid: deepSplit 2/3/4, minimum module size 20/30/50, merge cut 0.15/0.25. A solution must have >=3 non-grey modules, largest fraction <0.70, and grey fraction <0.50. Ranking then favors usable status, more genes, more non-grey modules, and smaller largest-module fraction. This adaptive choice did not use Diagnosis, but module counts are not directly comparable across the 15k/5k universes. CG25 did not reach target R2 within powers 1–20.

No module eigengene Diagnosis association survives within-region BH correction. The current eigengene formula is `ME ~ Diagnosis + Age + Sex + RIN + PMI`; all emphasized raw p-values are nominal.

Phase 2 records unchanged SHA-256 values for broad and frontal assignment inputs and exact stable-core shared columns (3,200 broad rows and 20,154 frontal rows). It did not rewrite eigengenes or connectivity. However, the network renderer recalculates sparse correlations/edges from unchanged expression and membership; the statement “edges preserved” is stronger than demonstrated. Use “deterministically reconstructed from unchanged inputs,” unless old/new edge lists are compared directly.

The corrected VS MEdarkred/global high+moderate test reproduces overlap 6, p=0.00217014, q=0.108507, after loss of PRELID2 and TMEM269-DT; it is nominal only. OFC MEtan retains FDR<0.10 for corrected global ML sets (high: overlap 2, p=0.00178064, q=0.0320514). Frontal corrected network inventory contains CG25, OFC, and dlPFC only—no pooled Frontal Cortex.

## 12. Structural core, network, and module-level ML

Structural-core membership uses module membership, expression availability, kME/hub rank—not ML confidence. Phase 2 refreshed ML flags on unchanged rows. “Structural core” and “ML-prioritized structural-core genes” must remain distinct concepts.

Step 13 retained outputs say features are fixed WGCNA module eigengenes/scores, scaling and PCA are trained within folds, and CV is `StratifiedGroupKFold` with 15 splits. The broad branches use a small set of pre-prioritized modules (VS four; Insula/HF one each), while frontal branches use many/all modules. Because the generator script is absent, these assertions cannot be checked from implementation or reproduced. Furthermore, modules/features and outcomes come from the same cohort, and some module choices were historically informed by diagnosis/Tier2/ML evidence. This is internally selected exploratory analysis, not independent validation. The saved lightweight permutation p-values operate on pooled CV predictions rather than repeating the full selection/training pipeline and should not be treated as confirmatory.

## 13. Cross-analysis independence

| Analysis | Same MDD data used for selection? | Same samples used for test? | Correct interpretation |
|---|---|---|---|
| Tier2 → pathway enrichment | yes for genes | annotation database external | functional annotation of selected genes |
| Tier2 → v1 ML | yes | yes | internal prioritization; performance optimistic |
| strict v2 selection | selection train-only within CV | held-out fold for metrics | internal validation only |
| ML shortlist → WGCNA enrichment | yes | WGCNA modules from same samples | convergence/context, not replication |
| ML shortlist → published AHBA overlap | yes | module membership externally fixed | external biological context, not MDD validation |
| WGCNA modules → ML2 | modules/some choices from same cohort | yes | internally selected exploratory prediction |
| all-32 AHBA scores → ML | modules externally fixed | MDD samples used for model/CV | internal validation with external feature definition |

## 14. Gene identifiers, universes, regions, and directions

- The corrected 671 rows are 671 unique region–Entrez pairs, representing 640 unique Entrez IDs and 640 symbols; no missing symbols or symbol-to-multiple-ID collisions occur in that table.
- AHBA published-module mapping reduces those 640 symbols to 276 present in the published-module background for the all-Tier2 test; mapping loss is explicit in the candidate audit.
- Project AHBA preprocessing reports 332/640 candidate symbols present in its 29,131-gene collapsed matrix.
- WGCNA biological annotation mapped only 1,796/4,267 prioritized-module rows (42.1%) through available msigdbr metadata; 2,471 were unmapped because `org.Hs.eg.db` was unavailable. This substantially limits GO/hub biological-theme coverage and must be disclosed. Region rates ranged from 25.2% for pooled frontal to 70.4% for VS.
- Tier2 universes are same-unit tested and database-mapped genes; published-AHBA uses the 12,042-symbol intersection background; WGCNA enrichment uses each network’s gene universe (15,000 broad, 5,000 frontal).
- Pooled Frontal Cortex is not equivalent to CG25, dlPFC, or OFC. Main-report figures should use the three original frontal subregions and must label pooled analyses as historical/supplementary.
- Diagnosis coding is Control reference. Positive beta/linear coefficient means MDD-up; negative means MDD-down. RF importance has no direction.

## 15. Figure/table and claim provenance

The Phase-2 manifest checksums match the corrected AHBA and WGCNA tables and corrected network figures. The preliminary manuscript package was intentionally not refreshed. Multiple active files still point to `10_ML_prioritized_Tier2_extended`, state 35/130/165, use old 0.30/0.20/0.15 score weights, or quote legacy M26/MEdarkred FDR results. Examples include the root README, main workflow/report documents, AHBA report, WGCNA module report, supervisor/Notion report, and some figure-copy inputs. These are active stale dependencies, not harmless historical mentions. Historical change logs and explicit “legacy” comparison documents are acceptable.

The main scientific wording must state: no validated biomarkers; no module Diagnosis BH-FDR result; corrected M26 and MEdarkred are nominal; OFC MEtan ML enrichment is an exploratory FDR-supported over-representation result; AHBA is context; module-level ML is exploratory; external replication is incomplete.

## 16. External validation status

No independent cohort validation of the gene ranking, MDD diagnosis coefficients, 10-module shortlist, or classifiers is complete. Cross-validation is internal validation. AHBA module overlap is annotation against an external brain atlas, not validation of an MDD signature. WGCNA enrichment and published-module mapping are not independent disease replication. A publication should request independent transcriptomic cohort replication with harmonized tissues and prospectively frozen thresholds/models; validation of direction, ranking, module preservation, and predictive calibration should be reported separately.

## 17. Environment, data, and fresh-clone assessment

`environment.yml` specifies Python 3.11 and an R/Bioconductor stack, but does not pin versions and omits several later dependencies (for example WGCNA/dynamicTreeCut, clusterProfiler, ReactomePA, msigdbr, lme4/lmerTest and test tooling). The present shell could not load its default PowerShell runtime due ACL; an ACL-safe temporary Python ran the checks. `pytest` is absent, although `unittest` runs the suite. R execution was not available for this audit.

The first clean-clone failure is Step 1: raw count and annotation files are expected at `D:/genes/genedata` and are neither tracked nor automatically downloaded. Subsequent blockers include absolute Windows paths, missing AHBA raw donor data/download automation, unpinned R/MSigDB/GO/KEGG/Reactome resources, event files whose provenance depends on Git blobs/output snapshots, missing Tier4/GSEA outputs, and missing Step 13 generator.

External resources requiring a Data Availability/licensing statement include GSE102556/GEO, the NCBI GRCh38.p13 annotation, Allen Human Brain Atlas data and terms, Hawrylycz et al. 2015 module membership, MSigDB (licensed/registered access and citation), GO, KEGG, Reactome, and R/Bioconductor annotations. Exact download dates/checksums and redistribution permissions are not comprehensively recorded.

## 18. Stale-reference search and tests

The repository-wide search distinguished valid historical comparison references from current-facing violations. Current-facing stale references occur in the root README and several report/manuscript documents; the legacy directory is correctly marked provenance-only in canonical manifests but remains wired into the historical numbered runner’s Step 10. The runner README warns about this, yet a fresh user could still invoke the wrong target.

Tests executed with:

```text
C:/Users/palom/AppData/Local/Temp/codex-python-ml-audit/python.exe -m unittest -v \
  tests.test_ml_priority_scoring \
  tests.test_canonical_ml_priority_result \
  tests.test_corrected_downstream_consistency
```

Result: 16 run, 16 passed, 0 failed, 0 skipped. The tests cover scorer rules, counts, denominators, manifest integrity, representative arithmetic, corrected flags/backgrounds, assignment membership, legacy-table comparison, and frontal exclusion. They do not validate normalization, gene models, Tier2/GSEA generation, WGCNA construction, database versions, module-level ML2 implementation, all figure values, or complete clean-clone execution; the audit script supplies spot checks but is not a replacement for pipeline integration tests.

## 19. Correction history and stability

- Legacy: 35 high, 130 moderate, 506 low.
- Corrected: 30 high, 119 moderate, 522 low.
- Transitions: five high→moderate, sixteen moderate→low, 650 unchanged; no upward transitions.
- Corrected-score denominator and weights changed as intended; RF sign and missing-sign gates caused no independent category changes in this dataset.
- M26 weakened from legacy FDR<0.10 to nominal; MEdarkred weakened from FDR<0.10 to nominal; OFC MEtan persisted at FDR<0.10.
- Tier2 membership/enrichment, the underlying Tier4 formula, WGCNA construction, and structural membership are upstream and unchanged. Missing Tier4/GSEA artifacts prevent numerical confirmation of their outputs.
- Conclusions most sensitive to one/few genes are M26 and MEdarkred. No correction creates a validated diagnosis biomarker.

## 20. Publication risk register

| ID | Severity | Category | Stage | Finding | Scientific impact | Code change? | Rerun? | Documentation? | Recommended action |
|---|---|---|---|---|---|---|---|---|---|
| PUB-01 | CRITICAL | G | Reports/manuscript | Current-facing files retain obsolete 35/130/165, old weights/paths and legacy M26/MEdarkred claims | submission would report superseded results | no | regenerate figures/tables only | yes | refresh manuscript/report package from corrected canonical outputs and re-audit claims |
| PUB-02 | MAJOR | F | Tier4/GSEA | Canonical ranking/GSEA outputs and complete run provenance not found | Tier4/GSEA claims cannot be verified or reproduced | maybe workflow packaging | rerun only if claims retained | yes | recover artifacts or regenerate from frozen upstream ranks with documented versions/settings |
| PUB-03 | MAJOR | F | Module-level ML2 | Generator script absent | implementation, leakage controls and reproduction cannot be independently checked | restore/add exact historical generator | possibly, after code recovery | yes | recover generator and pin manifest/checksums before publication |
| PUB-04 | MAJOR | E | AHBA module-score ML | post-enrichment sensitivity subsets reflect legacy FDR-selected modules | sensitivity labels/features no longer correspond to corrected enrichment | no for all32; yes for sensitivity | sensitivity only | yes | exclude stale subsets or regenerate them under a prospectively stated corrected rule |
| PUB-05 | MAJOR | F | End-to-end reproduction | raw counts/annotation/AHBA inputs local-only; hard-coded paths; incomplete unpinned environment | clean clone cannot regenerate results | yes, packaging/config only | no scientific rerun required to document | yes | add acquisition/checksum manifest, parameterized paths, complete lockfiles |
| PUB-06 | MODERATE | C | ML v1 | Tier2 feature set is defined on full dataset | CV performance is optimistic | no if treated as prioritization | no | yes | omit/soften predictive claims; distinguish v2 internal validation |
| PUB-07 | MODERATE | C | WGCNA ML2 | features/modules derive from same cohort; some module selection integrates same-cohort evidence | performance is not independent validation | no | no | yes | label internal/exploratory; do not call predictor validated |
| PUB-08 | MODERATE | D | Network provenance | Phase-2 narrative says edges preserved, renderer recomputes sparse edges | exact edge identity not established | no | optional lightweight comparison | yes | say “reconstructed from unchanged inputs” or compare edge lists |
| PUB-09 | MODERATE | F | Biological annotation | 57.9% of prioritized broad-module annotation rows unmapped because org.Hs.eg.db unavailable | broad theme/GO interpretation has incomplete coverage | environment fix | only if annotation claims are central | yes | disclose coverage; consider controlled reannotation as a future versioned analysis |
| PUB-10 | MODERATE | D | Gene models | most “LMM” outputs are fixed-effect LMs | method can be misreported | no | no | yes | name models correctly by unit |
| PUB-11 | MODERATE | B | WGCNA | adaptive grid, different 15k/5k universes, small CG25/no target R2 | module comparability/stability limited | no | no | yes | report as exploratory and avoid module-count comparison |
| PUB-12 | MODERATE | G | AHBA/WGCNA | M26 and MEdarkred are nominal after correction | overstatement would cross claimed FDR threshold | no | no | yes | use nominal/exploratory wording |
| PUB-13 | MINOR | F | Databases | package/resource versions and download dates incomplete | future enrichment may drift | no scientific change | no | yes | record versions/checksums/citations |
| PUB-14 | MINOR | D | Covariates | gene-wise model omits PMI while ME model includes it; rationale fragmented | readers may infer inconsistency | no | no | yes | explicitly distinguish model families and sensitivity rationale |
| PUB-15 | INFORMATIONAL | I | ML score | 671 rows, 30/119/522, 125/150 denominators and arithmetic verified | supports corrected publication record | no | no | no | retain canonical manifest/tests |
| PUB-16 | INFORMATIONAL | I | Donor CV | six gene-ML units have one sample/donor; pooled analyses use grouping/random effect | no detected donor leakage in audited splits | no | no | yes | report clearly |

# PUBLICATION BLOCKERS

1. **PUB-01:** submission-facing materials contain superseded numerical results and claims. Refresh them from the 30/119 corrected record; affected figures/tables and their captions need regeneration, not the upstream science.
2. **PUB-02:** Tier4/GSEA claims have no recoverable canonical result package. Recover exact outputs or rerun only Tier4/GSEA from frozen upstream effects with fully recorded databases/settings; dependent GSEA figures/tables and claims must then be regenerated.
3. **PUB-03:** Step 13 module-level ML2 lacks its generator. Recover the exact code and verify outputs, or remove the analysis from submission. If recovered code does not reproduce retained outputs, Step 13 figures/tables must be rerun.

# MUST FIX BEFORE SUBMISSION

- **PUB-04:** remove or version-regenerate legacy AHBA post-enrichment sensitivity feature sets; the all32 analysis can remain.
- **PUB-05:** provide a reproducible data-acquisition and environment record sufficient for reviewers/readers to understand how canonical files were generated.
- **PUB-12:** change corrected M26 and MEdarkred wording to nominal/exploratory everywhere and ensure OFC MEtan is described as an enrichment result, not diagnosis validation.
- State unambiguously that there are no validated biomarkers, no diagnosis BH-FDR-significant WGCNA modules, and no completed independent MDD cohort validation.

# SHOULD FIX BEFORE SUBMISSION

- **PUB-06/PUB-07:** limit ML claims to internal prioritization/internal validation; document same-cohort selection and correlated evidence components.
- **PUB-08:** soften “edge preservation” or verify exact old/new edge lists.
- **PUB-09:** disclose annotation mapping coverage and avoid strong theme claims from incompletely mapped modules.
- **PUB-10/PUB-11:** accurately name LM versus LMM and report WGCNA universe/grid/CG25 limitations.
- Add database/package versions, exact correction families, and region-specific universe definitions to methods/supplement.

# DOCUMENTATION ONLY

- Explain why gene-wise and eigengene formulas differ and that no result becomes diagnosis-FDR significant under tested ME formulas.
- Keep pooled Frontal Cortex separate from original frontal subregions and supplementary only.
- Distinguish structural core from its ML-prioritized annotation subset.
- Distinguish Tier2 over-representation from Tier4 ranked GSEA.

# OPTIONAL IMPROVEMENTS

- Add end-to-end fixture tests for normalization, LM/LMM extraction, tiers, enrichment universes, and figure-source consistency.
- Containerize or lock the complete Python/R stack and replace absolute paths with configuration.
- Add independent cohort replication and prospective module preservation/prediction analyses.

# VERIFIED / NO ACTION

- Requested Git ancestry and baseline are correct.
- Canonical ranking hash and counts are correct.
- Raw v1/v2 event checksums and denominators are correct.
- All 16 relevant tests pass.
- No within-unit donor leakage exists for the six gene-level ML units.
- Checked GeTMM, LM, BH, standardized-effect, Tier1–3, Tier2-enrichment, ML-score, AHBA, and WGCNA-enrichment arithmetic passes.
- Corrected published-AHBA shortlist analysis has no BH-FDR<0.10 result; M26 and MEdarkred are nominal; OFC MEtan retains an exploratory FDR<0.10 enrichment.
- Frontal corrected network inventory excludes pooled Frontal Cortex.
