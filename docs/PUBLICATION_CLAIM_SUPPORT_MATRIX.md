> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Publication claim-support matrix

This matrix applies to current report/manuscript material at baseline `ca25ec35...`. “Supported” means the strength of the proposed wording matches the canonical evidence; it does not mean biological validation.

| Claim | Canonical source/value | Evidence class | Independence | Supported? | Problem / appropriate wording |
|---|---|---|---|---|---|
| There are 671 Tier2 gene-region pairs | corrected ranking; 671 | descriptive | same discovery data | yes | say gene-region pairs, not necessarily unique genes |
| 35 high and 130 moderate candidates | legacy directory | unsupported/currently obsolete | — | no | current values are 30 and 119; 149 combined |
| ML confidence is a probability of truth | score implementation | unsupported | correlated components | no | call heuristic internal priority class |
| ML uses corrected weights | 0.40/0.15/0.10/0.15/0.10/0.05/0.05 | descriptive | — | yes | several reports still show old first three weights |
| v1 score uses Tier4 top-k evidence | scorer filters Tier2_fixed | unsupported/currently obsolete | — | no | Tier4 v1 runs are historical only |
| v1/v2 predictive performance validates biomarkers | internal CV | unsupported | no external cohort | no | performance is exploratory/internal; v1 feature pool is full-data Tier2 |
| High candidates are stable across 125 v1 and 150 v2 events | raw events and corrected score | internally ML-prioritized | same cohort/CV events | partly | “internally stable under the executed grids”; components are correlated |
| No gene-wise diagnosis finding is FDR significant | regional BH outputs | FDR-supported null conclusion | discovery cohort | yes | retain precise correction family |
| Tier2 represents significant DE genes | Tier2 p<.01 + abs effect>.8 | nominal/exploratory | discovery cohort | no if “significant” means FDR | call thresholded exploratory candidates |
| Tier2 pathway result GO:0017017 is enriched | p=.000241, BH=.02660 in HF down GO:MF family | FDR-supported enrichment | annotation database external; genes selected same data | yes | say over-represented among selected candidates, not disease mechanism |
| Tier4 ranked GSEA supports conclusions | canonical files absent | unsupported/provenance gap | unknown | no | remove until outputs/settings recovered or regenerated |
| M26 is enriched among corrected VS up ML genes at FDR<.10 | corrected p=.0280, q=.8823 | nominal/exploratory | external fixed module, internal gene selection | no | “nominal overlap/context”; not external MDD validation |
| M26 validates the MDD signature externally | same result | functional/context annotation | not an independent MDD cohort | no | AHBA is biological context only |
| VS MEdarkred retains FDR<.10 ML enrichment | p=.00217, q=.10851 | nominal/exploratory | same cohort genes/modules | no | state nominal only |
| OFC BA11 MEtan retains corrected ML enrichment | high set p=.001781, q=.03205 | FDR-supported enrichment | same cohort genes/modules | yes with caveat | exploratory within-study convergence, not diagnosis association |
| WGCNA modules are diagnosis-derived | construction code | descriptive | diagnosis not used in construction | no | modules are unsupervised with respect to diagnosis |
| WGCNA found diagnosis-significant modules | module-trait tables | unsupported | — | no | no module survives Diagnosis BH-FDR |
| Ten modules are validated disease modules | integrated shortlist | unsupported | multiple same-cohort evidence streams | no | exploratory prioritized modules/external-validation targets |
| Module ME p<.05 confirms MDD association | raw p only | nominal/exploratory | same cohort | no | explicitly label nominal and corrected-null |
| Broad Frontal Cortex equals frontal subregions | distinct assignments/samples | unsupported | — | no | pooled frontal is supplementary; main report uses CG25/dlPFC/OFC |
| Structural-core genes are ML-selected | stable-core method | descriptive mismatch | — | no | structural core is kME/hub/membership based; ML flags are annotation |
| Phase 2 preserved module assignments | assignment hashes | descriptive verified | — | yes | byte-preserved input assignments |
| Phase 2 preserved exact network edges | renderer recomputes correlations/edges | unsupported as written | — | no | say deterministically reconstructed from unchanged expression/membership unless compared |
| Module-level ML2 predicts diagnosis | internal grouped CV | internally validated/exploratory | modules and outcomes from same cohort | only with weak wording | “carried exploratory internal signal”; not diagnostic validation |
| Step 13 is reproducible | generator absent | unsupported | — | no | recover generator or omit from publication |
| AHBA all32 module-score ML is independent of corrected shortlist | fixed published modules | sensitivity/internal ML | direct feature definition independent | yes | still same MDD cohort and not external validation |
| AHBA FDR-selected module-score sensitivities remain current | legacy post-enrichment subsets | stale/legacy | selected from old enrichment | no | remove/regenerate; corrected enrichment has no FDR<.10 modules |
| Biological theme labels are formal mechanisms | interpretive annotation | functional/context annotation | incomplete mapping | no | keyword/manual interpretive summary only; mapping coverage 42.1% in audited table |
| Results are robust | mixed sensitivity evidence | sensitivity result | limited | only if scoped | specify what was stable: directions, counts, or threshold status; do not use globally |
| Results are reproducible from a fresh clone | dependency dry run | unsupported | — | no | local raw data, paths, environments, missing Tier4/Step13 block reproduction |
| Results are externally validated | no independent MDD cohort complete | unsupported | no | no | explicitly state external validation is incomplete |

## Required wording guardrails

- “Significant” must be paired with the exact adjusted threshold and family; raw p<.05 is “nominal.”
- “Enriched” is acceptable only with the named candidate set, universe, test, and adjusted status.
- “Predicts” should be replaced by “showed exploratory internal cross-validated signal” unless an independent cohort is evaluated.
- “Validated,” “biomarker,” “replicated,” and “independent evidence” are not supported for the MDD findings.
- AHBA and GO/MSigDB evidence are annotation/context, not validation of disease association.
