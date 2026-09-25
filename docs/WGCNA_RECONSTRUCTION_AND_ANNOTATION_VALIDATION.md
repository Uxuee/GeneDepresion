# WGCNA reconstruction and annotation validation

> **Status:** Phase-3B validation candidate. The existing WGCNA assignments remain canonical. This analysis does not establish biomarkers, causal modules, or external replication.

## Reconstruction inputs

The six analysis units used frozen `log2(GeTMM + 1)` matrices. The exact paths, SHA-256 hashes, dimensions, and sample IDs are in `inputs_manifest/wgcna_reconstruction_input_manifest.csv` and `wgcna_reconstruction_sample_ids.csv`. Metadata hashes are in `wgcna_reconstruction_metadata_manifest.csv`. The verified broad-region matrices are the retained Length>100 sensitivity outputs; an alternate `genedata` copy was rejected after it failed to reproduce Insula's historical soft power.

Diagnosis was not used for MAD filtering, adjacency, TOM construction, tree cutting, module merging, or solution selection.

## Exact recovered WGCNA method

The executable method is documented in `docs/WGCNA_CANONICAL_PARAMETER_AUDIT.md`. In brief: top-MAD filtering; signed bicor network with `maxPOutliers=0.05`; powers 1–20 and target signed R² 0.80 with a maximum-R² fallback; signed TOM; average-linkage clustering; an 18-combination dynamic-tree-cut/merge grid; and selection by usability, gene-universe size, non-grey module count, then largest-module fraction. CG25 correctly invoked the fallback and selected power 17.

## Reproducible R environment

Historical execution was under R 4.4.2 but had no complete lockfile. Phase 3B uses R 4.6.1 with an isolated library. Actual versions and sources are recorded in `environment/phase3b/renv.lock`, `package_versions.csv`, and `sessionInfo.txt`. Therefore the historical and validation environments are explicitly distinguished.

## Region-by-region reconstruction and canonical comparison

The equivalence thresholds were written before enrichment was evaluated: EXACT requires ARI=1 and complete mapped assignment agreement; EFFECTIVELY_EQUIVALENT requires ARI≥0.90, median best Jaccard≥0.90, mapped gene agreement≥0.90, and median absolute eigengene correlation≥0.95.

| Region | Status | Old modules | New modules | ARI | Median best Jaccard | Median |eigengene r| | ML conclusion stable? | Old mapping | New mapping | Recommended action |
|---|---|---:|---:|---:|---:|---:|---|---:|---:|---|
| Ventral Striatum | EXACT | 50 | 50 | 1e+00 | 1e+00 | 1e+00 | TRUE | 70.4% | 97.9% | retain existing canonical modules; independently verified |
| Insula | EXACT | 29 | 29 | 1e+00 | 1e+00 | 1e+00 | TRUE | 64.8% | 100.0% | retain existing canonical modules; independently verified |
| Hippocampal Formation | EXACT | 69 | 69 | 1e+00 | 1e+00 | 1e+00 | TRUE | 44.0% | 92.9% | retain existing canonical modules; independently verified |
| Cingulate gyrus 25 | EXACT | 27 | 27 | 1e+00 | 1e+00 | 1e+00 | TRUE | not available | 100.0% | retain existing canonical modules; independently verified |
| dlPFC BA8/9 | EXACT | 21 | 21 | 1e+00 | 1e+00 | 1e+00 | TRUE | not available | 87.3% | retain existing canonical modules; independently verified |
| OFC BA11 | EXACT | 18 | 18 | 1e+00 | 1e+00 | 1e+00 | TRUE | not available | 95.1% | retain existing canonical modules; independently verified |

Module colors were treated as arbitrary labels. Comparison used overlap matrices, Jaccard similarity, maximum-overlap alignment, Adjusted Rand Index, aligned eigengene correlations, and top-hub kME comparison. Full tables are under `module_comparison/` and `eigengene_validation/`.

## Eigengene and Diagnosis validation

The reconstructed eigengene models produced 0 module(s) with Diagnosis BH q<0.05 and 0 with q<0.10. The validation is a reproducibility check, not a new significance search. The executable historical formula and its duplicate alias fields are documented in the parameter audit; every old and reconstructed formula is retained in the comparison CSV.

## Corrected ML-gene enrichment validation

Only the canonical corrected ranking was used: 671 gene-region rows, with 30 high, 119 moderate, 522 low, and 149 high+moderate. Four sets were tested: global high, global high+moderate, same-region high, and same-region high+moderate. Tests were one-sided Fisher exact tests using the region-specific WGCNA network genes as the universe; BH correction was performed within region and gene set across non-grey modules.

OFC BA11 MEtan aligned to MEtan (Jaccard 1e+00). Global high+moderate enrichment: old p=3.81e-03, q=6.86e-02; reconstructed p=3.81e-03, q=6.86e-02.

Ventral Striatum MEdarkred aligned to MEdarkred (Jaccard 1e+00). Global high+moderate enrichment: old p=2.17e-03, q=1.09e-01; reconstructed p=2.17e-03, q=1.09e-01.

ML-gene overlap and module-eigengene Diagnosis association are different evidence layers and are not independent replication.

## Functional annotation

The ten pre-existing report-prioritized modules were aligned to their reconstructed counterparts and annotated with GO BP/MF/CC (`clusterProfiler` + `org.Hs.eg.db`), Reactome (`ReactomePA`), and MSigDB Hallmark (`msigdbr` + `clusterProfiler::enricher`). Entrez identifiers were converted explicitly with `org.Hs.eg.db`. Each region's mapped network genes—not all human genes—formed the enrichment background. Raw p, BH q, overlaps, mapped/unmapped module genes, and mapped universe sizes are retained in `functional_annotation/`.

## Mapping coverage

The historical incomplete layer mapped 1,796/4,267 genes (42.1%). The Phase-3B prioritized reconstructed modules mapped 2611/2783 (93.8%). Among identifiers present in both mapping collections, 1089 genes that were historically unmapped are newly mapped and 156 remain unmapped. The full transition table also labels genes unique to either module collection. Because the historical collection included pooled broad Frontal Cortex whereas Phase 3B uses original frontal subregions, the headline denominators differ and are not a paired biological comparison.

Unmapped identifiers are retained in `still_unmapped_prioritized_module_genes.csv`; they were not silently discarded.

## Evidence layers and interpretation

`prioritized_module_evidence_table.csv` separates Diagnosis beta/p/q, corrected ML overlap, network alignment, functional terms, and mapping coverage. Functional enrichment describes the genes in a module; it does not validate the module as disease-associated. Structural hub status, ML overlap, eigengene association, and pathway annotation answer different questions.

## Biological interpretation changes

The reconstruction does not require replacing or renaming the canonical modules. Clean annotation improves identifier coverage and provides versioned support for interpretation, but it does not convert exploratory modules into validated disease mechanisms.

## Reproducibility and provenance

Every reconstruction has a region manifest; comparison/annotation and frozen inputs have stage manifests; a recursive SHA-256 inventory is in `manifests/phase3b_output_file_manifest.csv`. `provenance/canonical_pipeline.json` registers Phase 3B as `validation_candidate`, not canonical.

## Recommended canonical decision

**A. retain existing canonical modules; reconstruction verifies them.**

This decision depends only on reconstruction similarity. It does not depend on which network yields more significant Diagnosis, ML-overlap, or pathway results.

## Validation tests

The Phase-3B R validator passed all 38 input, parameter, assignment, manifest, checksum, ML-list, enrichment-background, and mapping-accounting checks. The complete Python unittest discovery passed all 36 repository tests, including 9 Phase-3B tests. `audit/validate_canonical_pipeline.py` also passed with no errors or warnings and reverified the canonical ranking SHA-256. The test commands and runtime workaround are recorded in `manifests/phase3b_test_execution.json`.

## Limitations

- No WGCNA module is called a validated disease module, biomarker, causal module, or replicated module.
- Functional annotation and corrected ML overlap are exploratory support layers.
- External validation remains incomplete.
- Package-version differences between historical R 4.4.2 and Phase-3B R 4.6.1 are recorded explicitly.
- The pooled broad Frontal Cortex analysis is not part of the six-unit main validation.
