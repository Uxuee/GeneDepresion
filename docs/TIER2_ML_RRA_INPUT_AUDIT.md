> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Tier2 → ML → RRA input audit

Date: 2026-10-06  
Branch base: `3e41a8a4ad97781dfcd4a60628181bbf9c143074`  
Purpose: establish the candidate identities before implementing or running the primary Tier2 ML–RRA analysis.

## Canonical Tier2 source

The upstream source is:

`outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv`

- SHA-256: `1672ddf95a4f4c9e44ae42e7accc208c311e2e1da8533e1808dcf91bbf7b67a1`
- Git blob: `33b37fda8626c09783ab3cb7b7a394cc826e0d5b`
- Rows: 671 region–gene pairs
- Duplicate region–gene identities: 0
- Canonical registry stage: stage 06, “Evidence tiers and Tier2 candidate sets”
- Canonical generator: `code/length_gt_100_clean_pipeline/05_create_evidence_tiers.py`
- Packaging provenance commit: `a318f1a` (`Create clean Tier2-only pipeline from raw genes to ML-prioritized candidates`)

This artifact precedes historical ML scoring and is therefore the primary identity source. The new analysis reads only `region_or_subregion`, `gene_id`, and `gene_symbol`. Diagnosis coefficients, p-values, effect sizes, directions, historical ML scores, confidence labels, selection frequencies, and model-support fields are not imported into fitting or RRA.

## Audited regional counts

| Region/subregion | Tier2 pairs |
|---|---:|
| Ventral Striatum | 146 |
| Insula | 70 |
| Hippocampal Formation | 127 |
| CG25 / Cingulate gyrus 25 | 186 |
| dlPFC BA8/9 | 72 |
| OFC BA11 | 70 |
| **Total** | **671** |

The observed counts reproduce the validation expectations without using them to construct membership.

## Identity cross-check

The later corrected historical ML table contains the same 671 region–gene identities: 671 intersections, 0 upstream-only identities, and 0 later-only identities. This is a provenance cross-check only. The later table is not an input to the new analysis.

## Expression matching

All Tier2 IDs matched their corresponding frozen normalized-expression matrix exactly:

| Analysis unit | Tier2 IDs | Present in expression | Missing |
|---|---:|---:|---:|
| Ventral Striatum | 146 | 146 | 0 |
| Insula | 70 | 70 | 0 |
| Hippocampal Formation | 127 | 127 | 0 |
| CG25 | 186 | 186 | 0 |
| dlPFC BA8/9 | 72 | 72 | 0 |
| OFC BA11 | 70 | 70 | 0 |

Sample counts and diagnosis labels also reproduced the frozen all-gene run: 48, 48, 43, 27, 48, and 47 samples respectively. Pooled Frontal Cortex is excluded.

## Frozen boundaries

- The completed all-gene sensitivity output remains frozen at `outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/`.
- Its annotated consensus checksum is `0aefec5fcfa100deafa1525c2d31521ae8c00ca890ef4e1b4145b1b861c2c87e`.
- The old corrected canonical ML ranking remains frozen at SHA-256 `6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1`.
- No Tier4, WGCNA, AHBA, all-gene RRA rank, or historical ML quantity will define or alter the new candidate universe.

## Audit decision

The canonical Tier2 identities are established, all 671 pairs match expression safely, regional counts reproduce, and no stop condition is triggered. Implementation may proceed using this upstream identity file only.
