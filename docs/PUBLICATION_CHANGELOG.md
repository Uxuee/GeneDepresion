> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# Publication Changelog

## 2026-09-25 — canonicalisation pending approval

- Designated `10_ML_prioritized_Tier2_corrected/` as authoritative.
- Preserved `10_ML_prioritized_Tier2_extended/` as legacy provenance.
- Recorded corrected counts 671 total, 30 high, 119 moderate, 522 low.
- Recorded v1 denominator 125 and v2 denominator 150 non-baseline events.
- Documented corrected direction handling, baseline exclusion, missing-sign gate, and weights.
- Added canonical pipeline, provenance, result index, downstream audit, checklist, manifest, and integrity test.
- No upstream analysis, ML training, score, threshold, or scientific result was rerun or changed.
- No commit or push has been made.

## Retained correction audit

- Legacy: 35 / 130 / 506.
- Corrected: 30 / 119 / 522.
- Five high→moderate, sixteen moderate→low, 650 unchanged.
- Baseline removal caused the largest category effect.
- Verified weights caused one independent category change.
- RF sign correction and the missing-sign gate caused no independent category changes here.
- Fully corrected: 21 category changes.
