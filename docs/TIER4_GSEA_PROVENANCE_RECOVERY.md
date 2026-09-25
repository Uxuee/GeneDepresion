# Tier4 and GSEA provenance recovery

## Result

The seven per-region Tier4 ranked signatures are exactly recoverable from Git
commit `7888e981d0cf2c011e45632c3245f3333a2d8f84`. Their blob ids, row counts, and
SHA-256 values are frozen in `provenance/tier4_ranked_signature_recovery_manifest.json`.
No fold/event records were reconstructed from aggregate tables.

The recovered generator ranks all tested genes by the absolute value of
`sign(standardized_effect_size_diagnosis) × -log10(p_value_diagnosis)`, breaking
ties by ascending p-value. The stored signed score preserves direction; the
absolute score controls rank.

An additional four-group combined ranked file is recoverable at commit
`037becc4766df6f18a0ba1bbde93c0af989cf117` (blob
`08a938639f586fb0fa0d8ee379a21af17141d5a9`, 95,319 rows, SHA-256
`3865443d878042717a4eb28293763b5ef1d18598ad399c562a650a854e24727c`).

## GSEA search

The search covered the working tree, ignored/untracked files, all reachable
branches and commits, Git stashes, and unreachable objects reported by
`git fsck`. Searches included `fgsea`, `gseGO`, `gseKEGG`, `GSEA`, `NES`,
permutation settings, and plausible output schemas. No exact GSEA generator,
settings, or result artifact was recovered.

Therefore Tier4 ranking is provenance-recovered, but Tier4 GSEA is
`excluded_from_publication` until its method and inputs are independently
specified and rerun as a new versioned analysis.
