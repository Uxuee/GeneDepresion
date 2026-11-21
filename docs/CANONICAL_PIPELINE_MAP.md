# Canonical pipeline map

See [current descriptive workflow](PIPELINE.md) and [Step 3 reconciliation](STEP3_PROVENANCE_RECONCILIATION.md). The authoritative registry is `provenance/canonical_pipeline.json`; its prior definition is preserved in `provenance/consolidation/canonical_pipeline.before.json`.

Current primary order: Tier2-RRA. Corrected weighted v1 + strict-v2 scoring is the preserved historical canonical benchmark; it alone retains the 30/119/522 labels. The registry keeps old numeric identifiers for compatibility and adds descriptive identifiers for new analyses. No output path is renamed. All stages are frozen evidence without automatic execution commands.

Source data -> GeTMM -> frozen regional LM -> frozen Tier2 -> Tier2-RRA -> final synthesis. Weighted ML, all-gene RRA, statistical benchmarks and network annotations remain separate evidence streams. Completed local literature comparisons do not establish independent validation.
