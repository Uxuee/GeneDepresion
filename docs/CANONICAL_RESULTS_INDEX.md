# Current results index

See [PIPELINE.md](PIPELINE.md) for the stage table and `provenance/canonical_outputs_manifest.json` for hashes, provenance, availability and ranking roles.

- Current primary: `18_Tier2_ML_RRA_prioritization/Tier2_consensus_rra_rankings.csv`.
- Preserved weighted benchmark: `10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv`.
- Sensitivity: `17_ML_RRA_consensus/`.
- Dependent benchmarks: `19_Tier2_ML_benchmark/`, `20_Tier2_ML_conditional_null/`, `21_Statistical_method_benchmark/`.
- Evidence synthesis: `22_Final_evidence_synthesis/`.

These output paths are under `outputs/length_gt_100_clean_pipeline/`. Existing literature retrieval, mappings, statistical comparisons and intersections are under `external_validation/`; the newer primary comparisons are under `external_validation/rra_validation/tier2_branch_primary/`. They are completed local work, with independent-cohort validation still incomplete. Superseded weighted rankings, excluded Tier4 GSEA and unverified module-level ML2 do not become current by consolidation.
