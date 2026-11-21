# Final synthesis data and Git audit

## Scope and frozen starting point

- Repository: `D:/genes/GeneDepresion_RRA`
- Frozen input commit: `437b2b63d3337896cf71085160d1545df687b1ed`
- Synthesis branch: `codex/final-evidence-synthesis`
- The working tree was clean before this audit and synthesis work began.
- Outputs 17–21 are read-only evidence sources. Their checksum-ledger SHA-256 values at the start of this phase were:
  - `17_ML_RRA_consensus`: `117c65e6cc796a9c6757b9d9e92ed493cc8bb6300d0873ea3d5cd2801ba3634c`
  - `18_Tier2_ML_RRA_prioritization`: `3c45b0a92e4a5a3b8c9fc126097c312d492141734ea8e3fe72534a91152e950e`
  - `19_Tier2_ML_benchmark`: `4e17dc8234317efc9bce92aa03fb7ae63ff72a128392d202f2f31a105c09cf10`
  - `20_Tier2_ML_conditional_null`: `d95f66712ed0ba80b7acefee7c3ee992cc14923128a47a5ddcaf9cff12e0d545`
  - `21_Statistical_method_benchmark`: `5626960b44eb446c023fd621f2808916f591667666be631a89691237a1b025ed`

## Source-data locations

Raw RNA-seq counts remain in their existing read-only source path:

`D:/genes/GeneDepresion/outputs/step1_length_sensitivity/Length_gt_100/03_filtered_raw_counts/keep_all_samples/`

Normalized expression matrices remain in their existing upstream project/output paths, including the canonical GeTMM and log2(GeTMM + 1) products documented by the frozen statistical and ML manifests. They are not copied into Phase 22.

The final synthesis reads compact frozen products from:

- `17_ML_RRA_consensus/`
- `18_Tier2_ML_RRA_prioritization/`
- `19_Tier2_ML_benchmark/`
- `20_Tier2_ML_conditional_null/`
- `21_Statistical_method_benchmark/`
- canonical WGCNA assignment/context tables only where linkage is unambiguous

## What is tracked

The repository selectively tracks code, tests, documentation, manifests, checksum ledgers, compact result tables, and useful figures. Some prior milestone outputs under the globally ignored `outputs/` tree were intentionally force-added as publication/provenance records.

Phase 22 will add only:

- synthesis code and tests;
- the 671-row evidence matrix and compact top-10/top-25/FDR tables;
- six concise regional summaries;
- provenance/validation manifests and checksums;
- compact supervisor-facing figures;
- this audit and the final synthesis report.

## What is intentionally ignored or excluded

The root `.gitignore` excludes generated CSV/TSV data by default, the general `outputs/` tree, local configuration, caches, Python environments, local R libraries, and raw AHBA downloads. Phase 22 will not stage:

- raw-count matrices;
- full normalized-expression matrices;
- `.runtime/`, `.venv/`, R installations, or package libraries;
- `__pycache__/`, `.pytest_cache/`, temporary files, or editor/OS noise;
- permutation fold-level outputs under `20_Tier2_ML_conditional_null/permutations/`;
- duplicated copies of outputs 17–21;
- any source dataset solely to make the repository self-contained.

## Derived publication-level output

The only new output root is:

`outputs/length_gt_100_clean_pipeline/22_Final_evidence_synthesis/`

It contains compact derived tables, figures, region summaries, provenance, and checksums. Every master-matrix column is mapped to its frozen source in `manifest/evidence_provenance.csv`.

## Remote audit

At the start of this phase, `origin` resolved to the local path `D:/genes/GeneDepresion` for fetch and push, not to a GitHub URL. The synthesis may be committed locally, but this remote must not be treated as the requested GitHub destination. Per the task stop rule, no credentials or remotes will be changed automatically and no push will be attempted unless an expected GitHub origin is available.

## Audit conclusion

The planned synthesis can proceed without copying or staging raw matrices, normalized full matrices, runtimes, caches, or bulky permutation intermediates. File sizes will be checked again immediately before staging.
