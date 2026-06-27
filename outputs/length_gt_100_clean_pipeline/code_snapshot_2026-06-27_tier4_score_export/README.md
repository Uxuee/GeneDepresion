# Code Snapshot: Length > 100 Clean Pipeline

This folder is a reproducibility snapshot of the full D:/genes/GeneDepresion/code directory used with the Length > 100 clean pipeline.

Source project: D:/genes/GeneDepresion
Git branch: $branch
Git commit: $commit
Snapshot date: 2026-06-27

The snapshot includes the clean pipeline scripts and the upstream analysis scripts called by the workflow, including the Tier 4 ranking / ML-prioritized Tier 2 code. Result outputs are not duplicated here; this is code only.

Important Tier 4 rule now documented in the workflow:

- signed_score keeps the MDD direction.
- bs_score is used for ranking magnitude / top-k feature selection.
- Both up_in_MDD and down_in_MDD genes are allowed into Tier4 top-k feature sets.
- Final ML-prioritized tables include signed_score, bs_score, and direction.

See SNAPSHOT_MANIFEST.json for the full copied file list and git status at snapshot time.
