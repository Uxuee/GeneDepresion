# Upstream R Scripts

These R scripts are copied here for convenience so the clean pipeline snapshot is easier to inspect from inside code/length_gt_100_clean_pipeline.

Important: the current clean pipeline Step 2 does **not** rerun normalization. It reuses the completed Length > 100 normalized expression matrices from the main analysis output folder. The normalization method documented by the pipeline is:

filtered raw counts -> GeTMM -> log2(GeTMM + 1)

The original upstream normalization / sensitivity scripts are included here so the code package remains reproducible and auditable. The same scripts are also present one level up in the snapshot's code/ directory because this snapshot copied the full source D:/genes/GeneDepresion/code folder.

Most relevant normalization-related file:

- Step1_length_cutoff_sensitivity_GeTMM.R

Other R scripts are included because later pipeline stages and sensitivity analyses reused them.
