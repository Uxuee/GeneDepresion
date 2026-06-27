"""Step 3: Reuse broad-region LMM outputs.

Inputs: completed Length > 100 all-gene LMM table and region summary.
Outputs: clean LMM copies under outputs/length_gt_100_clean_pipeline/03_broad_region_LMM.
LMM means linear mixed model; this script does not rerun the model.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(3)

