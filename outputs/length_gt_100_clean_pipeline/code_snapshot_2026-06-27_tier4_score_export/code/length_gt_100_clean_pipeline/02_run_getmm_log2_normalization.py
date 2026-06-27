"""Step 2: Reuse GeTMM + log2 normalized expression matrices.

Inputs: completed normalized log2(GeTMM + 1) matrices by brain region.
Outputs: clean copies under outputs/length_gt_100_clean_pipeline/02_normalized_expression.
Does not recompute normalization unless future edits add that behavior.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(2)

