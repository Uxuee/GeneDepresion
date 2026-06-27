"""Step 1: Prepare Length > 100 inputs.

Inputs: completed Length > 100 gene-filter summary, metadata, and LMM table.
Outputs: clean input summaries under outputs/length_gt_100_clean_pipeline/01_inputs.
Reuses existing outputs and does not rerun filtering.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(1)

