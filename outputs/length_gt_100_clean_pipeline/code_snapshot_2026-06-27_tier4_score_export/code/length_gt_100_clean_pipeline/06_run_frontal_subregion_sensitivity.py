"""Step 6: Reuse frontal subregion fixed-effect sensitivity outputs.

Inputs: completed frontal subregion standardized-effect table.
Outputs: clean frontal Tier 2 genes and model feasibility documentation.
Does not rerun the frontal model.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(6)

