"""Step 5: Create clean evidence-tier outputs.

Inputs: standardized-effect-size full table.
Outputs: Tier 1, Tier 2, Tier 3, and all-gene evidence-tier tables.
Tier 4 is intentionally excluded from this clean Tier2-only pipeline.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(5)

