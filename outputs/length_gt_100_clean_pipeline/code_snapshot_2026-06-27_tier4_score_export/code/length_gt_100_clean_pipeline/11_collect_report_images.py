"""Step 11: Collect report-ready images.

Inputs: selected figures from the old analysis and clean Tier2FixedOnly ML outputs.
Outputs: images/tier2, images/enrichment, images/ml_v1, images/ml_prioritized, and inventory.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(11)

