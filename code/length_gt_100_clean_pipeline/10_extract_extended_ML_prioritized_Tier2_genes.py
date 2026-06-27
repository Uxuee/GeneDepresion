"""Step 10: Run the stable extended ML-prioritized Tier 2 procedure.

Inputs: Tier 2 gene sets plus existing ML v1 and strict-v2 feature-importance/stability outputs.
Outputs: clean extended ML-prioritized Tier 2 tables, gene lists, plots, and manifest.
Final genes remain ML-prioritized Tier 2 candidates; Tier 4 is support only.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(10)
