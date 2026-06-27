"""Step 10: Restore extended ML-prioritized Tier 2 genes.

Inputs: previous stable ML_prioritized_Tier2_genes output package.
Outputs: clean extended ML-prioritized Tier 2 tables, gene lists, plots, and manifest.
Final genes remain ML-prioritized Tier 2 candidates; Tier 4 is support only.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(10)
