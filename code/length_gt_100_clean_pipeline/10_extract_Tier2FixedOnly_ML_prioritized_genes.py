"""Step 10: Extract Tier2FixedOnly ML-prioritized genes.

Inputs: clean Tier2_fixed-only ML feature importance and clean Tier 2 candidate pools.
Outputs: ranked, high-confidence, moderate-confidence, gene-list, and plot files.
Does not use Tier 4 or strict_univariate feature sets.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(10)

