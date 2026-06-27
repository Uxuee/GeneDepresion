"""Archived sensitivity script: extract Tier2FixedOnly ML-prioritized genes.

Inputs: clean Tier2_fixed-only ML feature importance and clean Tier 2 candidate pools.
The main clean workflow now uses 10_extract_extended_ML_prioritized_Tier2_genes.py.
Tier2FixedOnly is retained as a sensitivity/archive result, not the main clean workflow result.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(10)
