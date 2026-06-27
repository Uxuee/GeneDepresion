"""Archived sensitivity script: create a clean ML v1 branch using Tier2_fixed only.

Inputs: existing ML_pilot_Tier2 performance, predictions, and feature-importance outputs.
The main clean workflow now uses 09_run_ML_v1_extended_prioritization.py.
Tier2FixedOnly is retained as a sensitivity/archive result, not the main clean workflow result.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(9)
