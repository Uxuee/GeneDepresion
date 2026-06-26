"""Step 4: Reuse standardized effect-size outputs.

Inputs: completed standardized-effect-size table.
Outputs: clean copies plus sign-check summary.
Definition: standardized_effect_size_diagnosis = beta_diagnosis / SD(expression of gene).
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(4)

