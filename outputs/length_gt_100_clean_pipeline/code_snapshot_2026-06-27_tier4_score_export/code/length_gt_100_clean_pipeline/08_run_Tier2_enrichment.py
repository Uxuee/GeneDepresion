"""Step 8: Reuse completed Tier 2 enrichment outputs.

Inputs: existing enrichment result tables and selected report-ready enrichment plots.
Outputs: clean enrichment folder with copied tables, summaries, and selected figures.
Does not rerun enrichment.
"""

from pipeline_lib import run_numbered_step


if __name__ == "__main__":
    run_numbered_step(8)

