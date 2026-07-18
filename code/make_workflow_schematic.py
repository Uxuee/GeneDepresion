"""Create a deterministic workflow schematic for the manuscript package.

This is not a generative image. It uses fixed text boxes and arrows so the
figure can be reproduced from repository code.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(r"D:/genes/GeneDepresion")
OUT = ROOT / "outputs/length_gt_100_clean_pipeline/14_Preliminary_Manuscript_Package/02_figures_main_original"

STEPS = [
    "GeTMM-log2\nexpression",
    "gene-wise\ndiagnosis models",
    "Tier2\ncandidate genes",
    "ML #1 gene\nprioritization",
    "AHBA published-\nmodule context",
    "MDD WGCNA\nmodules",
    "frontal subregion\nWGCNA",
    "biological annotation\n/ network plots",
    "module-level\nML #2",
    "robustness /\nexternal validation",
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(14, 4.8))
    ax.axis("off")
    coords = []
    for y in (0.68, 0.30):
        for x in (0.06, 0.28, 0.50, 0.72, 0.94):
            coords.append((x, y))
    for i, (x, y) in enumerate(coords):
        box = FancyBboxPatch(
            (x - 0.09, y - 0.09),
            0.18,
            0.18,
            boxstyle="round,pad=0.018,rounding_size=0.018",
            facecolor="#F7F8FA",
            edgecolor="#334E68",
            linewidth=1.4,
            transform=ax.transAxes,
        )
        ax.add_patch(box)
        ax.text(x, y, STEPS[i], ha="center", va="center", fontsize=9.5, transform=ax.transAxes)
        if i < len(coords) - 1:
            x2, y2 = coords[i + 1]
            style = "arc3,rad=-0.35" if i == 4 else "arc3,rad=0"
            arrow = FancyArrowPatch(
                (x + 0.095, y),
                (x2 - 0.095, y2),
                arrowstyle="-|>",
                mutation_scale=13,
                color="#334E68",
                linewidth=1.2,
                connectionstyle=style,
                transform=ax.transAxes,
            )
            ax.add_patch(arrow)
    ax.set_title(
        "Region-aware exploratory MDD transcriptomics workflow",
        fontsize=14,
        fontweight="bold",
        pad=16,
    )
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"Figure_1_workflow_schematic_REGENERATED_BY_SCRIPT.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
