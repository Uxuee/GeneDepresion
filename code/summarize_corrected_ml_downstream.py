"""Create deterministic Phase-2 comparison tables, figures, and reports."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from canonical_downstream_paths import CANONICAL_RANKING, LEGACY_RANKING, PHASE2_ROOT, PROJECT_ROOT

DOC = PROJECT_ROOT / "docs" / "DOWNSTREAM_CORRECTED_RESULT_COMPARISON.md"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    old = pd.read_csv(LEGACY_RANKING, dtype={"gene_id": str})
    new = pd.read_csv(CANONICAL_RANKING, dtype={"gene_id": str})
    high = "High-confidence ML-prioritized Tier 2"
    moderate = "Moderate-confidence ML-prioritized Tier 2"
    rows = []
    for region in sorted(set(old.region_or_subregion) | set(new.region_or_subregion)):
        for version, frame in [("legacy", old), ("corrected", new)]:
            sub = frame[frame.region_or_subregion.eq(region)]
            rows.append({"version": version, "region": region, "high": int(sub.ML_priority_class.eq(high).sum()), "moderate": int(sub.ML_priority_class.eq(moderate).sum()), "high_plus_moderate": int(sub.ML_priority_class.isin([high, moderate]).sum())})
    membership = pd.DataFrame(rows)
    membership.to_csv(PHASE2_ROOT / "ML_membership_counts_before_after_by_region.csv", index=False)

    ahba = pd.read_csv(PHASE2_ROOT / "AHBA_corrected/published_AHBA_ML_enrichment_before_after.csv")
    wgcna = pd.read_csv(PHASE2_ROOT / "WGCNA_ML_annotations_corrected/all_WGCNA_ML_enrichment_before_after.csv")
    changed_a = ahba[ahba.significance_changed].copy()
    changed_w = wgcna[wgcna.significance_changed].copy()

    fig_dir = PHASE2_ROOT / "figures"
    fig_dir.mkdir(exist_ok=True)
    pivot = membership.pivot(index="region", columns="version", values="high_plus_moderate")
    ax = pivot[["legacy", "corrected"]].plot(kind="bar", figsize=(10, 5), color=["#9E9E9E", "#3366AA"])
    ax.set_ylabel("High + moderate gene-region pairs")
    ax.set_title("ML-prioritised membership before and after correction")
    ax.legend(title="Result")
    plt.tight_layout()
    plt.savefig(fig_dir / "ML_membership_before_after_by_region.png", dpi=220)
    plt.close()

    corrected_w = pd.concat([
        pd.read_csv(PHASE2_ROOT / "WGCNA_ML_annotations_corrected/broad_regions_ML_module_enrichment_corrected.csv"),
        pd.read_csv(PHASE2_ROOT / "WGCNA_ML_annotations_corrected/frontal_subregions_ML_module_enrichment_corrected.csv"),
    ], ignore_index=True)
    show = corrected_w[(corrected_w.module_role == "coexpression_module") & (corrected_w.p_value < .05)].sort_values(["fdr_BH", "p_value"]).head(30)
    if len(show):
        labels = show.region + " | " + show.module + " | " + show.gene_set.str.replace("_ML_genes", "", regex=False)
        y = np.arange(len(show))
        fig, ax = plt.subplots(figsize=(11, max(5, .34 * len(show))))
        pts = ax.scatter(-np.log10(show.p_value.clip(lower=np.nextafter(0, 1))), y, s=35 + 20*show.overlap_count, c=show.fdr_BH, cmap="viridis_r")
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8); ax.invert_yaxis()
        ax.set_xlabel("-log10 one-sided Fisher p-value")
        ax.set_title("Corrected ML-gene enrichment in preserved WGCNA modules")
        fig.colorbar(pts, ax=ax, label="BH-FDR")
        fig.tight_layout(); fig.savefig(fig_dir / "corrected_WGCNA_ML_enrichment_dotplot.png", dpi=220); plt.close(fig)

    def md_table(frame: pd.DataFrame, cols: list[str]) -> str:
        if frame.empty: return "None."
        return frame[cols].to_markdown(index=False)

    a_cols = ["candidate_set", "scope", "published_module", "overlap_genes_old", "overlap_genes_corrected", "p_value_hypergeometric_old", "p_value_hypergeometric_corrected", "fdr_BH_old", "fdr_BH_corrected", "old_significance", "corrected_significance", "genes_lost", "genes_gained"]
    w_cols = ["analysis_unit", "region", "module", "gene_set", "overlap_count_old", "overlap_count_corrected", "p_value_old", "p_value_corrected", "fdr_BH_old", "fdr_BH_corrected", "significance_category_old", "significance_category_corrected", "genes_lost", "genes_gained"]
    lines = [
        "# Downstream Corrected Result Comparison", "",
        "## Scope", "",
        "This report isolates changes caused by replacing the legacy 35-high/130-moderate ranking with the canonical 30-high/119-moderate ranking. No expression normalisation, diagnosis model, tier, ML training/scoring, AHBA module, WGCNA module, eigengene, connectivity, or module-level ML model was rebuilt.", "",
        "## Membership", "", membership.to_markdown(index=False), "",
        "## AHBA", "",
        f"Across the directly comparable ML-dependent tests, {len(changed_a)} rows changed nominal/FDR significance category. The corrected analysis has no published-AHBA module enrichment surviving BH-FDR < 0.10. The former Ventral Striatum up-in-MDD M26 result moved from FDR<0.10 to nominal p<0.05 after TSGA10 left the overlap. These changes are annotation/context changes, not evidence for or against a validated disease mechanism.", "",
        md_table(changed_a, a_cols), "",
        "## WGCNA", "",
        f"Underlying module assignments were preserved byte-for-byte. {len(changed_w)} broad/frontal enrichment rows changed nominal/FDR category. Broad Ventral Striatum MEdarkred high+moderate enrichment moved from FDR<0.10 to nominal after PRELID2 and TMEM269-DT (gene IDs 153768 and 728975) left its overlap. In the original frontal subregions, no enrichment row changed category; corrected OFC BA11 MEtan retains FDR<0.10 enrichment for the global high set and global high+moderate set.", "",
        md_table(changed_w, w_cols), "",
        "## Stable core and networks", "",
        "The WGCNA structural core was not recalculated because its selection uses module membership, expression presence, and kME/hub rank—not ML priority. Versioned tables preserve every structural-core row and refresh only `is_high_confidence_ML` and `is_high_plus_moderate_ML`. Base network membership and edges are preserved; corrected node-annotation tables identify changed ML labels. Visual-only network rerenders were generated in the versioned Phase 2 folder. The main-report frontal rerender contains only the original frontal subregions (dlPFC BA8/9, OFC BA11, and CG25); pooled broad Frontal Cortex is excluded from that figure.", "",
        "## Other downstream analyses", "",
        "- Tier2 enrichment: unchanged input; not rerun.",
        "- Tier4 ranked GSEA: unchanged input; not rerun.",
        "- AHBA module-score ML: fixed published modules, independent of the shortlist; preserved.",
        "- MDD module/pathway-level ML2: fixed WGCNA module scores; preserved. Only ML membership annotations changed.",
        "- Preliminary manuscript package and dirty report copies: not overwritten. They remain queued for a separate publication refresh after human review of the changes above.", "",
        "## Interpretation safeguard", "",
        "Crossing a nominal or FDR threshold after changing the candidate list is reported as a sensitivity outcome, not a new biological discovery. The corrected results remain exploratory and require external validation.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")

    inventory_path = PHASE2_ROOT / "phase2_file_inventory.csv"
    manifest_path = PHASE2_ROOT / "phase2_manifest.json"
    inventory_paths = sorted(
        [p for p in PHASE2_ROOT.rglob("*") if p.is_file() and p not in {inventory_path, manifest_path}],
        key=lambda p: str(p).lower(),
    )
    inventory_rows = [{
        "relative_path": str(p.relative_to(PROJECT_ROOT)),
        "size_bytes": p.stat().st_size,
        "sha256": digest(p),
    } for p in inventory_paths]
    inventory_rows.append({
        "relative_path": str(inventory_path.relative_to(PROJECT_ROOT)),
        "size_bytes": "self",
        "sha256": "not_applicable_self_inventory",
    })
    inventory_rows.append({
        "relative_path": str(manifest_path.relative_to(PROJECT_ROOT)),
        "size_bytes": "circular_reference",
        "sha256": "not_applicable_manifest_inventory_cycle",
    })
    pd.DataFrame(inventory_rows).to_csv(inventory_path, index=False)

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    key_files = [
        PHASE2_ROOT / "AHBA_corrected/published_AHBA_ML_enrichment_corrected.csv",
        PHASE2_ROOT / "WGCNA_ML_annotations_corrected/broad_regions_ML_module_enrichment_corrected.csv",
        PHASE2_ROOT / "WGCNA_ML_annotations_corrected/frontal_subregions_ML_module_enrichment_corrected.csv",
        PHASE2_ROOT / "stable_core_annotations_corrected/structural_core_corrected_ML_annotation_summary.csv",
        PHASE2_ROOT / "WGCNA_ML_annotations_corrected/network_figures_broad/full_module_network_summary_top_modules.png",
        PHASE2_ROOT / "WGCNA_ML_annotations_corrected/network_figures_frontal_subregions/frontal_subregion_network_summary_corrected_ML.png",
        PHASE2_ROOT / "WGCNA_ML_annotations_corrected/network_figures_frontal_subregions/frontal_subregion_network_inventory.csv",
        inventory_path,
        DOC,
    ]
    manifest["key_output_sha256"] = {str(p.relative_to(PROJECT_ROOT)): digest(p) for p in key_files}
    manifest["commands"] = [
        "python code/regenerate_corrected_ml_downstream.py",
        "python code/create_mdd_wgcna_full_module_networks.py (with versioned output/enrichment environment variables)",
        "python code/render_corrected_frontal_wgcna_networks.py",
        "python code/summarize_corrected_ml_downstream.py",
    ]
    manifest["environment"] = {"python": "3.12 ACL-safe temporary runtime", "pandas": pd.__version__, "numpy": np.__version__, "scipy": "1.18.1", "matplotlib": matplotlib.__version__}
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(DOC)


if __name__ == "__main__":
    main()
