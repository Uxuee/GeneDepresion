"""Test enrichment of ML-prioritized MDD candidate genes in AHBA modules."""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CLEAN_DIR = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
STEP10_DIR = CLEAN_DIR / "10_ML_prioritized_Tier2_extended"
TIER2_DIR = CLEAN_DIR / "07_Tier2_candidate_gene_sets"
AHBA_DIR = CLEAN_DIR / "ahba_modules"
MODULE_DIR = AHBA_DIR / "modules"
OUT_DIR = AHBA_DIR / "enrichment"
FIG_DIR = AHBA_DIR / "figures"


def norm_symbol(x: object) -> str:
    return "" if pd.isna(x) else str(x).strip().upper()


def bh_fdr(pvals: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvals, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    n = len(pv)
    if n == 0:
        return pd.Series(out, index=pvals.index)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    tmp = np.empty(n)
    tmp[order] = np.minimum(adj, 1)
    out[valid] = tmp
    return pd.Series(out, index=pvals.index)


def load_candidates() -> tuple[pd.DataFrame, list[dict]]:
    ranked = pd.read_csv(STEP10_DIR / "final_ranked_ML_prioritized_Tier2_genes.csv", dtype={"gene_id": str})
    all_tier2 = pd.read_csv(TIER2_DIR / "Tier2_candidate_genes_all_regions.csv", dtype={"gene_id": str})
    ranked["gene_symbol_norm"] = ranked["gene_symbol"].map(norm_symbol)
    all_tier2["gene_symbol_norm"] = all_tier2["gene_symbol"].map(norm_symbol)

    high = ranked[ranked["ML_priority_class"].eq("High-confidence ML-prioritized Tier 2")]
    moderate = ranked[ranked["ML_priority_class"].eq("Moderate-confidence ML-prioritized Tier 2")]
    high_moderate = ranked[ranked["ML_priority_class"].isin([
        "High-confidence ML-prioritized Tier 2",
        "Moderate-confidence ML-prioritized Tier 2",
    ])]

    sets: list[dict] = [
        {"candidate_set": "high_confidence", "scope": "all_regions", "genes": set(high["gene_symbol_norm"])},
        {"candidate_set": "high_plus_moderate", "scope": "all_regions", "genes": set(high_moderate["gene_symbol_norm"])},
        {"candidate_set": "all_Tier2", "scope": "all_regions", "genes": set(all_tier2["gene_symbol_norm"])},
    ]
    for region, sub in ranked.groupby("region_or_subregion"):
        sets.append(
            {
                "candidate_set": "high_plus_moderate",
                "scope": f"region:{region}",
                "genes": set(sub[sub["ML_priority_class"].isin([
                    "High-confidence ML-prioritized Tier 2",
                    "Moderate-confidence ML-prioritized Tier 2",
                ])]["gene_symbol_norm"]),
            }
        )
        for direction, dsub in sub.groupby("statistical_direction"):
            genes = set(dsub[dsub["ML_priority_class"].isin([
                "High-confidence ML-prioritized Tier 2",
                "Moderate-confidence ML-prioritized Tier 2",
            ])]["gene_symbol_norm"])
            if len(genes) >= 5:
                sets.append(
                    {
                        "candidate_set": "high_plus_moderate",
                        "scope": f"region:{region}|direction:{direction}",
                        "genes": genes,
                    }
                )
    return ranked, sets


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    module_path = MODULE_DIR / "ahba_gene_to_module_assignments.csv"
    if not module_path.exists():
        msg = f"AHBA module assignment table not found: {module_path}. Run 02_build_ahba_modules.R first."
        (OUT_DIR / "AHBA_MODULE_ENRICHMENT_NOT_RUN.txt").write_text(msg + "\n", encoding="utf-8")
        print(msg)
        return

    modules = pd.read_csv(module_path)
    modules["gene_symbol_norm"] = modules["gene_symbol"].map(norm_symbol)
    modules = modules[modules["gene_symbol_norm"].ne("")]

    ranked, candidate_sets = load_candidates()
    mdd_universe = set(ranked["gene_symbol"].map(norm_symbol))
    ahba_universe = set(modules["gene_symbol_norm"])
    background = sorted(mdd_universe & ahba_universe)
    bg_set = set(background)
    module_groups = {m: set(g["gene_symbol_norm"]) & bg_set for m, g in modules.groupby("ahba_module")}

    rows = []
    M = len(bg_set)
    for cset in candidate_sets:
        genes = set(cset["genes"]) & bg_set
        N = len(genes)
        if M == 0 or N == 0:
            continue
        for module, module_genes in module_groups.items():
            K = len(module_genes)
            x = len(genes & module_genes)
            p = hypergeom.sf(x - 1, M, K, N) if x > 0 else 1.0
            rows.append(
                {
                    "candidate_set": cset["candidate_set"],
                    "scope": cset["scope"],
                    "ahba_module": module,
                    "background_genes": M,
                    "candidate_genes_in_background": N,
                    "module_genes_in_background": K,
                    "overlap_genes": x,
                    "overlap_gene_symbols": ";".join(sorted(genes & module_genes)),
                    "p_value_hypergeometric": p,
                    "odds_ratio_approx": ((x + 0.5) / (N - x + 0.5)) / ((K - x + 0.5) / (M - N - K + x + 0.5)),
                }
            )
    res = pd.DataFrame(rows)
    if res.empty:
        print("No enrichment rows generated after AHBA/MDD background intersection.")
        return
    res["fdr_BH"] = res.groupby(["candidate_set", "scope"], group_keys=False)["p_value_hypergeometric"].apply(bh_fdr)
    res = res.sort_values(["fdr_BH", "p_value_hypergeometric", "candidate_set", "scope"])
    res.to_csv(OUT_DIR / "ahba_module_enrichment_results.csv", index=False)
    pd.DataFrame({"background_gene_symbol": background}).to_csv(OUT_DIR / "ahba_mdd_background_genes.csv", index=False)

    counts = modules[modules["gene_symbol_norm"].isin(bg_set)].groupby("ahba_module")["gene_symbol_norm"].nunique().reset_index(name="background_genes")
    counts.to_csv(OUT_DIR / "ahba_module_background_gene_counts.csv", index=False)

    top = res.head(30).copy()
    top["minus_log10_fdr"] = -np.log10(top["fdr_BH"].clip(lower=np.nextafter(0, 1)))
    fig, ax = plt.subplots(figsize=(10, max(5, 0.28 * len(top))))
    labels = top["candidate_set"] + " | " + top["scope"] + " | " + top["ahba_module"].astype(str)
    ax.scatter(top["minus_log10_fdr"], range(len(top)), s=25 + 25 * top["overlap_genes"], c=top["odds_ratio_approx"], cmap="viridis")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("-log10(BH-FDR)")
    ax.set_title("AHBA module enrichment of MDD candidate gene sets")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "ahba_module_enrichment_dotplot.png", dpi=220)
    plt.close(fig)

    region_rows = res[res["scope"].str.startswith("region:") & ~res["scope"].str.contains("direction:", regex=False)]
    if not region_rows.empty:
        heat = region_rows.pivot_table(index="scope", columns="ahba_module", values="overlap_genes", fill_value=0)
        fig, ax = plt.subplots(figsize=(max(8, 0.35 * heat.shape[1]), max(4, 0.35 * heat.shape[0])))
        im = ax.imshow(heat.to_numpy(), aspect="auto", cmap="magma")
        fig.colorbar(im, ax=ax, label="Candidate gene overlap")
        ax.set_xticks(np.arange(heat.shape[1]))
        ax.set_xticklabels(heat.columns, rotation=45, ha="right")
        ax.set_yticks(np.arange(heat.shape[0]))
        ax.set_yticklabels([x.replace("region:", "") for x in heat.index])
        ax.set_title("Region-by-AHBA-module candidate gene overlap")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "region_by_ahba_module_overlap_heatmap.png", dpi=220)
        plt.close(fig)

    report = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "module_assignment_file": str(module_path),
        "background_definition": "genes present in both final MDD Tier2 tested universe and processed AHBA module universe",
        "background_gene_count": M,
        "n_candidate_sets_tested": len(candidate_sets),
        "n_enrichment_rows": len(res),
        "results_file": str(OUT_DIR / "ahba_module_enrichment_results.csv"),
        "interpretation": "Healthy-brain AHBA coexpression context only; not MDD biomarker validation.",
    }
    (OUT_DIR / "ahba_module_enrichment_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
