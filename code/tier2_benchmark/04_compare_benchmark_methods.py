#!/usr/bin/env python3
"""Compare frozen Tier2-RRA, LASSO alone, and Stability Selection transparently."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kendalltau, spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark"
PRIMARY = ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization"
ALLGENE = ROOT / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/consensus_rra_rankings.csv"
METHODS = ("Tier2-RRA", "LASSO", "Stability Selection")
EXPECTED_PRIMARY_SHA = "2e1ed5fc91aae752f57e46e61022a9077d7f1ab285cfc05a6973c3c6df23d358"
EXPECTED_ALLGENE_SHA = "0aefec5fcfa100deafa1525c2d31521ae8c00ca890ef4e1b4145b1b861c2c87e"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")


def save(fig: plt.Figure, stem: str) -> None:
    figures = OUTPUT / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures / f"{stem}.png", dpi=300, bbox_inches="tight")
    fixed_time = datetime(2026, 10, 10, tzinfo=timezone.utc)
    fig.savefig(
        figures / f"{stem}.pdf", bbox_inches="tight",
        metadata={"CreationDate": fixed_time, "ModDate": fixed_time},
    )
    plt.close(fig)


def pair_metrics(left: set[str], right: set[str]) -> tuple[int, float, float, float]:
    intersection = len(left & right)
    union = len(left | right)
    return (
        intersection,
        intersection / union if union else np.nan,
        intersection / len(left) if left else np.nan,
        intersection / len(right) if right else np.nan,
    )


def describe(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "mean": float(np.nanmean(array)), "median": float(np.nanmedian(array)),
        "sd": float(np.nanstd(array, ddof=1)), "minimum": float(np.nanmin(array)),
        "maximum": float(np.nanmax(array)), "n_pairs": int(np.isfinite(array).sum()),
    }


def method_tables() -> dict[str, pd.DataFrame]:
    rra = pd.read_csv(PRIMARY / "Tier2_consensus_rra_rankings.csv", dtype={"gene_id": str})
    rra = rra.assign(priority_rank=rra["consensus_rank"], correlation_rank=rra["consensus_rank"], support=True)
    lasso = pd.read_csv(OUTPUT / "LASSO_gene_rankings.csv", dtype={"gene_id": str})
    lasso = lasso.assign(priority_rank=lasso["lasso_export_rank"], correlation_rank=lasso["lasso_correlation_rank"], support=lasso["selection_count"].gt(0))
    stability = pd.read_csv(OUTPUT / "stability_selection_gene_rankings.csv", dtype={"gene_id": str})
    stability = stability.assign(priority_rank=stability["stability_export_rank"], correlation_rank=stability["stability_correlation_rank"], support=stability["selection_count"].gt(0))
    return {"Tier2-RRA": rra, "LASSO": lasso, "Stability Selection": stability}


def top_definition(n_universe: int, label: str) -> int:
    return {"top10": 10, "top25": 25, "top20percent": int(math.ceil(0.20 * n_universe))}[label]


def final_comparisons(tables: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    overlap_rows, correlation_rows = [], []
    regions = tables["Tier2-RRA"]["region"].drop_duplicates().tolist()
    for region in regions:
        region_tables = {method: table.loc[table["region"].eq(region)].copy() for method, table in tables.items()}
        universe = len(region_tables["Tier2-RRA"])
        for left_name, right_name in itertools.combinations(METHODS, 2):
            left = region_tables[left_name].set_index("gene_id")
            right = region_tables[right_name].set_index("gene_id")
            common = left.index.intersection(right.index)
            correlation_rows.append({
                "region": region, "method_left": left_name, "method_right": right_name, "n_common": len(common),
                "spearman_rho": float(spearmanr(left.loc[common, "correlation_rank"], right.loc[common, "correlation_rank"]).statistic),
                "kendall_tau": float(kendalltau(left.loc[common, "correlation_rank"], right.loc[common, "correlation_rank"]).statistic),
                "never_selected_rule": "shared worst-rank block",
            })
            for label in ("top10", "top25", "top20percent"):
                n_top = top_definition(universe, label)
                left_top = set(left.loc[left["support"]].nsmallest(n_top, "priority_rank").index)
                right_top = set(right.loc[right["support"]].nsmallest(n_top, "priority_rank").index)
                count, jaccard, frac_left, frac_right = pair_metrics(left_top, right_top)
                overlap_rows.append({
                    "region": region, "method_left": left_name, "method_right": right_name,
                    "top_definition": label, "requested_n": n_top, "left_list_size": len(left_top),
                    "right_list_size": len(right_top), "intersection_count": count, "jaccard": jaccard,
                    "overlap_fraction_left": frac_left, "overlap_fraction_right": frac_right,
                })
    return pd.DataFrame(overlap_rows), pd.DataFrame(correlation_rows)


def resample_maps() -> dict[str, dict[str, list[dict[str, float]]]]:
    maps: dict[str, dict[str, list[dict[str, float]]]] = {method: {} for method in METHODS}
    rra = pd.read_csv(OUTPUT / "RRA_resample_consensus_rankings.csv", dtype={"gene_id": str})
    for region, region_part in rra.groupby("region", sort=False):
        maps["Tier2-RRA"][region] = [
            dict(zip(part["gene_id"], part["resample_rank"]))
            for _, part in region_part.groupby(["repeat", "fold"], sort=False)
        ]
    lasso = pd.read_csv(OUTPUT / "LASSO_fold_rankings.csv", dtype={"gene_id": str})
    for region, region_part in lasso.groupby("region", sort=False):
        maps["LASSO"][region] = [
            dict(zip(part["gene_id"], part["within_fit_rank"]))
            for _, part in region_part.groupby("fit_id", sort=False)
        ]
    stability = pd.read_csv(OUTPUT / "stability_selection_run_details.csv", dtype={"gene_id": str}, low_memory=False)
    for region, region_part in stability.groupby("region", sort=False):
        region_maps = []
        for _, part in region_part.groupby("subsample_id", sort=False):
            selected = part.loc[part["selected"].astype(str).str.lower().eq("true")]
            region_maps.append(dict(zip(selected["gene_id"], selected["within_subsample_rank"])))
        maps["Stability Selection"][region] = region_maps
    return maps


def stability_summary(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    maps = resample_maps()
    rows = []
    for method in METHODS:
        for region, ranking_maps in maps[method].items():
            candidates = tables[method].loc[tables[method]["region"].eq(region), "gene_id"].tolist()
            n_universe = len(candidates)
            correlations_s, correlations_k = [], []
            jaccards = {label: [] for label in ("top10", "top25", "top20percent")}
            for left, right in itertools.combinations(ranking_maps, 2):
                left_vector = [left.get(gene, len(left) + 1) for gene in candidates]
                right_vector = [right.get(gene, len(right) + 1) for gene in candidates]
                correlations_s.append(float(spearmanr(left_vector, right_vector).statistic))
                correlations_k.append(float(kendalltau(left_vector, right_vector).statistic))
                for label in jaccards:
                    n_top = top_definition(n_universe, label)
                    left_top = set(sorted(left, key=lambda gene: (left[gene], gene))[:n_top])
                    right_top = set(sorted(right, key=lambda gene: (right[gene], gene))[:n_top])
                    jaccards[label].append(pair_metrics(left_top, right_top)[1])
            for metric, values in [("spearman_rank", correlations_s), ("kendall_rank", correlations_k)]:
                rows.append({"region": region, "method": method, "stability_metric": metric, "top_definition": "all_candidates", **describe(values)})
            for label, values in jaccards.items():
                rows.append({"region": region, "method": method, "stability_metric": "jaccard", "top_definition": label, **describe(values)})
    return pd.DataFrame(rows)


def coverage_summary(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for method, table in tables.items():
        for region, part in table.groupby("region", sort=False):
            if method == "Tier2-RRA":
                support_values = np.ones(len(part))
            elif method == "LASSO":
                support_values = part["selection_frequency"].to_numpy()
            else:
                support_values = part["selection_probability"].to_numpy()
            rows.append({
                "region": region, "method": method, "tier2_universe_size": len(part),
                "genes_with_any_support": int((support_values > 0).sum()),
                "proportion_with_any_support": float((support_values > 0).mean()),
                "support_mean": float(np.mean(support_values)), "support_median": float(np.median(support_values)),
                "support_minimum": float(np.min(support_values)), "support_maximum": float(np.max(support_values)),
            })
    return pd.DataFrame(rows)


def allgene_comparison(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    allgene = pd.read_csv(ALLGENE, dtype={"gene_id": str})
    rows = []
    for method, table in tables.items():
        for region, part in table.groupby("region", sort=False):
            universe = len(part)
            external = allgene.loc[allgene["region"].eq(region)].copy()
            external_ranks = external[["gene_id", "consensus_rank"]].rename(columns={"consensus_rank": "allgene_rank"})
            merged = part.merge(external_ranks, on="gene_id", how="inner")
            base = {
                "region": region, "method": method, "n_tier2_candidates": universe,
                "n_present_in_allgene_rra": len(merged),
                "spearman_rho_common": float(spearmanr(merged["correlation_rank"], merged["allgene_rank"]).statistic) if len(merged) >= 3 else np.nan,
            }
            for label in ("top10", "top25", "top20percent"):
                n_top = top_definition(universe, label)
                method_top = set(part.loc[part["support"]].nsmallest(n_top, "priority_rank")["gene_id"])
                allgene_top = set(external.nsmallest(n_top, "consensus_rank")["gene_id"])
                count, jaccard, frac_method, frac_allgene = pair_metrics(method_top, allgene_top)
                base[f"{label}_intersection"] = count
                base[f"{label}_jaccard"] = jaccard
                base[f"{label}_overlap_fraction_method"] = frac_method
                base[f"{label}_overlap_fraction_allgene"] = frac_allgene
            rows.append(base)
    return pd.DataFrame(rows)


def runtime_summary() -> pd.DataFrame:
    primary_manifest = json.loads((PRIMARY / "Tier2_run_manifest.json").read_text(encoding="utf-8"))
    lasso = pd.read_csv(OUTPUT / "LASSO_performance_summary.csv")
    stability = pd.read_csv(OUTPUT / "stability_selection_runtime.csv")
    return pd.DataFrame([
        {"method": "Tier2-RRA", "final_fits": 750, "inner_evaluation_fits": 6300,
         "runtime_seconds": primary_manifest["cumulative_fit_elapsed_seconds"], "models": 5,
         "tuning_burden": "five model families; 14 parameter settings total per split",
         "interpretability": "consensus ranking; model-level ranks retained"},
        {"method": "LASSO", "final_fits": 150, "inner_evaluation_fits": 1350,
         "runtime_seconds": float(lasso["runtime_seconds"].sum()), "models": 1,
         "tuning_burden": "three C values in 3-fold CV per outer split",
         "interpretability": "selection frequency with transparent tie-breakers; exact frozen events reused"},
        {"method": "Stability Selection", "final_fits": 600, "inner_evaluation_fits": 5400,
         "runtime_seconds": float(stability["runtime_seconds"].sum()), "models": 1,
         "tuning_burden": "three C values in 3-fold CV per half-sample",
         "interpretability": "continuous selection probability"},
    ])


def figures(tables: dict[str, pd.DataFrame], overlap: pd.DataFrame, stability: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    region_order = tables["Tier2-RRA"]["region"].drop_duplicates().tolist()
    for other, stem, title in [
        ("LASSO", "RRA_vs_LASSO_rank_comparison", "Tier2-RRA versus LASSO rank"),
        ("Stability Selection", "RRA_vs_stability_selection_rank_comparison", "Tier2-RRA versus Stability Selection rank"),
    ]:
        merged = tables["Tier2-RRA"][["region", "gene_id", "priority_rank"]].merge(
            tables[other][["region", "gene_id", "correlation_rank"]], on=["region", "gene_id"]
        )
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        for ax, region in zip(axes.flat, region_order):
            part = merged.loc[merged["region"].eq(region)]
            sns.scatterplot(data=part, x="priority_rank", y="correlation_rank", s=34, alpha=0.65, ax=ax)
            ax.set_title(region)
            ax.set_xlabel("")
            ax.set_ylabel("")
        fig.suptitle(title, y=0.995)
        fig.supxlabel("Tier2-RRA rank", y=0.02)
        fig.supylabel(f"{other} rank (unselected genes share worst rank)", x=0.015)
        fig.tight_layout(rect=(0.04, 0.04, 1, 0.96))
        save(fig, stem)
    plot_overlap = overlap.copy()
    plot_overlap["comparison"] = plot_overlap["method_left"] + " vs " + plot_overlap["method_right"]
    label_order = ["top10", "top25", "top20percent"]
    display_labels = ["Top 10", "Top 25", "Top 20%"]
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), sharey=True)
    handles = labels = None
    for ax, region in zip(axes.flat, region_order):
        part = plot_overlap.loc[plot_overlap["region"].eq(region)]
        sns.barplot(data=part, x="top_definition", y="jaccard", hue="comparison", order=label_order, ax=ax)
        ax.set_title(region)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks(range(len(display_labels)), labels=display_labels)
        handles, labels = ax.get_legend_handles_labels()
        if ax.legend_ is not None:
            ax.legend_.remove()
    fig.suptitle("Pairwise top-set agreement", y=0.995)
    fig.supxlabel("Top-set definition", y=0.02)
    fig.supylabel("Jaccard index", x=0.02)
    fig.legend(handles, labels, title="Comparison", loc="center left", bbox_to_anchor=(0.86, 0.5), fontsize=11)
    fig.tight_layout(rect=(0.04, 0.04, 0.86, 0.96))
    save(fig, "topset_Jaccard_by_region")
    plot_stability = stability.loc[stability["stability_metric"].eq("jaccard")].copy()
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), sharey=True)
    for ax, region in zip(axes.flat, region_order):
        part = plot_stability.loc[plot_stability["region"].eq(region)]
        sns.barplot(data=part, x="top_definition", y="mean", hue="method", order=label_order, ax=ax)
        ax.set_title(region)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks(range(len(display_labels)), labels=display_labels)
        handles, labels = ax.get_legend_handles_labels()
        if ax.legend_ is not None:
            ax.legend_.remove()
    fig.suptitle("Resampling stability by prioritisation method", y=0.995)
    fig.supxlabel("Top-set definition", y=0.02)
    fig.supylabel("Mean pairwise Jaccard index", x=0.02)
    fig.legend(handles, labels, title="Method", loc="center left", bbox_to_anchor=(0.86, 0.5), fontsize=11)
    fig.tight_layout(rect=(0.04, 0.04, 0.86, 0.96))
    save(fig, "resampling_stability_by_method")
    support_rows = []
    for method in ("LASSO", "Stability Selection"):
        column = "selection_frequency" if method == "LASSO" else "selection_probability"
        support_rows.append(tables[method][["region", column]].rename(columns={column: "selection_support"}).assign(method=method))
    support = pd.concat(support_rows, ignore_index=True)
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), sharex=True)
    for ax, region in zip(axes.flat, region_order):
        part = support.loc[support["region"].eq(region)]
        sns.histplot(data=part, x="selection_support", hue="method", bins=11, element="step", stat="density", common_norm=False, ax=ax)
        ax.set_title(region)
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.suptitle("Sparse-method selection support distributions", y=0.995)
    fig.supxlabel("Selection frequency / probability", y=0.02)
    fig.supylabel("Density", x=0.02)
    fig.tight_layout(rect=(0.04, 0.04, 1, 0.96))
    save(fig, "selection_support_distributions")


def main() -> int:
    if sha256(PRIMARY / "Tier2_consensus_rra_rankings.csv") != EXPECTED_PRIMARY_SHA:
        raise RuntimeError("Frozen primary Tier2-RRA output changed")
    if sha256(ALLGENE) != EXPECTED_ALLGENE_SHA:
        raise RuntimeError("Frozen all-gene RRA output changed")
    tables = method_tables()
    overlap, correlations = final_comparisons(tables)
    stability = stability_summary(tables)
    coverage = coverage_summary(tables)
    allgene = allgene_comparison(tables)
    runtime = runtime_summary()
    write_csv(overlap, OUTPUT / "method_pairwise_topset_overlap.csv")
    write_csv(correlations, OUTPUT / "method_rank_correlations.csv")
    write_csv(stability, OUTPUT / "method_stability_summary.csv")
    write_csv(coverage, OUTPUT / "method_coverage_summary.csv")
    write_csv(runtime, OUTPUT / "method_runtime_summary.csv")
    write_csv(allgene, OUTPUT / "benchmark_vs_allgene_RRA.csv")
    figures(tables, overlap, stability)
    print(f"Compared {len(correlations)} regional method pairs and {len(overlap)} top-set contrasts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
