#!/usr/bin/env python3
"""Build predeclared conditional-null summaries and diagnostic figures."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ConstantInputWarning, spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs/length_gt_100_clean_pipeline/20_Tier2_ML_conditional_null"
PERMS = OUTPUT / "permutations"
SUMMARY = OUTPUT / "summary"
FIGURES = OUTPUT / "figures"
OBSERVED = ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization"
BENCHMARK = ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark"
TIER2 = ROOT / "outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv"
MODELS = ["logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest"]
TOP_LABELS = ["top10", "top25", "top20percent"]
DISPLAY_MODELS = {
    "logistic_l1": "L1 logistic", "elastic_net": "Elastic Net", "logistic_l2": "L2 logistic",
    "linear_svm": "Linear SVM", "random_forest": "Random Forest", "dummy_most_frequent": "Dummy",
}
REGION_ORDER = ["Ventral Striatum", "Insula", "Hippocampal Formation", "CG25", "dlPFC BA8/9", "OFC BA11"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")


def top_n(universe: int, label: str) -> int:
    return {"top10": 10, "top25": 25, "top20percent": int(math.ceil(0.20 * universe))}[label]


def pairwise_stability(rank_maps: list[dict[str, float]], candidates: list[str]) -> dict[str, float]:
    jaccards = {label: [] for label in TOP_LABELS}
    correlations = []
    for left, right in itertools.combinations(rank_maps, 2):
        left_vector = [left.get(gene, len(left) + 1) for gene in candidates]
        right_vector = [right.get(gene, len(right) + 1) for gene in candidates]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConstantInputWarning)
            correlations.append(float(spearmanr(left_vector, right_vector).statistic))
        for label in TOP_LABELS:
            n = top_n(len(candidates), label)
            left_top = set(sorted(left, key=lambda gene: (left[gene], gene))[:n])
            right_top = set(sorted(right, key=lambda gene: (right[gene], gene))[:n])
            union = left_top | right_top
            jaccards[label].append(len(left_top & right_top) / len(union) if union else np.nan)
    result = {
        "spearman_mean": float(np.nanmean(correlations)),
        "spearman_median": float(np.nanmedian(correlations)),
    }
    result.update({f"jaccard_mean_{label}": float(np.nanmean(values)) for label, values in jaccards.items()})
    return result


def aggregate(rows: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    output_rows = []
    for keys, part in rows.groupby(group_columns, sort=False, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        identity = dict(zip(group_columns, keys))
        observed_values = part.loc[part["source"].eq("observed"), "value"]
        null = part.loc[part["source"].eq("permutation"), "value"].to_numpy(dtype=float)
        if len(observed_values) != 1 or len(null) != 20:
            raise RuntimeError(f"Incomplete observed/null grid: {identity}")
        observed = float(observed_values.iloc[0])
        direction = str(part["stronger_direction"].iloc[0])
        extreme = int(np.sum(null >= observed)) if direction == "larger" else int(np.sum(null <= observed))
        output_rows.append({
            **identity, "stronger_direction": direction, "observed_value": observed,
            "null_mean": float(np.mean(null)), "null_median": float(np.median(null)),
            "null_sd": float(np.std(null, ddof=1)), "null_q025": float(np.quantile(null, 0.025)),
            "null_q95": float(np.quantile(null, 0.95)), "null_q975": float(np.quantile(null, 0.975)),
            "n_permutations": len(null), "n_null_as_or_more_extreme": extreme,
            "empirical_p_pilot": (1 + extreme) / (1 + len(null)),
        })
    return pd.DataFrame(output_rows)


def load_candidates() -> tuple[pd.DataFrame, list[str]]:
    candidates = pd.read_csv(TIER2, dtype={"gene_id": str})
    candidates["region"] = candidates["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
    observed_regions = set(candidates["region"])
    if observed_regions != set(REGION_ORDER):
        raise RuntimeError(f"Unexpected Tier2 regions: {sorted(observed_regions)}")
    regions = REGION_ORDER
    return candidates, regions


def rra_stability_values(candidates: pd.DataFrame, regions: list[str], perm_dirs: list[Path]) -> pd.DataFrame:
    rows = []
    sources = [("observed", 0, BENCHMARK / "RRA_resample_consensus_rankings.csv")]
    sources += [("permutation", index, path / "RRA_resample_consensus_rankings.csv") for index, path in enumerate(perm_dirs, start=1)]
    for source, index, path in sources:
        table = pd.read_csv(path, dtype={"gene_id": str})
        for region in regions:
            genes = candidates.loc[candidates["region"].eq(region), "gene_id"].tolist()
            part = table.loc[table["region"].eq(region)]
            maps = [dict(zip(split["gene_id"], split["resample_rank"])) for _, split in part.groupby(["repeat", "fold"], sort=False)]
            metrics = pairwise_stability(maps, genes)
            for metric, value in metrics.items():
                rows.append({"family": "RRA resample stability", "region": region, "metric": metric,
                             "source": source, "permutation_index": index, "value": value, "stronger_direction": "larger"})
    return pd.DataFrame(rows)


def model_agreement_values(candidates: pd.DataFrame, regions: list[str], perm_dirs: list[Path]) -> pd.DataFrame:
    rows = []
    sources = [("observed", 0, OBSERVED / "Tier2_model_level_rra_rankings.csv")]
    sources += [("permutation", index, path / "Tier2_model_level_rra_rankings.csv") for index, path in enumerate(perm_dirs, start=1)]
    for source, index, path in sources:
        table = pd.read_csv(path, dtype={"gene_id": str})
        for region in regions:
            genes = candidates.loc[candidates["region"].eq(region), "gene_id"].tolist()
            model_maps = []
            for model in MODELS:
                part = table.loc[table["region"].eq(region) & table["model"].eq(model)]
                model_maps.append(dict(zip(part["gene_id"], part["model_rank"])))
            correlations = []
            for left, right in itertools.combinations(model_maps, 2):
                lv = [left.get(gene, len(left) + 1) for gene in genes]
                rv = [right.get(gene, len(right) + 1) for gene in genes]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", ConstantInputWarning)
                    correlations.append(float(spearmanr(lv, rv).statistic))
            rows.append({"family": "RRA model agreement", "region": region,
                         "metric": "model_to_model_spearman_mean", "source": source,
                         "permutation_index": index, "value": float(np.nanmean(correlations)),
                         "stronger_direction": "larger"})
    return pd.DataFrame(rows)


def rra_strength_values(candidates: pd.DataFrame, regions: list[str], perm_dirs: list[Path]) -> pd.DataFrame:
    rows = []
    sources = [("observed", 0, OBSERVED / "Tier2_consensus_rra_rankings.csv")]
    sources += [("permutation", index, path / "Tier2_consensus_rra_rankings.csv") for index, path in enumerate(perm_dirs, start=1)]
    for source, index, path in sources:
        table = pd.read_csv(path, dtype={"gene_id": str})
        for region in regions:
            part = table.loc[table["region"].eq(region)].sort_values(["consensus_rank", "gene_id"])
            universe = len(part)
            metrics: list[tuple[str, float, str]] = [("minimum_consensus_rra_score", float(part["consensus_rra_score"].min()), "smaller")]
            for label in TOP_LABELS:
                n = top_n(universe, label)
                metrics.append((f"median_consensus_rra_score_{label}", float(part.head(n)["consensus_rra_score"].median()), "smaller"))
            for threshold in (3, 4, 5):
                metrics.append((f"genes_supported_by_at_least_{threshold}_models", float(part["model_appearance_count"].ge(threshold).sum()), "larger"))
            for metric, value, direction in metrics:
                rows.append({"family": "RRA consensus strength", "region": region, "metric": metric,
                             "source": source, "permutation_index": index, "value": value, "stronger_direction": direction})
    return pd.DataFrame(rows)


def lasso_values(candidates: pd.DataFrame, regions: list[str], perm_dirs: list[Path]) -> pd.DataFrame:
    rows = []
    sources = [("observed", 0, BENCHMARK / "LASSO_gene_rankings.csv", BENCHMARK / "LASSO_fold_rankings.csv")]
    sources += [("permutation", index, path / "LASSO_gene_rankings.csv", path / "Tier2_fold_gene_rankings.csv") for index, path in enumerate(perm_dirs, start=1)]
    for source, index, ranking_path, fold_path in sources:
        ranking = pd.read_csv(ranking_path, dtype={"gene_id": str})
        folds = pd.read_csv(fold_path, dtype={"gene_id": str}, low_memory=False)
        if source == "permutation":
            folds = folds.loc[folds["model"].eq("logistic_l1")].copy()
        for region in regions:
            part = ranking.loc[ranking["region"].eq(region)].sort_values(["lasso_export_rank", "gene_id"])
            metrics: list[tuple[str, float]] = [
                ("maximum_selection_frequency", float(part["selection_frequency"].max())),
                ("mean_selection_frequency_top10", float(part.head(10)["selection_frequency"].mean())),
                ("mean_selection_frequency_top25", float(part.head(25)["selection_frequency"].mean())),
                ("genes_selected_at_least_20_of_25", float(part["selection_count"].ge(20).sum())),
                ("genes_selected_at_least_15_of_25", float(part["selection_count"].ge(15).sum())),
                ("genes_selected_at_least_10_of_25", float(part["selection_count"].ge(10).sum())),
            ]
            genes = candidates.loc[candidates["region"].eq(region), "gene_id"].tolist()
            fold_part = folds.loc[folds["region"].eq(region)]
            maps = [dict(zip(split["gene_id"], split["within_fit_rank"])) for _, split in fold_part.groupby("fit_id", sort=False)]
            stability = pairwise_stability(maps, genes)
            metrics.extend((metric, value) for metric, value in stability.items())
            for metric, value in metrics:
                rows.append({"family": "LASSO", "region": region, "metric": metric,
                             "source": source, "permutation_index": index, "value": value, "stronger_direction": "larger"})
    return pd.DataFrame(rows)


def topgene_behavior(candidates: pd.DataFrame, regions: list[str], perm_dirs: list[Path], method: str) -> pd.DataFrame:
    if method == "RRA":
        observed = pd.read_csv(OBSERVED / "Tier2_consensus_rra_rankings.csv", dtype={"gene_id": str})
        null_tables = [pd.read_csv(path / "Tier2_consensus_rra_rankings.csv", dtype={"gene_id": str}) for path in perm_dirs]
        rank_column, score_column, support_column = "consensus_rank", "consensus_rra_score", "model_appearance_count"
    else:
        observed = pd.read_csv(BENCHMARK / "LASSO_gene_rankings.csv", dtype={"gene_id": str})
        null_tables = [pd.read_csv(path / "LASSO_gene_rankings.csv", dtype={"gene_id": str}) for path in perm_dirs]
        rank_column, score_column, support_column = "lasso_export_rank", "selection_frequency", "selection_count"
    rows = []
    for region in regions:
        observed_part = observed.loc[observed["region"].eq(region)].sort_values([rank_column, "gene_id"]).head(25)
        universe = int(candidates["region"].eq(region).sum())
        n20 = top_n(universe, "top20percent")
        null_region = [table.loc[table["region"].eq(region)].set_index("gene_id") for table in null_tables]
        for row in observed_part.itertuples(index=False):
            gene = str(row.gene_id)
            ranks = np.array([float(table.loc[gene, rank_column]) for table in null_region])
            scores = np.array([float(table.loc[gene, score_column]) for table in null_region])
            observed_rank = float(getattr(row, rank_column))
            observed_score = float(getattr(row, score_column))
            score_extreme = int(np.sum(scores <= observed_score)) if method == "RRA" else int(np.sum(scores >= observed_score))
            rows.append({
                "region": region, "gene_id": gene, "gene_symbol": str(row.gene_symbol),
                "observed_rank": observed_rank, "observed_evidence": observed_score,
                "observed_model_or_fit_support": float(getattr(row, support_column)),
                "null_top10_count": int(np.sum(ranks <= 10)), "null_top25_count": int(np.sum(ranks <= 25)),
                "null_top20percent_count": int(np.sum(ranks <= n20)),
                "null_rank_median": float(np.median(ranks)), "null_rank_best": float(np.min(ranks)),
                "null_rank_worst": float(np.max(ranks)), "null_evidence_mean": float(np.mean(scores)),
                "null_evidence_median": float(np.median(scores)),
                "exploratory_rank_empirical_p_pilot": (1 + int(np.sum(ranks <= observed_rank))) / 21,
                "exploratory_evidence_empirical_p_pilot": (1 + score_extreme) / 21,
                "null_rank_distribution_json": json.dumps(ranks.astype(int).tolist()),
                "null_evidence_distribution_json": json.dumps(scores.tolist()),
            })
    return pd.DataFrame(rows)


def prediction_values(regions: list[str], perm_dirs: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed = pd.read_csv(OBSERVED / "Tier2_model_performance_summary.csv")
    null_tables = [pd.read_csv(path / "Tier2_model_performance_summary.csv") for path in perm_dirs]
    metric_columns = {
        "balanced_accuracy": "balanced_accuracy_mean", "f1": "f1_mean", "roc_auc": "roc_auc_mean",
        "precision": "precision_mean", "recall": "recall_mean",
    }
    rows = []
    best_rows = []
    for region in regions:
        for model in MODELS + ["dummy_most_frequent"]:
            obs_row = observed.loc[observed["region"].eq(region) & observed["model"].eq(model)].iloc[0]
            for metric, column in metric_columns.items():
                null = np.array([float(table.loc[table["region"].eq(region) & table["model"].eq(model), column].iloc[0]) for table in null_tables])
                observed_value = float(obs_row[column])
                rows.append({
                    "region": region, "model": model, "metric": metric, "observed_value": observed_value,
                    "null_mean": float(np.mean(null)), "null_median": float(np.median(null)),
                    "null_sd": float(np.std(null, ddof=1)), "null_q025": float(np.quantile(null, 0.025)),
                    "null_q95": float(np.quantile(null, 0.95)), "null_q975": float(np.quantile(null, 0.975)),
                    "n_permutations": 20, "n_null_as_or_more_extreme": int(np.sum(null >= observed_value)),
                    "empirical_p_pilot": (1 + int(np.sum(null >= observed_value))) / 21,
                })
        observed_best = float(observed.loc[observed["region"].eq(region) & observed["model"].isin(MODELS), "balanced_accuracy_mean"].max())
        null_best = np.array([float(table.loc[table["region"].eq(region) & table["model"].isin(MODELS), "balanced_accuracy_mean"].max()) for table in null_tables])
        best_rows.append({
            "region": region, "observed_best_of_five_balanced_accuracy": observed_best,
            "null_mean": float(np.mean(null_best)), "null_median": float(np.median(null_best)),
            "null_sd": float(np.std(null_best, ddof=1)), "null_q025": float(np.quantile(null_best, 0.025)),
            "null_q95": float(np.quantile(null_best, 0.95)), "null_q975": float(np.quantile(null_best, 0.975)),
            "n_permutations": 20, "n_null_as_or_more_extreme": int(np.sum(null_best >= observed_best)),
            "empirical_p_pilot": (1 + int(np.sum(null_best >= observed_best))) / 21,
            "null_distribution_json": json.dumps(null_best.tolist()),
        })
    return pd.DataFrame(rows), pd.DataFrame(best_rows)


def save_figure(fig: plt.Figure, name: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{name}.png", dpi=300, bbox_inches="tight")
    fixed = datetime(2026, 10, 11, tzinfo=timezone.utc)
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight", metadata={"CreationDate": fixed, "ModDate": fixed})
    plt.close(fig)


def facet_null_plot(values: pd.DataFrame, regions: list[str], family: str, metric: str, title: str, ylabel: str, name: str, log_x: bool = False) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    part = values.loc[values["family"].eq(family) & values["metric"].eq(metric)]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for ax, region in zip(axes.flat, regions):
        region_part = part.loc[part["region"].eq(region)]
        null = region_part.loc[region_part["source"].eq("permutation"), "value"]
        observed = float(region_part.loc[region_part["source"].eq("observed"), "value"].iloc[0])
        # Use deterministic vertical offsets instead of seaborn's random jitter so
        # repeated summary builds produce byte-identical figure/checksum ledgers.
        offsets = np.linspace(-0.12, 0.12, len(null))
        ax.scatter(null.to_numpy(), offsets, color="#4C72B0", s=45)
        ax.set_ylim(-0.25, 0.25)
        ax.set_yticks([])
        ax.axvline(observed, color="#C44E52", linewidth=2.5, label="Observed")
        if log_x:
            ax.set_xscale("log")
        ax.set_title(region); ax.set_xlabel(""); ax.set_ylabel("")
    fig.suptitle(title, y=0.995); fig.supxlabel(ylabel, y=0.02)
    handles = [plt.Line2D([0], [0], color="#C44E52", lw=2.5, label="Observed"),
               plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#4C72B0", markersize=8, label="Permuted labels")]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.88, 0.5))
    fig.tight_layout(rect=(0.03, 0.04, 0.88, 0.96)); save_figure(fig, name)


def build_figures(values: pd.DataFrame, regions: list[str], rra_top: pd.DataFrame, prediction: pd.DataFrame, best: pd.DataFrame) -> None:
    facet_null_plot(values, regions, "RRA resample stability", "spearman_mean",
                    "Observed RRA rank stability versus conditional null", "Mean pairwise Spearman correlation",
                    "observed_RRA_resampling_stability_vs_null")
    facet_null_plot(values, regions, "RRA resample stability", "jaccard_mean_top10",
                    "Observed top-10 RRA overlap versus conditional null", "Mean pairwise top-10 Jaccard index",
                    "observed_topset_Jaccard_vs_null")
    facet_null_plot(values, regions, "RRA consensus strength", "minimum_consensus_rra_score",
                    "Observed minimum consensus RRA score versus conditional null", "Minimum consensus RRA score (smaller is stronger)",
                    "observed_RRA_consensus_strength_vs_null", log_x=True)
    facet_null_plot(values, regions, "LASSO", "maximum_selection_frequency",
                    "Observed LASSO selection concentration versus conditional null", "Maximum selection frequency",
                    "observed_LASSO_selection_stability_vs_null")

    sns.set_theme(style="whitegrid", context="talk")
    selected = rra_top.loc[rra_top["observed_rank"].le(5)].copy()
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    for ax, region in zip(axes.flat, regions):
        part = selected.loc[selected["region"].eq(region)].sort_values("observed_rank")
        y = np.arange(len(part)); med = part["null_rank_median"].to_numpy(); low = part["null_rank_best"].to_numpy(); high = part["null_rank_worst"].to_numpy()
        ax.errorbar(med, y, xerr=[med - low, high - med], fmt="o", color="#4C72B0", label="Null median/range")
        ax.scatter(part["observed_rank"], y, marker="D", color="#C44E52", label="Observed rank", zorder=3)
        ax.set_yticks(y, labels=part["gene_symbol"].where(part["gene_symbol"].ne(""), part["gene_id"])); ax.invert_yaxis()
        ax.set_title(region); ax.set_xlabel(""); ax.set_ylabel("")
    fig.suptitle("Observed leading RRA genes versus permutation rank distributions", y=0.995)
    fig.supxlabel("Rank (lower is stronger; bars show null range)", y=0.02)
    handles, labels = axes.flat[0].get_legend_handles_labels(); fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.88, 0.5))
    fig.tight_layout(rect=(0.03, 0.04, 0.88, 0.96)); save_figure(fig, "observed_topgenes_vs_null_rank_distribution")

    pred = prediction.loc[prediction["metric"].eq("balanced_accuracy") & prediction["model"].isin(MODELS)].copy()
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), sharey=True)
    for ax, region in zip(axes.flat, regions):
        part = pred.loc[pred["region"].eq(region)].copy(); x = np.arange(len(part))
        ax.errorbar(x, part["null_mean"], yerr=[part["null_mean"] - part["null_q025"], part["null_q975"] - part["null_mean"]], fmt="o", color="#4C72B0", label="Null mean / 95% interval")
        ax.scatter(x, part["observed_value"], marker="D", color="#C44E52", label="Observed")
        ax.set_xticks(x, labels=[DISPLAY_MODELS[m] for m in part["model"]], rotation=35, ha="right"); ax.set_title(region); ax.set_xlabel(""); ax.set_ylabel("")
    fig.suptitle("Observed balanced accuracy versus conditional permutation null", y=0.995); fig.supylabel("Mean outer-fold balanced accuracy", x=0.02)
    handles, labels = axes.flat[0].get_legend_handles_labels(); fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.88, 0.5))
    fig.tight_layout(rect=(0.04, 0.04, 0.88, 0.96)); save_figure(fig, "observed_balanced_accuracy_vs_null")

    fig, ax = plt.subplots(figsize=(12, 7)); x = np.arange(len(best))
    ax.errorbar(x, best["null_mean"], yerr=[best["null_mean"] - best["null_q025"], best["null_q975"] - best["null_mean"]], fmt="o", color="#4C72B0", label="Null maximum mean / 95% interval")
    ax.scatter(x, best["observed_best_of_five_balanced_accuracy"], marker="D", s=70, color="#C44E52", label="Observed best of five")
    ax.set_xticks(x, labels=best["region"], rotation=25, ha="right"); ax.set_ylabel("Best-of-five mean balanced accuracy")
    ax.set_title("Best-of-five balanced accuracy versus conditional null"); ax.legend(); fig.tight_layout(); save_figure(fig, "best_of_five_balanced_accuracy_vs_null")


def main() -> int:
    candidates, regions = load_candidates()
    perm_dirs = [PERMS / f"perm_{index:04d}" for index in range(1, 21)]
    if not all((path / "permutation_completion.json").exists() for path in perm_dirs):
        raise RuntimeError("All 20 completed permutations are required")
    rra_stability = rra_stability_values(candidates, regions, perm_dirs)
    model_agreement = model_agreement_values(candidates, regions, perm_dirs)
    rra_strength = rra_strength_values(candidates, regions, perm_dirs)
    lasso = lasso_values(candidates, regions, perm_dirs)
    all_values = pd.concat([rra_stability, model_agreement, rra_strength, lasso], ignore_index=True)
    write_csv(all_values, SUMMARY / "null_metric_values.csv")
    rra_stability_summary = aggregate(pd.concat([rra_stability, model_agreement]), ["family", "region", "metric"])
    rra_strength_summary = aggregate(rra_strength, ["family", "region", "metric"])
    lasso_summary = aggregate(lasso, ["family", "region", "metric"])
    write_csv(rra_stability_summary, SUMMARY / "RRA_null_stability_summary.csv")
    write_csv(rra_strength_summary, SUMMARY / "RRA_null_consensus_strength.csv")
    write_csv(lasso_summary, SUMMARY / "LASSO_null_stability_summary.csv")
    rra_top = topgene_behavior(candidates, regions, perm_dirs, "RRA")
    lasso_top = topgene_behavior(candidates, regions, perm_dirs, "LASSO")
    write_csv(rra_top, SUMMARY / "RRA_observed_topgene_null_behavior.csv")
    write_csv(lasso_top, SUMMARY / "LASSO_observed_topgene_null_behavior.csv")
    prediction, best = prediction_values(regions, perm_dirs)
    write_csv(prediction, SUMMARY / "prediction_null_summary.csv")
    write_csv(best, SUMMARY / "best_of_five_prediction_null.csv")
    combined = pd.concat([rra_stability_summary, rra_strength_summary, lasso_summary], ignore_index=True)
    write_csv(combined, SUMMARY / "observed_vs_null_summary.csv")
    build_figures(all_values, regions, rra_top, prediction, best)

    runtimes = [json.loads((path / "ml_completion.json").read_text(encoding="utf-8"))["wall_runtime_seconds"] for path in perm_dirs]
    manifest = {
        "status": "20_permutation_pilot_complete", "permutation_indices": [1, 20],
        "permutation_seeds": [41001, 41020], "n_successful_permutations": 20, "n_failed_permutations": 0,
        "candidate_pairs": 671, "conditional_null": True,
        "minimum_empirical_pilot_probability": 1 / 21,
        "runtime_seconds_ml": {"mean": float(np.mean(runtimes)), "median": float(np.median(runtimes)),
                               "minimum": float(np.min(runtimes)), "maximum": float(np.max(runtimes)),
                               "sum_permutation_wall_times": float(np.sum(runtimes))},
        "frozen_observed_rra_sha256": sha256(OBSERVED / "Tier2_consensus_rra_rankings.csv"),
        "frozen_benchmark_checksum_ledger_sha256": sha256(BENCHMARK / "checksums.sha256"),
        "interpretation_warning": "Tier2 was selected using observed Diagnosis; this conditional null does not test the complete discovery pipeline.",
    }
    (OUTPUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print("Built conditional-null summaries and seven diagnostic figure pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
