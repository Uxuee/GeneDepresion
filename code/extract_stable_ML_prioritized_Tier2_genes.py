"""Extract stable ML-prioritized Tier 2 genes for Length > 100 MDD analysis.

This is a post-processing analysis only. It does not rerun the LMM, enrichment,
ML v1, or ML v2. It reads existing Tier 2 candidate gene lists and existing ML
feature-importance/stability outputs, then ranks Tier 2 genes by ML support.

Output:
  outputs/length_gt_100_elisa_style/ML_prioritized_Tier2_genes/
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
BASE_DIR = REPO_DIR / "outputs" / "length_gt_100_elisa_style"
OUT_DIR = BASE_DIR / "ML_prioritized_Tier2_genes"
PLOT_DIR = OUT_DIR / "plots"
GENE_LIST_DIR = OUT_DIR / "gene_lists"

TIER2_DIR = BASE_DIR / "DE_candidate_gene_sets" / "Tier2_working_exploratory"
V1_IMPORTANCE_DIR = BASE_DIR / "ML_pilot_Tier2" / "feature_importance"
V1_PERFORMANCE = BASE_DIR / "ML_pilot_Tier2" / "ML_performance_summary.csv"
V1_BEST = BASE_DIR / "ML_pilot_Tier2" / "best_model_by_region.csv"
V2_STABILITY_DIR = BASE_DIR / "ML_strict_v2" / "feature_stability"
V2_PERFORMANCE = BASE_DIR / "ML_strict_v2" / "ML_strict_v2_performance_summary.csv"
V2_BEST = BASE_DIR / "ML_strict_v2" / "best_strict_model_by_region.csv"

REGIONS = [
    ("four_group", "Ventral Striatum", "Ventral_Striatum"),
    ("four_group", "Insula", "Insula"),
    ("four_group", "Hippocampal Formation", "Hippocampal_Formation"),
    ("frontal_subregion", "Cingulate gyrus 25", "Cingulate_gyrus_25"),
    ("frontal_subregion", "dlPFC BA8/9", "dlPFC_BA8_9"),
    ("frontal_subregion", "OFC BA11", "OFC_BA11"),
]

MODEL_FAMILIES = [
    "logistic_l1",
    "elastic_net_logistic",
    "logistic_l2",
    "linear_svm",
    "random_forest_exploratory",
]
SPARSE_MODELS = {"logistic_l1", "elastic_net_logistic"}
LINEAR_MODELS = {"logistic_l1", "elastic_net_logistic", "logistic_l2", "linear_svm"}
TOP_MODELS = {"logistic_l2", "linear_svm", "random_forest_exploratory"}
RF_MODELS = {"random_forest_exploratory"}
SIGNED_MODELS = LINEAR_MODELS
BASELINE_MODELS = {"majority_class_baseline"}

SCORE_WEIGHTS = {
    "selection_frequency_combined": 0.40,
    "model_support_scaled": 0.15,
    "sparse_model_support_binary": 0.10,
    "importance_rank_score": 0.15,
    "sign_consistency_score": 0.10,
    "lmm_evidence_scaled": 0.05,
    "v2_support_score": 0.05,
}


def safe_name(value: str) -> str:
    text = "".join(ch if ch.isalnum() else "_" for ch in str(value))
    while "__" in text:
        text = text.replace("__", "_")
    return text.strip("_")


def reset_output() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for folder in [PLOT_DIR, GENE_LIST_DIR]:
        folder.mkdir(parents=True, exist_ok=True)
        for path in folder.glob("*"):
            if path.is_file():
                path.unlink()
    for path in OUT_DIR.glob("*"):
        if path.is_file():
            path.unlink()


def load_optional_csv(path: Path, used: list[str], missing: list[str], **kwargs) -> pd.DataFrame | None:
    if path.exists():
        print(f"Using: {path}")
        used.append(str(path))
        return pd.read_csv(path, **kwargs)
    print(f"Missing optional file: {path}")
    missing.append(str(path))
    return None


def load_tier2_tables(used: list[str], missing: list[str]) -> pd.DataFrame:
    rows = []
    for _, _, key in REGIONS:
        path = TIER2_DIR / f"{key}__Tier2_all.csv"
        if path.exists():
            print(f"Using Tier 2 list: {path}")
            used.append(str(path))
            df = pd.read_csv(path, dtype={"gene_id": str})
            df["gene_id"] = df["gene_id"].astype(str)
            rows.append(df)
        else:
            print(f"Missing Tier 2 list: {path}")
            missing.append(str(path))
    if not rows:
        raise FileNotFoundError("No Tier 2 gene lists were found.")
    tier2 = pd.concat(rows, ignore_index=True)
    tier2["gene_symbol"] = tier2["gene_symbol"].fillna("").astype(str)
    tier2["beta_diagnosis"] = pd.to_numeric(tier2["beta_diagnosis"], errors="coerce")
    tier2["p_value_diagnosis"] = pd.to_numeric(tier2["p_value_diagnosis"], errors="coerce")
    tier2["standardized_effect_size_diagnosis"] = pd.to_numeric(
        tier2["standardized_effect_size_diagnosis"], errors="coerce"
    )
    tier2["lmm_direction"] = np.select(
        [tier2["beta_diagnosis"] > 0, tier2["beta_diagnosis"] < 0],
        ["up_in_MDD", "down_in_MDD"],
        default="no_change",
    )
    tier2["lmm_evidence_score"] = (
        -np.log10(tier2["p_value_diagnosis"].clip(lower=1e-300))
        * tier2["standardized_effect_size_diagnosis"].abs()
    )
    return tier2


def annotate_importance_events(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Add selected/rank columns to fold-level v1/v2 importance events."""
    if df is None or df.empty:
        return pd.DataFrame()
    df = df.copy()
    # A majority-class baseline does not use gene features. Its feature rows are
    # bookkeeping from the upstream selector, not gene-level model evidence.
    if source == "v2" and "model" in df.columns:
        df = df[~df["model"].isin(BASELINE_MODELS)].copy()
    if df.empty:
        return pd.DataFrame()
    df["gene_id"] = df["gene_id"].astype(str)
    if "importance_value" not in df.columns and "coefficient" in df.columns:
        df["importance_value"] = pd.to_numeric(df["coefficient"], errors="coerce")
    df["importance_value"] = pd.to_numeric(df["importance_value"], errors="coerce")
    df["abs_importance"] = df["importance_value"].abs()
    # Only genuinely signed linear coefficients may inform direction. Random-
    # forest importance remains available for selection, support, and ranking.
    df["signed_importance_value"] = df["importance_value"].where(df["model"].isin(SIGNED_MODELS))
    if "feature_set" not in df.columns:
        df["feature_set"] = (
            df.get("feature_selection_method", "strict_v2").astype(str)
            + "_top"
            + df.get("top_k", "").astype(str)
        )
    combo_cols = [c for c in ["analysis_level", "region_or_subregion", "feature_set", "model", "repeat", "fold", "outer_repeat", "outer_fold"] if c in df.columns]
    df["rank_within_event"] = df.groupby(combo_cols)["abs_importance"].rank(method="min", ascending=False, na_option="bottom")
    df["event_feature_count"] = df.groupby(combo_cols)["gene_id"].transform("count")
    df["top20pct"] = df["rank_within_event"] <= np.ceil(df["event_feature_count"] * 0.20)
    df["top25"] = df["rank_within_event"] <= 25
    df["top50"] = df["rank_within_event"] <= 50
    if "selected" in df.columns:
        base_selected = df["selected"].fillna(False).astype(bool)
    elif "nonzero" in df.columns:
        base_selected = df["nonzero"].fillna(False).astype(bool)
    else:
        base_selected = df["abs_importance"] > 1e-10

    sparse_selected = df["model"].isin(SPARSE_MODELS) & (df["abs_importance"] > 1e-10)
    ranked_selected = df["model"].isin(TOP_MODELS) & df["top20pct"]
    other_selected = ~df["model"].isin(SPARSE_MODELS | TOP_MODELS) & base_selected
    df["selected_conservative"] = sparse_selected | ranked_selected | other_selected
    df["source"] = source
    return df


def summarize_events(events: pd.DataFrame, tier2: pd.DataFrame, source: str) -> pd.DataFrame:
    keys = tier2[["analysis_level", "region_or_subregion", "gene_id"]].drop_duplicates()
    if events is None or events.empty:
        out = keys.copy()
        out[f"selection_frequency_{source}"] = np.nan
        out[f"support_models_{source}"] = ""
        out[f"v2_support_flag"] = np.nan if source == "v2" else False
        return out

    event_cols = [c for c in ["feature_set", "model", "repeat", "fold", "outer_repeat", "outer_fold"] if c in events.columns]
    denom = (
        events[["analysis_level", "region_or_subregion"] + event_cols]
        .drop_duplicates()
        .groupby(["analysis_level", "region_or_subregion"])
        .size()
        .rename("n_events_region")
        .reset_index()
    )
    # Denominator is the complete non-baseline event grid for the region, not
    # the subset of events in which a particular Tier2 gene was eligible.
    events = events.merge(keys, on=["analysis_level", "region_or_subregion", "gene_id"], how="inner")
    selected = events[events["selected_conservative"]].copy()
    if selected.empty:
        out = keys.copy()
        out[f"selection_frequency_{source}"] = 0.0
        return out

    def sign_consistency(values: pd.Series) -> float:
        vals = pd.to_numeric(values, errors="coerce").dropna()
        vals = vals[vals.abs() > 1e-10]
        if vals.empty:
            return np.nan
        return max((vals > 0).mean(), (vals < 0).mean())

    def dominant_direction(values: pd.Series) -> str:
        vals = pd.to_numeric(values, errors="coerce").dropna()
        vals = vals[vals.abs() > 1e-10]
        if vals.empty:
            return "unavailable"
        pos = int((vals > 0).sum())
        neg = int((vals < 0).sum())
        if pos == neg:
            return "mixed"
        consistency = max(pos, neg) / (pos + neg)
        if consistency < 0.6:
            return "mixed"
        return "up_in_MDD" if pos > neg else "down_in_MDD"

    grouped = (
        selected.groupby(["analysis_level", "region_or_subregion", "gene_id"], as_index=False)
        .agg(
            selected_event_count=("selected_conservative", "sum"),
            model_support_count_source=("model", lambda x: len(set(x) & set(MODEL_FAMILIES))),
            support_models_source=("model", lambda x: ";".join(sorted(set(x) & set(MODEL_FAMILIES)))),
            sparse_model_support_source=("model", lambda x: bool(set(x) & SPARSE_MODELS)),
            linear_model_support_source=("model", lambda x: bool(set(x) & LINEAR_MODELS)),
            rf_support_source=("model", lambda x: bool(set(x) & RF_MODELS)),
            mean_abs_importance_source=("abs_importance", "mean"),
            median_abs_importance_source=("abs_importance", "median"),
            mean_rank_source=("rank_within_event", "mean"),
            top25_count_source=("top25", "sum"),
            top50_count_source=("top50", "sum"),
            near_top_count_source=("top20pct", "sum"),
            sign_consistency_source=("signed_importance_value", sign_consistency),
            ml_direction_source=("signed_importance_value", dominant_direction),
        )
    )
    grouped = grouped.merge(denom, on=["analysis_level", "region_or_subregion"], how="left")
    grouped[f"selection_frequency_{source}"] = grouped["selected_event_count"] / grouped["n_events_region"].replace(0, np.nan)
    grouped[f"rank_stability_{source}"] = grouped["near_top_count_source"] / grouped["selected_event_count"].replace(0, np.nan)
    rename = {
        "model_support_count_source": f"model_support_count_{source}",
        "support_models_source": f"support_models_{source}",
        "sparse_model_support_source": f"sparse_model_support_{source}",
        "linear_model_support_source": f"linear_model_support_{source}",
        "rf_support_source": f"rf_support_{source}",
        "mean_abs_importance_source": f"mean_abs_importance_{source}",
        "median_abs_importance_source": f"median_abs_importance_{source}",
        "mean_rank_source": f"mean_rank_{source}",
        "top25_count_source": f"top25_count_{source}",
        "top50_count_source": f"top50_count_{source}",
        "sign_consistency_source": f"sign_consistency_{source}",
        "ml_direction_source": f"ml_direction_{source}",
        "n_events_region": f"n_events_{source}",
    }
    grouped = grouped.rename(columns=rename)
    out = keys.merge(grouped, on=["analysis_level", "region_or_subregion", "gene_id"], how="left")
    out[f"selection_frequency_{source}"] = out[f"selection_frequency_{source}"].fillna(0.0)
    out[f"model_support_count_{source}"] = out[f"model_support_count_{source}"].fillna(0).astype(int)
    for col in [f"sparse_model_support_{source}", f"linear_model_support_{source}", f"rf_support_{source}"]:
        out[col] = out[col].fillna(False).astype(bool)
    if source == "v2":
        out["v2_support_flag"] = out[f"selection_frequency_{source}"] > 0
    return out


def minmax_by_region(series: pd.Series, groups: pd.Series) -> pd.Series:
    out = pd.Series(index=series.index, dtype=float)
    for _, idx in groups.groupby(groups).groups.items():
        vals = series.loc[idx].astype(float)
        min_val = vals.min(skipna=True)
        max_val = vals.max(skipna=True)
        if pd.isna(min_val) or pd.isna(max_val) or max_val == min_val:
            out.loc[idx] = 0.5
        else:
            out.loc[idx] = (vals - min_val) / (max_val - min_val)
    return out.fillna(0.0)


def build_priority_table(tier2: pd.DataFrame, v1_summary: pd.DataFrame, v2_summary: pd.DataFrame | None) -> pd.DataFrame:
    result = tier2.merge(v1_summary, on=["analysis_level", "region_or_subregion", "gene_id"], how="left")
    v2_available = v2_summary is not None and not v2_summary.empty and "selection_frequency_v2" in v2_summary.columns
    if v2_available:
        result = result.merge(v2_summary, on=["analysis_level", "region_or_subregion", "gene_id"], how="left")
    else:
        result["selection_frequency_v2"] = np.nan
        result["v2_support_flag"] = np.nan

    for source in ["v1", "v2"]:
        for col in [
            f"selection_frequency_{source}",
            f"model_support_count_{source}",
            f"mean_abs_importance_{source}",
            f"median_abs_importance_{source}",
            f"mean_rank_{source}",
            f"rank_stability_{source}",
            f"sign_consistency_{source}",
            f"n_events_{source}",
        ]:
            if col not in result.columns:
                result[col] = np.nan
        for col in [f"sparse_model_support_{source}", f"linear_model_support_{source}", f"rf_support_{source}"]:
            if col not in result.columns:
                result[col] = False

    result["selection_frequency_v1"] = result["selection_frequency_v1"].fillna(0.0)
    if v2_available:
        result["selection_frequency_v2"] = result["selection_frequency_v2"].fillna(0.0)
        result["selection_frequency_combined"] = 0.6 * result["selection_frequency_v2"] + 0.4 * result["selection_frequency_v1"]
        result["v2_support_score"] = result["v2_support_flag"].fillna(False).astype(bool).astype(float)
    else:
        result["selection_frequency_combined"] = result["selection_frequency_v1"]
        result["v2_support_score"] = 0.5

    result["model_support_count"] = (
        result["model_support_count_v1"].fillna(0).astype(int)
        + result["model_support_count_v2"].fillna(0).astype(int)
    )
    # Count a family once even if v1 and v2 both support it.
    for idx, row in result.iterrows():
        models = set()
        for col in ["support_models_v1", "support_models_v2"]:
            if col in result.columns and pd.notna(row.get(col)):
                models.update(str(row[col]).split(";"))
        models.discard("")
        result.at[idx, "model_support_count"] = len(models & set(MODEL_FAMILIES))

    result["sparse_model_support"] = result["sparse_model_support_v1"].fillna(False) | result["sparse_model_support_v2"].fillna(False)
    result["linear_model_support"] = result["linear_model_support_v1"].fillna(False) | result["linear_model_support_v2"].fillna(False)
    result["rf_support"] = result["rf_support_v1"].fillna(False) | result["rf_support_v2"].fillna(False)
    result["mean_abs_importance"] = result[["mean_abs_importance_v1", "mean_abs_importance_v2"]].mean(axis=1, skipna=True)
    result["median_abs_importance"] = result[["median_abs_importance_v1", "median_abs_importance_v2"]].median(axis=1, skipna=True)
    result["mean_rank"] = result[["mean_rank_v1", "mean_rank_v2"]].mean(axis=1, skipna=True)
    result["rank_stability"] = result[["rank_stability_v1", "rank_stability_v2"]].mean(axis=1, skipna=True).fillna(0.0)
    result["sign_consistency"] = result[["sign_consistency_v1", "sign_consistency_v2"]].mean(axis=1, skipna=True)

    def choose_direction(row: pd.Series) -> str:
        dirs = [row.get("ml_direction_v2"), row.get("ml_direction_v1")]
        dirs = [d for d in dirs if isinstance(d, str) and d not in {"", "nan", "unavailable"}]
        if not dirs:
            return "unavailable"
        if "mixed" in dirs:
            return "mixed"
        return dirs[0]

    result["ml_direction"] = result.apply(choose_direction, axis=1)
    result["lmm_ml_sign_agreement"] = np.where(
        result["ml_direction"].isin(["up_in_MDD", "down_in_MDD"]),
        result["ml_direction"].eq(result["lmm_direction"]),
        np.nan,
    )

    max_models_by_region = result.groupby("region_or_subregion")["model_support_count"].transform("max").replace(0, np.nan)
    result["model_support_scaled"] = (result["model_support_count"] / max_models_by_region).fillna(0.0)
    max_rank_by_region = result.groupby("region_or_subregion")["mean_rank"].transform("max")
    result["importance_rank_score"] = (1 - ((result["mean_rank"] - 1) / max_rank_by_region.replace(0, np.nan))).clip(0, 1).fillna(0.0)
    result["sign_consistency_score"] = result["sign_consistency"].fillna(0.5)
    result["lmm_evidence_scaled"] = minmax_by_region(result["lmm_evidence_score"], result["region_or_subregion"])
    result["sparse_model_support_binary"] = result["sparse_model_support"].astype(float)
    result["ML_priority_score"] = sum(
        weight * result[component] for component, weight in SCORE_WEIGHTS.items()
    ).clip(0, 1)

    agreement_allowed = result["lmm_ml_sign_agreement"].map(lambda x: True if pd.isna(x) else bool(x))
    agreement_not_false = result["lmm_ml_sign_agreement"].map(lambda x: True if pd.isna(x) else bool(x))
    high = (
        (result["ML_priority_score"] >= 0.70)
        & (result["selection_frequency_combined"] >= 0.50)
        & (result["model_support_count"] >= 2)
        & result["sparse_model_support"]
        & result["sign_consistency"].notna()
        & (result["sign_consistency"] >= 0.70)
        & agreement_allowed
    )
    moderate = (
        (result["ML_priority_score"] >= 0.45)
        & (result["selection_frequency_combined"] >= 0.25)
        & (result["model_support_count"] >= 1)
        & agreement_not_false
    )
    result["ML_priority_class"] = np.select(
        [high, moderate],
        ["High-confidence ML-prioritized Tier 2", "Moderate-confidence ML-prioritized Tier 2"],
        default="Exploratory ML-prioritized Tier 2",
    )
    result["enrichment_ready_flag"] = result["gene_symbol"].astype(str).str.len() > 0
    result = result.sort_values(["region_or_subregion", "ML_priority_score"], ascending=[True, False])
    return result


def make_summaries(result: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    region_rows = []
    for region, sub in result.groupby("region_or_subregion", sort=False):
        top10 = sub.sort_values("ML_priority_score", ascending=False).head(10)["gene_symbol"].tolist()
        region_rows.append(
            {
                "region_or_subregion": region,
                "n_Tier2_genes": len(sub),
                "n_high_confidence": int((sub["ML_priority_class"] == "High-confidence ML-prioritized Tier 2").sum()),
                "n_moderate_confidence": int((sub["ML_priority_class"] == "Moderate-confidence ML-prioritized Tier 2").sum()),
                "n_with_sparse_support": int(sub["sparse_model_support"].sum()),
                "n_with_v2_support": int(sub["v2_support_flag"].fillna(False).sum()),
                "median_ML_priority_score": float(sub["ML_priority_score"].median()),
                "top_10_gene_symbols": ";".join([x for x in top10 if isinstance(x, str) and x]),
            }
        )
    model_rows = []
    for region, sub in result.groupby("region_or_subregion", sort=False):
        for model in MODEL_FAMILIES:
            support = sub.apply(
                lambda row: model in str(row.get("support_models_v1", "")).split(";")
                or model in str(row.get("support_models_v2", "")).split(";"),
                axis=1,
            )
            model_rows.append({"region_or_subregion": region, "model_family": model, "n_Tier2_genes_supported": int(support.sum())})
    return pd.DataFrame(region_rows), pd.DataFrame(model_rows)


def write_gene_lists(result: pd.DataFrame) -> None:
    for _, _, key in REGIONS:
        region = key.replace("_", " ")
        sub = result[result["region_or_subregion"].map(safe_name) == key]
        high = sub[sub["ML_priority_class"] == "High-confidence ML-prioritized Tier 2"]
        moderate = sub[sub["ML_priority_class"] == "Moderate-confidence ML-prioritized Tier 2"]
        for label, table in [("high_confidence", high), ("moderate_confidence", moderate)]:
            symbols = table["gene_symbol"].dropna().astype(str)
            symbols = symbols[symbols.str.len() > 0].drop_duplicates().tolist()
            (GENE_LIST_DIR / f"{key}_{label}_symbols.txt").write_text("\n".join(symbols) + ("\n" if symbols else ""), encoding="utf-8")


def make_plots(result: pd.DataFrame, summary: pd.DataFrame) -> None:
    counts = summary.set_index("region_or_subregion")[["n_high_confidence", "n_moderate_confidence"]]
    counts.plot(kind="bar", figsize=(9, 4))
    plt.ylabel("Number of genes")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "barplot_high_moderate_counts_by_region.png", dpi=180)
    plt.close()

    for region, sub in result.groupby("region_or_subregion", sort=False):
        top = sub.sort_values("ML_priority_score", ascending=False).head(20)
        labels = top["gene_symbol"].fillna(top["gene_id"]).astype(str)
        plt.figure(figsize=(8, 5))
        plt.barh(labels[::-1], top["ML_priority_score"].to_numpy()[::-1])
        plt.xlabel("ML priority score")
        plt.title(f"{region}: top ML-prioritized Tier 2 genes")
        plt.tight_layout()
        plt.savefig(PLOT_DIR / f"{safe_name(region)}_top20_ML_priority_genes.png", dpi=180)
        plt.close()

        plt.figure(figsize=(5, 4))
        plt.scatter(sub["standardized_effect_size_diagnosis"], sub["ML_priority_score"], s=18, alpha=0.75)
        plt.xlabel("LMM standardized effect size")
        plt.ylabel("ML priority score")
        plt.title(region)
        plt.tight_layout()
        plt.savefig(PLOT_DIR / f"{safe_name(region)}_std_effect_vs_ML_priority.png", dpi=180)
        plt.close()

        plt.figure(figsize=(5, 4))
        neg_logp = -np.log10(sub["p_value_diagnosis"].clip(lower=1e-300))
        plt.scatter(neg_logp, sub["ML_priority_score"], s=18, alpha=0.75)
        plt.xlabel("-log10 diagnosis p-value")
        plt.ylabel("ML priority score")
        plt.title(region)
        plt.tight_layout()
        plt.savefig(PLOT_DIR / f"{safe_name(region)}_neglogp_vs_ML_priority.png", dpi=180)
        plt.close()

    heat_source = result.sort_values("ML_priority_score", ascending=False).groupby("region_or_subregion").head(10)
    heat = heat_source.pivot_table(index="gene_symbol", columns="region_or_subregion", values="ML_priority_score", aggfunc="max", fill_value=0)
    if not heat.empty and len(heat) <= 80:
        plt.figure(figsize=(9, max(5, 0.18 * len(heat))))
        plt.imshow(heat.to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=1)
        plt.colorbar(label="ML priority score")
        plt.yticks(np.arange(len(heat.index)), heat.index, fontsize=7)
        plt.xticks(np.arange(len(heat.columns)), heat.columns, rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(PLOT_DIR / "heatmap_top_stable_ML_prioritized_genes.png", dpi=180)
        plt.close()


def write_report(result: pd.DataFrame, summary: pd.DataFrame, used: list[str], missing: list[str]) -> None:
    top_tables = []
    for region, sub in result.groupby("region_or_subregion", sort=False):
        top = sub.sort_values("ML_priority_score", ascending=False).head(10)
        top_tables.append(f"### {region}\n\n" + top[["gene_symbol", "ML_priority_score", "ML_priority_class", "selection_frequency_combined", "model_support_count", "lmm_direction", "ml_direction"]].to_markdown(index=False))
    lines = [
        "# ML-prioritized Tier 2 genes for MDD/control separation",
        "",
        "## Purpose",
        "",
        "This analysis identifies Tier 2 exploratory DE candidate genes repeatedly prioritized by ML models for separating MDD vs Control.",
        "",
        "## Important Interpretation",
        "",
        "This is gene prioritization, not diagnostic classifier validation. ML v1 is used as Tier2-guided feature prioritization / exploratory classification screening. ML v2 is used as stricter support when a Tier 2 gene was selected in leakage-controlled folds.",
        "",
        "## Input Data",
        "",
        "Files used:",
        "",
        "\n".join(f"- `{path}`" for path in used),
        "",
        "Missing optional files:",
        "",
        "\n".join(f"- `{path}`" for path in missing) if missing else "None.",
        "",
        "## Methods",
        "",
        "The candidate pool was restricted to existing Tier 2 genes. Final v1 scoring used Tier2_fixed only; historical Tier4 v1 sets did not enter the final score. Strict v2 used nested train-only feature-selection/top-k evidence. Random forest was retained as non-directional selection/rank evidence, the majority-class baseline was excluded from gene-level evidence, and missing sign consistency could not satisfy the high-confidence sign gate. Predictive-performance metrics such as ROC-AUC did not enter the gene score. The component weights were 0.40/0.15/0.10/0.15/0.10/0.05/0.05.",
        "",
        "## Results By Region/Subregion",
        "",
        summary.to_markdown(index=False),
        "",
        "## Top Genes By Region/Subregion",
        "",
        "\n\n".join(top_tables),
        "",
        "## Interpretation",
        "",
        "Regions with more high/moderate genes have more internally supported ML-prioritized Tier 2 candidates. These genes are candidates for follow-up, not validated biomarkers.",
        "",
        "## Recommended Next Step",
        "",
        "Use these ML-prioritized Tier 2 genes for enrichment revisit, overlap with previous Tier 2 enrichment terms, AHBA/WGCNA module mapping, and pathway/module-level ML.",
        "",
        "## Warnings",
        "",
        "Do not claim these are validated biomarkers. Do not claim the classifier is externally validated. Use language such as ML-prioritized candidate genes, exploratory gene prioritization, internally supported, and requires validation.",
        "",
    ]
    (OUT_DIR / "ML_Prioritized_Tier2_Report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    used: list[str] = []
    missing: list[str] = []
    tier2 = load_tier2_tables(used, missing)

    v1_fold = load_optional_csv(V1_IMPORTANCE_DIR / "fold_level_feature_importance.csv", used, missing, dtype={"gene_id": str})
    _ = load_optional_csv(V1_IMPORTANCE_DIR / "feature_importance_summary.csv", used, missing, dtype={"gene_id": str})
    _ = load_optional_csv(V1_PERFORMANCE, used, missing)
    _ = load_optional_csv(V1_BEST, used, missing)
    v2_fold = load_optional_csv(V2_STABILITY_DIR / "selected_genes_per_fold_all_models.csv", used, missing, dtype={"gene_id": str})
    _ = load_optional_csv(V2_STABILITY_DIR / "selection_frequency_all_models.csv", used, missing, dtype={"gene_id": str})
    _ = load_optional_csv(V2_PERFORMANCE, used, missing)
    _ = load_optional_csv(V2_BEST, used, missing)

    required_missing = []
    if v1_fold is None:
        required_missing.append(str(V1_IMPORTANCE_DIR / "fold_level_feature_importance.csv"))
    if v2_fold is None:
        required_missing.append(str(V2_STABILITY_DIR / "selected_genes_per_fold_all_models.csv"))
    if required_missing:
        raise FileNotFoundError(
            "Cannot regenerate corrected ML-priority results without required event-level inputs: "
            + "; ".join(required_missing)
        )

    # Do not clear retained results until every required scoring input has been
    # loaded successfully.
    reset_output()

    if v1_fold is not None and "feature_set" in v1_fold.columns:
        v1_fold = v1_fold[v1_fold["feature_set"].eq("Tier2_fixed")].copy()
    v1_events = annotate_importance_events(v1_fold, "v1")
    v2_events = annotate_importance_events(v2_fold, "v2") if v2_fold is not None else pd.DataFrame()
    v1_summary = summarize_events(v1_events, tier2, "v1")
    v2_summary = summarize_events(v2_events, tier2, "v2") if not v2_events.empty else None

    result = build_priority_table(tier2, v1_summary, v2_summary)
    safe_p = pd.to_numeric(result["p_value_diagnosis"], errors="coerce").clip(lower=np.nextafter(0, 1))
    std_effect = pd.to_numeric(result["standardized_effect_size_diagnosis"], errors="coerce")
    # Tier 4-style score: keep the MDD direction in the sign, and use abs_score only for ranking magnitude.
    result["signed_score"] = np.sign(std_effect) * -np.log10(safe_p)
    result["abs_score"] = result["signed_score"].abs()
    result["direction"] = result["lmm_direction"]
    summary, model_support = make_summaries(result)
    write_gene_lists(result)
    make_plots(result, summary)
    write_report(result, summary, used, missing)

    output_cols = [
        "analysis_level",
        "region_or_subregion",
        "gene_id",
        "gene_symbol",
        "beta_diagnosis",
        "p_value_diagnosis",
        "standardized_effect_size_diagnosis",
        "signed_score",
        "abs_score",
        "direction",
        "lmm_direction",
        "selection_frequency_v1",
        "selection_frequency_v2",
        "selection_frequency_combined",
        "n_events_v1",
        "n_events_v2",
        "model_support_count",
        "sparse_model_support",
        "linear_model_support",
        "rf_support",
        "mean_abs_importance",
        "median_abs_importance",
        "mean_rank",
        "rank_stability",
        "sign_consistency",
        "sign_consistency_v1",
        "sign_consistency_v2",
        "ml_direction",
        "lmm_ml_sign_agreement",
        "lmm_evidence_score",
        "v2_support_flag",
        "model_support_scaled",
        "sparse_model_support_binary",
        "importance_rank_score",
        "sign_consistency_score",
        "lmm_evidence_scaled",
        "v2_support_score",
        "ML_priority_score",
        "ML_priority_class",
        "enrichment_ready_flag",
    ]
    result[output_cols].to_csv(OUT_DIR / "ML_prioritized_Tier2_genes_all_ranked.csv", index=False)
    result[result["ML_priority_class"].eq("High-confidence ML-prioritized Tier 2")][output_cols].to_csv(
        OUT_DIR / "ML_prioritized_Tier2_genes_high_confidence.csv", index=False
    )
    result[result["ML_priority_class"].eq("Moderate-confidence ML-prioritized Tier 2")][output_cols].to_csv(
        OUT_DIR / "ML_prioritized_Tier2_genes_moderate_confidence.csv", index=False
    )
    summary.to_csv(OUT_DIR / "ML_prioritized_Tier2_gene_summary_by_region.csv", index=False)
    model_support.to_csv(OUT_DIR / "model_support_summary.csv", index=False)

    manifest = {
        "script_path": str(Path(__file__).resolve()),
        "output_path": str(OUT_DIR),
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "input_files_used": used,
        "files_missing": missing,
        "number_of_rows_produced": int(len(result)),
        "regions_analyzed": [region for _, region, _ in REGIONS],
        "interpretation": "Gene prioritization only; not diagnostic classifier validation.",
    }
    (OUT_DIR / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("\nML-prioritized Tier 2 summary:")
    print(summary.to_string(index=False))
    print(f"\nRows produced: {len(result)}")
    print(f"Output folder: {OUT_DIR}")


if __name__ == "__main__":
    main()
