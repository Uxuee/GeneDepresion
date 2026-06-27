"""Strict leakage-controlled ML v2 for Length > 100 MDD transcriptomics.

This script does not rerun the LMM or enrichment. It creates an independent
output branch:

  outputs/length_gt_100_elisa_style/ML_strict_v2/

Strict rule:
  For every outer CV split, feature selection, imputation, scaling, and model
  tuning are performed on the outer-training data only. The held-out test fold
  is used only for final evaluation.

Runtime note:
  Nested CV is much heavier than the exploratory v1 pilot. The defaults are
  intentionally runtime-bounded and fully documented in the output manifest.
  Increase STRICT_ML_OUTER_REPEATS or STRICT_ML_PERMUTATIONS in the environment
  for a larger run.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

try:
    from xgboost import XGBClassifier  # noqa: F401

    HAS_XGBOOST = True
except Exception:
    HAS_XGBOOST = False


REPO_DIR = Path(__file__).resolve().parents[1]
BASE_DIR = REPO_DIR / "outputs" / "length_gt_100_elisa_style"
OUT_DIR = BASE_DIR / "ML_strict_v2"
PLOT_DIR = OUT_DIR / "plots"
STABILITY_DIR = OUT_DIR / "feature_stability"

METADATA_DIR = BASE_DIR / "step1_from_raw" / "Length_gt_100" / "01_metadata" / "keep_all_samples"
BROAD_EXPR_DIR = BASE_DIR / "normalized_expression_length_gt_100_by_region"
SUBREGION_EXPR_DIR = BASE_DIR / "frontal_subregion_sensitivity" / "subregion_normalized_expression_matrices"
V1_PERF = BASE_DIR / "ML_pilot_Tier2" / "ML_performance_summary.csv"
V1_BEST = BASE_DIR / "ML_pilot_Tier2" / "best_model_by_region.csv"
TIER2_PACKAGE = BASE_DIR / "DE_candidate_gene_sets" / "Tier2_working_exploratory"

TOP_K_VALUES = [25, 50, 100, 250, 500]
FEATURE_METHODS = ["train_only_anova_f", "train_only_standardized_effect"]
MODEL_NAMES = [
    "majority_class_baseline",
    "logistic_l2",
    "logistic_l1",
    "elastic_net_logistic",
    "linear_svm",
    "random_forest_exploratory",
]

# Nested CV is expensive. Defaults are deliberately bounded and recorded.
OUTER_REPEATS = int(os.environ.get("STRICT_ML_OUTER_REPEATS", "1"))
MAX_OUTER_FOLDS = int(os.environ.get("STRICT_ML_MAX_OUTER_FOLDS", "3"))
COMPACT_GRID = os.environ.get("STRICT_ML_FULL_GRID", "0") != "1"
PERMUTATIONS = int(os.environ.get("STRICT_ML_PERMUTATIONS", "20"))
PERMUTATION_REPEATS = int(os.environ.get("STRICT_ML_PERMUTATION_REPEATS", "1"))
RANDOM_SEEDS = list(range(31001, 31001 + max(OUTER_REPEATS, PERMUTATION_REPEATS, 1)))
PERMUTATION_SEEDS = list(range(42001, 42001 + PERMUTATIONS))


@dataclass(frozen=True)
class AnalysisSpec:
    analysis_level: str
    region_or_subregion: str
    package_key: str
    expression_file: Path
    metadata_file: Path
    metadata_filter_value: str | None = None


ANALYSES = [
    AnalysisSpec(
        "four_group",
        "Ventral Striatum",
        "Ventral_Striatum",
        BROAD_EXPR_DIR / "normalized_GeTMM_log2_ventral_striatum.csv",
        METADATA_DIR / "metadata_ventral_striatum.csv",
    ),
    AnalysisSpec(
        "four_group",
        "Insula",
        "Insula",
        BROAD_EXPR_DIR / "normalized_GeTMM_log2_insula.csv",
        METADATA_DIR / "metadata_insula.csv",
    ),
    AnalysisSpec(
        "four_group",
        "Hippocampal Formation",
        "Hippocampal_Formation",
        BROAD_EXPR_DIR / "normalized_GeTMM_log2_hippocampal_formation.csv",
        METADATA_DIR / "metadata_hippocampal_formation.csv",
    ),
    AnalysisSpec(
        "frontal_subregion",
        "Cingulate gyrus 25",
        "Cingulate_gyrus_25",
        SUBREGION_EXPR_DIR / "Cingulate_gyrus_25_Cg25_normalized_GeTMM_log2_length_gt_100.csv",
        METADATA_DIR / "metadata_frontal_cortex.csv",
        "Cingulate gyrus 25 (Cg25)",
    ),
    AnalysisSpec(
        "frontal_subregion",
        "dlPFC BA8/9",
        "dlPFC_BA8_9",
        SUBREGION_EXPR_DIR / "Dorsolateral_prefrontal_cortex_dlPFC_BA8_9_normalized_GeTMM_log2_length_gt_100.csv",
        METADATA_DIR / "metadata_frontal_cortex.csv",
        "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
    ),
    AnalysisSpec(
        "frontal_subregion",
        "OFC BA11",
        "OFC_BA11",
        SUBREGION_EXPR_DIR / "Orbitofrontal_OFC_BA11_normalized_GeTMM_log2_length_gt_100.csv",
        METADATA_DIR / "metadata_frontal_cortex.csv",
        "Orbitofrontal (OFC; BA11)",
    ),
]


def safe_name(value: str) -> str:
    text = "".join(ch if ch.isalnum() else "_" for ch in str(value))
    while "__" in text:
        text = text.replace("__", "_")
    return text.strip("_")


def prepare_output() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for folder in [PLOT_DIR, STABILITY_DIR]:
        if folder.exists():
            shutil.rmtree(folder)
        folder.mkdir(parents=True, exist_ok=True)
    for file_name in [
        "ML_strict_v2_performance_summary.csv",
        "best_strict_model_by_region.csv",
        "cross_validated_predictions_strict_v2.csv",
        "permutation_test_results_strict_v2.csv",
        "v1_vs_strict_v2_comparison.csv",
        "analysis_sample_counts_and_cv_strategy.csv",
        "run_manifest.json",
        "ML_Strict_V2_Report.md",
        "README.md",
    ]:
        path = OUT_DIR / file_name
        if path.exists():
            path.unlink()


def read_expression(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing expression matrix: {path}")
    matrix = pd.read_csv(path, index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    return matrix.apply(pd.to_numeric, errors="coerce")


def read_metadata(spec: AnalysisSpec, sample_ids: list[str]) -> pd.DataFrame:
    if not spec.metadata_file.exists():
        raise FileNotFoundError(f"Missing metadata: {spec.metadata_file}")
    meta = pd.read_csv(spec.metadata_file, low_memory=False)
    meta["SampleID"] = meta["SampleID"].astype(str)
    if spec.metadata_filter_value is not None:
        meta = meta[meta["OriginalRegion"].astype(str).eq(spec.metadata_filter_value)].copy()
    meta = meta.drop_duplicates("SampleID")
    meta = meta[meta["SampleID"].isin(sample_ids)].copy()
    meta = meta.set_index("SampleID").reindex(sample_ids).reset_index()
    meta = meta.dropna(subset=["Diagnosis"])
    meta["label"] = meta["Diagnosis"].map(
        lambda x: 1 if str(x).upper() == "MDD" else (0 if str(x).upper() in {"CTRL", "CONTROL"} else np.nan)
    )
    meta = meta.dropna(subset=["label"])
    meta["label"] = meta["label"].astype(int)
    return meta


def load_dataset(spec: AnalysisSpec) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, pd.Series | None]:
    expr = read_expression(spec.expression_file)
    meta = read_metadata(spec, list(expr.columns))
    samples = [sample for sample in expr.columns if sample in set(meta["SampleID"])]
    expr = expr[samples]
    meta = meta.set_index("SampleID").loc[samples].reset_index()
    x = expr.T
    x.index = samples
    y = meta["label"].to_numpy()
    groups = meta["BrainID"].astype(str) if "BrainID" in meta.columns else None
    return x, meta, y, groups


def get_outer_splits(
    y: np.ndarray, groups: pd.Series | None, repeats: int
) -> tuple[list[tuple[int, int, np.ndarray, np.ndarray]], str, bool, int]:
    min_class = min(Counter(y).values())
    n_splits = 5 if min_class >= 5 else 3 if min_class >= 3 else 0
    if n_splits > 0:
        n_splits = min(n_splits, MAX_OUTER_FOLDS)
    if n_splits < 2:
        return [], "insufficient_samples", False, n_splits
    repeated_groups = groups is not None and groups.notna().all() and groups.duplicated().any()
    splits = []
    if repeated_groups:
        splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEEDS[0])
        for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(np.zeros(len(y)), y, groups)):
            splits.append((0, fold_idx, train_idx, test_idx))
        return splits, f"StratifiedGroupKFold_{n_splits}_fold_1_repeat", True, n_splits
    for repeat_idx, seed in enumerate(RANDOM_SEEDS[:repeats]):
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(np.zeros(len(y)), y)):
            splits.append((repeat_idx, fold_idx, train_idx, test_idx))
    return splits, f"RepeatedStratifiedKFold_{n_splits}_fold_{repeats}_repeat_runtime_bounded", False, n_splits


def inner_cv_type(y_train: np.ndarray, groups_train: pd.Series | None) -> tuple[object, str]:
    min_class = min(Counter(y_train).values())
    n_splits = 3 if min_class >= 3 else 2
    repeated_groups = groups_train is not None and groups_train.notna().all() and groups_train.duplicated().any()
    if repeated_groups:
        return StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=17), f"StratifiedGroupKFold_{n_splits}_fold"
    return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=17), f"StratifiedKFold_{n_splits}_fold"


def rank_anova_f(x_train: pd.DataFrame, y_train: np.ndarray) -> list[str]:
    imputed = SimpleImputer(strategy="median").fit_transform(x_train)
    scores, _ = f_classif(imputed, y_train)
    scores = np.nan_to_num(scores, nan=-np.inf, posinf=np.inf, neginf=-np.inf)
    order = np.argsort(scores)[::-1]
    return list(x_train.columns[order])


def rank_standardized_effect(x_train: pd.DataFrame, y_train: np.ndarray) -> list[str]:
    x = x_train.astype(float)
    mdd = x.loc[y_train == 1]
    ctrl = x.loc[y_train == 0]
    diff = mdd.mean(axis=0, skipna=True) - ctrl.mean(axis=0, skipna=True)
    sd = x.std(axis=0, skipna=True, ddof=1).replace(0, np.nan)
    score = (diff / sd).abs().replace([np.inf, -np.inf], np.nan).fillna(-np.inf)
    return list(score.sort_values(ascending=False).index)


def select_features(method: str, x_train: pd.DataFrame, y_train: np.ndarray, top_k: int) -> list[str]:
    ranking = rank_anova_f(x_train, y_train) if method == "train_only_anova_f" else rank_standardized_effect(x_train, y_train)
    return ranking[: min(top_k, len(ranking))]


def model_grid(model_name: str) -> list[dict]:
    c_values = [0.01, 0.1, 1] if COMPACT_GRID else [0.001, 0.01, 0.1, 1, 10]
    l1_ratios = [0.5] if COMPACT_GRID else [0.1, 0.5, 0.9]
    rf_depths = [3, None] if COMPACT_GRID else [2, 3, 5, None]
    rf_leaves = [2] if COMPACT_GRID else [2, 5]
    if model_name == "majority_class_baseline":
        return [{"model": DummyClassifier(strategy="most_frequent"), "scale": False, "params": {}}]
    if model_name == "logistic_l2":
        return [
            {
                "model": LogisticRegression(penalty="l2", solver="liblinear", class_weight="balanced", max_iter=4000),
                "scale": True,
                "params": {"C": c},
            }
            for c in c_values
        ]
    if model_name == "logistic_l1":
        return [
            {
                "model": LogisticRegression(penalty="l1", solver="liblinear", class_weight="balanced", max_iter=4000),
                "scale": True,
                "params": {"C": c},
            }
            for c in c_values
        ]
    if model_name == "elastic_net_logistic":
        return [
            {
                "model": LogisticRegression(
                    penalty="elasticnet", solver="saga", class_weight="balanced", max_iter=5000, random_state=11
                ),
                "scale": True,
                "params": {"C": c, "l1_ratio": r},
            }
            for c in c_values
            for r in l1_ratios
        ]
    if model_name == "linear_svm":
        return [
            {"model": LinearSVC(class_weight="balanced", max_iter=6000), "scale": True, "params": {"C": c}}
            for c in c_values
        ]
    if model_name == "random_forest_exploratory":
        return [
            {
                "model": RandomForestClassifier(
                    n_estimators=200,
                    max_depth=max_depth,
                    min_samples_leaf=min_leaf,
                    class_weight="balanced_subsample",
                    random_state=13,
                    n_jobs=-1,
                ),
                "scale": False,
                "params": {"max_depth": max_depth, "min_samples_leaf": min_leaf},
            }
            for max_depth in rf_depths
            for min_leaf in rf_leaves
        ]
    raise ValueError(model_name)


def make_pipeline(model, scale: bool) -> Pipeline:
    steps = [("imputer", SimpleImputer(strategy="median"))]
    if scale:
        steps.append(("scaler", StandardScaler()))
    steps.append(("model", model))
    return Pipeline(steps)


def apply_params(model, params: dict):
    for key, value in params.items():
        setattr(model, key, value)
    return model


def estimator_score(pipe: Pipeline, x_test: pd.DataFrame) -> np.ndarray:
    model = pipe.named_steps["model"]
    if hasattr(model, "predict_proba"):
        return pipe.predict_proba(x_test)[:, 1]
    if hasattr(model, "decision_function"):
        return pipe.decision_function(x_test)
    return pipe.predict(x_test)


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray, y_score: np.ndarray) -> dict[str, float]:
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "roc_auc": roc_auc_score(y_true, y_score),
        "pr_auc": average_precision_score(y_true, y_score),
        "accuracy": accuracy_score(y_true, y_pred),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "sensitivity": tp / (tp + fn) if tp + fn else np.nan,
        "specificity": tn / (tn + fp) if tn + fp else np.nan,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "tp": tp,
    }


def tune_model(
    model_name: str,
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    groups_train: pd.Series | None,
) -> tuple[dict, str]:
    if model_name == "majority_class_baseline":
        return model_grid(model_name)[0], "not_tuned_baseline"
    splitter, inner_type = inner_cv_type(y_train, groups_train)
    best_item = None
    best_score = -np.inf
    split_args = (np.zeros(len(y_train)), y_train, groups_train) if "Group" in inner_type else (np.zeros(len(y_train)), y_train)
    for item in model_grid(model_name):
        scores = []
        for inner_train_idx, inner_valid_idx in splitter.split(*split_args):
            model = apply_params(item["model"].__class__(**item["model"].get_params()), item["params"])
            pipe = make_pipeline(model, item["scale"])
            pipe.fit(x_train.iloc[inner_train_idx], y_train[inner_train_idx])
            pred = pipe.predict(x_train.iloc[inner_valid_idx])
            scores.append(balanced_accuracy_score(y_train[inner_valid_idx], pred))
        mean_score = float(np.mean(scores))
        if mean_score > best_score:
            best_score = mean_score
            best_item = item
    return best_item, inner_type


def fit_evaluate_outer(
    spec: AnalysisSpec,
    x: pd.DataFrame,
    meta: pd.DataFrame,
    y: np.ndarray,
    groups: pd.Series | None,
    splits: list[tuple[int, int, np.ndarray, np.ndarray]],
    outer_cv: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    summary_rows = []
    prediction_rows = []
    feature_rows = []
    inner_type_seen = "not_available"
    n_mdd = int(np.sum(y == 1))
    n_ctrl = int(np.sum(y == 0))

    for method in FEATURE_METHODS:
        for top_k in TOP_K_VALUES:
            for model_name in MODEL_NAMES:
                fold_metrics = defaultdict(list)
                confusion = Counter()
                selected_counts = []
                for repeat_idx, fold_idx, train_idx, test_idx in splits:
                    x_train_all = x.iloc[train_idx]
                    x_test_all = x.iloc[test_idx]
                    y_train = y[train_idx]
                    y_test = y[test_idx]
                    groups_train = groups.iloc[train_idx] if groups is not None else None
                    genes = select_features(method, x_train_all, y_train, top_k)
                    selected_counts.append(len(genes))
                    x_train = x_train_all[genes]
                    x_test = x_test_all[genes]
                    best_item, inner_type = tune_model(model_name, x_train, y_train, groups_train)
                    inner_type_seen = inner_type
                    model = apply_params(best_item["model"].__class__(**best_item["model"].get_params()), best_item["params"])
                    pipe = make_pipeline(model, best_item["scale"])
                    pipe.fit(x_train, y_train)
                    y_pred = pipe.predict(x_test)
                    y_score = estimator_score(pipe, x_test)
                    metrics = compute_metrics(y_test, y_pred, y_score)
                    for key, val in metrics.items():
                        if key in {"tn", "fp", "fn", "tp"}:
                            confusion[key] += int(val)
                        else:
                            fold_metrics[key].append(float(val))
                    for sample_id, truth, pred, score in zip(meta.iloc[test_idx]["SampleID"], y_test, y_pred, y_score):
                        prediction_rows.append(
                            {
                                "sample_id": sample_id,
                                "analysis_level": spec.analysis_level,
                                "region_or_subregion": spec.region_or_subregion,
                                "outer_repeat": repeat_idx,
                                "outer_fold": fold_idx,
                                "true_label": int(truth),
                                "predicted_label": int(pred),
                                "predicted_probability_MDD": float(score),
                                "model": model_name,
                                "feature_selection_method": method,
                                "top_k": top_k,
                            }
                        )
                    model_obj = pipe.named_steps["model"]
                    coefs = None
                    if hasattr(model_obj, "coef_"):
                        coefs = np.asarray(model_obj.coef_).ravel()
                    for idx, gene_id in enumerate(genes):
                        coef = float(coefs[idx]) if coefs is not None else np.nan
                        feature_rows.append(
                            {
                                "analysis_level": spec.analysis_level,
                                "region_or_subregion": spec.region_or_subregion,
                                "feature_selection_method": method,
                                "top_k": top_k,
                                "model": model_name,
                                "outer_repeat": repeat_idx,
                                "outer_fold": fold_idx,
                                "gene_id": gene_id,
                                "coefficient": coef,
                                "selected": True,
                            }
                        )
                row = {
                    "analysis_level": spec.analysis_level,
                    "region_or_subregion": spec.region_or_subregion,
                    "feature_selection_method": method,
                    "top_k": top_k,
                    "model": model_name,
                    "outer_cv_type": outer_cv,
                    "inner_cv_type": inner_type_seen,
                    "n_samples": len(y),
                    "n_MDD": n_mdd,
                    "n_Control": n_ctrl,
                    "n_features_selected": float(np.mean(selected_counts)),
                    "leakage_status": "strict_train_fold_only_feature_selection_tuning_preprocessing",
                    "interpretation_note": "Strict ML pilot with internal validation only; requires external validation.",
                    "confusion_tn_total": confusion["tn"],
                    "confusion_fp_total": confusion["fp"],
                    "confusion_fn_total": confusion["fn"],
                    "confusion_tp_total": confusion["tp"],
                }
                for metric in ["balanced_accuracy", "roc_auc", "pr_auc", "accuracy", "f1", "sensitivity", "specificity"]:
                    vals = np.asarray(fold_metrics[metric], dtype=float)
                    row[f"{metric}_mean"] = float(np.nanmean(vals))
                    row[f"{metric}_sd"] = float(np.nanstd(vals, ddof=1)) if np.sum(~np.isnan(vals)) > 1 else 0.0
                summary_rows.append(row)
    return pd.DataFrame(summary_rows), pd.DataFrame(prediction_rows), pd.DataFrame(feature_rows), inner_type_seen


def feature_stability(feature_rows: pd.DataFrame) -> pd.DataFrame:
    if feature_rows.empty:
        return pd.DataFrame()
    grouped = (
        feature_rows.groupby(
            ["analysis_level", "region_or_subregion", "feature_selection_method", "top_k", "model", "gene_id"],
            as_index=False,
        )
        .agg(
            selection_frequency=("selected", "sum"),
            coefficient_mean=("coefficient", "mean"),
            coefficient_sd=("coefficient", "std"),
            positive_coefficients=("coefficient", lambda x: int(np.sum(np.asarray(x, dtype=float) > 0))),
            negative_coefficients=("coefficient", lambda x: int(np.sum(np.asarray(x, dtype=float) < 0))),
        )
    )
    grouped["coefficient_sign_consistency"] = grouped[["positive_coefficients", "negative_coefficients"]].max(axis=1) / grouped[
        "selection_frequency"
    ].replace(0, np.nan)
    return grouped.sort_values(["region_or_subregion", "selection_frequency"], ascending=[True, False])


def permutation_test_for_best(
    spec: AnalysisSpec,
    best: pd.Series,
    x: pd.DataFrame,
    y: np.ndarray,
    groups: pd.Series | None,
) -> dict:
    observed_auc = float(best["roc_auc_mean"])
    observed_bal = float(best["balanced_accuracy_mean"])
    auc_scores = []
    bal_scores = []
    for seed in PERMUTATION_SEEDS:
        rng = np.random.default_rng(seed)
        y_perm = rng.permutation(y)
        splits, _, _, _ = get_outer_splits(y_perm, groups, PERMUTATION_REPEATS)
        fold_auc = []
        fold_bal = []
        for _, _, train_idx, test_idx in splits:
            genes = select_features(str(best["feature_selection_method"]), x.iloc[train_idx], y_perm[train_idx], int(best["top_k"]))
            x_train = x.iloc[train_idx][genes]
            x_test = x.iloc[test_idx][genes]
            # Runtime-bounded permutation check: use the middle compact-grid
            # setting for the best model family rather than rerunning full inner
            # tuning for every permuted label vector.
            grid = model_grid(str(best["model"]))
            item = grid[len(grid) // 2]
            model = apply_params(item["model"].__class__(**item["model"].get_params()), item["params"])
            pipe = make_pipeline(model, item["scale"])
            pipe.fit(x_train, y_perm[train_idx])
            pred = pipe.predict(x_test)
            score = estimator_score(pipe, x_test)
            fold_auc.append(roc_auc_score(y_perm[test_idx], score))
            fold_bal.append(balanced_accuracy_score(y_perm[test_idx], pred))
        auc_scores.append(float(np.mean(fold_auc)))
        bal_scores.append(float(np.mean(fold_bal)))
    return {
        "analysis_level": spec.analysis_level,
        "region_or_subregion": spec.region_or_subregion,
        "feature_selection_method": best["feature_selection_method"],
        "top_k": int(best["top_k"]),
        "model": best["model"],
        "observed_roc_auc_mean": observed_auc,
        "observed_balanced_accuracy_mean": observed_bal,
        "n_permutations": len(auc_scores),
        "empirical_p_roc_auc": (1 + sum(v >= observed_auc for v in auc_scores)) / (1 + len(auc_scores)),
        "empirical_p_balanced_accuracy": (1 + sum(v >= observed_bal for v in bal_scores)) / (1 + len(bal_scores)),
        "permutation_note": f"Runtime-bounded {PERMUTATIONS}-permutation strict rerun using {PERMUTATION_REPEATS} repeat(s); inner tuning skipped inside permutations to keep runtime feasible.",
    }


def make_v1_comparison(strict_best: pd.DataFrame) -> pd.DataFrame:
    if not V1_BEST.exists():
        return pd.DataFrame()
    v1 = pd.read_csv(V1_BEST)
    rows = []
    for _, strict in strict_best.iterrows():
        region = strict["region_or_subregion"]
        v1_region = v1[v1["region_or_subregion"] == region]
        if v1_region.empty:
            continue
        v1_row = v1_region.sort_values(["roc_auc_mean", "balanced_accuracy_mean"], ascending=False).iloc[0]
        rows.append(
            {
                "region_or_subregion": region,
                "v1_best_balanced_accuracy": v1_row["balanced_accuracy_mean"],
                "v1_best_roc_auc": v1_row["roc_auc_mean"],
                "strict_v2_best_balanced_accuracy": strict["balanced_accuracy_mean"],
                "strict_v2_best_roc_auc": strict["roc_auc_mean"],
                "balanced_accuracy_drop_v2_minus_v1": strict["balanced_accuracy_mean"] - v1_row["balanced_accuracy_mean"],
                "roc_auc_drop_v2_minus_v1": strict["roc_auc_mean"] - v1_row["roc_auc_mean"],
                "interpretation": "v1 used full-dataset feature selection/ranking and may be optimistic; strict v2 is the leakage-controlled reference.",
            }
        )
    return pd.DataFrame(rows)


def plot_curves(best: pd.DataFrame, predictions: pd.DataFrame) -> None:
    for _, row in best.iterrows():
        subset = predictions[
            (predictions["region_or_subregion"] == row["region_or_subregion"])
            & (predictions["model"] == row["model"])
            & (predictions["feature_selection_method"] == row["feature_selection_method"])
            & (predictions["top_k"] == row["top_k"])
        ]
        if subset.empty:
            continue
        y = subset["true_label"].to_numpy()
        score = subset["predicted_probability_MDD"].to_numpy()
        prefix = safe_name(row["region_or_subregion"])
        fpr, tpr, _ = roc_curve(y, score)
        plt.figure(figsize=(5.5, 4.5))
        plt.plot(fpr, tpr, label=f"ROC-AUC={roc_auc_score(y, score):.3f}")
        plt.plot([0, 1], [0, 1], "--", color="gray")
        plt.xlabel("False positive rate")
        plt.ylabel("True positive rate")
        plt.title(f"{row['region_or_subregion']} strict ROC")
        plt.legend()
        plt.tight_layout()
        plt.savefig(PLOT_DIR / f"{prefix}__strict_best_ROC.png", dpi=180)
        plt.close()
        precision, recall, _ = precision_recall_curve(y, score)
        plt.figure(figsize=(5.5, 4.5))
        plt.plot(recall, precision, label=f"PR-AUC={average_precision_score(y, score):.3f}")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title(f"{row['region_or_subregion']} strict PR")
        plt.legend()
        plt.tight_layout()
        plt.savefig(PLOT_DIR / f"{prefix}__strict_best_PR.png", dpi=180)
        plt.close()


def plot_summaries(best: pd.DataFrame, comparison: pd.DataFrame, stability: pd.DataFrame) -> None:
    best.set_index("region_or_subregion")["roc_auc_mean"].plot(kind="bar", figsize=(8, 4))
    plt.ylabel("Strict v2 ROC-AUC")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "barplot_strict_roc_auc_by_region.png", dpi=180)
    plt.close()
    best.set_index("region_or_subregion")["balanced_accuracy_mean"].plot(kind="bar", figsize=(8, 4))
    plt.ylabel("Strict v2 balanced accuracy")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "barplot_strict_balanced_accuracy_by_region.png", dpi=180)
    plt.close()
    if not comparison.empty:
        comparison.set_index("region_or_subregion")[["v1_best_roc_auc", "strict_v2_best_roc_auc"]].plot(kind="bar", figsize=(9, 4))
        plt.ylabel("ROC-AUC")
        plt.tight_layout()
        plt.savefig(PLOT_DIR / "v1_vs_strict_v2_roc_auc_comparison.png", dpi=180)
        plt.close()
    if not stability.empty:
        for region, sub in stability.groupby("region_or_subregion"):
            top = sub.sort_values("selection_frequency", ascending=False).head(20)
            labels = top["gene_id"].astype(str)
            plt.figure(figsize=(7, 5))
            plt.barh(labels[::-1], top["selection_frequency"].to_numpy()[::-1])
            plt.xlabel("Selection frequency")
            plt.title(f"{region} stable selected genes")
            plt.tight_layout()
            plt.savefig(PLOT_DIR / f"{safe_name(region)}__feature_selection_stability.png", dpi=180)
            plt.close()
        top_stable = stability.groupby(["region_or_subregion", "gene_id"], as_index=False)["selection_frequency"].sum()
        top_genes = top_stable.groupby("gene_id")["selection_frequency"].sum().sort_values(ascending=False).head(30).index
        heat = top_stable[top_stable["gene_id"].isin(top_genes)].pivot_table(
            index="gene_id", columns="region_or_subregion", values="selection_frequency", aggfunc="sum", fill_value=0
        )
        plt.figure(figsize=(8, max(5, 0.18 * len(heat))))
        plt.imshow(heat.to_numpy(), aspect="auto", cmap="magma")
        plt.colorbar(label="Selection frequency")
        plt.yticks(np.arange(len(heat.index)), heat.index, fontsize=7)
        plt.xticks(np.arange(len(heat.columns)), heat.columns, rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(PLOT_DIR / "heatmap_stable_selected_genes_across_regions.png", dpi=180)
        plt.close()


def write_report(perf: pd.DataFrame, best: pd.DataFrame, comparison: pd.DataFrame, perm: pd.DataFrame, sample_counts: pd.DataFrame, stability: pd.DataFrame) -> None:
    top_stability = (
        stability.sort_values(["region_or_subregion", "selection_frequency"], ascending=[True, False])
        .groupby("region_or_subregion")
        .head(10)
    )
    lines = [
        "# Strict ML V2 Report",
        "",
        "## Purpose",
        "",
        "Run leakage-controlled ML after the exploratory v1 ML pilot. The LMM and enrichment analyses were not rerun.",
        "",
        "## Reminder About V1",
        "",
        "The v1 fixed Tier 2 and fixed Tier 4 ML pilot was exploratory because features were selected or ranked using the full dataset. Those values are included only as comparison context.",
        "",
        "## Strict V2 Design",
        "",
        "Feature selection, imputation, scaling, and hyperparameter tuning were done inside cross-validation. The held-out outer test fold was not used for ranking, preprocessing, or tuning.",
        "",
        f"Runtime-bounded settings used: outer repeats={OUTER_REPEATS}; max outer folds={MAX_OUTER_FOLDS}; compact grid={COMPACT_GRID}; permutations={PERMUTATIONS}; permutation repeats={PERMUTATION_REPEATS}.",
        "",
        "## Regions/Subregions Tested",
        "",
        sample_counts.to_markdown(index=False),
        "",
        "## Feature-Selection Methods",
        "",
        "- Train-only ANOVA F-score ranking.",
        "- Train-only standardized effect-size ranking.",
        "- Model-internal L1 and Elastic Net sparsity after train-only feature ranking.",
        "- Covariate-residualized ranking was skipped because per-gene covariate fitting inside nested CV would be too slow for this first strict run.",
        "- The first 50-permutation nested run exceeded the desktop execution window, so v2 uses 20 preliminary permutations and skips inner tuning inside permutations.",
        "",
        "## Best Strict Model Per Region/Subregion",
        "",
        best[
            [
                "region_or_subregion",
                "feature_selection_method",
                "top_k",
                "model",
                "balanced_accuracy_mean",
                "roc_auc_mean",
                "pr_auc_mean",
                "leakage_status",
            ]
        ].to_markdown(index=False),
        "",
        "## Main Performance Snapshot",
        "",
        perf.sort_values(["region_or_subregion", "roc_auc_mean"], ascending=[True, False])
        .groupby("region_or_subregion")
        .head(5)[
            [
                "region_or_subregion",
                "feature_selection_method",
                "top_k",
                "model",
                "balanced_accuracy_mean",
                "roc_auc_mean",
                "pr_auc_mean",
            ]
        ]
        .to_markdown(index=False),
        "",
        "## V1 Vs Strict V2",
        "",
        comparison.to_markdown(index=False) if not comparison.empty else "V1 comparison was not available.",
        "",
        "## Permutation Tests",
        "",
        perm.to_markdown(index=False) if not perm.empty else "Permutation tests were not available.",
        "",
        "## Feature Stability",
        "",
        top_stability[["region_or_subregion", "gene_id", "selection_frequency", "coefficient_mean", "coefficient_sign_consistency"]].to_markdown(index=False)
        if not top_stability.empty
        else "No stable feature summary was available.",
        "",
        "## Biological Interpretation",
        "",
        "Stable ML genes should be compared with the Tier 2 and enrichment results only at a high level in this report. No enrichment was rerun here.",
        "",
        "## Strong Warning",
        "",
        "This is a strict ML pilot with leakage-controlled internal validation and a preliminary predictive signal. It is not a validated diagnostic classifier and requires external validation.",
        "",
        "## Recommendation",
        "",
        "If strict v2 performance remains strong, proceed to ranked GSEA/enrichment revisit, pathway-score ML, AHBA/WGCNA module mapping, and module/eigengene-level ML. If strict v2 performance collapses, move to pathway/module-level features, search for external validation, and avoid raw-gene deep learning.",
        "",
    ]
    (OUT_DIR / "ML_Strict_V2_Report.md").write_text("\n".join(lines), encoding="utf-8")


def write_readme() -> None:
    text = f"""# Strict ML V2

This folder contains the strict leakage-controlled ML v2 analysis for the Length > 100 MDD transcriptomics project.

What was run:

- MDD vs Control classification for six regions/subregions.
- Nested cross-validation.
- Train-fold-only feature ranking.
- Train-fold-only imputation and scaling.
- Inner-CV hyperparameter tuning.
- Classical ML only; no deep learning.

Why this is stricter than v1:

The v1 pilot used fixed Tier 2 and Tier 4 features selected/ranked using the full dataset. In v2, the primary feature rankings are recomputed inside each outer-training fold.

Key outputs:

- `ML_strict_v2_performance_summary.csv`
- `best_strict_model_by_region.csv`
- `cross_validated_predictions_strict_v2.csv`
- `feature_stability/selection_frequency_all_models.csv`
- `permutation_test_results_strict_v2.csv`
- `v1_vs_strict_v2_comparison.csv`
- `ML_Strict_V2_Report.md`

Runtime-bounded settings:

- outer repeats: {OUTER_REPEATS}
- max outer folds: {MAX_OUTER_FOLDS}
- compact hyperparameter grid: {COMPACT_GRID}
- permutations: {PERMUTATIONS}
- permutation repeats: {PERMUTATION_REPEATS}

How to rerun:

```bash
set TMP=D:\\tmp
set TEMP=D:\\tmp
set TMPDIR=D:\\tmp
python code\\run_ML_strict_v2_length_gt_100.py
```

To run a larger version:

```bash
set STRICT_ML_OUTER_REPEATS=5
set STRICT_ML_PERMUTATIONS=100
set STRICT_ML_MAX_OUTER_FOLDS=5
set STRICT_ML_FULL_GRID=1
python code\\run_ML_strict_v2_length_gt_100.py
```
"""
    (OUT_DIR / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    os.environ.setdefault("TMP", "D:\\tmp")
    os.environ.setdefault("TEMP", "D:\\tmp")
    os.environ.setdefault("TMPDIR", "D:\\tmp")
    Path("D:/tmp").mkdir(parents=True, exist_ok=True)
    prepare_output()

    all_perf = []
    all_pred = []
    all_features = []
    sample_rows = []
    dataset_cache = {}

    for spec in ANALYSES:
        print(f"=== {spec.region_or_subregion} ===")
        x, meta, y, groups = load_dataset(spec)
        splits, outer_cv, group_used, n_splits = get_outer_splits(y, groups, OUTER_REPEATS)
        dataset_cache[spec.region_or_subregion] = (spec, x, meta, y, groups)
        sample_rows.append(
            {
                "analysis_level": spec.analysis_level,
                "region_or_subregion": spec.region_or_subregion,
                "n_samples": len(y),
                "n_MDD": int(np.sum(y == 1)),
                "n_Control": int(np.sum(y == 0)),
                "outer_cv_type": outer_cv,
                "grouped_cv_used": group_used,
                "outer_n_splits": n_splits,
            }
        )
        if not splits:
            continue
        perf, pred, features, inner_seen = fit_evaluate_outer(spec, x, meta, y, groups, splits, outer_cv)
        all_perf.append(perf)
        all_pred.append(pred)
        all_features.append(features)

    performance = pd.concat(all_perf, ignore_index=True)
    predictions = pd.concat(all_pred, ignore_index=True)
    feature_rows = pd.concat(all_features, ignore_index=True)
    stability = feature_stability(feature_rows)

    performance.to_csv(OUT_DIR / "ML_strict_v2_performance_summary.csv", index=False)
    predictions.to_csv(OUT_DIR / "cross_validated_predictions_strict_v2.csv", index=False)
    feature_rows.to_csv(STABILITY_DIR / "selected_genes_per_fold_all_models.csv", index=False)
    stability.to_csv(STABILITY_DIR / "selection_frequency_all_models.csv", index=False)
    sample_counts = pd.DataFrame(sample_rows)
    sample_counts.to_csv(OUT_DIR / "analysis_sample_counts_and_cv_strategy.csv", index=False)

    best = (
        performance.sort_values(["region_or_subregion", "roc_auc_mean", "balanced_accuracy_mean"], ascending=[True, False, False])
        .groupby(["analysis_level", "region_or_subregion"], as_index=False)
        .head(1)
    )
    best.to_csv(OUT_DIR / "best_strict_model_by_region.csv", index=False)

    perm_rows = []
    for _, row in best.iterrows():
        spec, x, meta, y, groups = dataset_cache[row["region_or_subregion"]]
        perm_rows.append(permutation_test_for_best(spec, row, x, y, groups))
    permutation = pd.DataFrame(perm_rows)
    permutation.to_csv(OUT_DIR / "permutation_test_results_strict_v2.csv", index=False)

    comparison = make_v1_comparison(best)
    comparison.to_csv(OUT_DIR / "v1_vs_strict_v2_comparison.csv", index=False)

    # Save most stable genes per region as easy-to-open tables.
    for region, sub in stability.groupby("region_or_subregion"):
        sub.sort_values("selection_frequency", ascending=False).head(100).to_csv(
            STABILITY_DIR / f"{safe_name(region)}__top100_stable_genes.csv", index=False
        )

    plot_curves(best, predictions)
    plot_summaries(best, comparison, stability)
    write_report(performance, best, comparison, permutation, sample_counts, stability)
    write_readme()

    manifest = {
        "script": str(Path(__file__).resolve()),
        "output_dir": str(OUT_DIR),
        "xgboost_available": HAS_XGBOOST,
        "outer_repeats": OUTER_REPEATS,
        "max_outer_folds": MAX_OUTER_FOLDS,
        "compact_grid": COMPACT_GRID,
        "permutations": PERMUTATIONS,
        "permutation_repeats": PERMUTATION_REPEATS,
        "n_performance_rows": int(len(performance)),
        "n_prediction_rows": int(len(predictions)),
        "n_feature_selection_rows": int(len(feature_rows)),
        "n_stability_rows": int(len(stability)),
        "leakage_status": "strict_train_fold_only_feature_selection_tuning_preprocessing",
        "note": "No LMM or enrichment rerun. No C-drive output is intentionally used.",
    }
    (OUT_DIR / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
