"""Run a first classical ML pilot for Length > 100 MDD transcriptomics.

This script does not rerun the main LMM. It uses existing normalized
log2(GeTMM + 1) expression matrices plus the packaged Tier 2 and Tier 4 gene
sets to create an exploratory ML branch:

  outputs/length_gt_100_elisa_style/ML_pilot_Tier2/

Important interpretation:
  Fixed Tier 2 genes and pre-ranked Tier 4 genes were derived from the full
  dataset, so their classifier performance is exploratory and may be optimistic.
  The strict_univariate_topk mode recomputes simple feature selection inside
  each training fold and is included as a leakage-controlled pilot reference.
"""

from __future__ import annotations

import json
import math
import shutil
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import clone
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
    from xgboost import XGBClassifier

    HAS_XGBOOST = True
except Exception:
    XGBClassifier = None
    HAS_XGBOOST = False


REPO_DIR = Path(__file__).resolve().parents[1]
BASE_DIR = REPO_DIR / "outputs" / "length_gt_100_elisa_style"
OUT_DIR = BASE_DIR / "ML_pilot_Tier2"
MATRIX_OUT_DIR = OUT_DIR / "ml_matrices"
IMPORTANCE_OUT_DIR = OUT_DIR / "feature_importance"
PLOT_DIR = OUT_DIR / "plots"

DE_SET_DIR = BASE_DIR / "DE_candidate_gene_sets"
TIER2_DIR = DE_SET_DIR / "Tier2_working_exploratory"
RANKED_DIR = DE_SET_DIR / "ranked_signatures"
METADATA_DIR = BASE_DIR / "step1_from_raw" / "Length_gt_100" / "01_metadata" / "keep_all_samples"
BROAD_EXPR_DIR = BASE_DIR / "normalized_expression_length_gt_100_by_region"
SUBREGION_EXPR_DIR = BASE_DIR / "frontal_subregion_sensitivity" / "subregion_normalized_expression_matrices"

TOP_K_VALUES = [25, 50, 100, 250, 500]
# Runtime-bounded pilot setting. A first full 20-repeat attempt exceeded the
# desktop execution window, so this pilot uses five repeats and records that
# explicitly in cv_type and the report.
RANDOM_SEEDS = list(range(2026, 2031))
PERMUTATION_SEEDS = list(range(9001, 9021))
STRICT_TOP_K_VALUES = [100]
MIN_FEATURES = 2


@dataclass(frozen=True)
class AnalysisSpec:
    analysis_level: str
    region_or_subregion: str
    package_key: str
    expression_file: Path
    metadata_file: Path
    metadata_filter: Callable[[pd.DataFrame], pd.Series] | None = None
    documentation_only: bool = False


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
        lambda m: m["OriginalRegion"].astype(str).eq("Cingulate gyrus 25 (Cg25)"),
    ),
    AnalysisSpec(
        "frontal_subregion",
        "dlPFC BA8/9",
        "dlPFC_BA8_9",
        SUBREGION_EXPR_DIR / "Dorsolateral_prefrontal_cortex_dlPFC_BA8_9_normalized_GeTMM_log2_length_gt_100.csv",
        METADATA_DIR / "metadata_frontal_cortex.csv",
        lambda m: m["OriginalRegion"].astype(str).eq("Dorsolateral prefrontal cortex (dlPFC; BA8/9)"),
    ),
    AnalysisSpec(
        "frontal_subregion",
        "OFC BA11",
        "OFC_BA11",
        SUBREGION_EXPR_DIR / "Orbitofrontal_OFC_BA11_normalized_GeTMM_log2_length_gt_100.csv",
        METADATA_DIR / "metadata_frontal_cortex.csv",
        lambda m: m["OriginalRegion"].astype(str).eq("Orbitofrontal (OFC; BA11)"),
    ),
    AnalysisSpec(
        "four_group",
        "broad Frontal Cortex",
        "broad_Frontal_Cortex",
        BROAD_EXPR_DIR / "normalized_GeTMM_log2_frontal_cortex.csv",
        METADATA_DIR / "metadata_frontal_cortex.csv",
        documentation_only=True,
    ),
]


def safe_name(value: str) -> str:
    text = "".join(ch if ch.isalnum() else "_" for ch in str(value))
    while "__" in text:
        text = text.replace("__", "_")
    return text.strip("_")


def ensure_clean_output() -> None:
    """Create a fresh ML output branch without touching other result folders."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for subdir in [MATRIX_OUT_DIR, IMPORTANCE_OUT_DIR, PLOT_DIR]:
        if subdir.exists():
            shutil.rmtree(subdir)
        subdir.mkdir(parents=True, exist_ok=True)


def read_expression_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Expression matrix not found: {path}")
    matrix = pd.read_csv(path, index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    return matrix


def read_metadata(spec: AnalysisSpec, sample_ids: list[str]) -> pd.DataFrame:
    if not spec.metadata_file.exists():
        raise FileNotFoundError(f"Metadata file not found: {spec.metadata_file}")
    metadata = pd.read_csv(spec.metadata_file, low_memory=False)
    if "SampleID" not in metadata.columns:
        raise ValueError(f"Metadata is missing SampleID: {spec.metadata_file}")
    metadata["SampleID"] = metadata["SampleID"].astype(str)
    if spec.metadata_filter is not None:
        metadata = metadata.loc[spec.metadata_filter(metadata)].copy()
    metadata = metadata.loc[metadata["SampleID"].isin(sample_ids)].copy()
    metadata = metadata.drop_duplicates(subset=["SampleID"])
    metadata = metadata.set_index("SampleID").reindex(sample_ids).reset_index()
    metadata = metadata.dropna(subset=["Diagnosis"])
    metadata["Diagnosis_binary"] = metadata["Diagnosis"].map(
        lambda x: 1 if str(x).upper() == "MDD" else (0 if str(x).upper() in {"CTRL", "CONTROL"} else np.nan)
    )
    metadata = metadata.dropna(subset=["Diagnosis_binary"])
    metadata["Diagnosis_binary"] = metadata["Diagnosis_binary"].astype(int)
    return metadata


def load_tier2_genes(package_key: str) -> pd.DataFrame:
    path = TIER2_DIR / f"{package_key}__Tier2_all.csv"
    if not path.exists():
        raise FileNotFoundError(f"Tier 2 gene table not found: {path}")
    genes = pd.read_csv(path, dtype={"gene_id": str})
    genes["gene_id"] = genes["gene_id"].astype(str)
    return genes


def load_ranked_genes(package_key: str) -> pd.DataFrame:
    path = RANKED_DIR / f"{package_key}__Tier4_ranked_signature.csv"
    if not path.exists():
        raise FileNotFoundError(f"Tier 4 ranked signature table not found: {path}")
    ranked = pd.read_csv(path, dtype={"gene_id": str})
    ranked["gene_id"] = ranked["gene_id"].astype(str)
    if "signed_standardized_logp_score" not in ranked.columns:
        raise ValueError(f"Ranked signature is missing signed_standardized_logp_score: {path}")
    ranked["abs_rank_score_for_ml"] = pd.to_numeric(ranked["signed_standardized_logp_score"], errors="coerce").abs()
    ranked = ranked.sort_values("abs_rank_score_for_ml", ascending=False)
    return ranked


def make_models(random_state: int) -> dict[str, object]:
    models = {
        "majority_class_baseline": DummyClassifier(strategy="most_frequent"),
        "logistic_l2": LogisticRegression(
            penalty="l2", solver="liblinear", class_weight="balanced", max_iter=5000, random_state=random_state
        ),
        "logistic_l1": LogisticRegression(
            penalty="l1", solver="liblinear", class_weight="balanced", max_iter=5000, random_state=random_state
        ),
        "elastic_net_logistic": LogisticRegression(
            penalty="elasticnet",
            solver="saga",
            l1_ratio=0.5,
            class_weight="balanced",
            max_iter=8000,
            random_state=random_state,
        ),
        "linear_svm": LinearSVC(class_weight="balanced", max_iter=8000, random_state=random_state),
        "random_forest_exploratory": RandomForestClassifier(
            n_estimators=100,
            max_features="sqrt",
            class_weight="balanced_subsample",
            min_samples_leaf=2,
            random_state=random_state,
            n_jobs=-1,
        ),
    }
    if HAS_XGBOOST:
        models["xgboost_exploratory"] = XGBClassifier(
            n_estimators=150,
            max_depth=2,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            eval_metric="logloss",
            random_state=random_state,
            n_jobs=1,
        )
    return models


def make_pipeline(model: object, needs_scaling: bool) -> Pipeline:
    steps = [("imputer", SimpleImputer(strategy="median"))]
    if needs_scaling:
        steps.append(("scaler", StandardScaler()))
    steps.append(("model", model))
    return Pipeline(steps)


def score_from_estimator(estimator: Pipeline, x_test: pd.DataFrame) -> np.ndarray:
    model = estimator.named_steps["model"]
    if hasattr(model, "predict_proba"):
        return estimator.predict_proba(x_test)[:, 1]
    if hasattr(model, "decision_function"):
        return estimator.decision_function(x_test)
    return estimator.predict(x_test)


def fold_metrics(y_true: np.ndarray, y_pred: np.ndarray, y_score: np.ndarray) -> dict[str, float]:
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "roc_auc": roc_auc_score(y_true, y_score) if len(np.unique(y_true)) == 2 else np.nan,
        "pr_auc": average_precision_score(y_true, y_score) if len(np.unique(y_true)) == 2 else np.nan,
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "sensitivity": tp / (tp + fn) if (tp + fn) else np.nan,
        "specificity": tn / (tn + fp) if (tn + fp) else np.nan,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "tp": tp,
    }


def build_cv_splits(y: np.ndarray, groups: pd.Series | None) -> tuple[list[tuple[int, int, np.ndarray, np.ndarray]], str, bool]:
    counts = Counter(y)
    min_class = min(counts.values())
    n_splits = 5 if min_class >= 5 else 3 if min_class >= 3 else 0
    if n_splits < 2:
        return [], "insufficient_samples_for_cv", False

    use_groups = False
    if groups is not None and groups.notna().all() and groups.astype(str).duplicated().any():
        use_groups = True

    splits: list[tuple[int, int, np.ndarray, np.ndarray]] = []
    if use_groups:
        splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEEDS[0])
        for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(np.zeros(len(y)), y, groups.astype(str))):
            splits.append((0, fold_idx, train_idx, test_idx))
        cv_type = f"StratifiedGroupKFold_{n_splits}_fold"
    else:
        for repeat_idx, seed in enumerate(RANDOM_SEEDS):
            splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
            for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(np.zeros(len(y)), y)):
                splits.append((repeat_idx, fold_idx, train_idx, test_idx))
        cv_type = f"RepeatedStratifiedKFold_{n_splits}_fold_5_repeats_runtime_bounded"
    return splits, cv_type, use_groups


def select_univariate_features(x_train: pd.DataFrame, y_train: np.ndarray, candidate_genes: list[str], k: int) -> list[str]:
    """Select top-k genes using only the current training fold."""
    train_block = x_train[candidate_genes]
    imputed = SimpleImputer(strategy="median").fit_transform(train_block)
    scores, _ = f_classif(imputed, y_train)
    scores = np.nan_to_num(scores, nan=-np.inf, posinf=np.inf, neginf=-np.inf)
    order = np.argsort(scores)[::-1]
    return [candidate_genes[i] for i in order[: min(k, len(candidate_genes))]]


def save_matrix_package(
    spec: AnalysisSpec,
    feature_set: str,
    x: pd.DataFrame,
    metadata: pd.DataFrame,
    feature_table: pd.DataFrame,
) -> None:
    prefix = f"{safe_name(spec.analysis_level)}__{safe_name(spec.region_or_subregion)}__{safe_name(feature_set)}"
    x.to_csv(MATRIX_OUT_DIR / f"{prefix}__expression_matrix.csv")
    metadata[["SampleID", "Diagnosis", "Diagnosis_binary"]].to_csv(MATRIX_OUT_DIR / f"{prefix}__labels.csv", index=False)
    metadata.to_csv(MATRIX_OUT_DIR / f"{prefix}__sample_metadata.csv", index=False)
    feature_table.to_csv(MATRIX_OUT_DIR / f"{prefix}__feature_list.csv", index=False)


def summarize_values(values: list[float]) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    if len(arr) == 0 or np.all(np.isnan(arr)):
        return np.nan, np.nan
    return float(np.nanmean(arr)), float(np.nanstd(arr, ddof=1)) if np.sum(~np.isnan(arr)) > 1 else 0.0


def evaluate_fixed_feature_set(
    spec: AnalysisSpec,
    feature_set: str,
    x: pd.DataFrame,
    y: np.ndarray,
    metadata: pd.DataFrame,
    splits: list[tuple[int, int, np.ndarray, np.ndarray]],
    cv_type: str,
    group_used: bool,
    leakage_status: str,
    interpretation_note: str,
) -> tuple[list[dict], pd.DataFrame, pd.DataFrame]:
    rows = []
    prediction_rows = []
    importance_rows = []
    n_mdd = int(np.sum(y == 1))
    n_control = int(np.sum(y == 0))

    for model_name in make_models(RANDOM_SEEDS[0]):
        fold_scores = defaultdict(list)
        confusion_totals = Counter()
        all_y_true = []
        all_y_score = []
        all_sample_ids = []

        for repeat_idx, fold_idx, train_idx, test_idx in splits:
            seed = RANDOM_SEEDS[repeat_idx % len(RANDOM_SEEDS)]
            base_model = make_models(seed)[model_name]
            needs_scaling = model_name not in {"majority_class_baseline", "random_forest_exploratory", "xgboost_exploratory"}
            estimator = make_pipeline(clone(base_model), needs_scaling)
            x_train = x.iloc[train_idx]
            x_test = x.iloc[test_idx]
            y_train = y[train_idx]
            y_test = y[test_idx]
            estimator.fit(x_train, y_train)
            y_pred = estimator.predict(x_test)
            y_score = score_from_estimator(estimator, x_test)
            metrics = fold_metrics(y_test, y_pred, y_score)

            for key, value in metrics.items():
                if key in {"tn", "fp", "fn", "tp"}:
                    confusion_totals[key] += int(value)
                else:
                    fold_scores[key].append(float(value))

            all_y_true.extend(y_test.tolist())
            all_y_score.extend(np.asarray(y_score).ravel().tolist())
            all_sample_ids.extend(metadata.iloc[test_idx]["SampleID"].tolist())

            model = estimator.named_steps["model"]
            if hasattr(model, "coef_"):
                coefs = np.asarray(model.coef_).ravel()
                for gene, coef in zip(x.columns, coefs):
                    if abs(coef) > 1e-12:
                        importance_rows.append(
                            {
                                "analysis_level": spec.analysis_level,
                                "region_or_subregion": spec.region_or_subregion,
                                "feature_set": feature_set,
                                "model": model_name,
                                "repeat": repeat_idx,
                                "fold": fold_idx,
                                "gene_id": gene,
                                "importance_value": coef,
                                "importance_type": "linear_coefficient",
                                "nonzero": True,
                                "direction": "positive_MDD" if coef > 0 else "negative_MDD",
                            }
                        )
            elif hasattr(model, "feature_importances_"):
                for gene, val in zip(x.columns, model.feature_importances_):
                    if val > 0:
                        importance_rows.append(
                            {
                                "analysis_level": spec.analysis_level,
                                "region_or_subregion": spec.region_or_subregion,
                                "feature_set": feature_set,
                                "model": model_name,
                                "repeat": repeat_idx,
                                "fold": fold_idx,
                                "gene_id": gene,
                                "importance_value": val,
                                "importance_type": "tree_importance_exploratory",
                                "nonzero": val > 0,
                                "direction": "not_signed",
                            }
                        )

        summary = {
            "analysis_level": spec.analysis_level,
            "region_or_subregion": spec.region_or_subregion,
            "feature_set": feature_set,
            "n_features": x.shape[1],
            "model": model_name,
            "cv_type": cv_type,
            "n_samples": len(y),
            "n_MDD": n_mdd,
            "n_Control": n_control,
            "group_splitting_used": group_used,
            "leakage_status": leakage_status,
            "interpretation_note": interpretation_note,
            "confusion_tn_total": confusion_totals["tn"],
            "confusion_fp_total": confusion_totals["fp"],
            "confusion_fn_total": confusion_totals["fn"],
            "confusion_tp_total": confusion_totals["tp"],
        }
        for metric in ["accuracy", "balanced_accuracy", "roc_auc", "pr_auc", "f1", "sensitivity", "specificity"]:
            mean, sd = summarize_values(fold_scores[metric])
            summary[f"{metric}_mean"] = mean
            summary[f"{metric}_sd"] = sd
        rows.append(summary)

        prediction_rows.append(
            pd.DataFrame(
                {
                    "analysis_level": spec.analysis_level,
                    "region_or_subregion": spec.region_or_subregion,
                    "feature_set": feature_set,
                    "model": model_name,
                    "SampleID": all_sample_ids,
                    "y_true": all_y_true,
                    "y_score": all_y_score,
                }
            )
        )

    return rows, pd.concat(prediction_rows, ignore_index=True), pd.DataFrame(importance_rows)


def evaluate_strict_topk(
    spec: AnalysisSpec,
    feature_set: str,
    all_gene_x: pd.DataFrame,
    candidate_genes: list[str],
    k: int,
    y: np.ndarray,
    metadata: pd.DataFrame,
    splits: list[tuple[int, int, np.ndarray, np.ndarray]],
    cv_type: str,
    group_used: bool,
) -> tuple[list[dict], pd.DataFrame, pd.DataFrame]:
    rows = []
    prediction_rows = []
    importance_rows = []
    n_mdd = int(np.sum(y == 1))
    n_control = int(np.sum(y == 0))

    for model_name in make_models(RANDOM_SEEDS[0]):
        fold_scores = defaultdict(list)
        confusion_totals = Counter()
        all_y_true = []
        all_y_score = []
        all_sample_ids = []

        for repeat_idx, fold_idx, train_idx, test_idx in splits:
            seed = RANDOM_SEEDS[repeat_idx % len(RANDOM_SEEDS)]
            selected = select_univariate_features(all_gene_x.iloc[train_idx], y[train_idx], candidate_genes, k)
            base_model = make_models(seed)[model_name]
            needs_scaling = model_name not in {"majority_class_baseline", "random_forest_exploratory", "xgboost_exploratory"}
            estimator = make_pipeline(clone(base_model), needs_scaling)
            x_train = all_gene_x.iloc[train_idx][selected]
            x_test = all_gene_x.iloc[test_idx][selected]
            y_train = y[train_idx]
            y_test = y[test_idx]
            estimator.fit(x_train, y_train)
            y_pred = estimator.predict(x_test)
            y_score = score_from_estimator(estimator, x_test)
            metrics = fold_metrics(y_test, y_pred, y_score)

            for key, value in metrics.items():
                if key in {"tn", "fp", "fn", "tp"}:
                    confusion_totals[key] += int(value)
                else:
                    fold_scores[key].append(float(value))

            all_y_true.extend(y_test.tolist())
            all_y_score.extend(np.asarray(y_score).ravel().tolist())
            all_sample_ids.extend(metadata.iloc[test_idx]["SampleID"].tolist())

            model = estimator.named_steps["model"]
            if hasattr(model, "coef_"):
                for gene, coef in zip(selected, np.asarray(model.coef_).ravel()):
                    if abs(coef) > 1e-12:
                        importance_rows.append(
                            {
                                "analysis_level": spec.analysis_level,
                                "region_or_subregion": spec.region_or_subregion,
                                "feature_set": feature_set,
                                "model": model_name,
                                "repeat": repeat_idx,
                                "fold": fold_idx,
                                "gene_id": gene,
                                "importance_value": coef,
                                "importance_type": "strict_fold_linear_coefficient",
                                "nonzero": True,
                                "direction": "positive_MDD" if coef > 0 else "negative_MDD",
                            }
                        )
            elif hasattr(model, "feature_importances_"):
                for gene, val in zip(selected, model.feature_importances_):
                    if val > 0:
                        importance_rows.append(
                            {
                                "analysis_level": spec.analysis_level,
                                "region_or_subregion": spec.region_or_subregion,
                                "feature_set": feature_set,
                                "model": model_name,
                                "repeat": repeat_idx,
                                "fold": fold_idx,
                                "gene_id": gene,
                                "importance_value": val,
                                "importance_type": "strict_fold_tree_importance_exploratory",
                                "nonzero": val > 0,
                                "direction": "not_signed",
                            }
                        )

        summary = {
            "analysis_level": spec.analysis_level,
            "region_or_subregion": spec.region_or_subregion,
            "feature_set": feature_set,
            "n_features": k,
            "model": model_name,
            "cv_type": cv_type,
            "n_samples": len(y),
            "n_MDD": n_mdd,
            "n_Control": n_control,
            "group_splitting_used": group_used,
            "leakage_status": "leakage_controlled_training_fold_univariate_selection",
            "interpretation_note": "Feature selection was recomputed inside each training fold; still a first pilot without external validation.",
            "confusion_tn_total": confusion_totals["tn"],
            "confusion_fp_total": confusion_totals["fp"],
            "confusion_fn_total": confusion_totals["fn"],
            "confusion_tp_total": confusion_totals["tp"],
        }
        for metric in ["accuracy", "balanced_accuracy", "roc_auc", "pr_auc", "f1", "sensitivity", "specificity"]:
            mean, sd = summarize_values(fold_scores[metric])
            summary[f"{metric}_mean"] = mean
            summary[f"{metric}_sd"] = sd
        rows.append(summary)

        prediction_rows.append(
            pd.DataFrame(
                {
                    "analysis_level": spec.analysis_level,
                    "region_or_subregion": spec.region_or_subregion,
                    "feature_set": feature_set,
                    "model": model_name,
                    "SampleID": all_sample_ids,
                    "y_true": all_y_true,
                    "y_score": all_y_score,
                }
            )
        )

    return rows, pd.concat(prediction_rows, ignore_index=True), pd.DataFrame(importance_rows)


def run_permutation_test(
    spec: AnalysisSpec,
    best_row: pd.Series,
    feature_data: dict[str, pd.DataFrame | tuple[pd.DataFrame, list[str], int]],
    y: np.ndarray,
    metadata: pd.DataFrame,
    groups: pd.Series | None,
) -> dict:
    observed_auc = float(best_row["roc_auc_mean"])
    observed_bal = float(best_row["balanced_accuracy_mean"])
    feature_set = str(best_row["feature_set"])
    model_name = str(best_row["model"])
    perm_auc = []
    perm_bal = []

    for seed in PERMUTATION_SEEDS:
        rng = np.random.default_rng(seed)
        y_perm = rng.permutation(y)
        splits, _, _ = build_cv_splits(y_perm, groups)
        if not splits:
            continue
        # Keep permutation time bounded by using the first three CV splits.
        splits = splits[: min(3, len(splits))]
        model_scores_auc = []
        model_scores_bal = []
        for repeat_idx, fold_idx, train_idx, test_idx in splits:
            base_model = make_models(seed)[model_name]
            needs_scaling = model_name not in {"majority_class_baseline", "random_forest_exploratory", "xgboost_exploratory"}
            estimator = make_pipeline(clone(base_model), needs_scaling)
            if feature_set.startswith("strict_univariate_top"):
                all_gene_x, candidate_genes, k = feature_data[feature_set]  # type: ignore[misc]
                selected = select_univariate_features(all_gene_x.iloc[train_idx], y_perm[train_idx], candidate_genes, k)
                x_train = all_gene_x.iloc[train_idx][selected]
                x_test = all_gene_x.iloc[test_idx][selected]
            else:
                x = feature_data[feature_set]  # type: ignore[assignment]
                x_train = x.iloc[train_idx]
                x_test = x.iloc[test_idx]
            estimator.fit(x_train, y_perm[train_idx])
            y_score = score_from_estimator(estimator, x_test)
            y_pred = estimator.predict(x_test)
            model_scores_auc.append(roc_auc_score(y_perm[test_idx], y_score))
            model_scores_bal.append(balanced_accuracy_score(y_perm[test_idx], y_pred))
        perm_auc.append(float(np.mean(model_scores_auc)))
        perm_bal.append(float(np.mean(model_scores_bal)))

    auc_p = (1 + sum(v >= observed_auc for v in perm_auc)) / (1 + len(perm_auc)) if perm_auc else np.nan
    bal_p = (1 + sum(v >= observed_bal for v in perm_bal)) / (1 + len(perm_bal)) if perm_bal else np.nan
    return {
        "analysis_level": spec.analysis_level,
        "region_or_subregion": spec.region_or_subregion,
        "feature_set": feature_set,
        "model": model_name,
        "observed_roc_auc_mean": observed_auc,
        "observed_balanced_accuracy_mean": observed_bal,
        "n_permutations": len(perm_auc),
        "permutation_cv_splits_per_permutation": 3,
        "empirical_p_roc_auc": auc_p,
        "empirical_p_balanced_accuracy": bal_p,
        "permutation_note": "Preliminary 20-label permutation test using three CV splits per permutation to bound runtime.",
    }


def plot_best_curves(predictions: pd.DataFrame, best_rows: pd.DataFrame) -> None:
    for _, row in best_rows.iterrows():
        subset = predictions[
            (predictions["analysis_level"] == row["analysis_level"])
            & (predictions["region_or_subregion"] == row["region_or_subregion"])
            & (predictions["feature_set"] == row["feature_set"])
            & (predictions["model"] == row["model"])
        ]
        if subset.empty or subset["y_true"].nunique() < 2:
            continue
        y_true = subset["y_true"].to_numpy()
        y_score = subset["y_score"].to_numpy()
        prefix = f"{safe_name(row['analysis_level'])}__{safe_name(row['region_or_subregion'])}__best_{safe_name(row['model'])}"

        fpr, tpr, _ = roc_curve(y_true, y_score)
        auc_val = roc_auc_score(y_true, y_score)
        plt.figure(figsize=(5.5, 4.5))
        plt.plot(fpr, tpr, label=f"ROC-AUC={auc_val:.3f}")
        plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
        plt.xlabel("False positive rate")
        plt.ylabel("True positive rate")
        plt.title(f"{row['region_or_subregion']} best ROC")
        plt.legend(loc="lower right")
        plt.tight_layout()
        plt.savefig(PLOT_DIR / f"{prefix}__ROC_curve.png", dpi=180)
        plt.close()

        precision, recall, _ = precision_recall_curve(y_true, y_score)
        pr_auc = average_precision_score(y_true, y_score)
        plt.figure(figsize=(5.5, 4.5))
        plt.plot(recall, precision, label=f"PR-AUC={pr_auc:.3f}")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title(f"{row['region_or_subregion']} best PR")
        plt.legend(loc="lower left")
        plt.tight_layout()
        plt.savefig(PLOT_DIR / f"{prefix}__PR_curve.png", dpi=180)
        plt.close()


def plot_summary(performance: pd.DataFrame, importance_summary: pd.DataFrame) -> None:
    primary = performance[performance["model"] != "majority_class_baseline"].copy()
    primary = primary.sort_values(["region_or_subregion", "roc_auc_mean"], ascending=[True, False])

    plt.figure(figsize=(12, 6))
    top = primary.groupby(["region_or_subregion", "model"], as_index=False)["roc_auc_mean"].max()
    pivot = top.pivot(index="region_or_subregion", columns="model", values="roc_auc_mean")
    pivot.plot(kind="bar", figsize=(12, 5))
    plt.ylabel("Best ROC-AUC across feature sets")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "barplot_best_roc_auc_by_region_model.png", dpi=180)
    plt.close()

    top = primary.groupby(["region_or_subregion", "model"], as_index=False)["balanced_accuracy_mean"].max()
    pivot = top.pivot(index="region_or_subregion", columns="model", values="balanced_accuracy_mean")
    pivot.plot(kind="bar", figsize=(12, 5))
    plt.ylabel("Best balanced accuracy across feature sets")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "barplot_best_balanced_accuracy_by_region_model.png", dpi=180)
    plt.close()

    heat = primary.pivot_table(
        index=["region_or_subregion", "feature_set"], columns="model", values="roc_auc_mean", aggfunc="max"
    )
    plt.figure(figsize=(12, max(6, 0.22 * len(heat))))
    plt.imshow(heat.fillna(np.nan).to_numpy(), aspect="auto", vmin=0, vmax=1, cmap="viridis")
    plt.colorbar(label="ROC-AUC")
    plt.yticks(np.arange(len(heat.index)), [f"{a} | {b}" for a, b in heat.index], fontsize=7)
    plt.xticks(np.arange(len(heat.columns)), heat.columns, rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "heatmap_roc_auc_feature_sets_models.png", dpi=180)
    plt.close()

    coef = importance_summary[
        importance_summary["model"].isin(["elastic_net_logistic", "logistic_l1"])
        & importance_summary["importance_type"].str.contains("coefficient", na=False)
    ].copy()
    if not coef.empty:
        top_coef = coef.sort_values("selection_count", ascending=False).head(30)
        plt.figure(figsize=(9, 8))
        labels = top_coef["gene_symbol"].fillna(top_coef["gene_id"]).astype(str)
        plt.barh(labels[::-1], top_coef["selection_count"].to_numpy()[::-1])
        plt.xlabel("Nonzero coefficient count across folds")
        plt.title("Most stable LASSO / Elastic Net features")
        plt.tight_layout()
        plt.savefig(PLOT_DIR / "top_lasso_elastic_net_feature_stability.png", dpi=180)
        plt.close()


def make_importance_summary(importance: pd.DataFrame, gene_annotations: pd.DataFrame) -> pd.DataFrame:
    if importance.empty:
        return pd.DataFrame()
    summary = (
        importance.groupby(["analysis_level", "region_or_subregion", "feature_set", "model", "gene_id", "importance_type"], as_index=False)
        .agg(
            selection_count=("nonzero", "sum"),
            mean_importance=("importance_value", "mean"),
            sd_importance=("importance_value", "std"),
            positive_count=("importance_value", lambda x: int(np.sum(np.asarray(x) > 0))),
            negative_count=("importance_value", lambda x: int(np.sum(np.asarray(x) < 0))),
        )
    )
    summary["sign_consistency"] = summary[["positive_count", "negative_count"]].max(axis=1) / summary["selection_count"].replace(0, np.nan)
    gene_annotations = gene_annotations.drop_duplicates(subset=["gene_id"])
    return summary.merge(gene_annotations, on="gene_id", how="left")


def write_report(performance: pd.DataFrame, best_rows: pd.DataFrame, permutation: pd.DataFrame, documentation_rows: list[dict]) -> None:
    report_path = OUT_DIR / "ML_Pilot_Tier2_Report.md"
    main_cols = [
        "analysis_level",
        "region_or_subregion",
        "feature_set",
        "n_features",
        "model",
        "balanced_accuracy_mean",
        "roc_auc_mean",
        "pr_auc_mean",
        "leakage_status",
    ]
    best_cols = [
        "region_or_subregion",
        "feature_set",
        "model",
        "balanced_accuracy_mean",
        "roc_auc_mean",
        "pr_auc_mean",
        "leakage_status",
    ]

    lines = [
        "# ML Pilot: Tier 2 and Ranked Length > 100 MDD Features",
        "",
        "## Purpose",
        "",
        "This first ML pilot trains baseline classical models to classify MDD vs Control from normalized log2(GeTMM + 1) expression features. The main LMM was not rerun.",
        "",
        "## Important Interpretation",
        "",
        "Fixed Tier 2 genes were selected using the full dataset, so fixed Tier 2 ML is exploratory and may be optimistic. Pre-ranked Tier 4 top-k genes also use a full-dataset ranking in this pilot. These results are an exploratory classification screen, not a validated diagnostic classifier.",
        "",
        "The stricter `strict_univariate_topk` feature sets recompute a simple univariate diagnosis ranking inside each training fold. This is more leakage-controlled, but it is still only a first pilot without external validation.",
        "",
        "## Regions and Subregions Tested",
        "",
        "- Broad regions: Ventral Striatum, Insula, Hippocampal Formation",
        "- Frontal subregions: Cingulate gyrus 25, dlPFC BA8/9, OFC BA11",
        "- Broad Frontal Cortex was documented only and not used as a primary frontal ML set because it has only one Tier 2 gene.",
        "",
        "## Feature Sets",
        "",
        "- `Tier2_fixed`: Tier 2 exploratory DE candidate genes.",
        "- `Tier4_top25`, `Tier4_top50`, `Tier4_top100`, `Tier4_top250`, `Tier4_top500`: genes ranked by `abs(signed_standardized_logp_score)`.",
        "- `strict_univariate_top100`: top-100 selected inside each training fold as a runtime-bounded leakage-control reference.",
        "",
        "## Models",
        "",
        "Majority-class baseline, L2 logistic regression, L1 logistic regression, Elastic Net logistic regression, linear SVM, and exploratory Random Forest were run. XGBoost was skipped because it was not installed." if not HAS_XGBOOST else "Majority-class baseline, L2 logistic regression, L1 logistic regression, Elastic Net logistic regression, linear SVM, exploratory Random Forest, and exploratory XGBoost were run.",
        "",
        "Scaling and imputation were fit inside each CV fold using training data only.",
        "",
        "A first full 20-repeat CV attempt exceeded the desktop execution window, so this committed pilot uses five repeated stratified folds for non-grouped datasets. Treat this as a first ML screen and expand repeats in the next strict ML stage.",
        "",
        "## Best Model Per Region/Subregion",
        "",
        best_rows[best_cols].to_markdown(index=False),
        "",
        "## Main Performance Table",
        "",
        performance.sort_values(["region_or_subregion", "roc_auc_mean"], ascending=[True, False])
        .head(40)[main_cols]
        .to_markdown(index=False),
        "",
        "## Permutation Tests",
        "",
        permutation.to_markdown(index=False) if not permutation.empty else "Permutation tests were not available.",
        "",
        "## Broad Frontal Cortex Documentation",
        "",
        pd.DataFrame(documentation_rows).to_markdown(index=False) if documentation_rows else "No documentation-only rows were generated.",
        "",
        "## Recommended Next Stage",
        "",
        "- If ML signal is weak, move to pathway/module-level features.",
        "- If ML signal is moderate, run strict nested CV with fold-internal feature selection and hyperparameter tuning.",
        "- If gene-level ML is unstable, use AHBA/WGCNA modules or pathway scores instead of raw genes.",
        "- Do not start deep learning until stable lower-dimensional features or external validation data are available.",
        "",
    ]
    report_path.write_text("\n".join(lines), encoding="utf-8")


def write_readme() -> None:
    text = """# Length > 100 Tier 2 ML Pilot

This folder contains a first classical ML pilot for the Length > 100 MDD analysis.

What was run:

- MDD vs Control classification using normalized log2(GeTMM + 1) expression.
- Fixed Tier 2 exploratory DE features.
- Fixed Tier 4 ranked top-k features.
- A stricter fold-internal univariate top-k feature selection mode.
- Classical models only: majority baseline, logistic regression, linear SVM, Random Forest, and XGBoost only if installed.

Most important outputs:

- `ML_performance_summary.csv`: main performance table.
- `best_model_by_region.csv`: best preliminary model per region/subregion.
- `permutation_test_results.csv`: preliminary label-permutation checks for the best model in each region/subregion.
- `ml_matrices/`: saved expression matrices, labels, metadata, and feature lists.
- `feature_importance/`: coefficient and feature-importance summaries.
- `plots/`: ROC, PR, performance, and feature-stability plots.
- `ML_Pilot_Tier2_Report.md`: human-readable summary.

Exploratory warning:

Tier 2 and Tier 4 rankings were computed from the full dataset. Performance from fixed Tier 2 or fixed Tier 4 feature sets may be optimistic and must not be reported as unbiased predictive performance. The strict univariate mode is a better leakage-control reference, but external validation is still required.

How to rerun:

```bash
python code/run_ML_pilot_Tier2_length_gt_100.py
```
"""
    (OUT_DIR / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    ensure_clean_output()
    all_summary_rows = []
    all_predictions = []
    all_importance = []
    all_gene_annotations = []
    all_permutation_rows = []
    documentation_rows = []

    for spec in ANALYSES:
        print(f"\n=== {spec.analysis_level}: {spec.region_or_subregion} ===")
        expression = read_expression_matrix(spec.expression_file)
        metadata = read_metadata(spec, list(expression.columns))
        common_samples = [sample for sample in expression.columns if sample in set(metadata["SampleID"])]
        expression = expression[common_samples]
        metadata = metadata.set_index("SampleID").loc[common_samples].reset_index()
        x_all = expression.T
        x_all.index = common_samples
        y = metadata["Diagnosis_binary"].to_numpy()
        if len(y) == 0:
            print(f"Skipping {spec.region_or_subregion}: no aligned MDD/Control samples after metadata matching.")
            documentation_rows.append(
                {
                    "analysis_level": spec.analysis_level,
                    "region_or_subregion": spec.region_or_subregion,
                    "n_samples": 0,
                    "n_MDD": 0,
                    "n_Control": 0,
                    "tier2_genes_packaged": np.nan,
                    "tier2_genes_available_in_expression": np.nan,
                    "documentation_only": spec.documentation_only,
                    "note": "No aligned MDD/Control samples after metadata matching.",
                }
            )
            continue
        groups = metadata["BrainID"] if "BrainID" in metadata.columns else None
        n_mdd = int(np.sum(y == 1))
        n_control = int(np.sum(y == 0))

        tier2 = load_tier2_genes(spec.package_key)
        ranked = load_ranked_genes(spec.package_key)
        gene_annotations = pd.concat(
            [
                tier2[["gene_id", "gene_symbol"]].drop_duplicates(),
                ranked[["gene_id", "gene_symbol"]].drop_duplicates(),
            ],
            ignore_index=True,
        ).drop_duplicates(subset=["gene_id"])
        all_gene_annotations.append(gene_annotations)

        tier2_available = [g for g in tier2["gene_id"].tolist() if g in x_all.columns]
        documentation_rows.append(
            {
                "analysis_level": spec.analysis_level,
                "region_or_subregion": spec.region_or_subregion,
                "n_samples": len(y),
                "n_MDD": n_mdd,
                "n_Control": n_control,
                "tier2_genes_packaged": len(tier2),
                "tier2_genes_available_in_expression": len(tier2_available),
                "documentation_only": spec.documentation_only,
                "note": "Broad Frontal Cortex was not used as a primary frontal ML set." if spec.documentation_only else "",
            }
        )

        if spec.documentation_only:
            continue

        splits, cv_type, group_used = build_cv_splits(y, groups)
        if not splits:
            print("Skipping ML: insufficient samples for CV.")
            continue

        feature_data: dict[str, pd.DataFrame | tuple[pd.DataFrame, list[str], int]] = {}

        fixed_feature_sets: list[tuple[str, list[str], pd.DataFrame, str, str]] = []
        fixed_feature_sets.append(
            (
                "Tier2_fixed",
                tier2_available,
                tier2[tier2["gene_id"].isin(tier2_available)].copy(),
                "exploratory_full_dataset_Tier2_selection",
                "Tier 2 was selected using the full dataset; performance may be optimistic.",
            )
        )
        for k in TOP_K_VALUES:
            genes = [g for g in ranked["gene_id"].head(k).tolist() if g in x_all.columns]
            feature_table = ranked[ranked["gene_id"].isin(genes)].copy()
            fixed_feature_sets.append(
                (
                    f"Tier4_top{k}",
                    genes,
                    feature_table,
                    "precomputed_full_dataset_ranking_preliminary",
                    "Tier 4 ranking was precomputed using the full dataset; performance is preliminary.",
                )
            )

        for feature_set, genes, feature_table, leakage_status, note in fixed_feature_sets:
            if len(genes) < MIN_FEATURES:
                print(f"Skipping {feature_set}: only {len(genes)} available features.")
                continue
            x = x_all[genes].copy()
            feature_data[feature_set] = x
            save_matrix_package(spec, feature_set, x, metadata, feature_table)
            rows, pred, imp = evaluate_fixed_feature_set(
                spec, feature_set, x, y, metadata, splits, cv_type, group_used, leakage_status, note
            )
            all_summary_rows.extend(rows)
            all_predictions.append(pred)
            if not imp.empty:
                all_importance.append(imp)

        candidate_genes = [g for g in x_all.columns if g in set(ranked["gene_id"])]
        for k in STRICT_TOP_K_VALUES:
            if len(candidate_genes) < MIN_FEATURES:
                continue
            feature_set = f"strict_univariate_top{k}"
            feature_data[feature_set] = (x_all, candidate_genes, k)
            # Save a placeholder feature list because exact genes are selected fold-by-fold.
            pd.DataFrame(
                {
                    "feature_selection_mode": [feature_set],
                    "candidate_gene_universe": [len(candidate_genes)],
                    "n_features_selected_inside_each_training_fold": [min(k, len(candidate_genes))],
                }
            ).to_csv(
                MATRIX_OUT_DIR
                / f"{safe_name(spec.analysis_level)}__{safe_name(spec.region_or_subregion)}__{safe_name(feature_set)}__feature_list.csv",
                index=False,
            )
            rows, pred, imp = evaluate_strict_topk(spec, feature_set, x_all, candidate_genes, k, y, metadata, splits, cv_type, group_used)
            all_summary_rows.extend(rows)
            all_predictions.append(pred)
            if not imp.empty:
                all_importance.append(imp)

        region_perf = pd.DataFrame([r for r in all_summary_rows if r["region_or_subregion"] == spec.region_or_subregion])
        eligible = region_perf[region_perf["model"] != "majority_class_baseline"].copy()
        if not eligible.empty:
            best = eligible.sort_values(["roc_auc_mean", "balanced_accuracy_mean"], ascending=False).iloc[0]
            all_permutation_rows.append(run_permutation_test(spec, best, feature_data, y, metadata, groups))

    performance = pd.DataFrame(all_summary_rows)
    predictions = pd.concat(all_predictions, ignore_index=True) if all_predictions else pd.DataFrame()
    importance = pd.concat(all_importance, ignore_index=True) if all_importance else pd.DataFrame()
    gene_annotations = pd.concat(all_gene_annotations, ignore_index=True).drop_duplicates(subset=["gene_id"])
    importance_summary = make_importance_summary(importance, gene_annotations)
    permutation = pd.DataFrame(all_permutation_rows)

    if performance.empty:
        raise RuntimeError("No ML performance rows were generated.")

    ordered_cols = [
        "analysis_level",
        "region_or_subregion",
        "feature_set",
        "n_features",
        "model",
        "cv_type",
        "n_samples",
        "n_MDD",
        "n_Control",
        "accuracy_mean",
        "accuracy_sd",
        "balanced_accuracy_mean",
        "balanced_accuracy_sd",
        "roc_auc_mean",
        "roc_auc_sd",
        "pr_auc_mean",
        "pr_auc_sd",
        "f1_mean",
        "f1_sd",
        "sensitivity_mean",
        "specificity_mean",
        "leakage_status",
        "interpretation_note",
        "group_splitting_used",
        "confusion_tn_total",
        "confusion_fp_total",
        "confusion_fn_total",
        "confusion_tp_total",
    ]
    performance = performance[[c for c in ordered_cols if c in performance.columns]]
    performance.to_csv(OUT_DIR / "ML_performance_summary.csv", index=False)
    predictions.to_csv(OUT_DIR / "cross_validated_predictions.csv", index=False)
    importance.to_csv(IMPORTANCE_OUT_DIR / "fold_level_feature_importance.csv", index=False)
    importance_summary.to_csv(IMPORTANCE_OUT_DIR / "feature_importance_summary.csv", index=False)
    pd.DataFrame(documentation_rows).to_csv(OUT_DIR / "analysis_documentation_and_sample_counts.csv", index=False)
    permutation.to_csv(OUT_DIR / "permutation_test_results.csv", index=False)

    best_rows = (
        performance[performance["model"] != "majority_class_baseline"]
        .sort_values(["region_or_subregion", "roc_auc_mean", "balanced_accuracy_mean"], ascending=[True, False, False])
        .groupby(["analysis_level", "region_or_subregion"], as_index=False)
        .head(1)
    )
    best_rows.to_csv(OUT_DIR / "best_model_by_region.csv", index=False)

    plot_best_curves(predictions, best_rows)
    plot_summary(performance, importance_summary)
    write_report(performance, best_rows, permutation, documentation_rows)
    write_readme()

    manifest = {
        "script": str(Path(__file__).resolve()),
        "output_dir": str(OUT_DIR),
        "xgboost_available": HAS_XGBOOST,
        "n_performance_rows": int(len(performance)),
        "n_prediction_rows": int(len(predictions)),
        "n_feature_importance_rows": int(len(importance)),
        "n_permutation_rows": int(len(permutation)),
        "note": "No LMM was rerun. Fixed Tier 2 and Tier 4 ML is exploratory; strict_univariate_topk is fold-internal.",
    }
    (OUT_DIR / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
