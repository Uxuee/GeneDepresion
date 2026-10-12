from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC


OUTER_SEEDS = [2026, 2027, 2028, 2029, 2030]
SELECTORS = ["anova_f", "absolute_standardized_mean_difference"]
K_GRID = [25, 50, 100, 250, 500]
MODEL_ORDER = ["logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest"]
ZERO_TOLERANCE = 1e-12


@dataclass(frozen=True)
class AnalysisSpec:
    region: str
    expression_relpath: str
    original_region: str


ANALYSES = [
    AnalysisSpec(
        "Ventral Striatum",
        "outputs/length_gt_100_elisa_style/normalized_expression_length_gt_100_by_region/normalized_GeTMM_log2_ventral_striatum.csv",
        "Nucleus Accumbens (Nac)",
    ),
    AnalysisSpec(
        "Insula",
        "outputs/length_gt_100_elisa_style/normalized_expression_length_gt_100_by_region/normalized_GeTMM_log2_insula.csv",
        "Anterior Insula (aINS)",
    ),
    AnalysisSpec(
        "Hippocampal Formation",
        "outputs/length_gt_100_elisa_style/normalized_expression_length_gt_100_by_region/normalized_GeTMM_log2_hippocampal_formation.csv",
        "Subiculum (Sub)",
    ),
    AnalysisSpec(
        "CG25",
        "outputs/length_gt_100_elisa_style/frontal_subregion_sensitivity/subregion_normalized_expression_matrices/Cingulate_gyrus_25_Cg25_normalized_GeTMM_log2_length_gt_100.csv",
        "Cingulate gyrus 25 (Cg25)",
    ),
    AnalysisSpec(
        "dlPFC BA8/9",
        "outputs/length_gt_100_elisa_style/frontal_subregion_sensitivity/subregion_normalized_expression_matrices/Dorsolateral_prefrontal_cortex_dlPFC_BA8_9_normalized_GeTMM_log2_length_gt_100.csv",
        "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
    ),
    AnalysisSpec(
        "OFC BA11",
        "outputs/length_gt_100_elisa_style/frontal_subregion_sensitivity/subregion_normalized_expression_matrices/Orbitofrontal_OFC_BA11_normalized_GeTMM_log2_length_gt_100.csv",
        "Orbitofrontal (OFC; BA11)",
    ),
]


def safe_name(value: str) -> str:
    return "_".join(filter(None, "".join(c if c.isalnum() else "_" for c in value).split("_")))


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def stable_seed(*parts: Any) -> int:
    text = "|".join(map(str, parts)).encode("utf-8")
    return int.from_bytes(hashlib.sha256(text).digest()[:4], "big") % (2**31 - 1)


def read_gene_symbols(input_root: Path) -> dict[str, str]:
    path = input_root / "outputs/length_gt_100_elisa_style/standardized_effect_size_4groups/LMM_all_genes_with_standardized_effect_size_4groups.csv"
    if not path.exists():
        return {}
    columns = pd.read_csv(path, nrows=0).columns
    wanted = [c for c in ["gene_id", "gene_symbol"] if c in columns]
    if len(wanted) != 2:
        return {}
    table = pd.read_csv(path, usecols=wanted, dtype={"gene_id": str}, low_memory=False)
    table = table.dropna(subset=["gene_id"]).drop_duplicates("gene_id")
    return dict(zip(table["gene_id"].astype(str), table["gene_symbol"].fillna("").astype(str)))


def load_dataset(input_root: Path, spec: AnalysisSpec) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    expression_path = input_root / spec.expression_relpath
    metadata_path = input_root / "outputs/length_gt_100_clean_pipeline/01_inputs/sample_metadata_used.csv"
    if not expression_path.exists():
        raise FileNotFoundError(f"Missing normalized expression matrix: {expression_path}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Missing metadata: {metadata_path}")

    expression = pd.read_csv(expression_path, index_col=0, low_memory=False)
    expression.index = expression.index.astype(str)
    expression.columns = expression.columns.astype(str)
    expression = expression.apply(pd.to_numeric, errors="coerce")
    if expression.index.duplicated().any():
        raise ValueError(f"Duplicate gene IDs in {expression_path}")

    metadata = pd.read_csv(metadata_path, low_memory=False)
    metadata["SampleID"] = metadata["SampleID"].astype(str)
    metadata = metadata[metadata["OriginalRegion"].astype(str).eq(spec.original_region)].copy()
    metadata = metadata.drop_duplicates("SampleID")
    metadata["label"] = metadata["Diagnosis"].astype(str).str.upper().map({"MDD": 1, "CTRL": 0, "CONTROL": 0})
    metadata = metadata.dropna(subset=["label"])
    metadata["label"] = metadata["label"].astype(int)

    meta_ids = set(metadata["SampleID"])
    sample_ids = [sample for sample in expression.columns if sample in meta_ids]
    if not sample_ids:
        raise ValueError(f"No expression/metadata sample overlap for {spec.region}")
    metadata = metadata.set_index("SampleID").loc[sample_ids].reset_index()
    expression = expression[sample_ids].T
    expression.index = sample_ids
    expression.columns = expression.columns.astype(str)
    labels = metadata["label"].to_numpy(dtype=int)
    if set(np.unique(labels)) != {0, 1}:
        raise ValueError(f"Both classes are required for {spec.region}")
    return expression, metadata, labels


def selector_scores(x_train: np.ndarray, y_train: np.ndarray, method: str) -> np.ndarray:
    if method == "anova_f":
        scores, _ = f_classif(x_train, y_train)
        return np.nan_to_num(scores, nan=-np.inf, posinf=np.inf, neginf=-np.inf)
    if method == "absolute_standardized_mean_difference":
        case = x_train[y_train == 1]
        control = x_train[y_train == 0]
        difference = np.nanmean(case, axis=0) - np.nanmean(control, axis=0)
        sd = np.nanstd(x_train, axis=0, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            scores = np.abs(difference / sd)
        return np.nan_to_num(scores, nan=-np.inf, posinf=np.inf, neginf=-np.inf)
    raise ValueError(f"Unknown selector: {method}")


def ranked_feature_indices(x_train: np.ndarray, y_train: np.ndarray, method: str) -> np.ndarray:
    scores = selector_scores(x_train, y_train, method)
    # Stable secondary key is the original feature position.
    return np.lexsort((np.arange(len(scores)), -scores))


def model_parameter_grid(model_name: str) -> list[dict[str, Any]]:
    c_values = [0.01, 0.1, 1.0]
    if model_name in {"logistic_l1", "logistic_l2", "linear_svm"}:
        return [{"C": c} for c in c_values]
    if model_name == "elastic_net":
        return [{"C": c, "l1_ratio": 0.5} for c in c_values]
    if model_name == "random_forest":
        return [{"max_depth": depth, "min_samples_leaf": 2} for depth in [3, None]]
    raise ValueError(model_name)


def build_model(model_name: str, params: dict[str, Any], seed: int):
    if model_name == "logistic_l1":
        return LogisticRegression(
            penalty="l1", solver="liblinear", class_weight="balanced", max_iter=4000, C=params["C"], random_state=seed
        )
    if model_name == "elastic_net":
        return LogisticRegression(
            penalty="elasticnet",
            solver="saga",
            class_weight="balanced",
            max_iter=5000,
            C=params["C"],
            l1_ratio=params["l1_ratio"],
            random_state=seed,
        )
    if model_name == "logistic_l2":
        return LogisticRegression(
            penalty="l2", solver="liblinear", class_weight="balanced", max_iter=4000, C=params["C"], random_state=seed
        )
    if model_name == "linear_svm":
        return LinearSVC(class_weight="balanced", max_iter=6000, C=params["C"], random_state=seed)
    if model_name == "random_forest":
        return RandomForestClassifier(
            n_estimators=200,
            max_depth=params["max_depth"],
            min_samples_leaf=params["min_samples_leaf"],
            class_weight="balanced_subsample",
            random_state=seed,
            n_jobs=-1,
        )
    raise ValueError(model_name)


def needs_scaling(model_name: str) -> bool:
    return model_name != "random_forest"


@dataclass
class FittedTrainingPipeline:
    imputer: SimpleImputer
    selector_method: str
    selected_indices: np.ndarray
    scaler: StandardScaler | None
    model: Any

    def transform(self, x: np.ndarray) -> np.ndarray:
        values = self.imputer.transform(x)[:, self.selected_indices]
        return self.scaler.transform(values) if self.scaler is not None else values

    def predict(self, x: np.ndarray) -> np.ndarray:
        return self.model.predict(self.transform(x))

    def decision_score(self, x: np.ndarray) -> np.ndarray:
        values = self.transform(x)
        if hasattr(self.model, "predict_proba"):
            return self.model.predict_proba(values)[:, 1]
        return self.model.decision_function(values)


def fit_training_pipeline(
    x_train: np.ndarray,
    y_train: np.ndarray,
    selector_method: str,
    k: int,
    model_name: str,
    model_params: dict[str, Any],
    seed: int,
) -> FittedTrainingPipeline:
    imputer = SimpleImputer(strategy="median", keep_empty_features=True)
    x_imputed = imputer.fit_transform(x_train)
    order = ranked_feature_indices(x_imputed, y_train, selector_method)
    selected = order[: min(int(k), x_imputed.shape[1])]
    x_selected = x_imputed[:, selected]
    scaler = StandardScaler().fit(x_selected) if needs_scaling(model_name) else None
    x_model = scaler.transform(x_selected) if scaler is not None else x_selected
    model = build_model(model_name, model_params, seed)
    model.fit(x_model, y_train)
    return FittedTrainingPipeline(imputer, selector_method, selected, scaler, model)


def inner_tune(
    x_train: np.ndarray,
    y_train: np.ndarray,
    model_name: str,
    inner_splits: list[tuple[np.ndarray, np.ndarray]],
    seed: int,
    selectors: Iterable[str] = SELECTORS,
    k_grid: Iterable[int] = K_GRID,
) -> dict[str, Any]:
    candidates = []
    for selector_method in selectors:
        for k in k_grid:
            k_effective = min(int(k), x_train.shape[1])
            for params in model_parameter_grid(model_name):
                fold_scores = []
                for inner_fold, (fit_idx, valid_idx) in enumerate(inner_splits):
                    pipeline = fit_training_pipeline(
                        x_train[fit_idx],
                        y_train[fit_idx],
                        selector_method,
                        k_effective,
                        model_name,
                        params,
                        stable_seed(seed, "inner", inner_fold, selector_method, k_effective, json.dumps(params, sort_keys=True)),
                    )
                    predictions = pipeline.predict(x_train[valid_idx])
                    fold_scores.append(balanced_accuracy_score(y_train[valid_idx], predictions))
                candidates.append(
                    {
                        "selector": selector_method,
                        "k": k_effective,
                        "params": params,
                        "mean_balanced_accuracy": float(np.mean(fold_scores)),
                        "fold_scores": fold_scores,
                    }
                )
    # Deterministic tie-breaking prefers smaller k, selector order, and simpler serialized parameters.
    selector_order = {name: idx for idx, name in enumerate(SELECTORS)}
    candidates.sort(
        key=lambda row: (
            -row["mean_balanced_accuracy"],
            row["k"],
            selector_order[row["selector"]],
            json.dumps(row["params"], sort_keys=True),
        )
    )
    return candidates[0]


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray, y_score: np.ndarray) -> dict[str, float]:
    roc = float("nan")
    if len(np.unique(y_true)) == 2:
        roc = float(roc_auc_score(y_true, y_score))
    return {
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "roc_auc": roc,
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
    }


def extract_importance(
    pipeline: FittedTrainingPipeline,
    model_name: str,
    gene_ids: np.ndarray,
    symbol_map: dict[str, str],
) -> pd.DataFrame:
    selected_gene_ids = gene_ids[pipeline.selected_indices].astype(str)
    if model_name == "random_forest":
        raw = np.asarray(pipeline.model.feature_importances_, dtype=float)
        sign = np.full(len(raw), np.nan)
        keep = raw > 0
    else:
        raw = np.asarray(pipeline.model.coef_).reshape(-1).astype(float)
        sign = np.sign(raw)
        keep = np.abs(raw) > ZERO_TOLERANCE if model_name in {"logistic_l1", "elastic_net"} else np.ones(len(raw), dtype=bool)
    frame = pd.DataFrame(
        {
            "gene_id": selected_gene_ids[keep],
            "gene_symbol": [symbol_map.get(gene, "") for gene in selected_gene_ids[keep]],
            "raw_importance": raw[keep],
            "absolute_importance": np.abs(raw[keep]),
            "sign": sign[keep],
        }
    )
    frame = frame.sort_values(["absolute_importance", "gene_id"], ascending=[False, True], kind="mergesort").reset_index(drop=True)
    frame["within_fit_rank"] = np.arange(1, len(frame) + 1)
    return frame


def make_outer_splits(y: np.ndarray, n_splits: int = 5, seeds: Iterable[int] = OUTER_SEEDS):
    if min(np.bincount(y)) < n_splits:
        raise ValueError(f"At least {n_splits} samples per class are required")
    for repeat, seed in enumerate(seeds):
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        for fold, (train_idx, test_idx) in enumerate(splitter.split(np.zeros(len(y)), y)):
            yield repeat, fold, seed, train_idx, test_idx


def make_inner_splits(y_train: np.ndarray, seed: int) -> tuple[list[tuple[np.ndarray, np.ndarray]], int]:
    min_class = int(min(np.bincount(y_train)))
    n_splits = 3 if min_class >= 3 else 2
    if min_class < 2:
        raise ValueError("At least two samples per class are required for inner CV")
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(splitter.split(np.zeros(len(y_train)), y_train)), n_splits


def atomic_write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def atomic_write_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    temporary.replace(path)
