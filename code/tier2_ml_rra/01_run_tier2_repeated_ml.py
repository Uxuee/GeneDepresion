#!/usr/bin/env python3
"""Fit repeated Tier2-restricted ML models and export fold rankings."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import sklearn
from sklearn.dummy import DummyClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import StandardScaler

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "code" / "ml_rra"))

from common import (  # noqa: E402
    ANALYSES,
    MODEL_ORDER,
    OUTER_SEEDS,
    ZERO_TOLERANCE,
    atomic_write_csv,
    atomic_write_json,
    build_model,
    compute_metrics,
    load_dataset,
    make_inner_splits,
    make_outer_splits,
    model_parameter_grid,
    needs_scaling,
    safe_name,
    sha256_file,
    stable_seed,
)


TIER2_SOURCE = "outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv"
EXPECTED_SOURCE_SHA256 = "1672ddf95a4f4c9e44ae42e7accc208c311e2e1da8533e1808dcf91bbf7b67a1"
EXPECTED_COUNTS = {
    "Ventral Striatum": 146,
    "Insula": 70,
    "Hippocampal Formation": 127,
    "CG25": 186,
    "dlPFC BA8/9": 72,
    "OFC BA11": 70,
}
TIER2_REGION_NAMES = {"CG25": "Cingulate gyrus 25"}
DUMMY_MODEL = "dummy_most_frequent"


@dataclass
class FixedFeaturePipeline:
    imputer: SimpleImputer
    scaler: StandardScaler | None
    model: Any

    def transform(self, values: np.ndarray) -> np.ndarray:
        transformed = self.imputer.transform(values)
        return self.scaler.transform(transformed) if self.scaler is not None else transformed

    def predict(self, values: np.ndarray) -> np.ndarray:
        return self.model.predict(self.transform(values))

    def decision_score(self, values: np.ndarray) -> np.ndarray:
        transformed = self.transform(values)
        if hasattr(self.model, "predict_proba"):
            probabilities = self.model.predict_proba(transformed)
            classes = list(self.model.classes_)
            return probabilities[:, classes.index(1)] if 1 in classes else np.zeros(len(values))
        return self.model.decision_function(transformed)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=REPO_ROOT)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization",
    )
    parser.add_argument("--regions", nargs="*", help="Optional exact analysis-unit names")
    parser.add_argument("--outer-folds", type=int, default=5)
    parser.add_argument("--outer-repeats", type=int, default=5)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def load_tier2_identities(root: Path) -> pd.DataFrame:
    path = root / TIER2_SOURCE
    if not path.exists():
        raise FileNotFoundError(f"Missing canonical Tier2 source: {path}")
    if sha256_file(path) != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("Canonical Tier2 source checksum differs from the audited value")
    table = pd.read_csv(
        path,
        usecols=["region_or_subregion", "gene_id", "gene_symbol"],
        dtype={"region_or_subregion": str, "gene_id": str, "gene_symbol": str},
    ).rename(columns={"region_or_subregion": "tier2_region"})
    table["gene_symbol"] = table["gene_symbol"].fillna("")
    if table.duplicated(["tier2_region", "gene_id"]).any():
        raise RuntimeError("Duplicate region–gene identities in the canonical Tier2 source")
    if len(table) != 671:
        raise RuntimeError(f"Expected 671 Tier2 region–gene pairs, found {len(table)}")
    observed = table["tier2_region"].value_counts().to_dict()
    remapped = {TIER2_REGION_NAMES.get(region, region): count for region, count in EXPECTED_COUNTS.items()}
    if observed != remapped:
        raise RuntimeError(f"Tier2 regional counts differ: {observed}")
    return table


def fit_fixed_pipeline(
    x_train: np.ndarray,
    y_train: np.ndarray,
    model_name: str,
    params: dict[str, Any],
    seed: int,
) -> FixedFeaturePipeline:
    imputer = SimpleImputer(strategy="median", keep_empty_features=True)
    imputed = imputer.fit_transform(x_train)
    scaler = StandardScaler().fit(imputed) if needs_scaling(model_name) else None
    model_input = scaler.transform(imputed) if scaler is not None else imputed
    model = build_model(model_name, params, seed)
    model.fit(model_input, y_train)
    return FixedFeaturePipeline(imputer, scaler, model)


def tune_model(
    x_train: np.ndarray,
    y_train: np.ndarray,
    model_name: str,
    inner_splits: list[tuple[np.ndarray, np.ndarray]],
    seed: int,
) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    for params in model_parameter_grid(model_name):
        scores = []
        for inner_fold, (fit_idx, valid_idx) in enumerate(inner_splits):
            pipeline = fit_fixed_pipeline(
                x_train[fit_idx],
                y_train[fit_idx],
                model_name,
                params,
                stable_seed(seed, "inner", inner_fold, json.dumps(params, sort_keys=True)),
            )
            scores.append(balanced_accuracy_score(y_train[valid_idx], pipeline.predict(x_train[valid_idx])))
        candidates.append(
            {
                "params": params,
                "mean_balanced_accuracy": float(np.mean(scores)),
                "fold_scores": scores,
            }
        )
    candidates.sort(key=lambda row: (-row["mean_balanced_accuracy"], json.dumps(row["params"], sort_keys=True)))
    return candidates[0]


def importance_table(
    pipeline: FixedFeaturePipeline,
    model_name: str,
    gene_ids: np.ndarray,
    symbol_map: dict[str, str],
) -> pd.DataFrame:
    if model_name == "random_forest":
        raw = np.asarray(pipeline.model.feature_importances_, dtype=float)
        sign = np.full(len(raw), np.nan)
        keep = raw > 0
    else:
        raw = np.asarray(pipeline.model.coef_).reshape(-1).astype(float)
        sign = np.sign(raw)
        keep = np.abs(raw) > ZERO_TOLERANCE if model_name in {"logistic_l1", "elastic_net"} else np.ones(len(raw), dtype=bool)
    result = pd.DataFrame(
        {
            "gene_id": gene_ids[keep].astype(str),
            "gene_symbol": [symbol_map.get(gene, "") for gene in gene_ids[keep].astype(str)],
            "raw_importance": raw[keep],
            "absolute_importance": np.abs(raw[keep]),
            "sign": sign[keep],
        }
    )
    result = result.sort_values(["absolute_importance", "gene_id"], ascending=[False, True], kind="mergesort").reset_index(drop=True)
    result["within_fit_rank"] = np.arange(1, len(result) + 1)
    return result


def read_checkpoint(path: Path, dtype: dict[str, Any] | None = None) -> pd.DataFrame:
    return pd.read_csv(path, dtype=dtype, low_memory=False) if path.exists() else pd.DataFrame()


def concatenate(paths: list[Path], dtype: dict[str, Any] | None = None) -> pd.DataFrame:
    frames = [pd.read_csv(path, dtype=dtype, low_memory=False) for path in paths if path.exists()]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main() -> int:
    args = parse_arguments()
    root = args.input_root.resolve()
    output = args.output_dir.resolve()
    checkpoint_dir = output / "logs/checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    tier2 = load_tier2_identities(root)
    specs = [spec for spec in ANALYSES if not args.regions or spec.region in args.regions]
    unknown = set(args.regions or []) - {spec.region for spec in specs}
    if unknown:
        raise ValueError(f"Unknown regions: {sorted(unknown)}")
    if not specs:
        raise ValueError("No regions selected")
    if not 1 <= args.outer_repeats <= len(OUTER_SEEDS):
        raise ValueError("outer-repeats must be between 1 and 5")
    seeds = OUTER_SEEDS[: args.outer_repeats]
    all_rank_paths: list[Path] = []
    all_prediction_paths: list[Path] = []
    all_performance_paths: list[Path] = []
    unit_records: list[dict[str, Any]] = []
    assembly_start = time.time()

    for spec in specs:
        full_x, metadata, y = load_dataset(root, spec)
        tier2_name = TIER2_REGION_NAMES.get(spec.region, spec.region)
        identities = tier2[tier2["tier2_region"].eq(tier2_name)].copy()
        tier2_ids = identities["gene_id"].astype(str).tolist()
        missing = sorted(set(tier2_ids) - set(full_x.columns.astype(str)))
        if missing:
            raise RuntimeError(f"{spec.region}: {len(missing)} Tier2 genes are absent from expression: {missing[:10]}")
        x_frame = full_x.loc[:, tier2_ids].copy()
        x = x_frame.to_numpy(dtype=float)
        gene_ids = x_frame.columns.to_numpy(dtype=str)
        symbol_map = dict(zip(identities["gene_id"].astype(str), identities["gene_symbol"].fillna("")))
        safe = safe_name(spec.region)
        rank_path = checkpoint_dir / f"{safe}_Tier2_fold_gene_rankings.csv"
        prediction_path = checkpoint_dir / f"{safe}_Tier2_outer_predictions.csv"
        performance_path = checkpoint_dir / f"{safe}_Tier2_outer_performance.csv"
        all_rank_paths.append(rank_path)
        all_prediction_paths.append(prediction_path)
        all_performance_paths.append(performance_path)
        rank_rows = (read_checkpoint(rank_path, {"gene_id": str}) if args.resume else pd.DataFrame()).to_dict("records")
        prediction_rows = (read_checkpoint(prediction_path) if args.resume else pd.DataFrame()).to_dict("records")
        performance_rows = (read_checkpoint(performance_path) if args.resume else pd.DataFrame()).to_dict("records")
        completed = {str(row["fit_id"]) for row in performance_rows}
        splits = list(make_outer_splits(y, n_splits=args.outer_folds, seeds=seeds))
        print(
            f"[{spec.region}] {len(y)} samples, {len(gene_ids)} fixed Tier2 genes, "
            f"{len(splits) * (len(MODEL_ORDER) + 1)} model/reference fits",
            flush=True,
        )
        for split_number, (repeat, fold, outer_seed, train_idx, test_idx) in enumerate(splits, start=1):
            inner_seed = stable_seed(spec.region, outer_seed, fold, "Tier2_inner_split")
            inner_splits, n_inner = make_inner_splits(y[train_idx], inner_seed)
            for model_name in MODEL_ORDER:
                fit_id = f"{safe}__r{repeat + 1}__f{fold + 1}__{model_name}"
                if fit_id in completed:
                    continue
                fit_start = time.time()
                model_seed = stable_seed(spec.region, outer_seed, fold, "Tier2", model_name)
                best = tune_model(x[train_idx], y[train_idx], model_name, inner_splits, model_seed)
                pipeline = fit_fixed_pipeline(x[train_idx], y[train_idx], model_name, best["params"], model_seed)
                y_pred = pipeline.predict(x[test_idx])
                y_score = pipeline.decision_score(x[test_idx])
                metrics = compute_metrics(y[test_idx], y_pred, y_score)
                importance = importance_table(pipeline, model_name, gene_ids, symbol_map)
                shared = {
                    "fit_id": fit_id,
                    "region": spec.region,
                    "repeat": repeat + 1,
                    "fold": fold + 1,
                    "outer_seed": outer_seed,
                    "model": model_name,
                    "candidate_universe": "canonical_Tier2_fixed",
                    "n_tier2_features": len(gene_ids),
                }
                rank_rows.extend({**shared, **row} for row in importance.to_dict("records"))
                prediction_rows.extend(
                    {
                        **shared,
                        "sample_id": x_frame.index[sample_index],
                        "y_true": int(y[sample_index]),
                        "y_pred": int(y_pred[local_index]),
                        "y_score": float(y_score[local_index]),
                    }
                    for local_index, sample_index in enumerate(test_idx)
                )
                performance_rows.append(
                    {
                        **shared,
                        "model_params_json": json.dumps(best["params"], sort_keys=True),
                        "inner_folds": n_inner,
                        "inner_mean_balanced_accuracy": best["mean_balanced_accuracy"],
                        "n_train": len(train_idx),
                        "n_test": len(test_idx),
                        "n_ranked_genes": len(importance),
                        "elapsed_seconds": time.time() - fit_start,
                        "contributes_gene_ranking": True,
                        **metrics,
                    }
                )
                completed.add(fit_id)
                if rank_rows:
                    atomic_write_csv(pd.DataFrame(rank_rows), rank_path)
                atomic_write_csv(pd.DataFrame(prediction_rows), prediction_path)
                atomic_write_csv(pd.DataFrame(performance_rows), performance_path)
                print(
                    f"  split {split_number}/{len(splits)} {model_name}: "
                    f"BA={metrics['balanced_accuracy']:.3f} ranked={len(importance)}",
                    flush=True,
                )

            dummy_fit_id = f"{safe}__r{repeat + 1}__f{fold + 1}__{DUMMY_MODEL}"
            if dummy_fit_id not in completed:
                fit_start = time.time()
                dummy = DummyClassifier(strategy="most_frequent")
                dummy.fit(x[train_idx], y[train_idx])
                y_pred = dummy.predict(x[test_idx])
                probabilities = dummy.predict_proba(x[test_idx])
                classes = list(dummy.classes_)
                y_score = probabilities[:, classes.index(1)] if 1 in classes else np.zeros(len(test_idx))
                metrics = compute_metrics(y[test_idx], y_pred, y_score)
                shared = {
                    "fit_id": dummy_fit_id,
                    "region": spec.region,
                    "repeat": repeat + 1,
                    "fold": fold + 1,
                    "outer_seed": outer_seed,
                    "model": DUMMY_MODEL,
                    "candidate_universe": "performance_reference_only",
                    "n_tier2_features": 0,
                }
                prediction_rows.extend(
                    {
                        **shared,
                        "sample_id": x_frame.index[sample_index],
                        "y_true": int(y[sample_index]),
                        "y_pred": int(y_pred[local_index]),
                        "y_score": float(y_score[local_index]),
                    }
                    for local_index, sample_index in enumerate(test_idx)
                )
                performance_rows.append(
                    {
                        **shared,
                        "model_params_json": json.dumps({"strategy": "most_frequent"}),
                        "inner_folds": 0,
                        "inner_mean_balanced_accuracy": np.nan,
                        "n_train": len(train_idx),
                        "n_test": len(test_idx),
                        "n_ranked_genes": 0,
                        "elapsed_seconds": time.time() - fit_start,
                        "contributes_gene_ranking": False,
                        **metrics,
                    }
                )
                completed.add(dummy_fit_id)
                atomic_write_csv(pd.DataFrame(prediction_rows), prediction_path)
                atomic_write_csv(pd.DataFrame(performance_rows), performance_path)

        region_performance = pd.DataFrame(performance_rows)
        expression_path = root / spec.expression_relpath
        unit_records.append(
            {
                "region": spec.region,
                "tier2_source_region": tier2_name,
                "n_samples": len(y),
                "n_mdd": int(y.sum()),
                "n_control": int((y == 0).sum()),
                "n_tier2_candidates": len(gene_ids),
                "n_tier2_matched_expression": len(gene_ids),
                "n_tier2_missing_expression": 0,
                "expression_path": spec.expression_relpath,
                "expression_sha256": sha256_file(expression_path),
                "cumulative_fit_elapsed_seconds": float(region_performance["elapsed_seconds"].sum()),
            }
        )

    rankings = concatenate(all_rank_paths, {"gene_id": str})
    predictions = concatenate(all_prediction_paths)
    performance = concatenate(all_performance_paths)
    feature_fits_expected = len(specs) * len(seeds) * args.outer_folds * len(MODEL_ORDER)
    dummy_fits_expected = len(specs) * len(seeds) * args.outer_folds
    feature_performance = performance[performance["model"].isin(MODEL_ORDER)]
    dummy_performance = performance[performance["model"].eq(DUMMY_MODEL)]
    if feature_performance["fit_id"].nunique() != feature_fits_expected:
        raise RuntimeError("Feature-model fit count is incomplete")
    if dummy_performance["fit_id"].nunique() != dummy_fits_expected:
        raise RuntimeError("Dummy-reference fit count is incomplete")
    if rankings.duplicated(["fit_id", "gene_id"]).any():
        raise RuntimeError("Duplicate gene within a fold ranking")
    valid_pairs = set(zip(tier2["tier2_region"].replace({"Cingulate gyrus 25": "CG25"}), tier2["gene_id"]))
    observed_pairs = set(zip(rankings["region"], rankings["gene_id"].astype(str)))
    if not observed_pairs.issubset(valid_pairs):
        raise RuntimeError("A non-Tier2 region–gene pair entered a fold ranking")

    atomic_write_csv(rankings, output / "Tier2_fold_gene_rankings.csv")
    atomic_write_csv(predictions, output / "Tier2_outer_predictions.csv")
    atomic_write_csv(performance, output / "Tier2_outer_performance.csv")
    summary = (
        performance.groupby(["region", "model"], as_index=False)
        .agg(
            n_outer_fits=("fit_id", "nunique"),
            contributes_gene_ranking=("contributes_gene_ranking", "first"),
            balanced_accuracy_mean=("balanced_accuracy", "mean"),
            balanced_accuracy_sd=("balanced_accuracy", "std"),
            f1_mean=("f1", "mean"),
            roc_auc_mean=("roc_auc", "mean"),
            precision_mean=("precision", "mean"),
            recall_mean=("recall", "mean"),
        )
    )
    atomic_write_csv(summary, output / "Tier2_model_performance_summary.csv")
    atomic_write_csv(pd.DataFrame(unit_records), output / "Tier2_analysis_unit_manifest.csv")
    manifest = {
        "status": "repeated_ml_complete",
        "primary_prioritization_candidate": True,
        "source_commit": "3e41a8a4ad97781dfcd4a60628181bbf9c143074",
        "branch": "codex/tier2-rra-prioritization",
        "tier2_source": TIER2_SOURCE,
        "tier2_source_sha256": EXPECTED_SOURCE_SHA256,
        "tier2_identity_columns_used": ["region_or_subregion", "gene_id", "gene_symbol"],
        "tier2_candidate_count": len(tier2),
        "analysis_units": unit_records,
        "outer_cv": {"folds": args.outer_folds, "seeds": seeds},
        "inner_cv": {"folds": 3, "fallback_folds": 2, "metric": "balanced_accuracy"},
        "feature_selector": None,
        "top_k_grid": None,
        "models_contributing_gene_rankings": MODEL_ORDER,
        "performance_reference_model": DUMMY_MODEL,
        "n_completed_feature_model_fits": feature_fits_expected,
        "n_completed_dummy_reference_fits": dummy_fits_expected,
        "cumulative_fit_elapsed_seconds": float(performance["elapsed_seconds"].sum()),
        "assembly_runtime_seconds": time.time() - assembly_start,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit_learn": sklearn.__version__,
        "frozen_allgene_consensus_sha256": "0aefec5fcfa100deafa1525c2d31521ae8c00ca890ef4e1b4145b1b861c2c87e",
        "frozen_corrected_ml_sha256": "6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1",
        "interpretation": "within-study ML prioritization of prespecified exploratory Tier2 candidates",
        "performance_warning": "Tier2 was defined using the full cohort; metrics are descriptive/internal and not unbiased diagnostic validation.",
    }
    atomic_write_json(manifest, output / "Tier2_run_manifest.json")
    print(f"Completed {feature_fits_expected} feature-model fits plus {dummy_fits_expected} dummy references")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
