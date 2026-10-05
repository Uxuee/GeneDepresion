#!/usr/bin/env python3
"""Run repeated nested CV and export one gene ranking per outer-fit/model."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn

from common import (
    ANALYSES,
    K_GRID,
    MODEL_ORDER,
    OUTER_SEEDS,
    SELECTORS,
    atomic_write_csv,
    atomic_write_json,
    compute_metrics,
    extract_importance,
    fit_training_pipeline,
    inner_tune,
    load_dataset,
    make_inner_splits,
    make_outer_splits,
    read_gene_symbols,
    safe_name,
    sha256_file,
    stable_seed,
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus",
    )
    parser.add_argument("--regions", nargs="*", help="Optional exact region names")
    parser.add_argument("--outer-folds", type=int, default=5)
    parser.add_argument("--outer-repeats", type=int, default=5)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def read_checkpoint(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False) if path.exists() else pd.DataFrame()


def concat_checkpoints(paths: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(path, low_memory=False) for path in paths if path.exists()]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main() -> int:
    args = arguments()
    root = args.input_root.resolve()
    output = args.output_dir.resolve()
    checkpoint_dir = output / "logs/checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    selected_specs = [spec for spec in ANALYSES if not args.regions or spec.region in args.regions]
    missing = set(args.regions or []) - {spec.region for spec in selected_specs}
    if missing:
        raise ValueError(f"Unknown regions: {sorted(missing)}")
    if not selected_specs:
        raise ValueError("No analysis units selected")
    if not 1 <= args.outer_repeats <= len(OUTER_SEEDS):
        raise ValueError("outer-repeats must be between 1 and 5")

    seeds = OUTER_SEEDS[: args.outer_repeats]
    symbols = read_gene_symbols(root)
    input_records: list[dict] = []
    all_rank_paths: list[Path] = []
    all_prediction_paths: list[Path] = []
    all_performance_paths: list[Path] = []
    start_all = time.time()

    for spec in selected_specs:
        region_start = time.time()
        x_frame, metadata, y = load_dataset(root, spec)
        x = x_frame.to_numpy(dtype=float)
        gene_ids = x_frame.columns.to_numpy(dtype=str)
        expression_path = root / spec.expression_relpath
        metadata_path = root / "outputs/length_gt_100_clean_pipeline/01_inputs/sample_metadata_used.csv"
        safe = safe_name(spec.region)
        rank_path = checkpoint_dir / f"{safe}_fold_gene_rankings.csv"
        prediction_path = checkpoint_dir / f"{safe}_outer_fold_predictions.csv"
        performance_path = checkpoint_dir / f"{safe}_outer_fold_performance.csv"
        all_rank_paths.append(rank_path)
        all_prediction_paths.append(prediction_path)
        all_performance_paths.append(performance_path)

        ranks = read_checkpoint(rank_path) if args.resume else pd.DataFrame()
        predictions = read_checkpoint(prediction_path) if args.resume else pd.DataFrame()
        performances = read_checkpoint(performance_path) if args.resume else pd.DataFrame()
        completed = set(performances.get("fit_id", pd.Series(dtype=str)).astype(str))
        rank_rows = ranks.to_dict("records")
        prediction_rows = predictions.to_dict("records")
        performance_rows = performances.to_dict("records")

        splits = list(make_outer_splits(y, n_splits=args.outer_folds, seeds=seeds))
        print(
            f"[{spec.region}] {len(y)} samples ({int(y.sum())} MDD/{int((y == 0).sum())} control), "
            f"{len(gene_ids)} genes, {len(splits) * len(MODEL_ORDER)} fits",
            flush=True,
        )
        for split_number, (repeat, fold, outer_seed, train_idx, test_idx) in enumerate(splits, start=1):
            inner_seed = stable_seed(spec.region, outer_seed, fold, "inner_split")
            inner_splits, n_inner = make_inner_splits(y[train_idx], inner_seed)
            for model_name in MODEL_ORDER:
                fit_id = f"{safe}__r{repeat + 1}__f{fold + 1}__{model_name}"
                if fit_id in completed:
                    continue
                fit_start = time.time()
                model_seed = stable_seed(spec.region, outer_seed, fold, model_name)
                best = inner_tune(x[train_idx], y[train_idx], model_name, inner_splits, model_seed)
                pipeline = fit_training_pipeline(
                    x[train_idx],
                    y[train_idx],
                    best["selector"],
                    best["k"],
                    model_name,
                    best["params"],
                    model_seed,
                )
                y_pred = pipeline.predict(x[test_idx])
                y_score = pipeline.decision_score(x[test_idx])
                metrics = compute_metrics(y[test_idx], y_pred, y_score)
                importance = extract_importance(pipeline, model_name, gene_ids, symbols)
                shared = {
                    "fit_id": fit_id,
                    "region": spec.region,
                    "repeat": repeat + 1,
                    "fold": fold + 1,
                    "outer_seed": outer_seed,
                    "model": model_name,
                    "selected_selector": best["selector"],
                    "selected_k": best["k"],
                }
                for row in importance.to_dict("records"):
                    rank_rows.append({**shared, **row})
                for local_index, sample_index in enumerate(test_idx):
                    prediction_rows.append(
                        {
                            **shared,
                            "sample_id": x_frame.index[sample_index],
                            "y_true": int(y[sample_index]),
                            "y_pred": int(y_pred[local_index]),
                            "y_score": float(y_score[local_index]),
                        }
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
                        **metrics,
                    }
                )
                atomic_write_csv(pd.DataFrame(rank_rows), rank_path)
                atomic_write_csv(pd.DataFrame(prediction_rows), prediction_path)
                atomic_write_csv(pd.DataFrame(performance_rows), performance_path)
                completed.add(fit_id)
                print(
                    f"  split {split_number}/{len(splits)} {model_name}: "
                    f"{best['selector']} k={best['k']} BA={metrics['balanced_accuracy']:.3f} "
                    f"ranked={len(importance)} ({time.time() - fit_start:.1f}s)",
                    flush=True,
                )

        input_records.append(
            {
                "region": spec.region,
                "original_region": spec.original_region,
                "expression_path": spec.expression_relpath,
                "expression_sha256": sha256_file(expression_path),
                "metadata_path": str(metadata_path.relative_to(root)).replace("\\", "/"),
                "metadata_sha256": sha256_file(metadata_path),
                "n_samples": len(y),
                "n_mdd": int(y.sum()),
                "n_control": int((y == 0).sum()),
                "n_genes": len(gene_ids),
                "n_unique_brain_ids": int(metadata["BrainID"].nunique()) if "BrainID" in metadata else None,
                "elapsed_seconds": time.time() - region_start,
            }
        )

    rankings = concat_checkpoints(all_rank_paths)
    predictions = concat_checkpoints(all_prediction_paths)
    performance = concat_checkpoints(all_performance_paths)
    required_fits = len(selected_specs) * len(seeds) * args.outer_folds * len(MODEL_ORDER)
    if performance["fit_id"].nunique() != required_fits:
        raise RuntimeError(f"Expected {required_fits} completed fits, found {performance['fit_id'].nunique()}")
    if rankings.duplicated(["fit_id", "gene_id"]).any():
        raise RuntimeError("Duplicate genes found within at least one fold ranking")

    atomic_write_csv(rankings, output / "fold_gene_rankings.csv")
    atomic_write_csv(predictions, output / "outer_fold_predictions.csv")
    atomic_write_csv(performance, output / "outer_fold_performance.csv")
    summary = (
        performance.groupby(["region", "model"], as_index=False)
        .agg(
            n_outer_fits=("fit_id", "nunique"),
            balanced_accuracy_mean=("balanced_accuracy", "mean"),
            balanced_accuracy_sd=("balanced_accuracy", "std"),
            f1_mean=("f1", "mean"),
            roc_auc_mean=("roc_auc", "mean"),
            precision_mean=("precision", "mean"),
            recall_mean=("recall", "mean"),
        )
    )
    atomic_write_csv(summary, output / "model_performance_summary.csv")
    atomic_write_csv(pd.DataFrame(input_records), output / "analysis_unit_manifest.csv")
    manifest = {
        "status": "nested_ml_complete",
        "experimental": True,
        "source_commit": "11d0c63167090aef2c6aa014859bca8e27128eed",
        "branch": "codex/rra-prioritization-redesign",
        "analysis_units": input_records,
        "outer_cv": {"folds": args.outer_folds, "seeds": seeds},
        "inner_cv": {"folds": 3, "fallback_folds": 2, "metric": "balanced_accuracy"},
        "selectors": SELECTORS,
        "k_grid": K_GRID,
        "models": MODEL_ORDER,
        "n_completed_outer_model_fits": required_fits,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit_learn": sklearn.__version__,
        "runtime_seconds": time.time() - start_all,
        "canonical_outputs_modified": False,
        "notes": [
            "Tier2, Tier4, legacy candidate lists, and old confidence labels were not used for training or aggregation.",
            "All selection, imputation, scaling, and tuning were confined to training data.",
            "Outer-fold predictive metrics are descriptive and do not weight gene rankings.",
        ],
    }
    atomic_write_json(manifest, output / "run_manifest.json")
    print(f"Completed {required_fits} fits in {(time.time() - start_all) / 60:.1f} minutes", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
