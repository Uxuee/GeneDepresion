#!/usr/bin/env python3
"""Run one deterministic, resumable conditional-null Tier2 ML permutation."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.dummy import DummyClassifier

ROOT = Path(__file__).resolve().parents[2]
PRIMARY_SCRIPT = ROOT / "code/tier2_ml_rra/01_run_tier2_repeated_ml.py"
spec = importlib.util.spec_from_file_location("frozen_tier2_ml", PRIMARY_SCRIPT)
if spec is None or spec.loader is None:
    raise RuntimeError("Cannot load frozen Tier2 ML implementation")
frozen = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = frozen
spec.loader.exec_module(frozen)

sys.path.insert(0, str(ROOT / "code/ml_rra"))
from common import (  # noqa: E402
    ANALYSES,
    MODEL_ORDER,
    OUTER_SEEDS,
    atomic_write_csv,
    atomic_write_json,
    compute_metrics,
    load_dataset,
    make_inner_splits,
    make_outer_splits,
    safe_name,
    sha256_file,
    stable_seed,
)

EXPECTED_OBSERVED_SHA = "2e1ed5fc91aae752f57e46e61022a9077d7f1ab285cfc05a6973c3c6df23d358"
EXPECTED_BENCHMARK_LEDGER_SHA = "4e17dc8234317efc9bce92aa03fb7ae63ff72a128392d202f2f31a105c09cf10"
EXPECTED_ALLGENE_SHA = "0aefec5fcfa100deafa1525c2d31521ae8c00ca890ef4e1b4145b1b861c2c87e"
DUMMY_MODEL = "dummy_most_frequent"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutation-index", type=int, required=True, choices=range(1, 1000))
    parser.add_argument("--output-root", type=Path, default=ROOT / "outputs/length_gt_100_clean_pipeline/20_Tier2_ML_conditional_null")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def permute_labels(y: np.ndarray, region: str, permutation_seed: int) -> tuple[np.ndarray, int]:
    region_seed = stable_seed("Tier2_conditional_null", permutation_seed, region)
    rng = np.random.default_rng(region_seed)
    permuted = rng.permutation(y)
    # An unchanged shuffle is not informative. A deterministic cross-class swap
    # preserves exact counts while ensuring that labels actually change.
    if np.array_equal(permuted, y):
        zero = int(np.flatnonzero(permuted == 0)[0])
        one = int(np.flatnonzero(permuted == 1)[0])
        permuted[[zero, one]] = permuted[[one, zero]]
    if sorted(permuted.tolist()) != sorted(y.tolist()):
        raise RuntimeError(f"{region}: permutation changed class counts")
    return permuted.astype(int), region_seed


def checksum_guard() -> None:
    observed = ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization/Tier2_consensus_rra_rankings.csv"
    benchmark = ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark/checksums.sha256"
    allgene = ROOT / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/consensus_rra_rankings.csv"
    expected = [(observed, EXPECTED_OBSERVED_SHA), (benchmark, EXPECTED_BENCHMARK_LEDGER_SHA), (allgene, EXPECTED_ALLGENE_SHA)]
    for path, digest in expected:
        if not path.exists() or sha256_file(path) != digest:
            raise RuntimeError(f"Frozen-input checksum differs: {path}")


def complete_hashes(output: Path) -> dict[str, str]:
    names = [
        "permuted_label_mapping.csv", "outer_fold_assignments.csv",
        "Tier2_fold_gene_rankings.csv", "Tier2_outer_predictions.csv",
        "Tier2_outer_performance.csv", "Tier2_model_performance_summary.csv",
        "Tier2_analysis_unit_manifest.csv",
    ]
    return {name: sha256_file(output / name) for name in names if (output / name).exists()}


def main() -> int:
    args = arguments()
    checksum_guard()
    permutation_seed = 41000 + args.permutation_index
    output = args.output_root.resolve() / "permutations" / f"perm_{args.permutation_index:04d}"
    completion_path = output / "ml_completion.json"
    if args.resume and completion_path.exists():
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
        current = complete_hashes(output)
        if completion.get("raw_output_sha256") == current and completion.get("n_feature_model_fits") == 750:
            print(f"perm_{args.permutation_index:04d}: complete checkpoint verified; no fits recomputed")
            return 0
        raise RuntimeError("Completion marker exists but raw output checksums do not match")

    checkpoint_dir = output / "logs/checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    tier2 = frozen.load_tier2_identities(ROOT)
    label_rows: list[dict] = []
    fold_rows: list[dict] = []
    rank_paths: list[Path] = []
    prediction_paths: list[Path] = []
    performance_paths: list[Path] = []
    unit_records: list[dict] = []
    run_start = time.time()

    for analysis in ANALYSES:
        full_x, metadata, observed_y = load_dataset(ROOT, analysis)
        y, region_permutation_seed = permute_labels(observed_y, analysis.region, permutation_seed)
        if np.array_equal(y, observed_y):
            raise RuntimeError(f"{analysis.region}: labels did not change")
        label_rows.extend(
            {
                "permutation_index": args.permutation_index,
                "permutation_seed": permutation_seed,
                "region_seed": region_permutation_seed,
                "region": analysis.region,
                "sample_id": str(sample_id),
                "observed_label": int(observed),
                "permuted_label": int(permuted),
            }
            for sample_id, observed, permuted in zip(full_x.index, observed_y, y)
        )
        tier2_name = frozen.TIER2_REGION_NAMES.get(analysis.region, analysis.region)
        identities = tier2.loc[tier2["tier2_region"].eq(tier2_name)].copy()
        gene_ids_list = identities["gene_id"].astype(str).tolist()
        missing = sorted(set(gene_ids_list) - set(full_x.columns.astype(str)))
        if missing:
            raise RuntimeError(f"{analysis.region}: missing Tier2 genes: {missing[:10]}")
        x_frame = full_x.loc[:, gene_ids_list].copy()
        x = x_frame.to_numpy(dtype=float)
        gene_ids = x_frame.columns.to_numpy(dtype=str)
        symbol_map = dict(zip(identities["gene_id"].astype(str), identities["gene_symbol"].fillna("")))
        region_safe = safe_name(analysis.region)
        rank_path = checkpoint_dir / f"{region_safe}_Tier2_fold_gene_rankings.csv"
        prediction_path = checkpoint_dir / f"{region_safe}_Tier2_outer_predictions.csv"
        performance_path = checkpoint_dir / f"{region_safe}_Tier2_outer_performance.csv"
        rank_paths.append(rank_path)
        prediction_paths.append(prediction_path)
        performance_paths.append(performance_path)
        rank_rows = (frozen.read_checkpoint(rank_path, {"gene_id": str}) if args.resume else pd.DataFrame()).to_dict("records")
        prediction_rows = (frozen.read_checkpoint(prediction_path) if args.resume else pd.DataFrame()).to_dict("records")
        performance_rows = (frozen.read_checkpoint(performance_path) if args.resume else pd.DataFrame()).to_dict("records")
        completed = {str(row["fit_id"]) for row in performance_rows}
        splits = list(make_outer_splits(y, n_splits=5, seeds=OUTER_SEEDS))
        for repeat, fold, outer_seed, train_idx, test_idx in splits:
            for role, indices in (("train", train_idx), ("test", test_idx)):
                fold_rows.extend(
                    {
                        "permutation_index": args.permutation_index,
                        "region": analysis.region,
                        "repeat": repeat + 1,
                        "fold": fold + 1,
                        "outer_seed": outer_seed,
                        "sample_id": str(x_frame.index[index]),
                        "permuted_label": int(y[index]),
                        "role": role,
                    }
                    for index in indices
                )
            inner_seed = stable_seed(analysis.region, permutation_seed, outer_seed, fold, "Tier2_inner_split")
            inner_splits, n_inner = make_inner_splits(y[train_idx], inner_seed)
            for model_name in MODEL_ORDER:
                fit_id = f"p{args.permutation_index:04d}__{region_safe}__r{repeat + 1}__f{fold + 1}__{model_name}"
                if fit_id in completed:
                    continue
                fit_start = time.time()
                model_seed = stable_seed(analysis.region, permutation_seed, outer_seed, fold, "Tier2", model_name)
                best = frozen.tune_model(x[train_idx], y[train_idx], model_name, inner_splits, model_seed)
                pipeline = frozen.fit_fixed_pipeline(x[train_idx], y[train_idx], model_name, best["params"], model_seed)
                y_pred = pipeline.predict(x[test_idx])
                y_score = pipeline.decision_score(x[test_idx])
                metrics = compute_metrics(y[test_idx], y_pred, y_score)
                importance = frozen.importance_table(pipeline, model_name, gene_ids, symbol_map)
                shared = {
                    "permutation_index": args.permutation_index, "permutation_seed": permutation_seed,
                    "fit_id": fit_id, "region": analysis.region, "repeat": repeat + 1,
                    "fold": fold + 1, "outer_seed": outer_seed, "model": model_name,
                    "candidate_universe": "canonical_Tier2_fixed", "n_tier2_features": len(gene_ids),
                }
                rank_rows.extend({**shared, **row} for row in importance.to_dict("records"))
                prediction_rows.extend(
                    {**shared, "sample_id": str(x_frame.index[index]), "y_true": int(y[index]),
                     "y_pred": int(y_pred[local]), "y_score": float(y_score[local])}
                    for local, index in enumerate(test_idx)
                )
                performance_rows.append(
                    {**shared, "model_params_json": json.dumps(best["params"], sort_keys=True),
                     "inner_folds": n_inner, "inner_mean_balanced_accuracy": best["mean_balanced_accuracy"],
                     "n_train": len(train_idx), "n_test": len(test_idx), "n_ranked_genes": len(importance),
                     "elapsed_seconds": time.time() - fit_start, "contributes_gene_ranking": True, **metrics}
                )
                completed.add(fit_id)
                if rank_rows:
                    atomic_write_csv(pd.DataFrame(rank_rows), rank_path)
                atomic_write_csv(pd.DataFrame(prediction_rows), prediction_path)
                atomic_write_csv(pd.DataFrame(performance_rows), performance_path)

            dummy_id = f"p{args.permutation_index:04d}__{region_safe}__r{repeat + 1}__f{fold + 1}__{DUMMY_MODEL}"
            if dummy_id not in completed:
                dummy = DummyClassifier(strategy="most_frequent")
                dummy.fit(x[train_idx], y[train_idx])
                y_pred = dummy.predict(x[test_idx])
                probabilities = dummy.predict_proba(x[test_idx])
                classes = list(dummy.classes_)
                y_score = probabilities[:, classes.index(1)] if 1 in classes else np.zeros(len(test_idx))
                metrics = compute_metrics(y[test_idx], y_pred, y_score)
                shared = {
                    "permutation_index": args.permutation_index, "permutation_seed": permutation_seed,
                    "fit_id": dummy_id, "region": analysis.region, "repeat": repeat + 1, "fold": fold + 1,
                    "outer_seed": outer_seed, "model": DUMMY_MODEL,
                    "candidate_universe": "performance_reference_only", "n_tier2_features": 0,
                }
                prediction_rows.extend(
                    {**shared, "sample_id": str(x_frame.index[index]), "y_true": int(y[index]),
                     "y_pred": int(y_pred[local]), "y_score": float(y_score[local])}
                    for local, index in enumerate(test_idx)
                )
                performance_rows.append(
                    {**shared, "model_params_json": json.dumps({"strategy": "most_frequent"}),
                     "inner_folds": 0, "inner_mean_balanced_accuracy": np.nan,
                     "n_train": len(train_idx), "n_test": len(test_idx), "n_ranked_genes": 0,
                     "elapsed_seconds": 0.0, "contributes_gene_ranking": False, **metrics}
                )
                completed.add(dummy_id)
                atomic_write_csv(pd.DataFrame(prediction_rows), prediction_path)
                atomic_write_csv(pd.DataFrame(performance_rows), performance_path)

        unit_performance = pd.DataFrame(performance_rows)
        unit_records.append(
            {"permutation_index": args.permutation_index, "permutation_seed": permutation_seed,
             "region": analysis.region, "tier2_source_region": tier2_name, "n_samples": len(y),
             "n_mdd": int(y.sum()), "n_control": int((y == 0).sum()),
             "n_tier2_candidates": len(gene_ids), "n_tier2_matched_expression": len(gene_ids),
             "n_tier2_missing_expression": 0, "expression_path": analysis.expression_relpath,
             "expression_sha256": sha256_file(ROOT / analysis.expression_relpath),
             "cumulative_fit_elapsed_seconds": float(unit_performance["elapsed_seconds"].sum())}
        )
        print(f"perm_{args.permutation_index:04d} {analysis.region}: 125 feature fits complete", flush=True)

    rankings = frozen.concatenate(rank_paths, {"gene_id": str})
    predictions = frozen.concatenate(prediction_paths)
    performance = frozen.concatenate(performance_paths)
    if performance.loc[performance["model"].isin(MODEL_ORDER), "fit_id"].nunique() != 750:
        raise RuntimeError("Permutation feature-model fits are incomplete")
    if performance.loc[performance["model"].eq(DUMMY_MODEL), "fit_id"].nunique() != 150:
        raise RuntimeError("Permutation dummy fits are incomplete")
    if rankings.duplicated(["fit_id", "gene_id"]).any():
        raise RuntimeError("Duplicate gene within permutation fit ranking")

    atomic_write_csv(pd.DataFrame(label_rows), output / "permuted_label_mapping.csv")
    atomic_write_csv(pd.DataFrame(fold_rows), output / "outer_fold_assignments.csv")
    atomic_write_csv(rankings, output / "Tier2_fold_gene_rankings.csv")
    atomic_write_csv(predictions, output / "Tier2_outer_predictions.csv")
    atomic_write_csv(performance, output / "Tier2_outer_performance.csv")
    summary = performance.groupby(["region", "model"], as_index=False).agg(
        n_outer_fits=("fit_id", "nunique"), contributes_gene_ranking=("contributes_gene_ranking", "first"),
        balanced_accuracy_mean=("balanced_accuracy", "mean"), balanced_accuracy_sd=("balanced_accuracy", "std"),
        f1_mean=("f1", "mean"), roc_auc_mean=("roc_auc", "mean"),
        precision_mean=("precision", "mean"), recall_mean=("recall", "mean"),
    )
    atomic_write_csv(summary, output / "Tier2_model_performance_summary.csv")
    atomic_write_csv(pd.DataFrame(unit_records), output / "Tier2_analysis_unit_manifest.csv")
    completion = {
        "status": "ml_complete", "permutation_index": args.permutation_index,
        "permutation_seed": permutation_seed, "candidate_pairs": len(tier2),
        "n_feature_model_fits": 750, "n_dummy_fits": 150,
        "wall_runtime_seconds": time.time() - run_start,
        "cumulative_fit_elapsed_seconds": float(performance["elapsed_seconds"].sum()),
        "raw_output_sha256": complete_hashes(output),
        "python": platform.python_version(), "scikit_learn": sklearn.__version__,
        "conditional_null_warning": "Tier2 was fixed after selection with the observed labels; this is not an end-to-end global-null test.",
    }
    atomic_write_json(completion, completion_path)
    print(f"perm_{args.permutation_index:04d}: 750 feature fits complete in {completion['wall_runtime_seconds']:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
