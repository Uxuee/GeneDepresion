#!/usr/bin/env python3
"""Deterministic classical half-sampling Stability Selection with L1 logistic regression."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code/ml_rra"))
from common import ANALYSES, ZERO_TOLERANCE, load_dataset, sha256_file, stable_seed  # noqa: E402

TIER2 = ROOT / "outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv"
PRIMARY = ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization/Tier2_consensus_rra_rankings.csv"
EXPECTED_TIER2_SHA = "1672ddf95a4f4c9e44ae42e7accc208c311e2e1da8533e1808dcf91bbf7b67a1"
EXPECTED_PRIMARY_SHA = "2e1ed5fc91aae752f57e46e61022a9077d7f1ab285cfc05a6973c3c6df23d358"
C_GRID = (0.01, 0.1, 1.0)
SUBSAMPLE_SEEDS = tuple(range(31001, 31101))


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--subsamples", type=int, default=100)
    return parser.parse_args()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False, lineterminator="\n")
    temporary.replace(path)


def fit_l1(x_train: np.ndarray, y_train: np.ndarray, c_value: float, seed: int):
    imputer = SimpleImputer(strategy="median", keep_empty_features=True)
    imputed = imputer.fit_transform(x_train)
    scaler = StandardScaler().fit(imputed)
    scaled = scaler.transform(imputed)
    model = LogisticRegression(
        penalty="l1", solver="liblinear", class_weight="balanced", max_iter=4000,
        C=c_value, random_state=seed,
    )
    model.fit(scaled, y_train)
    return imputer, scaler, model


def tune_l1(x_sub: np.ndarray, y_sub: np.ndarray, seed: int) -> tuple[float, float]:
    splitter = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    candidates = []
    for c_value in C_GRID:
        scores = []
        for fold, (fit_idx, valid_idx) in enumerate(splitter.split(x_sub, y_sub)):
            imputer, scaler, model = fit_l1(x_sub[fit_idx], y_sub[fit_idx], c_value, stable_seed(seed, "inner", fold, c_value))
            valid = scaler.transform(imputer.transform(x_sub[valid_idx]))
            scores.append(balanced_accuracy_score(y_sub[valid_idx], model.predict(valid)))
        candidates.append((float(np.mean(scores)), c_value))
    candidates.sort(key=lambda item: (-item[0], item[1]))
    return candidates[0][1], candidates[0][0]


def stratified_half_sample(y: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    pieces = []
    for label in (0, 1):
        available = np.flatnonzero(y == label)
        pieces.append(np.sort(rng.choice(available, size=len(available) // 2, replace=False)))
    return np.sort(np.concatenate(pieces))


def main() -> int:
    args = arguments()
    if args.subsamples != 100:
        raise ValueError("The frozen prospective design requires exactly 100 subsamples")
    if sha256_file(TIER2) != EXPECTED_TIER2_SHA or sha256_file(PRIMARY) != EXPECTED_PRIMARY_SHA:
        raise RuntimeError("Frozen input checksum mismatch")
    output = args.output_dir.resolve()
    checkpoints = output / "logs/checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    runtime_path = output / "stability_selection_runtime.csv"
    prior_runtime = pd.read_csv(runtime_path) if runtime_path.exists() else pd.DataFrame()
    tier2 = pd.read_csv(TIER2, dtype={"gene_id": str}, usecols=["region_or_subregion", "gene_id", "gene_symbol"])
    tier2["region"] = tier2["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
    all_paths: list[Path] = []
    runtime_rows = []
    for spec in ANALYSES:
        checkpoint = checkpoints / f"{spec.region.replace(' ', '_').replace('/', '_')}_stability_events.csv"
        all_paths.append(checkpoint)
        if args.resume and checkpoint.exists():
            saved = pd.read_csv(checkpoint)
            if saved["subsample_id"].nunique() == 100:
                previous = prior_runtime.loc[prior_runtime.get("region", pd.Series(dtype=str)).eq(spec.region)]
                elapsed = float(previous["runtime_seconds"].iloc[0]) if len(previous) else np.nan
                prior_status = str(previous["status"].iloc[0]) if len(previous) else "resumed"
                runtime_rows.append({"region": spec.region, "n_subsamples": 100, "runtime_seconds": elapsed, "status": prior_status})
                print(f"[{spec.region}] resumed 100 complete subsamples")
                continue
        expression, metadata, y = load_dataset(ROOT, spec)
        candidates = tier2.loc[tier2["region"].eq(spec.region)].copy()
        genes = candidates["gene_id"].to_numpy(dtype=str)
        missing = sorted(set(genes) - set(expression.columns.astype(str)))
        if missing:
            raise RuntimeError(f"{spec.region}: Tier2 genes missing from expression: {missing[:5]}")
        x = expression.loc[:, genes].to_numpy(dtype=float)
        symbols = dict(zip(candidates["gene_id"], candidates["gene_symbol"].fillna("")))
        started = time.time()
        events = []
        for subsample_id, base_seed in enumerate(SUBSAMPLE_SEEDS, start=1):
            sample_seed = stable_seed("stability_selection", spec.region, base_seed)
            indices = stratified_half_sample(y, sample_seed)
            x_sub, y_sub = x[indices], y[indices]
            chosen_c, inner_score = tune_l1(x_sub, y_sub, stable_seed(sample_seed, "tuning"))
            imputer, scaler, model = fit_l1(x_sub, y_sub, chosen_c, stable_seed(sample_seed, "final"))
            coefficients = model.coef_.reshape(-1).astype(float)
            absolute = np.abs(coefficients)
            selected = absolute > ZERO_TOLERANCE
            selected_order = np.flatnonzero(selected)[np.argsort(-absolute[selected], kind="stable")]
            ranks = np.full(len(genes), np.nan)
            ranks[selected_order] = np.arange(1, len(selected_order) + 1)
            selected_samples = "|".join(metadata.iloc[indices]["SampleID"].astype(str))
            for index, gene_id in enumerate(genes):
                events.append({
                    "region": spec.region, "subsample_id": subsample_id, "base_seed": base_seed,
                    "derived_seed": sample_seed, "n_subsample": len(indices),
                    "n_control_subsample": int((y_sub == 0).sum()), "n_mdd_subsample": int((y_sub == 1).sum()),
                    "selected_sample_ids_sha256": hashlib.sha256(selected_samples.encode("utf-8")).hexdigest(),
                    "chosen_C": chosen_c, "inner_mean_balanced_accuracy": inner_score,
                    "gene_id": gene_id, "gene_symbol": symbols.get(gene_id, ""),
                    "coefficient": coefficients[index], "absolute_coefficient": absolute[index],
                    "sign": np.sign(coefficients[index]) if selected[index] else 0,
                    "selected": bool(selected[index]), "within_subsample_rank": ranks[index],
                })
            if subsample_id % 10 == 0:
                print(f"[{spec.region}] {subsample_id}/100 subsamples")
        write_csv(pd.DataFrame(events), checkpoint)
        runtime_rows.append({"region": spec.region, "n_subsamples": 100, "runtime_seconds": time.time() - started, "status": "computed"})

    details = pd.concat([pd.read_csv(path, dtype={"gene_id": str}) for path in all_paths], ignore_index=True)
    write_csv(details, output / "stability_selection_run_details.csv")
    summaries = []
    for (region, gene_id), part in details.groupby(["region", "gene_id"], sort=False):
        selected = part.loc[part["selected"].astype(str).str.lower().eq("true")]
        signs = selected["sign"].loc[selected["sign"].ne(0)]
        summaries.append({
            "region": region, "gene_id": gene_id, "gene_symbol": part["gene_symbol"].iloc[0],
            "valid_subsamples": int(part["subsample_id"].nunique()),
            "selection_count": int(len(selected)), "selection_probability": float(len(selected) / part["subsample_id"].nunique()),
            "mean_abs_coefficient_when_selected": float(selected["absolute_coefficient"].mean()) if len(selected) else 0.0,
            "mean_rank_when_selected": float(selected["within_subsample_rank"].mean()) if len(selected) else np.nan,
            "sign_consistency_when_selected": float(signs.value_counts().max() / len(signs)) if len(signs) else np.nan,
        })
    summary = pd.DataFrame(summaries)
    final_tables = []
    for region, part in summary.groupby("region", sort=False):
        part = part.sort_values(
            ["selection_probability", "mean_abs_coefficient_when_selected", "gene_id"],
            ascending=[False, False, True], kind="mergesort",
        ).reset_index(drop=True)
        part["stability_export_rank"] = np.arange(1, len(part) + 1)
        supported = part["selection_count"].gt(0)
        part["stability_correlation_rank"] = part["stability_export_rank"].where(supported, int(supported.sum()) + 1)
        part["above_predeclared_0_60_threshold"] = part["selection_probability"].ge(0.60)
        final_tables.append(part)
    final = pd.concat(final_tables, ignore_index=True)
    write_csv(final, output / "stability_selection_gene_rankings.csv")
    write_csv(pd.DataFrame(runtime_rows), runtime_path)
    metadata = {
        "implementation": "transparent classical stratified half-sampling with L1 logistic regression",
        "package": "scikit-learn", "package_version": __import__("sklearn").__version__,
        "subsample_fraction": "floor(n per diagnosis class / 2)", "n_subsamples_per_region": 100,
        "seeds": [31001, 31100], "c_grid": list(C_GRID), "inner_folds": 3,
        "tuning_metric": "balanced_accuracy", "selection_tolerance": ZERO_TOLERANCE,
        "descriptive_threshold": 0.60, "threshold_role": "descriptive only; continuous selection probability is primary",
        "tier2_source_sha256": EXPECTED_TIER2_SHA, "frozen_primary_sha256": EXPECTED_PRIMARY_SHA,
    }
    (output / "stability_selection_metadata.json").write_bytes(json.dumps(metadata, indent=2, sort_keys=True).encode("utf-8"))
    print(f"Stability Selection complete: {len(final)} genes, {details[['region','subsample_id']].drop_duplicates().shape[0]} subsamples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
