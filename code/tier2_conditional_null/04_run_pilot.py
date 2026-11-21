#!/usr/bin/env python3
"""Run/resume a range of conditional-null permutations."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs/length_gt_100_clean_pipeline/20_Tier2_ML_conditional_null"
PYTHON = ROOT / ".runtime/python/python.exe"
RSCRIPT = Path("D:/genes/GeneDepresion/environment/phase3b/R-4.6.1/bin/Rscript.exe")
R_LIBRARY = ROOT / "environment/ml_rra/library"
FROZEN_RRA = ROOT / "code/tier2_ml_rra/02_run_tier2_rra.R"
RESAMPLE_RRA = ROOT / "code/tier2_conditional_null/02_build_permutation_resample_rra.R"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--end", type=int, default=20)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_completion(directory: Path) -> dict | None:
    marker = directory / "permutation_completion.json"
    if not marker.exists():
        return None
    payload = json.loads(marker.read_text(encoding="utf-8"))
    for name, expected in payload["files_sha256"].items():
        path = directory / name
        if not path.exists() or sha256(path) != expected:
            raise RuntimeError(f"Completed permutation checksum mismatch: {path}")
    return payload


def run_command(command: list[str], env: dict[str, str] | None = None) -> None:
    subprocess.run(command, cwd=ROOT, env=env, check=True)


def run_one(index: int, resume: bool) -> dict:
    directory = OUTPUT / "permutations" / f"perm_{index:04d}"
    existing = validate_completion(directory)
    if existing is not None and resume:
        print(f"perm_{index:04d}: complete checksums verified; skipped", flush=True)
        return existing
    start = time.time()
    command = [str(PYTHON), str(ROOT / "code/tier2_conditional_null/01_run_permutation_ml.py"), "--permutation-index", str(index)]
    if resume:
        command.append("--resume")
    run_command(command)
    env = os.environ.copy()
    env["R_LIBS_USER"] = str(R_LIBRARY)
    run_command([str(RSCRIPT), str(FROZEN_RRA), str(directory)], env)
    run_command([str(RSCRIPT), str(RESAMPLE_RRA), str(directory)], env)
    run_command([str(PYTHON), str(ROOT / "code/tier2_conditional_null/03_finalize_permutation.py"), str(directory)])
    payload = validate_completion(directory)
    if payload is None:
        raise RuntimeError(f"perm_{index:04d}: final completion record absent")
    payload["orchestration_wall_seconds"] = time.time() - start
    return payload


def write_manifest() -> None:
    rows = []
    permutations = OUTPUT / "permutations"
    if permutations.exists():
        for directory in sorted(permutations.glob("perm_*")):
            payload = validate_completion(directory)
            if payload:
                rows.append({
                    "permutation_index": payload["permutation_index"],
                    "permutation_seed": payload["permutation_seed"],
                    "status": payload["status"],
                    "n_feature_model_fits": payload["n_feature_model_fits"],
                    "n_lasso_outer_fits": payload["n_lasso_outer_fits"],
                    "wall_runtime_seconds_ml": payload["wall_runtime_seconds_ml"],
                    "directory": str(directory.relative_to(ROOT)).replace("\\", "/"),
                    "completion_sha256": sha256(directory / "permutation_completion.json"),
                })
    OUTPUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).sort_values("permutation_index").to_csv(OUTPUT / "permutation_manifest.csv", index=False, lineterminator="\n")


def main() -> int:
    args = arguments()
    if not (1 <= args.start <= args.end <= 99):
        raise ValueError("Require 1 <= start <= end <= 99")
    if not (1 <= args.workers <= 4):
        raise ValueError("workers must be 1-4")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    indices = list(range(args.start, args.end + 1))
    if args.workers == 1:
        for index in indices:
            run_one(index, args.resume)
            write_manifest()
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_one, index, args.resume): index for index in indices}
            for future in concurrent.futures.as_completed(futures):
                future.result()
                write_manifest()
    write_manifest()
    print(f"Pilot orchestration complete for permutations {args.start}-{args.end}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
