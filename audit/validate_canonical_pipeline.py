"""Validate the canonical pipeline registry and its publication guardrails."""

from __future__ import annotations

import hashlib
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "provenance" / "canonical_pipeline.json"
REPORT = ROOT / "audit" / "canonical_pipeline_validation.json"
REQUIRED = {"stage_id", "stage_name", "status", "canonical_script", "canonical_inputs", "canonical_outputs", "input_sha256", "output_sha256", "scientific_dependency", "git_milestone", "random_seed", "environment_reference", "rerun_required", "historical_replacements", "notes"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_spec_exists(spec: str) -> bool:
    _, commit, path = spec.split(":", 2)
    return subprocess.run(["git", "cat-file", "-e", f"{commit}:{path}"], cwd=ROOT, capture_output=True).returncode == 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-missing-local-artifacts",
        action="store_true",
        help="Fresh-clone mode: record intentionally uncommitted artifacts as notices while still checking every available checksum.",
    )
    args = parser.parse_args()
    data = json.loads(REGISTRY.read_text(encoding="utf-8"))
    errors: list[str] = []
    warnings: list[str] = []
    notices: list[str] = []
    ids: list[str] = []
    allowed = set(data["allowed_statuses"])
    for stage in data["stages"]:
        missing = REQUIRED - set(stage)
        if missing:
            errors.append(f"{stage.get('stage_id')}: missing fields {sorted(missing)}")
        ids.append(stage["stage_id"])
        if stage["status"] not in allowed:
            errors.append(f"{stage['stage_id']}: invalid status {stage['status']}")
        script = stage.get("canonical_script")
        scripts = script if isinstance(script, list) else ([script] if script else [])
        for item in scripts:
            if item.startswith("git:"):
                if not git_spec_exists(item):
                    errors.append(f"{stage['stage_id']}: missing Git script object {item}")
            elif not (ROOT / item).exists():
                warnings.append(f"{stage['stage_id']}: registered script is absent from working tree: {item}")
        if stage["status"] == "canonical":
            for output in stage["canonical_outputs"]:
                if output.startswith("${"):
                    continue
                output_path = ROOT / output
                if not output_path.exists():
                    message = f"{stage['stage_id']}: local canonical artifact is unavailable in this checkout: {output}"
                    if args.allow_missing_local_artifacts:
                        notices.append(message)
                    else:
                        errors.append(message)
                elif output_path.is_file():
                    recorded = stage["output_sha256"].get(output_path.name) or stage["output_sha256"].get(output)
                    if recorded and sha256(output_path) != recorded:
                        errors.append(f"{stage['stage_id']}: canonical output checksum mismatch: {output}")
                if "10_ML_prioritized_Tier2_extended" in output:
                    errors.append(f"{stage['stage_id']}: canonical output points to superseded ML directory")
    if len(ids) != len(set(ids)):
        errors.append("duplicate stage_id values")
    ranking = ROOT / "outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv"
    expected = data["canonical_ml_result"]["sha256"]
    if not ranking.exists() or sha256(ranking) != expected:
        errors.append("canonical corrected ranking is absent or its SHA-256 differs")
    if ranking.exists():
        with ranking.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        classes = {name: sum(r["ML_priority_class"] == name for r in rows) for name in ("High-confidence ML-prioritized Tier 2", "Moderate-confidence ML-prioritized Tier 2", "Exploratory ML-prioritized Tier 2")}
        if (len(rows), classes["High-confidence ML-prioritized Tier 2"], classes["Moderate-confidence ML-prioritized Tier 2"], classes["Exploratory ML-prioritized Tier 2"]) != (671, 30, 119, 522):
            errors.append(f"canonical corrected ranking counts differ: rows={len(rows)}, classes={classes}")
    tier4 = json.loads((ROOT / "provenance/tier4_ranked_signature_recovery_manifest.json").read_text(encoding="utf-8"))
    for record in tier4["files"]:
        if subprocess.run(["git", "cat-file", "-e", record["blob"]], cwd=ROOT, capture_output=True).returncode:
            errors.append(f"missing Tier4 Git blob {record['blob']}")
    phase2 = json.loads((ROOT / "outputs/length_gt_100_clean_pipeline/16_corrected_ML_downstream/phase2_manifest.json").read_text(encoding="utf-8"))
    if phase2.get("canonical_ranking_sha256") != expected or phase2.get("counts") != {"total": 671, "high": 30, "moderate": 119, "high_plus_moderate": 149}:
        errors.append("Phase-2 manifest canonical ranking identity/counts differ")
    for rel, expected_hash in phase2.get("key_output_sha256", {}).items():
        path = ROOT / Path(rel.replace("\\", "/"))
        if not path.exists():
            errors.append(f"Phase-2 manifest output missing: {rel}")
        elif sha256(path) != expected_hash:
            errors.append(f"Phase-2 manifest checksum mismatch: {rel}")
    stages = {s["stage_id"]: s for s in data["stages"]}
    if stages["17"]["status"] != "excluded_from_publication":
        errors.append("Stage 17 must remain excluded_from_publication until its generator is recovered")
    active_scripts = [
        ROOT / "code/run_canonical_pipeline.py",
        ROOT / "code/regenerate_corrected_ml_downstream.py",
        ROOT / "code/summarize_corrected_ml_downstream.py",
        ROOT / "code/canonical_downstream_paths.py",
    ]
    for path in active_scripts:
        if "10_ML_prioritized_Tier2_extended" in path.read_text(encoding="utf-8"):
            errors.append(f"active corrected script references superseded ranking: {path.relative_to(ROOT)}")
    current_readme = (ROOT / "README.md").read_text(encoding="utf-8").split("<!-- HISTORICAL_CONTENT_START -->", 1)[0]
    for stale in ("35 high", "130 moderate", "165 high+moderate", "v2 denominator = 180"):
        if stale.lower() in current_readme.lower():
            errors.append(f"current README contains stale value: {stale}")
    result = {"ok": not errors, "errors": errors, "warnings": warnings, "notices": notices, "stage_count": len(ids), "canonical_ranking_sha256": expected}
    if not args.allow_missing_local_artifacts:
        REPORT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
