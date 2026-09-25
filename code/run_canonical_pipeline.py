"""Registry-driven entry point for the canonical GeneDepresion workflow.

This runner does not duplicate scientific algorithms. It validates the registry
and can execute only stages that explicitly declare a command. Historical,
excluded, superseded, or unresolved stages are never executed here.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from project_paths import load_project_paths

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "provenance" / "canonical_pipeline.json"
ALLOWED_TO_RUN = {"canonical", "canonical_external_input"}


def main() -> int:
    parser = argparse.ArgumentParser(description="Inspect or run registered canonical stages.")
    parser.add_argument("--stage", action="append", help="Stage id to select; repeat as needed.")
    parser.add_argument("--status", action="store_true", help="Print the canonical stage table.")
    parser.add_argument("--validate", action="store_true", help="Run the canonical registry validator.")
    parser.add_argument("--dry-run", action="store_true", help="Show commands without executing them.")
    parser.add_argument("--config", help="Optional local project-path YAML file.")
    for key in ("project_root", "raw_data_root", "output_root", "external_resource_root", "r_library_path"):
        parser.add_argument(f"--{key.replace('_', '-')}", dest=key)
    args = parser.parse_args()
    overrides = {key: getattr(args, key) for key in ("project_root", "raw_data_root", "output_root", "external_resource_root", "r_library_path") if getattr(args, key)}
    paths = load_project_paths(args.config, overrides)
    data = json.loads(REGISTRY.read_text(encoding="utf-8"))
    stages = data["stages"]
    selected = stages if not args.stage else [s for s in stages if s["stage_id"] in set(args.stage)]
    if args.stage and len(selected) != len(set(args.stage)):
        found = {s["stage_id"] for s in selected}
        parser.error(f"unknown stage(s): {sorted(set(args.stage) - found)}")
    if args.status or (not args.validate and not args.dry_run and not args.stage):
        print(f"project_root={paths.get('project_root', ROOT)}")
        for stage in stages:
            print(f"{stage['stage_id']:>3}  {stage['status']:<37} {stage['stage_name']}")
    if args.validate:
        return subprocess.call([sys.executable, str(ROOT / "audit" / "validate_canonical_pipeline.py")])
    for stage in selected:
        status = stage["status"]
        command = stage.get("execution", {}).get("command")
        if status not in ALLOWED_TO_RUN:
            print(f"REFUSE {stage['stage_id']}: status={status}")
            continue
        if not command:
            print(f"NO-OP  {stage['stage_id']}: no canonical execution command registered")
            continue
        print(f"RUN   {stage['stage_id']}: {' '.join(command)}")
        if not args.dry_run:
            subprocess.run(command, cwd=paths.get("project_root", ROOT), check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
