"""Check external-data placement without downloading or exposing local data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from project_paths import load_project_paths  # noqa: E402

MANIFEST = ROOT / "provenance" / "external_data_manifest.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_template(template: str, paths: dict[str, Path]) -> Path | None:
    value = template.strip()
    if not value or value.startswith("package:"):
        return None
    for key, root in paths.items():
        value = value.replace("${" + key + "}", str(root))
    if "${" in value:
        return None
    return Path(value).resolve()


def check(config: str | None = None) -> tuple[int, list[str]]:
    paths = load_project_paths(config)
    rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8", newline="")))
    messages: list[str] = []
    failures = 0
    for row in rows:
        if row.get("check_kind") == "package":
            messages.append(f"PACKAGE {row['resource_name']}: restore through the recorded stage environment")
            continue
        destination = resolve_template(row.get("local_destination", ""), paths)
        required = row.get("required", "no").lower() == "yes"
        if destination is None:
            if required:
                failures += 1
                messages.append(
                    f"MISSING CONFIG {row['resource_name']}: set path roots in "
                    "config/project_paths.local.yaml (copy config/project_paths.example.yaml)"
                )
            continue
        exists = destination.is_dir() if row.get("check_kind") == "directory" else destination.is_file()
        if not exists:
            prefix = "MISSING" if required else "OPTIONAL MISSING"
            failures += int(required)
            messages.append(f"{prefix} {row['resource_name']}: expected {destination}")
            continue
        expected = row.get("sha256", "").strip()
        if expected and destination.is_file():
            observed = sha256(destination)
            if observed != expected:
                failures += 1
                messages.append(f"CHECKSUM MISMATCH {row['resource_name']}: {destination}")
                continue
        messages.append(f"OK {row['resource_name']}: {destination}")
    if failures:
        messages.append(
            f"External-data preflight found {failures} required issue(s). "
            "See data/README.md and provenance/external_data_manifest.csv."
        )
        return 2, messages
    messages.append("External-data preflight passed for all required configured files.")
    return 0, messages


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", help="Path to project_paths YAML file")
    args = parser.parse_args()
    status, messages = check(args.config)
    print("\n".join(messages))
    return status


if __name__ == "__main__":
    raise SystemExit(main())
