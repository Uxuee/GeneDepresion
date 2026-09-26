"""Small, dependency-free helpers for deterministic provenance manifests."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path, root: Path | None = None) -> dict[str, object]:
    label = path.relative_to(root).as_posix() if root else path.as_posix()
    return {"path": label, "sha256": sha256(path), "bytes": path.stat().st_size}


def write_manifest(path: Path, payload: dict[str, object]) -> None:
    body = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        **payload,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(body, indent=2) + "\n", encoding="utf-8")


def new_stage_manifest(
    *,
    stage: str,
    command: list[str],
    script: Path,
    inputs: list[Path],
    outputs: list[Path],
    parameters: dict[str, object],
    random_seeds: dict[str, object],
    environment: dict[str, object],
    resource_versions: dict[str, object],
    upstream_stage_ids: list[str],
    repository_root: Path,
) -> dict[str, object]:
    """Build the required provenance record for a newly generated stage."""
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repository_root, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    return {
        "stage": stage,
        "git_commit": commit,
        "command": command,
        "script": file_record(script, repository_root),
        "inputs": [file_record(path, repository_root) for path in inputs],
        "outputs": [file_record(path, repository_root) for path in outputs],
        "parameters": parameters,
        "random_seeds": random_seeds,
        "environment": environment,
        "database_or_resource_versions": resource_versions,
        "upstream_canonical_stage_ids": upstream_stage_ids,
    }
