"""Resolve machine-local project paths without embedding them in analysis code."""

from __future__ import annotations

import os
from pathlib import Path

KEYS = ("project_root", "raw_data_root", "output_root", "external_resource_root", "r_library_path")
ENV = {key: f"GENEDEPRESSION_{key.upper()}" for key in KEYS}


def _read_simple_yaml(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        key, value = line.split(":", 1)
        if key.strip() in KEYS:
            values[key.strip()] = value.strip().strip("'\"")
    return values


def load_project_paths(
    config_path: str | Path | None = None,
    overrides: dict[str, str | Path] | None = None,
) -> dict[str, Path]:
    """Resolve paths in precedence order: overrides, environment, config, defaults."""
    repo = Path(__file__).resolve().parents[1]
    local = Path(config_path) if config_path else repo / "config" / "project_paths.local.yaml"
    raw = _read_simple_yaml(local)
    raw.setdefault("project_root", str(repo))
    raw.setdefault("output_root", str(repo / "outputs" / "length_gt_100_clean_pipeline"))
    overrides = overrides or {}
    result: dict[str, Path] = {}
    for key in KEYS:
        value = overrides.get(key) or os.environ.get(ENV[key], raw.get(key, ""))
        if value:
            result[key] = Path(value).expanduser().resolve()
    return result
