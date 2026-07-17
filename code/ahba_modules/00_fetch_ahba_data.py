"""Locate or download AHBA microarray data for the module-analysis stage.

Run from the project root:

    python code/ahba_modules/00_fetch_ahba_data.py

The script first checks data/ahba for existing donor-level AHBA files. If they
are missing, it attempts abagen.fetch_microarray(donors="all",
data_dir="data/ahba"). Raw downloads are intentionally ignored by git.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "ahba"
OUT_DIR = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline" / "11_AHBA_modules_enrichment"
LOG_PATH = OUT_DIR / "ahba_data_inventory.json"

EXPECTED_NAMES = {"MicroarrayExpression.csv", "Probes.csv", "SampleAnnot.csv"}


def inventory_files(data_dir: Path) -> list[dict]:
    rows = []
    for path in sorted(data_dir.rglob("*")):
        if path.is_file() and path.name not in {".gitkeep", "README.md"}:
            rows.append(
                {
                    "path": str(path),
                    "relative_path": str(path.relative_to(data_dir)),
                    "file_name": path.name,
                    "size_bytes": path.stat().st_size,
                    "parent": path.parent.name,
                }
            )
    return rows


def has_ahba_core_files(rows: list[dict]) -> bool:
    names = {row["file_name"] for row in rows}
    return bool(EXPECTED_NAMES & names) and any(row["file_name"] == "MicroarrayExpression.csv" for row in rows)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = inventory_files(DATA_DIR)
    attempted_download = False
    download_error = None

    if not has_ahba_core_files(rows):
        attempted_download = True
        try:
            import abagen  # type: ignore

            print("No complete local AHBA microarray inventory found. Trying abagen.fetch_microarray...")
            fetched = abagen.fetch_microarray(donors="all", data_dir=str(DATA_DIR))
            print(f"abagen.fetch_microarray returned: {fetched}")
            rows = inventory_files(DATA_DIR)
        except Exception as exc:  # noqa: BLE001
            download_error = repr(exc)
            print("AHBA automatic download failed.")
            print("Install abagen and retry, or manually download donor microarray files into data/ahba.")
            print("Expected donor files include MicroarrayExpression.csv, Probes.csv, and SampleAnnot.csv.")
            print(f"Error: {download_error}")

    donors = sorted({row["parent"] for row in rows if row["file_name"] in EXPECTED_NAMES})
    manifest = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "data_dir": str(DATA_DIR),
        "ahba_core_files_present": has_ahba_core_files(rows),
        "attempted_download": attempted_download,
        "download_error": download_error,
        "donor_or_parent_ids_detected": donors,
        "n_files": len(rows),
        "files": rows,
        "manual_download_note": (
            "If ahba_core_files_present is false, place AHBA donor folders under data/ahba. "
            "Each donor folder should contain MicroarrayExpression.csv, Probes.csv, and SampleAnnot.csv."
        ),
    }
    LOG_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote AHBA data inventory: {LOG_PATH}")
    print(f"Core AHBA files present: {manifest['ahba_core_files_present']}")


if __name__ == "__main__":
    main()
