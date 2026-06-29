"""Fetch or load published AHBA canonical module assignments.

This script standardizes Hawrylycz et al. 2015 Supplementary Table 6 into a
small gene-to-module table for downstream enrichment. Automatic download is
kept conservative because publisher supplement URLs can change; if no local
file is found, the script writes manual instructions and exits gracefully.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MANUAL_DIR = PROJECT_ROOT / "data" / "ahba" / "published_modules"
OUT_DIR = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline" / "ahba_modules" / "published_modules"
OUT_TABLE = OUT_DIR / "ahba_published_gene_to_module_assignments.csv"
SOURCE_REFERENCE = "Hawrylycz et al. 2015, Canonical genetic signatures of the adult human brain, Supplementary Table 6"
DEFAULT_DOWNLOAD_URLS = [
    "https://static-content.springer.com/esm/art%3A10.1038%2Fnn.4171/MediaObjects/41593_2015_BFnn4171_MOESM97_ESM.csv",
]

GENE_COLUMN_CANDIDATES = [
    "gene_symbol", "gene symbol", "genesymbol", "symbol", "gene", "gene_name", "gene name", "hgnc_symbol"
]
MODULE_COLUMN_CANDIDATES = [
    "published_module", "module label", "module_label", "module", "module assignment", "module_assignment",
    "consensus module", "consensus_module", "signature", "signature module", "eigengene module",
    "module color", "module_color"
]
LABEL_COLUMN_CANDIDATES = ["module_label", "module label", "module color", "module_color", "label", "module name", "module_name", "description"]


def clean_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).strip().lower()).strip()


def normalize_symbol(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def find_column(columns: Iterable[str], candidates: list[str]) -> str | None:
    cleaned = {clean_name(c): c for c in columns}
    for candidate in candidates:
        key = clean_name(candidate)
        if key in cleaned:
            return cleaned[key]
    for col in columns:
        col_clean = clean_name(col)
        if any(clean_name(candidate) in col_clean for candidate in candidates):
            return col
    return None


def find_eigengene_corr_column(columns: Iterable[str]) -> str | None:
    for col in columns:
        key = clean_name(col)
        if key in {"eigengene correlation", "eigengene_correlation", "module membership", "kme", "correlation"}:
            return col
    for col in columns:
        key = clean_name(col)
        if ("eigengene" in key or "module" in key or "kme" in key) and ("cor" in key or "membership" in key):
            return col
    return None


def read_candidate_file(path: Path) -> list[tuple[str, pd.DataFrame]]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return [(path.name, pd.read_csv(path))]
    if suffix in {".tsv", ".txt"}:
        return [(path.name, pd.read_csv(path, sep=None, engine="python"))]
    if suffix in {".xlsx", ".xls"}:
        sheets = pd.read_excel(path, sheet_name=None)
        return [(f"{path.name}:{sheet}", df) for sheet, df in sheets.items()]
    return []


def standardize_table(source_name: str, df: pd.DataFrame) -> pd.DataFrame | None:
    if df.empty:
        return None
    gene_col = find_column(df.columns, GENE_COLUMN_CANDIDATES)
    module_col = find_column(df.columns, MODULE_COLUMN_CANDIDATES)
    if gene_col is None or module_col is None:
        return None
    label_col = find_column(df.columns, LABEL_COLUMN_CANDIDATES)
    corr_col = find_eigengene_corr_column(df.columns)

    out = pd.DataFrame({
        "gene_symbol": df[gene_col].map(normalize_symbol),
        "published_module": df[module_col].astype(str).str.strip(),
        "module_label": df[label_col].astype(str).str.strip() if label_col else df[module_col].astype(str).str.strip(),
        "eigengene_correlation": pd.to_numeric(df[corr_col], errors="coerce") if corr_col else np.nan,
        "source_table": source_name,
        "source_reference": SOURCE_REFERENCE,
    })
    out = out[out["gene_symbol"].ne("") & out["published_module"].ne("")]
    out = out[~out["published_module"].str.lower().isin({"nan", "none", "na"})]
    out = out.drop_duplicates(["gene_symbol", "published_module"])
    if out.empty:
        return None
    return out


def candidate_input_files(explicit: str | None) -> list[Path]:
    if explicit:
        return [Path(explicit)]
    files: list[Path] = []
    for pattern in ("*.csv", "*.tsv", "*.txt", "*.xlsx", "*.xls"):
        files.extend(sorted(MANUAL_DIR.glob(pattern)))
    return [p for p in files if p.name.lower() != "readme.md"]


def maybe_download_module_table() -> Path | None:
    urls = []
    env_url = os.environ.get("AHBA_PUBLISHED_MODULE_URL", "").strip()
    if env_url:
        urls.append(env_url)
    urls.extend(DEFAULT_DOWNLOAD_URLS)

    MANUAL_DIR.mkdir(parents=True, exist_ok=True)
    for url in urls:
        filename = url.split("/")[-1] or "Hawrylycz_2015_Supplementary_Table_6.csv"
        if "MOESM97" in filename:
            filename = "Hawrylycz_2015_Supplementary_Table_6.csv"
        target = MANUAL_DIR / filename
        if target.exists():
            return target
        try:
            print(f"Attempting download of published AHBA module table: {url}")
            urllib.request.urlretrieve(url, target)
            return target
        except Exception as exc:  # noqa: BLE001 - keep fallback graceful for publisher/network issues.
            print(f"Download failed for {url}: {exc}")
    return None


def write_manual_notice(reason: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    notice = OUT_DIR / "PUBLISHED_AHBA_MODULES_NOT_LOADED.txt"
    notice.write_text(
        "Published AHBA canonical module assignments were not standardized yet.\n\n"
        f"Reason: {reason}\n\n"
        "Manual fallback:\n"
        "1. Download Hawrylycz et al. 2015 Supplementary Table 6.\n"
        f"2. Place the CSV/XLSX/TSV file in: {MANUAL_DIR}\n"
        "3. Run: python code/ahba_modules/04_fetch_or_load_published_ahba_modules.py\n"
        "4. Run: python code/ahba_modules/05_test_published_ahba_module_enrichment.py\n",
        encoding="utf-8",
    )
    print(notice.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="Optional explicit path to Hawrylycz 2015 Supplementary Table 6 CSV/XLSX/TSV file.")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MANUAL_DIR.mkdir(parents=True, exist_ok=True)

    files = candidate_input_files(args.input)
    if not files:
        downloaded = maybe_download_module_table()
        files = [downloaded] if downloaded else []
    if not files:
        write_manual_notice("No local Supplementary Table 6 file was found and automatic download did not succeed.")
        return

    standardized_parts: list[pd.DataFrame] = []
    inspected = []
    for path in files:
        if not path.exists():
            inspected.append({"file": str(path), "status": "missing"})
            continue
        try:
            tables = read_candidate_file(path)
        except Exception as exc:  # noqa: BLE001 - audit all failed readers.
            inspected.append({"file": str(path), "status": "read_failed", "error": str(exc)})
            continue
        matched = 0
        for source_name, df in tables:
            standardized = standardize_table(source_name, df)
            inspected.append({
                "file": str(path),
                "source_table": source_name,
                "rows": int(len(df)),
                "columns": list(map(str, df.columns)),
                "status": "matched" if standardized is not None else "no_gene_module_columns_detected",
            })
            if standardized is not None:
                standardized_parts.append(standardized)
                matched += 1
        print(f"Inspected {path}: matched {matched} table(s)")

    if not standardized_parts:
        (OUT_DIR / "published_module_load_attempt_audit.json").write_text(json.dumps(inspected, indent=2), encoding="utf-8")
        write_manual_notice("Files were present, but no table with detectable gene-symbol and module-assignment columns was found.")
        return

    out = pd.concat(standardized_parts, ignore_index=True)
    out = out.drop_duplicates(["gene_symbol", "published_module"])
    out = out.sort_values(["published_module", "gene_symbol"])
    out.to_csv(OUT_TABLE, index=False)

    module_sizes = out.groupby("published_module")["gene_symbol"].nunique().reset_index(name="module_gene_count")
    module_sizes = module_sizes.sort_values("module_gene_count", ascending=False)
    module_sizes.to_csv(OUT_DIR / "ahba_published_module_size_table.csv", index=False)

    report = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "source_reference": SOURCE_REFERENCE,
        "input_files": [str(p) for p in files],
        "standardized_table": str(OUT_TABLE),
        "rows": int(len(out)),
        "unique_gene_symbols": int(out["gene_symbol"].nunique()),
        "modules": int(out["published_module"].nunique()),
        "module_size_table": str(OUT_DIR / "ahba_published_module_size_table.csv"),
        "inspected_tables": inspected,
    }
    (OUT_DIR / "published_module_load_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
