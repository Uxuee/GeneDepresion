"""Preprocess AHBA microarray data for gene-module analysis.

This script expects donor-level AHBA files under data/ahba, typically:
MicroarrayExpression.csv, Probes.csv, and SampleAnnot.csv. It collapses probes
to one row per gene symbol by selecting the probe with the highest mean
expression across samples, then writes a gene x sample matrix.
"""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "ahba"
CLEAN_DIR = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
OUT_DIR = CLEAN_DIR / "11_AHBA_modules_enrichment" / "processed"
CANONICAL_ML_RANKING = CLEAN_DIR / "10_ML_prioritized_Tier2_corrected" / "ML_prioritized_Tier2_genes_all_ranked.csv"


def normalize_gene_symbol(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", "", str(value).strip().upper())


def find_donor_sets() -> list[dict[str, Path]]:
    expression_files = sorted(DATA_DIR.rglob("MicroarrayExpression.csv"))
    donor_sets: list[dict[str, Path]] = []
    for expr in expression_files:
        folder = expr.parent
        probes = folder / "Probes.csv"
        samples = folder / "SampleAnnot.csv"
        if probes.exists() and samples.exists():
            donor_sets.append({"donor": folder.name, "expression": expr, "probes": probes, "samples": samples})
    return donor_sets


def read_expression(path: Path) -> pd.DataFrame:
    """Read AHBA expression as probe x sample, tolerating header/no-header variants."""
    df = pd.read_csv(path, header=None)
    # Allen files are commonly no-header numeric probe x sample matrices.
    if df.shape[1] > 1 and pd.to_numeric(df.iloc[:, 0], errors="coerce").notna().mean() < 0.5:
        df = pd.read_csv(path)
    return df


def pick_col(columns: list[str], candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    for c in columns:
        lc = c.lower()
        if any(cand.lower() in lc for cand in candidates):
            return c
    return None


def process_donor(donor: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expr = read_expression(donor["expression"])
    probes = pd.read_csv(donor["probes"])
    samples = pd.read_csv(donor["samples"])

    gene_col = pick_col(list(probes.columns), ["gene_symbol", "gene symbol", "symbol", "gene_name"])
    probe_col = pick_col(list(probes.columns), ["probe_id", "probe id", "probe_name", "id"])
    if gene_col is None:
        raise ValueError(f"Could not identify gene-symbol column in {donor['probes']}")

    if len(expr) == len(probes):
        expr.index = probes[probe_col].astype(str).to_list() if probe_col else probes.index.astype(str)
    elif expr.shape[1] == len(probes):
        expr = expr.T
        expr.index = probes[probe_col].astype(str).to_list() if probe_col else probes.index.astype(str)
    else:
        raise ValueError(
            f"Expression dimensions {expr.shape} do not align with probe rows {len(probes)} for {donor['donor']}"
        )

    expr = expr.apply(pd.to_numeric, errors="coerce")
    sample_names = [f"{donor['donor']}_S{i+1}" for i in range(expr.shape[1])]
    expr.columns = sample_names

    mapping = pd.DataFrame(
        {
            "donor": donor["donor"],
            "probe_id": expr.index.astype(str),
            "gene_symbol_raw": probes[gene_col].astype(str).to_list(),
        }
    )
    mapping["gene_symbol"] = mapping["gene_symbol_raw"].map(normalize_gene_symbol)
    keep = mapping["gene_symbol"].ne("") & mapping["gene_symbol"].ne("NAN")
    expr = expr.loc[keep.to_numpy()]
    mapping = mapping.loc[keep].reset_index(drop=True)
    expr.index = mapping["gene_symbol"]

    sample_meta = samples.copy()
    sample_meta.insert(0, "ahba_sample_id", sample_names[: len(sample_meta)])
    sample_meta.insert(0, "donor", donor["donor"])
    if len(sample_meta) != len(sample_names):
        sample_meta = pd.DataFrame({"donor": donor["donor"], "ahba_sample_id": sample_names})

    return expr, mapping, sample_meta


def collapse_probes(combined: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    combined = combined.copy()
    combined["missing_fraction"] = combined.isna().mean(axis=1)
    expr_only = combined.drop(columns=["missing_fraction"])
    probe_table = expr_only.reset_index(drop=True)
    probe_table.insert(0, "gene_symbol", expr_only.index.astype(str))
    probe_table.insert(1, "probe_id", expr_only.index.astype(str))
    sample_columns = list(expr_only.columns)

    report_rows = []
    rows = []
    labels = []
    for gene, gene_block in probe_table.groupby("gene_symbol", sort=True):
        gene_block_expr = gene_block[sample_columns]
        probe_means = gene_block_expr.mean(axis=1, skipna=True)
        # Gene symbols are duplicated across probes; select by position within the block.
        best_position = int(probe_means.to_numpy().argmax())
        selected_expr = gene_block_expr.iloc[best_position]
        selected_probe = gene_block.iloc[best_position]["probe_id"]
        report_rows.append(
            {
                "gene_symbol": gene,
                "n_probes": int(len(gene_block)),
                "selected_probe_id": str(selected_probe),
                "selected_probe_mean_expression": float(probe_means.iloc[best_position]),
            }
        )
        rows.append(selected_expr)
        labels.append(gene)
    collapsed = pd.DataFrame(rows)
    collapsed.index = labels
    return collapsed, pd.DataFrame(report_rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    donor_sets = find_donor_sets()
    if not donor_sets:
        msg = (
            "No donor-level AHBA MicroarrayExpression.csv/Probes.csv/SampleAnnot.csv sets were found. "
            "Run code/ahba_modules/00_fetch_ahba_data.py or manually place AHBA files under data/ahba."
        )
        (OUT_DIR / "AHBA_PREPROCESSING_NOT_RUN.txt").write_text(msg + "\n", encoding="utf-8")
        print(msg)
        return

    exprs = []
    mappings = []
    metas = []
    for donor in donor_sets:
        expr, mapping, meta = process_donor(donor)
        exprs.append(expr)
        mappings.append(mapping)
        metas.append(meta)

    combined = pd.concat(exprs, axis=1, join="outer")
    missing_fraction = combined.isna().mean(axis=1)
    filtered = combined.loc[missing_fraction <= 0.50].copy()
    collapsed, collapse_report = collapse_probes(filtered)
    collapsed = collapsed.sort_index()

    final_ranked = pd.read_csv(CANONICAL_ML_RANKING, dtype={"gene_id": str})
    mdd_symbols = set(final_ranked["gene_symbol"].dropna().map(normalize_gene_symbol))
    ahba_symbols = set(collapsed.index)
    overlap = sorted(mdd_symbols & ahba_symbols)

    collapsed.to_csv(OUT_DIR / "ahba_gene_expression_matrix_gene_by_sample.csv")
    pd.concat(mappings, ignore_index=True).to_csv(OUT_DIR / "ahba_probe_gene_mapping_long.csv", index=False)
    pd.concat(metas, ignore_index=True).to_csv(OUT_DIR / "ahba_sample_metadata.csv", index=False)
    collapse_report.to_csv(OUT_DIR / "ahba_probe_collapse_report.csv", index=False)

    report = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "n_donor_sets": len(donor_sets),
        "donors": [d["donor"] for d in donor_sets],
        "n_genes_after_missingness_and_probe_collapse": int(collapsed.shape[0]),
        "n_samples": int(collapsed.shape[1]),
        "missingness_filter": "genes/probes with >50% missing values removed before probe collapse",
        "probe_collapse_method": "select one probe per gene symbol by highest mean expression across AHBA samples",
        "n_unique_mdd_candidate_symbols": len(mdd_symbols),
        "n_mdd_candidate_symbols_present_in_ahba": len(overlap),
        "mdd_candidate_symbol_overlap_file": str(OUT_DIR / "mdd_candidate_symbols_present_in_ahba.txt"),
    }
    (OUT_DIR / "ahba_preprocessing_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    (OUT_DIR / "mdd_candidate_symbols_present_in_ahba.txt").write_text("\n".join(overlap) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
