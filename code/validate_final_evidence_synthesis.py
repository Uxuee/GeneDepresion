from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
BASE = REPO / "outputs" / "length_gt_100_clean_pipeline"
OUT = BASE / "22_Final_evidence_synthesis"
MATRIX = OUT / "tables" / "Tier2_final_evidence_matrix.csv"
EXPECTED = {"Ventral Striatum": 146, "Insula": 70, "Hippocampal Formation": 127, "CG25": 186, "dlPFC BA8/9": 72, "OFC BA11": 70}
LEDGERS = {
    BASE / "17_ML_RRA_consensus" / "checksums.sha256": "117c65e6cc796a9c6757b9d9e92ed493cc8bb6300d0873ea3d5cd2801ba3634c",
    BASE / "18_Tier2_ML_RRA_prioritization" / "checksums.sha256": "3c45b0a92e4a5a3b8c9fc126097c312d492141734ea8e3fe72534a91152e950e",
    BASE / "19_Tier2_ML_benchmark" / "checksums.sha256": "4e17dc8234317efc9bce92aa03fb7ae63ff72a128392d202f2f31a105c09cf10",
    BASE / "20_Tier2_ML_conditional_null" / "checksums.sha256": "d95f66712ed0ba80b7acefee7c3ee992cc14923128a47a5ddcaf9cff12e0d545",
    BASE / "21_Statistical_method_benchmark" / "checksums.sha256": "5626960b44eb446c023fd621f2808916f591667666be631a89691237a1b025ed",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype={"gene_id": str})
    if "gene_id" in frame:
        frame["gene_id"] = frame.gene_id.str.replace(r"\.0$", "", regex=True)
    return frame


def run_checks(check_determinism: bool = False) -> dict[str, bool]:
    matrix = read(MATRIX)
    source = read(BASE / "21_Statistical_method_benchmark" / "tier2_statistical_method_audit.csv")
    provenance = pd.read_csv(OUT / "manifest" / "evidence_provenance.csv")
    fdr = read(OUT / "tables" / "method_specific_FDR_findings.csv")
    rra_null = read(BASE / "20_Tier2_ML_conditional_null" / "summary" / "RRA_observed_topgene_null_behavior.csv")
    source_keys = set(zip(source.region, source.gene_id))
    matrix_keys = set(zip(matrix.region, matrix.gene_id))
    null_keys = set(zip(rra_null.region, rra_null.gene_id))
    filled_null_keys = set(zip(matrix.loc[matrix.RRA_null_rank_median.notna(), "region"], matrix.loc[matrix.RRA_null_rank_median.notna(), "gene_id"]))

    checks = {
        "01_exact_671_rows": len(matrix) == 671 and not matrix.duplicated(["region", "gene_id"]).any(),
        "02_regional_counts": matrix.groupby("region").size().to_dict() == EXPECTED,
        "03_tier2_identities_unchanged": matrix_keys == source_keys and matrix.Tier2_membership.astype(bool).all(),
        "04_frozen_upstream_ledgers_unchanged": all(sha256(path) == expected for path, expected in LEDGERS.items()),
        "05_every_column_has_provenance": set(matrix.columns) == set(provenance.column_name) and provenance.column_name.is_unique,
        "06_no_weighted_score": not any("weighted" in column.lower() or "composite" in column.lower() for column in matrix.columns),
        "07_no_new_confidence_category": not any("confidence" in column.lower() or "validated" in column.lower() for column in matrix.columns),
        "08_rank_direction_consistent": all(sorted(group.Tier2_RRA_consensus_rank.astype(int)) == list(range(1, len(group) + 1)) for _, group in matrix.groupby("region")),
        "09_missing_null_evidence_not_invented": filled_null_keys == null_keys and matrix.AHBA_gene_level_context.isna().all(),
        "10_allgene_labeled_sensitivity": provenance.loc[provenance.column_name.str.startswith("allgene_RRA"), "interpretation"].str.contains("sensitivity", case=False).all(),
        "11_conditional_null_labeled_conditional": provenance.loc[provenance.column_name.str.contains("null", case=False), "source_analysis"].str.contains("conditional", case=False).all(),
        "12_wgcna_ahba_context_only": matrix.WGCNA_context_caveat.dropna().str.contains("context only").all() and matrix.AHBA_context_caveat.str.contains("excluded from prioritisation").all(),
        "13_method_specific_fdr_separate": len(fdr) == 10 and set(fdr.method) == {"limma-voom", "DESeq2 same-universe BH"} and fdr.finding_label.eq("method-specific sensitivity finding").all(),
        "14_required_compact_outputs_exist": all(path.exists() for path in [OUT / "tables" / "Tier2_top10_by_region_evidence.csv", OUT / "tables" / "Tier2_top25_by_region_evidence.csv", OUT / "checksums.sha256", REPO / "docs" / "FINAL_CROSS_METHOD_EVIDENCE_SYNTHESIS.md"]),
        "15_upstream_sources_read_only_by_design": all(str(path).find("22_Final_evidence_synthesis") < 0 for path in LEDGERS),
    }
    if check_determinism:
        before = sha256(OUT / "checksums.sha256")
        subprocess.run([sys.executable, str(REPO / "code" / "final_evidence_synthesis.py")], cwd=REPO, check=True, capture_output=True, text=True)
        checks["16_deterministic_regeneration"] = sha256(OUT / "checksums.sha256") == before
    return checks


def main() -> None:
    checks = run_checks(check_determinism=True)
    for name, passed in checks.items():
        print(f"{'PASS' if passed else 'FAIL'} {name}")
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise SystemExit(f"Validation failed: {failed}")
    print(f"Validation passed: {len(checks)} checks")


if __name__ == "__main__":
    main()
