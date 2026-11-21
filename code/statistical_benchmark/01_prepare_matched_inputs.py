from __future__ import annotations

import io
import json
import subprocess

import numpy as np
import pandas as pd

from benchmark_common import (
    ANNOTATION,
    BROAD_LM,
    EXPECTED_GENES,
    EXPECTED_SAMPLES,
    EXPECTED_TIER2,
    FRONTAL_BLOB,
    INPUTS,
    LEGACY_FILTERED_STEP1,
    LEGACY_REPO,
    OUT,
    REGIONS,
    SLUGS,
    STEP1,
    TIER2,
    ensure_dirs,
    sha256,
)

BROAD_SOURCE = {
    "Ventral Striatum": ("Ventral_Striatum", "ventral_striatum"),
    "Insula": ("Insula", "insula"),
    "Hippocampal Formation": ("Hippocampal_Formation", "hippocampal_formation"),
}

FRONTAL_SOURCE = {
    "CG25": "Cingulate gyrus 25 (Cg25)",
    "dlPFC BA8/9": "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
    "OFC BA11": "Orbitofrontal (OFC; BA11)",
}


def read_matrix(path):
    frame = pd.read_csv(path, dtype={0: str}, low_memory=False)
    gene_col = frame.columns[0]
    frame = frame.rename(columns={gene_col: "GeneID"})
    frame["GeneID"] = frame["GeneID"].astype(str)
    if frame["GeneID"].duplicated().any():
        raise ValueError(f"Duplicate GeneID in {path}")
    return frame.set_index("GeneID")


def read_metadata(path):
    meta = pd.read_csv(path, dtype={"SampleID": str, "BrainID": str}, low_memory=False)
    required = ["SampleID", "Diagnosis", "Age", "Sex", "RIN"]
    missing = [column for column in required if column not in meta]
    if missing:
        raise ValueError(f"Missing metadata columns {missing} in {path}")
    meta["SampleID"] = meta["SampleID"].astype(str)
    meta["Diagnosis"] = meta["Diagnosis"].replace({"Control": "CTRL", "CON": "CTRL"}).astype(str)
    meta["Age"] = pd.to_numeric(meta["Age"], errors="coerce")
    meta["RIN"] = pd.to_numeric(meta["RIN"], errors="coerce")
    meta = meta.dropna(subset=["Diagnosis", "Age", "Sex", "RIN"]).copy()
    if set(meta["Diagnosis"]) != {"CTRL", "MDD"}:
        raise ValueError(f"Unexpected Diagnosis values in {path}: {sorted(meta['Diagnosis'].unique())}")
    if meta["SampleID"].duplicated().any():
        raise ValueError(f"Duplicate SampleID in {path}")
    return meta


def load_frontal_lm():
    completed = subprocess.run(
        ["git", "-C", str(LEGACY_REPO), "cat-file", "blob", FRONTAL_BLOB],
        check=True,
        capture_output=True,
    )
    frame = pd.read_csv(io.BytesIO(completed.stdout), dtype={"gene_id": str, "GeneID": str}, low_memory=False)
    return frame


def normalize_lm(frame, region, region_column=None, region_value=None):
    data = frame.copy()
    if region_column is not None:
        data = data[data[region_column].astype(str).eq(region_value)].copy()
    rename = {
        "GeneID": "gene_id",
        "Symbol": "gene_symbol",
        "diagnosis_effect": "beta_diagnosis",
        "diagnosis_p_value": "p_value_diagnosis",
        "diagnosis_FDR_BH": "fdr_BH",
    }
    data = data.rename(columns={key: value for key, value in rename.items() if key in data and value not in data})
    required = ["gene_id", "gene_symbol", "beta_diagnosis", "p_value_diagnosis", "fdr_BH", "standardized_effect_size_diagnosis"]
    missing = [column for column in required if column not in data]
    if missing:
        raise ValueError(f"Frozen LM for {region} lacks {missing}; columns={list(data.columns)}")
    data["gene_id"] = data["gene_id"].astype(str)
    data.insert(0, "region", region)
    ordered = data.sort_values(["p_value_diagnosis", "gene_id"], kind="mergesort")
    data["lm_primary_rank"] = pd.Series(np.arange(1, len(data) + 1), index=ordered.index).reindex(data.index)
    return data[["region", "gene_id", "gene_symbol", "beta_diagnosis", "p_value_diagnosis", "fdr_BH", "standardized_effect_size_diagnosis", "lm_primary_rank"]]


def main():
    ensure_dirs()
    tier2 = pd.read_csv(TIER2, dtype={"gene_id": str}, low_memory=False)
    tier2["region"] = tier2["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
    tier2 = tier2[tier2["region"].isin(REGIONS)].copy()
    if len(tier2) != 671 or tier2.groupby("region").size().to_dict() != EXPECTED_TIER2:
        raise ValueError("Frozen Tier2 universe does not match 671-pair validation target")

    broad_lm = pd.read_csv(BROAD_LM, dtype={"gene_id": str}, low_memory=False)
    frontal_lm = load_frontal_lm()
    all_lm = []
    sample_rows = []
    gene_rows = []

    frontal_counts = read_matrix(LEGACY_FILTERED_STEP1 / "03_filtered_raw_counts" / "keep_all_samples" / "filtered_raw_counts_frontal_cortex.csv")
    frontal_meta_path = STEP1 / "01_metadata" / "metadata_clean_matched_all_samples_with_brain_area.csv"
    frontal_meta_raw = pd.read_csv(frontal_meta_path, dtype={"SampleID": str, "BrainID": str}, low_memory=False)
    frontal_meta = read_metadata(frontal_meta_path)

    for region in REGIONS:
        slug = SLUGS[region]
        if region in BROAD_SOURCE:
            lm_name, file_stub = BROAD_SOURCE[region]
            counts = read_matrix(LEGACY_FILTERED_STEP1 / "03_filtered_raw_counts" / "keep_all_samples" / f"filtered_raw_counts_{file_stub}.csv")
            metadata_path = STEP1 / "01_metadata" / "keep_all_samples" / f"metadata_{file_stub}.csv"
            meta_raw = pd.read_csv(metadata_path, dtype={"SampleID": str}, low_memory=False)
            meta = read_metadata(metadata_path)
            lm = normalize_lm(broad_lm[broad_lm["brain_region"].eq(lm_name)], region)
        else:
            original = FRONTAL_SOURCE[region]
            meta_raw = frontal_meta_raw[frontal_meta_raw["OriginalRegion"].astype(str).eq(original)].copy()
            meta = frontal_meta[frontal_meta["OriginalRegion"].astype(str).eq(original)].copy()
            counts = frontal_counts.copy()
            region_col = "original_frontal_subregion" if "original_frontal_subregion" in frontal_lm else "OriginalRegion"
            lm = normalize_lm(frontal_lm, region, region_col, original)

        samples = meta["SampleID"].tolist()
        missing_samples = sorted(set(samples) - set(counts.columns.astype(str)))
        extra_lm_genes = sorted(set(lm["gene_id"]) - set(counts.index))
        missing_lm_genes = sorted(set(counts.index) - set(lm["gene_id"]))
        if missing_samples or extra_lm_genes or missing_lm_genes:
            raise ValueError(
                f"Unsafe match for {region}: missing_samples={len(missing_samples)}, "
                f"LM_not_counts={len(extra_lm_genes)}, counts_not_LM={len(missing_lm_genes)}"
            )
        counts = counts.loc[lm["gene_id"], samples]
        numeric = counts.apply(pd.to_numeric, errors="raise")
        if not np.equal(numeric.to_numpy(), np.floor(numeric.to_numpy())).all() or (numeric.to_numpy() < 0).any():
            raise ValueError(f"Non-integer or negative counts for {region}")
        if len(meta) != EXPECTED_SAMPLES[region] or len(counts) != EXPECTED_GENES[region]:
            raise ValueError(f"Expected sample/gene count failed for {region}: {len(meta)}, {len(counts)}")

        numeric.to_csv(INPUTS / f"{slug}_raw_counts.csv.gz", compression="gzip")
        meta[["SampleID", "Diagnosis", "Age", "Sex", "RIN"]].to_csv(INPUTS / f"{slug}_metadata.csv", index=False)
        lm.to_csv(INPUTS / f"{slug}_frozen_lm.csv.gz", index=False, compression="gzip")
        all_lm.append(lm)

        sample_rows.append({
            "region": region,
            "total_raw_count_samples": len(counts.columns),
            "samples_in_frozen_LM": len(meta),
            "samples_available_limma": len(meta),
            "samples_available_DESeq2": len(meta),
            "excluded_samples": int(len(meta_raw) - len(meta)),
            "exclusion_reason": "missing Diagnosis/Age/Sex/RIN" if len(meta_raw) != len(meta) else "none",
            "exact_matched_set": True,
        })
        gene_rows.append({
            "region": region,
            "frozen_LM_eligible_genes": len(lm),
            "matched_raw_count_genes": len(counts),
            "genes_tested_limma": len(counts),
            "genes_tested_DESeq2": len(counts),
            "unmatched_genes": 0,
            "unmatched_reason": "none",
            "exact_matched_universe": True,
        })

    pd.concat(all_lm, ignore_index=True).to_csv(INPUTS / "frozen_lm_reference_all_regions.csv.gz", index=False, compression="gzip")
    pd.DataFrame(sample_rows).to_csv(OUT / "sample_accounting.csv", index=False)
    pd.DataFrame(gene_rows).to_csv(OUT / "gene_universe_accounting.csv", index=False)
    source_manifest = {
        "raw_counts_source": str(LEGACY_FILTERED_STEP1 / "03_filtered_raw_counts" / "keep_all_samples"),
        "raw_master_source": str(STEP1.parent / "GSE102556_raw_counts_GRCh38.p13_NCBI.tsv"),
        "annotation_source": str(ANNOTATION),
        "broad_lm_source": str(BROAD_LM),
        "broad_lm_sha256": sha256(BROAD_LM),
        "frontal_lm_git_blob": FRONTAL_BLOB,
        "tier2_source": str(TIER2),
        "tier2_sha256": sha256(TIER2),
        "matched_sample_sets": True,
        "matched_gene_universes": True,
        "count_input_integer": True,
        "n_tier2_pairs": 671,
        "frozen_checksum_ledgers": {
            str(path.relative_to(OUT.parents[2])): sha256(path)
            for path in [
                OUT.parent / "18_Tier2_ML_RRA_prioritization" / "checksums.sha256",
                OUT.parent / "19_Tier2_ML_benchmark" / "checksums.sha256",
                OUT.parent / "20_Tier2_ML_conditional_null" / "checksums.sha256",
            ]
            if path.exists()
        },
    }
    (INPUTS / "input_provenance.json").write_text(json.dumps(source_manifest, indent=2) + "\n", encoding="utf-8")
    print("Prepared exact matched inputs for six regions")


if __name__ == "__main__":
    main()
