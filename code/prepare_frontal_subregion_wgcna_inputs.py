#!/usr/bin/env python
"""Create frontal subregion GeTMM-log2 matrices for WGCNA.

The clean pipeline currently stores one broad Frontal Cortex normalized matrix.
This helper splits that matrix by the sample metadata OriginalRegion labels so
OFC BA11, dlPFC BA8/9, and Cingulate gyrus 25 can be analyzed separately without
overwriting the existing broad Frontal Cortex WGCNA outputs.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path("D:/genes/GeneDepresion")
PIPELINE_ROOT = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
EXPR_DIR = PIPELINE_ROOT / "02_normalized_expression"
OUT_DIR = EXPR_DIR / "frontal_subregion_wgcna_inputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FRONTAL_MATRIX = EXPR_DIR / "normalized_GeTMM_log2_frontal_cortex.csv"
METADATA_FILE = PIPELINE_ROOT / "01_inputs" / "sample_metadata_used.csv"

SUBREGIONS = {
    "OFC_BA11": {
        "label": "OFC BA11",
        "patterns": ("Orbitofrontal", "OFC", "BA11"),
        "output": "normalized_GeTMM_log2_OFC_BA11.csv",
    },
    "dlPFC_BA8_9": {
        "label": "dlPFC BA8/9",
        "patterns": ("Dorsolateral prefrontal", "dlPFC", "BA8", "BA9"),
        "output": "normalized_GeTMM_log2_dlPFC_BA8_9.csv",
    },
    "CG25": {
        "label": "Cingulate gyrus 25 / CG25",
        "patterns": ("Cingulate gyrus 25", "Cg25", "CG25", "BA25"),
        "output": "normalized_GeTMM_log2_CG25.csv",
    },
}


def contains_any(value: object, patterns: tuple[str, ...]) -> bool:
    text = str(value).lower()
    return any(pattern.lower() in text for pattern in patterns)


def main() -> None:
    if not FRONTAL_MATRIX.exists():
        raise FileNotFoundError(FRONTAL_MATRIX)
    if not METADATA_FILE.exists():
        raise FileNotFoundError(METADATA_FILE)

    metadata = pd.read_csv(METADATA_FILE)
    sample_col = "SampleID" if "SampleID" in metadata.columns else "geo_accession"
    region_col = "OriginalRegion" if "OriginalRegion" in metadata.columns else "tissue"

    frontal = pd.read_csv(FRONTAL_MATRIX)
    gene_col = frontal.columns[0]
    matrix_samples = set(frontal.columns[1:])

    manifest_rows: list[dict[str, object]] = []
    for key, spec in SUBREGIONS.items():
        mask = metadata[region_col].apply(lambda x: contains_any(x, spec["patterns"]))
        meta_sub = metadata.loc[mask].copy()
        samples = [s for s in meta_sub[sample_col].astype(str).tolist() if s in matrix_samples]
        out_path = OUT_DIR / spec["output"]
        meta_path = OUT_DIR / spec["output"].replace(".csv", "_sample_metadata.csv")

        if samples:
            frontal[[gene_col] + samples].to_csv(out_path, index=False)
            meta_sub[meta_sub[sample_col].astype(str).isin(samples)].to_csv(meta_path, index=False)

        manifest_rows.append(
            {
                "subregion_key": key,
                "subregion_label": spec["label"],
                "patterns": ";".join(spec["patterns"]),
                "found": bool(samples),
                "n_samples": len(samples),
                "n_genes": len(frontal),
                "matrix": str(out_path) if samples else "",
                "sample_metadata": str(meta_path) if samples else "",
                "note": (
                    "Derived from broad Frontal Cortex normalized GeTMM-log2 matrix using metadata OriginalRegion."
                    if samples
                    else "No matching samples found in metadata/matrix."
                ),
            }
        )

    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "source_matrix": str(FRONTAL_MATRIX),
        "metadata_file": str(METADATA_FILE),
        "output_dir": str(OUT_DIR),
        "subregions": manifest_rows,
    }
    pd.DataFrame(manifest_rows).to_csv(OUT_DIR / "frontal_subregion_wgcna_input_manifest.csv", index=False)
    (OUT_DIR / "frontal_subregion_wgcna_input_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print(pd.DataFrame(manifest_rows).to_string(index=False))


if __name__ == "__main__":
    main()
