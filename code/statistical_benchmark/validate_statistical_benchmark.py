from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from benchmark_common import EXPECTED_GENES, EXPECTED_SAMPLES, EXPECTED_TIER2, INPUTS, OUT, REGIONS, REPO, SLUGS, TIER2, sha256


def rank_primary(frame, p_column):
    ordered = frame.sort_values([p_column, "gene_id"], na_position="last", kind="mergesort")
    ranks = pd.Series(np.arange(1, len(frame) + 1), index=ordered.index)
    return ranks.reindex(frame.index).astype(int)


def collect_checks():
    checks = {}
    sample = pd.read_csv(OUT / "sample_accounting.csv")
    genes = pd.read_csv(OUT / "gene_universe_accounting.csv")
    tier2 = pd.read_csv(TIER2, dtype={"gene_id": str})
    tier2["region"] = tier2.region_or_subregion.replace({"Cingulate gyrus 25": "CG25"})
    tier2 = tier2[tier2.region.isin(REGIONS)]
    provenance = json.loads((INPUTS / "input_provenance.json").read_text(encoding="utf-8"))
    manifest = json.loads((OUT / "statistical_benchmark_manifest.json").read_text(encoding="utf-8"))
    script = (REPO / "code" / "statistical_benchmark" / "02_run_limma_deseq2.R").read_text(encoding="utf-8")
    limma_all = pd.read_csv(OUT / "limma_voom_results_all_regions.csv", dtype={"gene_id": str}, low_memory=False)
    deseq_all = pd.read_csv(OUT / "deseq2_results_all_regions.csv", dtype={"gene_id": str}, low_memory=False)
    correlations = pd.read_csv(OUT / "method_rank_correlations.csv")
    directions = pd.read_csv(OUT / "method_direction_agreement.csv")

    checks["01_exact_six_regions"] = list(sample.region) == REGIONS and set(limma_all.region) == set(REGIONS) and set(deseq_all.region) == set(REGIONS)
    checks["02_tier2_671_exact_counts"] = len(tier2) == 671 and tier2.groupby("region").size().to_dict() == EXPECTED_TIER2
    checks["03_frozen_tier2_hash_unchanged"] = sha256(TIER2) == provenance["tier2_sha256"] == manifest["tier2_sha256"]
    checks["04_frozen_LM_hash_unchanged"] = sha256(Path(provenance["broad_lm_source"])) == provenance["broad_lm_sha256"]

    integer_ok = True
    sample_ok = True
    gene_ok = True
    ids_ok = True
    ranks_ok = True
    for region in REGIONS:
        slug = SLUGS[region]
        counts = pd.read_csv(INPUTS / f"{slug}_raw_counts.csv.gz", index_col=0)
        metadata = pd.read_csv(INPUTS / f"{slug}_metadata.csv", dtype={"SampleID": str})
        lm = pd.read_csv(INPUTS / f"{slug}_frozen_lm.csv.gz", dtype={"gene_id": str})
        limma = limma_all[limma_all.region.eq(region)].copy()
        deseq = deseq_all[deseq_all.region.eq(region)].copy()
        integer_ok &= bool((counts.to_numpy() >= 0).all() and np.equal(counts.to_numpy(), np.floor(counts.to_numpy())).all())
        sample_ok &= list(counts.columns.astype(str)) == metadata.SampleID.tolist() and len(metadata) == EXPECTED_SAMPLES[region]
        gene_ok &= len(counts) == EXPECTED_GENES[region] and set(counts.index.astype(str)) == set(lm.gene_id) == set(limma.gene_id) == set(deseq.gene_id)
        ids_ok &= not lm.gene_id.duplicated().any() and not limma.gene_id.duplicated().any() and not deseq.gene_id.duplicated().any()
        ranks_ok &= np.array_equal(rank_primary(limma, "raw_p").to_numpy(), limma.primary_rank.astype(int).to_numpy())
        ranks_ok &= np.array_equal(rank_primary(deseq, "raw_p").to_numpy(), deseq.primary_rank.astype(int).to_numpy())

    checks["05_DESeq2_raw_integer_counts"] = integer_ok
    checks["06_limma_DGEList_from_raw_counts"] = "DGEList(counts = counts)" in script and "raw_counts.csv.gz" in script
    checks["07_no_GeTMM_for_DESeq2"] = "GeTMM" not in script and "DESeqDataSetFromMatrix" in script
    checks["08_no_GeTMM_for_voom"] = "GeTMM" not in script and "limma::voom" in script
    checks["09_same_matched_samples"] = sample_ok and sample.exact_matched_set.all() and dict(zip(sample.region, sample.samples_in_frozen_LM)) == EXPECTED_SAMPLES
    checks["10_same_eligible_genes"] = gene_ok and genes.exact_matched_universe.all() and dict(zip(genes.region, genes.frozen_LM_eligible_genes)) == EXPECTED_GENES
    checks["11_contrast_MDD_minus_CTRL"] = "levels = c(\"CTRL\", \"MDD\")" in script and "c(\"Diagnosis\", \"MDD\", \"CTRL\")" in script
    checks["12_positive_effect_is_higher_MDD"] = checks["11_contrast_MDD_minus_CTRL"] and manifest["design"].endswith("MDD minus CTRL")
    checks["13_age_sex_RIN_in_design"] = "~ Diagnosis + Age + Sex + RIN" in script and "~ Age + Sex + RIN + Diagnosis" in script
    checks["14_no_Tier2_in_fitting"] = "tier2" not in script.lower()
    checks["15_no_ML_rank_in_fitting"] = "ml_priority" not in script.lower() and "rra" not in script.lower()
    checks["16_BH_denominator_documented"] = "same eligible tested gene universe" in (REPO / "docs" / "UPSTREAM_STATISTICAL_BENCHMARK_DESIGN.md").read_text(encoding="utf-8")
    checks["17_DESeq2_native_and_same_BH_separate"] = {"native_padj", "BH_same_universe"}.issubset(deseq_all.columns)
    checks["18_rank_definitions_reproduce"] = ranks_ok and manifest["primary_rank"].startswith("two-sided raw p")
    checks["19_gene_ids_one_to_one"] = ids_ok and (genes.unmatched_genes == 0).all()
    checks["20_frozen_checksum_ledgers_unchanged"] = all(sha256(REPO / rel) == digest for rel, digest in provenance["frozen_checksum_ledgers"].items())

    reproduced = True
    for _, row in correlations.iterrows():
        region = row.region
        lm = pd.read_csv(INPUTS / f"{SLUGS[region]}_frozen_lm.csv.gz", dtype={"gene_id": str})[["gene_id", "lm_primary_rank"]]
        l = limma_all[limma_all.region.eq(region)][["gene_id", "primary_rank"]]
        d = deseq_all[deseq_all.region.eq(region)][["gene_id", "primary_rank"]]
        mapping = {"Frozen LM": lm.rename(columns={"lm_primary_rank": "rank"}), "limma-voom": l.rename(columns={"primary_rank": "rank"}), "DESeq2": d.rename(columns={"primary_rank": "rank"})}
        merged = mapping[row.method_A].merge(mapping[row.method_B], on="gene_id", suffixes=("_a", "_b"))
        reproduced &= np.isclose(stats.spearmanr(merged.rank_a, merged.rank_b).statistic, row.spearman_rho, atol=1e-12)
    checks["21_summaries_reproduce_from_saved_results"] = reproduced and len(directions) == 54
    checks["22_manifest_complete"] = manifest["status"] == "complete" and manifest["candidate_tier_changed"] is False
    return checks


def main():
    checks = collect_checks()
    for name, passed in checks.items():
        print(("PASS" if passed else "FAIL"), name)
    if not all(checks.values()):
        raise SystemExit(1)
    print(f"Validation passed: {len(checks)} checks")


if __name__ == "__main__":
    main()
