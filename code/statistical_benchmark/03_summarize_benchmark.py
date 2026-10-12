from __future__ import annotations

import json
import math
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from benchmark_common import FIGURES, INPUTS, OUT, REGIONS, SLUGS, SUMMARY, TIER2, ensure_dirs, sha256

METHODS = ["Frozen LM", "limma-voom", "DESeq2"]
PAIRS = list(combinations(METHODS, 2))
TOP_DEFINITIONS = [("top10", 10), ("top25", 25), ("top50", 50), ("top1pct", 0.01), ("top5pct", 0.05)]
LEADERS = {
    "Ventral Striatum": "MORN5",
    "Insula": "TMEM161B-DT",
    "Hippocampal Formation": "TOR2A",
    "CG25": "LOC107985795",
    "dlPFC BA8/9": "LINC00269",
    "OFC BA11": "LOC105373386",
}


def bh(values):
    values = np.asarray(values, dtype=float)
    result = np.full(len(values), np.nan)
    valid = np.isfinite(values)
    pv = values[valid]
    order = np.argsort(pv, kind="mergesort")
    adjusted = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty(len(pv))
    restored[order] = np.minimum(adjusted, 1)
    result[valid] = restored
    return result


def save_figure(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


def load_data():
    lm = pd.read_csv(INPUTS / "frozen_lm_reference_all_regions.csv.gz", dtype={"gene_id": str}, low_memory=False)
    limma = pd.read_csv(OUT / "limma_voom_results_all_regions.csv", dtype={"gene_id": str}, low_memory=False)
    deseq = pd.read_csv(OUT / "deseq2_results_all_regions.csv", dtype={"gene_id": str}, low_memory=False)
    symbols = lm[["region", "gene_id", "gene_symbol"]].drop_duplicates()
    limma = limma.drop(columns=["gene_symbol"], errors="ignore").merge(symbols, on=["region", "gene_id"], how="left")
    deseq = deseq.drop(columns=["gene_symbol"], errors="ignore").merge(symbols, on=["region", "gene_id"], how="left")
    limma.to_csv(OUT / "limma_voom_results_all_regions.csv", index=False)
    deseq.to_csv(OUT / "deseq2_results_all_regions.csv", index=False)
    tier2 = pd.read_csv(TIER2, dtype={"gene_id": str}, low_memory=False)
    tier2["region"] = tier2["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
    tier2 = tier2[tier2["region"].isin(REGIONS)].copy()
    return lm, limma, deseq, tier2


def method_frame(lm, limma, deseq, region):
    a = lm[lm.region.eq(region)].copy().rename(columns={
        "beta_diagnosis": "Frozen LM_effect", "p_value_diagnosis": "Frozen LM_p",
        "fdr_BH": "Frozen LM_BH", "lm_primary_rank": "Frozen LM_rank"})
    b = limma[limma.region.eq(region)].copy().rename(columns={
        "logFC": "limma-voom_effect", "raw_p": "limma-voom_p", "BH": "limma-voom_BH",
        "primary_rank": "limma-voom_rank", "signed_evidence": "limma-voom_signed"})
    c = deseq[deseq.region.eq(region)].copy().rename(columns={
        "log2FoldChange": "DESeq2_effect", "raw_p": "DESeq2_p", "BH_same_universe": "DESeq2_BH",
        "native_padj": "DESeq2_native_padj", "primary_rank": "DESeq2_rank", "signed_evidence": "DESeq2_signed"})
    merged = a[["region", "gene_id", "gene_symbol", "Frozen LM_effect", "Frozen LM_p", "Frozen LM_BH", "Frozen LM_rank", "standardized_effect_size_diagnosis"]]
    merged = merged.merge(b[["gene_id", "limma-voom_effect", "limma-voom_p", "limma-voom_BH", "limma-voom_rank", "limma-voom_signed"]], on="gene_id", validate="one_to_one")
    merged = merged.merge(c[["gene_id", "DESeq2_effect", "DESeq2_p", "DESeq2_BH", "DESeq2_native_padj", "DESeq2_rank", "DESeq2_signed"]], on="gene_id", validate="one_to_one")
    merged["Frozen LM_signed"] = np.sign(merged["Frozen LM_effect"]) * -np.log10(merged["Frozen LM_p"].clip(lower=np.nextafter(0, 1)))
    return merged


def direction_and_rank(all_merged, tier_keys):
    direction_rows, rank_rows, overlap_rows = [], [], []
    for region, merged in all_merged.items():
        region_tier = {gene for reg, gene in tier_keys if reg == region}
        n = len(merged)
        for first, second in PAIRS:
            valid = merged[[f"{first}_effect", f"{second}_effect"]].notna().all(axis=1)
            same = np.sign(merged.loc[valid, f"{first}_effect"]) == np.sign(merged.loc[valid, f"{second}_effect"])
            scopes = {
                "full_shared": valid,
                "frozen_Tier2": valid & merged.gene_id.isin(region_tier),
                "top50_union": valid & ((merged[f"{first}_rank"] <= 50) | (merged[f"{second}_rank"] <= 50)),
            }
            for scope, mask in scopes.items():
                signs = np.sign(merged.loc[mask, f"{first}_effect"]) == np.sign(merged.loc[mask, f"{second}_effect"])
                direction_rows.append({
                    "region": region, "method_A": first, "method_B": second, "scope": scope,
                    "n_compared": int(mask.sum()), "same_sign_n": int(signs.sum()),
                    "opposite_sign_n": int((~signs).sum()), "same_sign_proportion": float(signs.mean()) if len(signs) else np.nan,
                })
            ranks = merged[[f"{first}_rank", f"{second}_rank"]].dropna()
            rank_rows.append({
                "region": region, "method_A": first, "method_B": second, "rank_definition": "two_sided_raw_p_GeneID_tiebreak",
                "n_genes": len(ranks), "spearman_rho": stats.spearmanr(ranks.iloc[:, 0], ranks.iloc[:, 1]).statistic,
                "kendall_tau": stats.kendalltau(ranks.iloc[:, 0], ranks.iloc[:, 1]).statistic,
            })
            for label, definition in TOP_DEFINITIONS:
                size = definition if isinstance(definition, int) else max(1, math.ceil(n * definition))
                set_a = set(merged.loc[merged[f"{first}_rank"] <= size, "gene_id"])
                set_b = set(merged.loc[merged[f"{second}_rank"] <= size, "gene_id"])
                intersection = len(set_a & set_b)
                overlap_rows.append({
                    "region": region, "method_A": first, "method_B": second, "top_set": label, "set_size": size,
                    "intersection_n": intersection, "jaccard": intersection / len(set_a | set_b),
                    "overlap_proportion_smaller_set": intersection / min(len(set_a), len(set_b)),
                })
    return pd.DataFrame(direction_rows), pd.DataFrame(rank_rows), pd.DataFrame(overlap_rows)


def tier2_audit(all_merged, tier2):
    records = []
    for region, merged in all_merged.items():
        selected = tier2[tier2.region.eq(region)][["gene_id"]].merge(merged, on="gene_id", validate="one_to_one")
        n = len(merged)
        selected["LM_limma_direction_agreement"] = np.sign(selected["Frozen LM_effect"]) == np.sign(selected["limma-voom_effect"])
        selected["LM_DESeq2_direction_agreement"] = np.sign(selected["Frozen LM_effect"]) == np.sign(selected["DESeq2_effect"])
        selected["three_method_direction_agreement"] = selected[["Frozen LM_effect", "limma-voom_effect", "DESeq2_effect"]].apply(lambda row: len(set(np.sign(row.dropna()))) == 1, axis=1)
        selected["limma_percentile_rank"] = 100 * (selected["limma-voom_rank"] - 1) / max(n - 1, 1)
        selected["DESeq2_percentile_rank"] = 100 * (selected["DESeq2_rank"] - 1) / max(n - 1, 1)
        records.append(selected)
    audit = pd.concat(records, ignore_index=True)
    audit = audit.rename(columns={"standardized_effect_size_diagnosis": "LM_standardized_effect"})
    return audit


def enrichment(all_merged, tier_keys):
    rows = []
    for region, merged in all_merged.items():
        is_tier = merged.gene_id.isin({gene for reg, gene in tier_keys if reg == region})
        for method in ["limma-voom", "DESeq2"]:
            ranks = merged[f"{method}_rank"].astype(float)
            u = stats.mannwhitneyu(ranks[is_tier], ranks[~is_tier], alternative="two-sided")
            n1, n0 = int(is_tier.sum()), int((~is_tier).sum())
            base = {
                "region": region, "method": method, "threshold": "rank_distribution",
                "tier2_n": n1, "eligible_n": len(merged), "tier2_median_percentile": float((100 * (ranks[is_tier] - 1) / (len(merged) - 1)).median()),
                "non_tier2_median_percentile": float((100 * (ranks[~is_tier] - 1) / (len(merged) - 1)).median()),
                "rank_biserial_tier2_higher": float(1 - 2 * u.statistic / (n1 * n0)), "raw_p": float(u.pvalue),
                "observed_tier2": np.nan, "expected_tier2": np.nan, "enrichment_ratio": np.nan,
            }
            rows.append(base)
            for fraction in [0.01, 0.05, 0.10]:
                cutoff = max(1, math.ceil(len(merged) * fraction))
                top = ranks <= cutoff
                observed = int((is_tier & top).sum())
                expected = n1 * cutoff / len(merged)
                table = [[observed, cutoff - observed], [n1 - observed, len(merged) - cutoff - n1 + observed]]
                fisher = stats.fisher_exact(table, alternative="greater")
                rows.append({
                    "region": region, "method": method, "threshold": f"top_{int(fraction*100)}pct", "tier2_n": n1,
                    "eligible_n": len(merged), "tier2_median_percentile": base["tier2_median_percentile"],
                    "non_tier2_median_percentile": base["non_tier2_median_percentile"],
                    "rank_biserial_tier2_higher": base["rank_biserial_tier2_higher"], "raw_p": float(fisher.pvalue),
                    "observed_tier2": observed, "expected_tier2": expected,
                    "enrichment_ratio": observed / expected if expected else np.nan,
                })
    output = pd.DataFrame(rows)
    mask = output.threshold.str.startswith("top_")
    output["BH_predeclared_enrichment_tests"] = np.nan
    output.loc[mask, "BH_predeclared_enrichment_tests"] = bh(output.loc[mask, "raw_p"])
    return output


def fdr_and_top(all_merged, deseq):
    fdr_rows = []
    for region, merged in all_merged.items():
        method_specs = [
            ("Frozen LM", "Frozen LM_p", "Frozen LM_BH"),
            ("limma-voom", "limma-voom_p", "limma-voom_BH"),
            ("DESeq2_same_universe_BH", "DESeq2_p", "DESeq2_BH"),
            ("DESeq2_native_padj", "DESeq2_p", "DESeq2_native_padj"),
        ]
        for method, pcol, qcol in method_specs:
            fdr_rows.append({
                "region": region, "method": method, "eligible_genes": len(merged), "nonmissing_raw_p": int(merged[pcol].notna().sum()),
                "raw_p_lt_0_01": int((merged[pcol] < 0.01).sum()), "BH_or_padj_lt_0_05": int((merged[qcol] < 0.05).sum()),
                "BH_or_padj_lt_0_10": int((merged[qcol] < 0.10).sum()), "minimum_raw_p": merged[pcol].min(), "minimum_adjusted": merged[qcol].min(),
            })
    fixed = []
    for region, merged in all_merged.items():
        for method in METHODS:
            for _, row in merged.sort_values([f"{method}_rank", "gene_id"], kind="mergesort").head(20).iterrows():
                fixed.append({"region": region, "method": method, "rank": int(row[f"{method}_rank"]), "gene_id": row.gene_id,
                              "gene_symbol": row.gene_symbol, "effect": row[f"{method}_effect"], "raw_p": row[f"{method}_p"],
                              "adjusted_p": row[f"{method}_BH"]})
    return pd.DataFrame(fdr_rows), pd.DataFrame(fixed)


def leader_audit(tier_audit):
    rows = []
    for region, symbol in LEADERS.items():
        match = tier_audit[(tier_audit.region.eq(region)) & (tier_audit.gene_symbol.eq(symbol))]
        if len(match) != 1:
            raise ValueError(f"Expected one leader match for {region} {symbol}; got {len(match)}")
        rows.append(match.iloc[0])
    return pd.DataFrame(rows)


def fdr_supported_audit(all_merged, tier_keys):
    rows = []
    for region, merged in all_merged.items():
        selected = merged[(merged["limma-voom_BH"] < 0.05) | (merged["DESeq2_BH"] < 0.05)].copy()
        for _, row in selected.iterrows():
            rows.append({
                "region": region, "gene_id": row.gene_id, "gene_symbol": row.gene_symbol,
                "frozen_Tier2": (region, row.gene_id) in tier_keys,
                "LM_effect": row["Frozen LM_effect"], "LM_raw_p": row["Frozen LM_p"], "LM_BH": row["Frozen LM_BH"], "LM_rank": row["Frozen LM_rank"],
                "limma_effect": row["limma-voom_effect"], "limma_raw_p": row["limma-voom_p"], "limma_BH": row["limma-voom_BH"], "limma_rank": row["limma-voom_rank"],
                "DESeq2_effect": row["DESeq2_effect"], "DESeq2_raw_p": row["DESeq2_p"], "DESeq2_same_universe_BH": row["DESeq2_BH"],
                "DESeq2_native_padj": row["DESeq2_native_padj"], "DESeq2_rank": row["DESeq2_rank"],
                "three_method_direction_agreement": len(set(np.sign([row["Frozen LM_effect"], row["limma-voom_effect"], row["DESeq2_effect"]]))) == 1,
            })
    return pd.DataFrame(rows)


def build_figures(all_merged, direction, correlations, overlap, tier_audit, fdr):
    sns.set_theme(style="whitegrid")
    for first, second, name in [
        ("Frozen LM", "limma-voom", "LM_vs_limma_signed_rank_scatter"),
        ("Frozen LM", "DESeq2", "LM_vs_DESeq2_signed_rank_scatter"),
        ("limma-voom", "DESeq2", "limma_vs_DESeq2_signed_rank_scatter"),
    ]:
        fig, axes = plt.subplots(2, 3, figsize=(15, 9))
        for ax, region in zip(axes.flat, REGIONS):
            data = all_merged[region]
            ax.scatter(data[f"{first}_signed"], data[f"{second}_signed"], s=3, alpha=0.2, rasterized=True)
            ax.axhline(0, color="grey", lw=0.7); ax.axvline(0, color="grey", lw=0.7)
            ax.set_title(region); ax.set_xlabel(first); ax.set_ylabel(second)
        fig.suptitle(f"Signed evidence: {first} versus {second}")
        fig.tight_layout(); save_figure(fig, name)

    full = direction[direction.scope.eq("full_shared")].copy()
    fig, ax = plt.subplots(figsize=(12, 6))
    full["pair"] = full.method_A + " vs " + full.method_B
    sns.barplot(data=full, x="region", y="same_sign_proportion", hue="pair", ax=ax)
    ax.set_ylim(0, 1); ax.tick_params(axis="x", rotation=30); ax.set_ylabel("Same-sign proportion")
    ax.set_title("Genome-wide diagnosis-effect direction agreement")
    fig.tight_layout(); save_figure(fig, "effect_direction_agreement_summary")

    pivot = correlations.pivot(index="region", columns=["method_A", "method_B"], values="spearman_rho").reindex(REGIONS)
    fig, ax = plt.subplots(figsize=(9, 5)); sns.heatmap(pivot, annot=True, fmt=".2f", cmap="vlag", vmin=-1, vmax=1, ax=ax)
    ax.set_title("Spearman correlation of two-sided raw-p ranks")
    fig.tight_layout(); save_figure(fig, "spearman_rank_correlation_by_region")

    for method, column, name in [("limma-voom", "limma_percentile_rank", "Tier2_percentile_rank_limma"), ("DESeq2", "DESeq2_percentile_rank", "Tier2_percentile_rank_DESeq2")]:
        fig, ax = plt.subplots(figsize=(11, 6)); sns.boxplot(data=tier_audit, x="region", y=column, order=REGIONS, ax=ax, color="#4C72B0")
        ax.set_ylim(100, 0); ax.tick_params(axis="x", rotation=30); ax.set_ylabel("Percentile rank (0 = strongest)")
        ax.set_title(f"Frozen Tier2 positions in {method} genome-wide ranking")
        fig.tight_layout(); save_figure(fig, name)

    plot_overlap = overlap[overlap.top_set.isin(["top10", "top50", "top5pct"])].copy()
    plot_overlap["comparison"] = plot_overlap.method_A + " vs " + plot_overlap.method_B + " | " + plot_overlap.top_set
    fig, ax = plt.subplots(figsize=(13, 7)); sns.barplot(data=plot_overlap, x="region", y="jaccard", hue="comparison", errorbar=None, ax=ax)
    ax.tick_params(axis="x", rotation=30); ax.set_title("Top-set overlap across statistical methods")
    fig.tight_layout(); save_figure(fig, "topset_Jaccard_summary")

    plot_fdr = fdr[fdr.method.isin(["Frozen LM", "limma-voom", "DESeq2_same_universe_BH"])].copy()
    melted = plot_fdr.melt(id_vars=["region", "method"], value_vars=["raw_p_lt_0_01", "BH_or_padj_lt_0_05"], var_name="threshold", value_name="genes")
    melted["series"] = melted.method + " | " + melted.threshold
    fig, ax = plt.subplots(figsize=(13, 7)); sns.barplot(data=melted, x="region", y="genes", hue="series", errorbar=None, ax=ax)
    ax.set_yscale("symlog", linthresh=1); ax.tick_params(axis="x", rotation=30); ax.set_title("Diagnosis findings by method (raw p<0.01 and BH<0.05 rows combined)")
    fig.tight_layout(); save_figure(fig, "finding_counts_by_method_region")


def main():
    ensure_dirs()
    lm, limma, deseq, tier2 = load_data()
    all_merged = {region: method_frame(lm, limma, deseq, region) for region in REGIONS}
    tier_keys = set(zip(tier2.region, tier2.gene_id))
    direction, correlations, overlap = direction_and_rank(all_merged, tier_keys)
    tier_audit = tier2_audit(all_merged, tier2)
    enrich = enrichment(all_merged, tier_keys)
    fdr, top = fdr_and_top(all_merged, deseq)
    leaders = leader_audit(tier_audit)
    fdr_genes = fdr_supported_audit(all_merged, tier_keys)

    direction.to_csv(OUT / "method_direction_agreement.csv", index=False)
    correlations.to_csv(OUT / "method_rank_correlations.csv", index=False)
    overlap.to_csv(OUT / "method_topset_overlap.csv", index=False)
    tier_audit.to_csv(OUT / "tier2_statistical_method_audit.csv", index=False)
    enrich.to_csv(OUT / "tier2_rank_enrichment.csv", index=False)
    fdr.to_csv(OUT / "fdr_summary_by_region_method.csv", index=False)
    top.to_csv(OUT / "top_genes_by_region_method.csv", index=False)
    leaders.to_csv(OUT / "tier2_rra_leading_gene_statistical_audit.csv", index=False)
    fdr_genes.to_csv(OUT / "fdr_supported_genes_method_comparison.csv", index=False)
    pd.concat(all_merged.values(), ignore_index=True).to_csv(
        SUMMARY / "three_method_merged_all_genes.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    build_figures(all_merged, direction, correlations, overlap, tier_audit, fdr)

    versions = pd.read_csv(OUT / "software_versions.csv").set_index("component")["value"].to_dict()
    manifest = {
        "status": "complete",
        "regions": REGIONS,
        "n_regions": 6,
        "n_frozen_tier2_pairs": 671,
        "candidate_tier_changed": False,
        "same_samples_all_methods": True,
        "same_eligible_genes_all_methods": True,
        "limma_input": "raw integer counts; DGEList; TMM; voom; lmFit; eBayes",
        "deseq2_input": "raw integer counts; standard DESeq2 median-of-ratios and DESeq defaults",
        "design": "Diagnosis + Age + Sex + RIN; explicit MDD minus CTRL",
        "primary_rank": "two-sided raw p ascending; GeneID deterministic tie-break; missing p last",
        "software": versions,
        "tier2_sha256": sha256(TIER2),
    }
    (OUT / "statistical_benchmark_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print("Built statistical-method concordance summaries and nine figure pairs")


if __name__ == "__main__":
    main()
