from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from benchmark_common import OUT, REPO, REGIONS, sha256

REPORT = REPO / "docs" / "UPSTREAM_STATISTICAL_METHOD_BENCHMARK_REPORT.md"


def fmt(value, digits=3):
    if pd.isna(value):
        return "NA"
    if isinstance(value, (bool,)): return "yes" if value else "no"
    if abs(float(value)) < 0.001 and float(value) != 0: return f"{float(value):.2e}"
    return f"{float(value):.{digits}f}"


def markdown_table(frame: pd.DataFrame) -> str:
    """Render a small DataFrame without pandas' optional tabulate dependency."""
    columns = [str(column) for column in frame.columns]
    rows = [[str(value) for value in row] for row in frame.itertuples(index=False, name=None)]
    widths = [
        max(len(columns[index]), *(len(row[index]) for row in rows))
        for index in range(len(columns))
    ]

    def render(row):
        return "| " + " | ".join(value.ljust(widths[index]) for index, value in enumerate(row)) + " |"

    return "\n".join([
        render(columns),
        "| " + " | ".join("-" * width for width in widths) + " |",
        *(render(row) for row in rows),
    ])


def main():
    sample = pd.read_csv(OUT / "sample_accounting.csv")
    genes = pd.read_csv(OUT / "gene_universe_accounting.csv")
    versions = pd.read_csv(OUT / "software_versions.csv")
    direction = pd.read_csv(OUT / "method_direction_agreement.csv")
    correlations = pd.read_csv(OUT / "method_rank_correlations.csv")
    overlap = pd.read_csv(OUT / "method_topset_overlap.csv")
    enrichment = pd.read_csv(OUT / "tier2_rank_enrichment.csv")
    fdr = pd.read_csv(OUT / "fdr_summary_by_region_method.csv")
    leaders = pd.read_csv(OUT / "tier2_rra_leading_gene_statistical_audit.csv", dtype={"gene_id": str})
    fdr_genes = pd.read_csv(OUT / "fdr_supported_genes_method_comparison.csv", dtype={"gene_id": str})

    account = sample[["region", "samples_in_frozen_LM"]].merge(genes[["region", "frozen_LM_eligible_genes"]], on="region")
    account.columns = ["Region", "Matched samples", "Matched eligible genes"]

    corr = correlations[correlations.method_A.eq("Frozen LM")].copy()
    corr["comparison"] = corr.method_A + " vs " + corr.method_B
    corr_table = corr.pivot(index="region", columns="comparison", values="spearman_rho").reindex(REGIONS).reset_index()
    corr_table.columns.name = None
    for column in corr_table.columns[1:]: corr_table[column] = corr_table[column].map(lambda value: fmt(value))

    full_dir = direction[(direction.scope.eq("full_shared")) & (direction.method_A.eq("Frozen LM"))].copy()
    full_dir["comparison"] = full_dir.method_A + " vs " + full_dir.method_B
    dir_table = full_dir.pivot(index="region", columns="comparison", values="same_sign_proportion").reindex(REGIONS).reset_index()
    dir_table.columns.name = None
    for column in dir_table.columns[1:]: dir_table[column] = dir_table[column].map(lambda value: f"{100*value:.1f}%")

    enrich_table = enrichment[(enrichment.threshold.eq("top_1pct"))].copy()
    enrich_table = enrich_table[["region", "method", "observed_tier2", "expected_tier2", "enrichment_ratio", "tier2_median_percentile"]]
    enrich_table.columns = ["Region", "Method", "Tier2 in top 1%", "Expected", "Enrichment ratio", "Tier2 median percentile"]
    enrich_table["Expected"] = enrich_table["Expected"].map(lambda value: fmt(value, 2))
    enrich_table["Enrichment ratio"] = enrich_table["Enrichment ratio"].map(lambda value: fmt(value, 1))
    enrich_table["Tier2 median percentile"] = enrich_table["Tier2 median percentile"].map(lambda value: fmt(value, 2))

    fdr_table = fdr[["region", "method", "raw_p_lt_0_01", "BH_or_padj_lt_0_05", "BH_or_padj_lt_0_10", "minimum_adjusted"]].copy()
    fdr_table.columns = ["Region", "Method", "Raw p<0.01", "Adjusted<0.05", "Adjusted<0.10", "Minimum adjusted"]
    fdr_table["Minimum adjusted"] = fdr_table["Minimum adjusted"].map(lambda value: fmt(value))

    leader_table = leaders[["region", "gene_symbol", "Frozen LM_effect", "Frozen LM_p", "Frozen LM_rank", "limma-voom_effect", "limma-voom_p", "limma-voom_rank", "DESeq2_effect", "DESeq2_p", "DESeq2_rank", "three_method_direction_agreement"]].copy()
    leader_table.columns = ["Region", "Tier2-RRA #1 gene", "LM effect", "LM p", "LM rank", "limma effect", "limma p", "limma rank", "DESeq2 effect", "DESeq2 p", "DESeq2 rank", "Three-way direction"]
    for column in ["LM effect", "LM p", "limma effect", "limma p", "DESeq2 effect", "DESeq2 p"]:
        leader_table[column] = leader_table[column].map(fmt)
    leader_table["Three-way direction"] = leader_table["Three-way direction"].map({True: "yes", False: "no"})

    fdr_gene_table = fdr_genes[["region", "gene_symbol", "frozen_Tier2", "LM_BH", "limma_BH", "DESeq2_same_universe_BH", "DESeq2_native_padj", "three_method_direction_agreement"]].copy()
    fdr_gene_table.columns = ["Region", "Gene", "Frozen Tier2", "LM BH", "limma BH", "DESeq2 same-universe BH", "DESeq2 native padj", "Three-way direction"]
    for column in ["LM BH", "limma BH", "DESeq2 same-universe BH", "DESeq2 native padj"]:
        fdr_gene_table[column] = fdr_gene_table[column].map(fmt)

    text = f"""# Upstream statistical-method benchmark report

## 1. Why Tier2 needed benchmarking

Tier2 was defined using a gene-wise LM on `log2(GeTMM + 1)`. This analysis asks whether its directions and genome-wide positions remain similar when the same samples and eligible genes are analysed using limma-voom and DESeq2. This is **within-study statistical-method concordance**, not independent replication.

## 2. Frozen current LM and Tier2 definition

**SOURCE-DERIVED FACT:** The frozen model was `expression ~ Diagnosis + Age + Sex + RIN`, with Control as reference and positive diagnosis effects meaning higher expression in MDD. Tier2 remains exactly 671 region–gene pairs selected by raw diagnosis p < 0.01 and absolute standardized diagnosis effect > 0.8. No membership was changed.

The LM used `log2(GeTMM + 1)`. Its standardized effect equals the diagnosis beta divided by the gene-wise expression SD. BH correction was regional across the full eligible LM universe.

## 3. Matched samples and genes

**SOURCE-DERIVED FACT:** Matching was exact. No benchmark sample or gene was added or removed relative to its regional frozen LM universe.

{markdown_table(account)}

The exact filtered raw-count source is `D:/genes/GeneDepresion/outputs/step1_length_sensitivity/Length_gt_100/03_filtered_raw_counts/keep_all_samples/`. A same-sized but non-identical later Step1 export was rejected before fitting.

## 4. limma-voom method

Raw integer counts were passed through `edgeR::DGEList`, TMM normalization, standard `limma::voom`, `lmFit`, and `eBayes`, using `~ Diagnosis + Age + Sex + RIN`. No quality weights or result-driven tuning were used. The diagnosis coefficient is MDD minus Control.

## 5. DESeq2 method

Raw integer counts were fitted with standard DESeq2 using `~ Age + Sex + RIN + Diagnosis` and the explicit MDD-versus-Control contrast. Native DESeq2 `padj` and a separate BH correction across the same eligible regional universe are both retained.

{markdown_table(versions)}

## 6. Genome-wide direction agreement

**SOURCE-DERIVED FACT:** Genome-wide direction agreement between the LM and benchmark methods was 92.9%–96.9%. All 671 frozen Tier2 pairs had the same direction in all three methods.

{markdown_table(dir_table)}

**INTERPRETATION:** Direction is highly robust to statistical-model choice in this cohort.

## 7. Genome-wide rank agreement

The primary rank is two-sided raw p ascending, with GeneID as deterministic tie-breaker. Exact coefficient magnitudes were not compared as if they were interchangeable.

{markdown_table(corr_table)}

**SOURCE-DERIVED FACT:** LM-versus-limma and LM-versus-DESeq2 Spearman correlations were 0.901–0.966 across regions. Kendall correlations are retained in `method_rank_correlations.csv`.

**INTERPRETATION:** Genome-wide ordering is strongly concordant, although the top few positions are not identical.

## 8. Top-set overlap

Top 10, 25, 50, 1%, and 5% pairwise intersections and Jaccard indices are retained in `method_topset_overlap.csv`. They show partial rather than identical leading-gene lists, which is expected from different count and variance models.

## 9. Direct audit of Tier2

The complete 671-row table `tier2_statistical_method_audit.csv` preserves LM, limma, and DESeq2 effects, p-values, adjusted values, ranks, percentiles, and direction-agreement flags separately. No consensus score or new tier was constructed.

## 10. Tier2 enrichment in benchmark rankings

{markdown_table(enrich_table)}

**SOURCE-DERIVED FACT:** Regional Tier2 median positions were within the strongest 0.18%–0.60% of both benchmark rankings. Top-1% enrichment ratios ranged from 67.6-fold to 98.5-fold, and every predeclared top-set enrichment remained below BH 0.05.

**INTERPRETATION:** Tier2 is strongly concentrated near the top of both benchmarks. Because Tier2 and all benchmark methods use the same cohort, this is expected to some degree and must not be called external validation.

## 11. FDR results

{markdown_table(fdr_table)}

**SOURCE-DERIVED FACT:** The frozen LM produced no BH<0.05 gene. limma-voom produced one in OFC BA11 (`QPRT`). The explicit same-universe DESeq2 BH produced eight in Ventral Striatum and one in OFC BA11; native DESeq2 `padj` produced eleven and one, respectively. No other region produced BH/padj<0.05 genes.

The genes reaching BH<0.05 in limma or the explicit DESeq2 comparison were:

{markdown_table(fdr_gene_table)}

All have the same effect direction across the three methods. These are method-specific same-cohort findings. They do not replace Tier2, and they should be independently checked before any manuscript claim.

## 12. Region-specific interpretation

- **Ventral Striatum:** very high directional/rank concordance; DESeq2 alone produces eight same-universe BH<0.05 genes. This is supportive but method-specific.
- **Insula:** high concordance and strong Tier2 concentration; no method produces an FDR-supported gene.
- **Hippocampal Formation:** the lowest, but still strong, LM rank concordance (about 0.90); no FDR-supported gene.
- **CG25:** very high LM/limma and LM/DESeq2 concordance despite the smaller sample size; no BH<0.05 gene.
- **dlPFC BA8/9:** high direction/rank concordance and strong Tier2 ranking concentration; no BH<0.05 gene.
- **OFC BA11:** high concordance; `QPRT` is rank 1 in all three methods and reaches BH<0.05 in limma and DESeq2, but not the frozen LM.

## 13. Current Tier2-RRA top-gene audit

{markdown_table(leader_table)}

**SOURCE-DERIVED FACT:** All six current Tier2-RRA rank-one genes retain the same direction across methods. Their benchmark ranks range from 3 to 49 under limma and 5 to 349 under DESeq2; CG25 `LOC107985795` and OFC `LOC105373386` show the largest DESeq2 rank weakening.

## 14. Limitations

1. All methods use the same cohort; this is not replication.
2. Tier2 was selected using the frozen LM, so enrichment in correlated benchmark rankings is not independent.
3. Regional sample sizes are small, particularly CG25.
4. LM, limma-voom, and DESeq2 use different normalization and distributional assumptions.
5. Effect magnitudes are on related but non-identical scales.
6. Method agreement does not establish a biological mechanism or biomarker.
7. No external cohort was included.
8. Cross-method agreement cannot convert an absent LM FDR result into LM significance.

## 15. Recommendation

**INTERPRETATION:** The result is closest to prospective **Scenario A**, with the important same-cohort limitation: effect directions, genome-wide rankings, and Tier2 concentration are strongly concordant across methods. Tier2 is robust enough to retain unchanged as an exploratory candidate screen.

Do not create a consensus tier or modify Tier2. Report the limma/DESeq2 FDR findings as method-specific sensitivity results, not replacements for the frozen analysis. Continue to describe Tier2 and downstream ML as exploratory prioritisation requiring external validation.
"""
    REPORT.write_text(text, encoding="utf-8")

    manifest_path = OUT / "statistical_benchmark_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.update({
        "report": str(REPORT),
        "overall_interpretation": "Scenario A with same-cohort limitation",
        "recommendation": "retain frozen Tier2 unchanged",
        "frozen_LM_BH005_total": int(fdr[(fdr.method.eq("Frozen LM"))].BH_or_padj_lt_0_05.sum()),
        "limma_BH005_total": int(fdr[(fdr.method.eq("limma-voom"))].BH_or_padj_lt_0_05.sum()),
        "DESeq2_same_universe_BH005_total": int(fdr[(fdr.method.eq("DESeq2_same_universe_BH"))].BH_or_padj_lt_0_05.sum()),
    })
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    checksum_path = OUT / "checksums.sha256"
    files = sorted(path for path in OUT.rglob("*") if path.is_file() and path != checksum_path)
    lines = [f"{sha256(path)}  {path.relative_to(OUT).as_posix()}" for path in files]
    checksum_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote supervisor report and checksum ledger for {len(files)} output files")


if __name__ == "__main__":
    main()
