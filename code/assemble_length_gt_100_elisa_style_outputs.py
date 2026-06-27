"""Assemble the Length > 100 Elisa-style output branch from completed pipeline outputs."""

from pathlib import Path
import math
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "outputs" / "length_gt_100_elisa_style"
STEP1 = OUT / "step1_from_raw" / "Length_gt_100"
STEP3 = REPO / "outputs" / "step3_length_sensitivity"
ANNOT = Path("D:/genes/genedata/Human.GRCh38.p13.annot.tsv")
RAW = Path("D:/genes/genedata/GSE102556_raw_counts_GRCh38.p13_NCBI.tsv")

REGIONS = {
    "Frontal_Cortex": {
        "stub": "frontal_cortex",
        "method": "lmm_getmm_log2",
        "result": STEP3 / "Length_gt_100" / "lmm_getmm_log2" / "lmm_getmm_log2_Frontal_Cortex_results.csv",
        "previous": STEP3 / "Length_gt_20" / "lmm_getmm_log2" / "lmm_getmm_log2_Frontal_Cortex_results.csv",
    },
    "Insula": {
        "stub": "insula",
        "method": "lm_getmm_log2_no_repeated_donors",
        "result": STEP3 / "Length_gt_100" / "lm_getmm_log2" / "lm_getmm_log2_Insula_results.csv",
        "previous": STEP3 / "Length_gt_20" / "lm_getmm_log2" / "lm_getmm_log2_Insula_results.csv",
    },
    "Ventral_Striatum": {
        "stub": "ventral_striatum",
        "method": "lm_getmm_log2_no_repeated_donors",
        "result": STEP3 / "Length_gt_100" / "lm_getmm_log2" / "lm_getmm_log2_Ventral_Striatum_results.csv",
        "previous": STEP3 / "Length_gt_20" / "lm_getmm_log2" / "lm_getmm_log2_Ventral_Striatum_results.csv",
    },
    "Hippocampal_Formation": {
        "stub": "hippocampal_formation",
        "method": "lm_getmm_log2_no_repeated_donors",
        "result": STEP3 / "Length_gt_100" / "lm_getmm_log2" / "lm_getmm_log2_Hippocampal_Formation_results.csv",
        "previous": STEP3 / "Length_gt_20" / "lm_getmm_log2" / "lm_getmm_log2_Hippocampal_Formation_results.csv",
    },
}


def mkdirs():
    for rel in [
        "normalized_expression_length_gt_100_by_region",
        "filtered_raw_counts_length_gt_100_by_region",
        "normalization_factors_length_gt_100_by_region",
        "qc_pca_length_gt_100",
        "volcano_plots_nominal_length_gt_100",
        "source_completed_pipeline_files",
    ]:
        (OUT / rel).mkdir(parents=True, exist_ok=True)


def p_adjust_bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    tmp = np.empty(n)
    tmp[order] = np.minimum(adj, 1)
    out[valid] = tmp
    return out


def load_annotation():
    annot = pd.read_csv(ANNOT, sep="\t", dtype={"GeneID": str}, low_memory=False)
    annot["Length"] = pd.to_numeric(annot["Length"], errors="coerce")
    return annot.drop_duplicates("GeneID").set_index("GeneID")


def write_filter_summary(annot):
    raw_first_col = pd.read_csv(RAW, sep="\t", usecols=[0], dtype=str)
    raw_genes = raw_first_col.iloc[:, 0].astype(str).drop_duplicates()
    common = raw_genes[raw_genes.isin(annot.index)]
    lengths = annot.loc[common, "Length"]
    after = lengths[lengths > 100]
    step1_summary = pd.read_csv(STEP1 / "07_qc" / "Step1_QC_summary_by_analysis_group_keep_all_samples.csv")

    def stats_for(x, prefix):
        return {
            f"{prefix}_min": x.min(),
            f"{prefix}_q1": x.quantile(0.25),
            f"{prefix}_median": x.median(),
            f"{prefix}_mean": x.mean(),
            f"{prefix}_q3": x.quantile(0.75),
            f"{prefix}_max": x.max(),
        }

    rows = [{
        "scope": "all_raw_count_genes_with_annotation",
        "genes_before_filtering": len(common),
        "genes_after_Length_gt_100_filtering": int((lengths > 100).sum()),
        "genes_removed_by_length_filter": int((lengths <= 100).sum()),
        "genes_after_group_expression_filter": "",
        **stats_for(lengths.dropna(), "length_before"),
        **stats_for(after.dropna(), "length_after"),
    }]
    for row in step1_summary.itertuples():
        rows.append({
            "scope": row.brain_area,
            "genes_before_filtering": len(common),
            "genes_after_Length_gt_100_filtering": int(row.genes_after_length_filter),
            "genes_removed_by_length_filter": len(common) - int(row.genes_after_length_filter),
            "genes_after_group_expression_filter": int(row.genes_after_group_expression_filter),
            **stats_for(lengths.dropna(), "length_before"),
            **stats_for(after.dropna(), "length_after"),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "gene_filter_length_gt_100_summary.csv", index=False)
    return out


def write_qc(region, stub):
    expr_path = STEP1 / "02_normalized_expression" / "keep_all_samples" / f"normalized_GeTMM_log2_{stub}.csv"
    counts_path = STEP1 / "03_filtered_raw_counts" / "keep_all_samples" / f"filtered_raw_counts_{stub}.csv"
    factors_path = STEP1 / "04_normalization_factors" / "keep_all_samples" / f"GeTMM_normalization_factors_{stub}.csv"
    meta_path = STEP1 / "01_metadata" / "keep_all_samples" / f"metadata_{stub}.csv"
    expr = pd.read_csv(expr_path, index_col=0)
    meta = pd.read_csv(meta_path, dtype={"SampleID": str, "BrainID": str})
    shutil.copy2(expr_path, OUT / "normalized_expression_length_gt_100_by_region" / expr_path.name)
    shutil.copy2(counts_path, OUT / "filtered_raw_counts_length_gt_100_by_region" / counts_path.name)
    shutil.copy2(factors_path, OUT / "normalization_factors_length_gt_100_by_region" / factors_path.name)

    sample_counts = meta.groupby(["Diagnosis"], observed=False).size().reset_index(name="n_samples")
    sample_counts["brain_region"] = region
    sample_counts.to_csv(OUT / "qc_pca_length_gt_100" / f"sample_counts_{stub}.csv", index=False)

    x = expr.T
    pcs = PCA(n_components=5).fit_transform(StandardScaler().fit_transform(x))
    pca = PCA(n_components=5).fit(StandardScaler().fit_transform(x))
    scores = pd.DataFrame(pcs, columns=[f"PC{i}" for i in range(1, 6)])
    scores.insert(0, "SampleID", x.index)
    scores = scores.merge(meta, on="SampleID", how="left")
    scores.to_csv(OUT / "qc_pca_length_gt_100" / f"pca_scores_{stub}.csv", index=False)
    pd.DataFrame({"PC": [f"PC{i}" for i in range(1, 6)], "variance_explained": pca.explained_variance_ratio_}).to_csv(
        OUT / "qc_pca_length_gt_100" / f"pca_variance_{stub}.csv", index=False
    )

    assoc = []
    for pc in ["PC1", "PC2", "PC3"]:
        for cov in ["Diagnosis", "Sex", "OriginalRegion"]:
            groups = [g[pc].dropna() for _, g in scores.groupby(cov, observed=False) if g[pc].notna().sum() > 1]
            assoc.append({"brain_region": region, "PC": pc, "covariate": cov, "p_value": stats.kruskal(*groups).pvalue if len(groups) > 1 else np.nan})
        for cov in ["Age", "RIN"]:
            valid = scores[[pc, cov]].dropna()
            assoc.append({"brain_region": region, "PC": pc, "covariate": cov, "p_value": stats.pearsonr(valid[pc], valid[cov]).pvalue if len(valid) > 2 else np.nan})
    pd.DataFrame(assoc).to_csv(OUT / "qc_pca_length_gt_100" / f"pc_covariate_associations_{stub}.csv", index=False)

    plt.figure(figsize=(7, 6))
    for diagnosis, group in scores.groupby("Diagnosis", observed=False):
        plt.scatter(group["PC1"], group["PC2"], label=diagnosis, alpha=0.8)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    plt.title(f"PCA: {region}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "qc_pca_length_gt_100" / f"pca_{stub}.png", dpi=160)
    plt.close()
    return expr, meta


def convert_results(region, cfg, annot, expr, meta):
    raw = pd.read_csv(cfg["result"], dtype={"GeneID": str})
    raw["gene_id"] = raw["GeneID"].astype(str)
    raw = raw.merge(annot[["Symbol", "Length"]], left_on="gene_id", right_index=True, how="left")
    mdd = meta.loc[meta["Diagnosis"] == "MDD", "SampleID"]
    ctrl = meta.loc[meta["Diagnosis"] == "CTRL", "SampleID"]
    result = pd.DataFrame({
        "gene_id": raw["gene_id"],
        "gene_symbol": raw["Symbol"],
        "brain_region": region,
        "gene_length": raw["Length"],
        "beta_diagnosis": raw["diagnosis_effect"],
        "p_value_diagnosis": raw["diagnosis_p_value"],
        "standard_error_diagnosis": raw.get("standard_error", np.nan),
        "mean_expression_MDD": expr.reindex(raw["gene_id"])[mdd].mean(axis=1).values,
        "mean_expression_Control": expr.reindex(raw["gene_id"])[ctrl].mean(axis=1).values,
        "model_source": cfg["method"],
    })
    result["direction"] = np.where(result["beta_diagnosis"] > 0, "up_in_MDD", np.where(result["beta_diagnosis"] < 0, "down_in_MDD", ""))
    result["signed_logp_score"] = np.sign(result["beta_diagnosis"]) * -np.log10(result["p_value_diagnosis"].clip(lower=np.finfo(float).tiny))
    result["fdr_BH"] = p_adjust_bh(result["p_value_diagnosis"])
    return result


def write_volcano(region, df):
    plot = df.dropna(subset=["beta_diagnosis", "p_value_diagnosis"])
    if plot.empty:
        return
    plt.figure(figsize=(7, 6))
    color = np.where(plot["p_value_diagnosis"] < 0.05, "#b23a48", "#667085")
    plt.scatter(plot["beta_diagnosis"], -np.log10(plot["p_value_diagnosis"].clip(lower=np.finfo(float).tiny)), s=8, c=color, alpha=0.65)
    plt.axhline(-math.log10(0.05), color="#444", linestyle="--", linewidth=1)
    plt.axvline(0, color="#444", linewidth=1)
    plt.xlabel("Diagnosis beta (MDD vs Control)")
    plt.ylabel("-log10(raw p-value)")
    plt.title(f"Nominal volcano: {region}")
    plt.tight_layout()
    plt.savefig(OUT / "volcano_plots_nominal_length_gt_100" / f"volcano_nominal_{region}.png", dpi=180)
    plt.close()


def summarize(all_results):
    rows = []
    for region, g in all_results.groupby("brain_region", sort=False):
        rows.append({
            "brain_region": region,
            "genes_tested": len(g),
            "nominal_p_lt_0.05": int((g["p_value_diagnosis"] < 0.05).sum()),
            "nominal_p_lt_0.01": int((g["p_value_diagnosis"] < 0.01).sum()),
            "nominal_p_lt_0.001": int((g["p_value_diagnosis"] < 0.001).sum()),
            "upregulated_in_MDD_p_lt_0.05": int(((g["p_value_diagnosis"] < 0.05) & (g["beta_diagnosis"] > 0)).sum()),
            "downregulated_in_MDD_p_lt_0.05": int(((g["p_value_diagnosis"] < 0.05) & (g["beta_diagnosis"] < 0)).sum()),
            "FDR_BH_lt_0.05": int((g["fdr_BH"] < 0.05).sum()),
            "minimum_raw_p_value": g["p_value_diagnosis"].min(),
            "minimum_fdr_BH": g["fdr_BH"].min(),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "summary_by_region_length_gt_100.csv", index=False)
    return out


def previous_summary():
    rows = []
    for region, cfg in REGIONS.items():
        old = pd.read_csv(cfg["previous"])
        rows.append({
            "brain_region": region,
            "previous_genes_tested": len(old),
            "previous_nominal_p_lt_0.05": int((old["diagnosis_p_value"] < 0.05).sum()),
            "previous_nominal_p_lt_0.01": int((old["diagnosis_p_value"] < 0.01).sum()),
            "previous_nominal_p_lt_0.001": int((old["diagnosis_p_value"] < 0.001).sum()),
            "previous_FDR_BH_lt_0.05": int((old["diagnosis_FDR_BH"] < 0.05).sum()),
            "previous_source_file": str(cfg["previous"].relative_to(REPO)),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "previous_pipeline_comparison_reference.csv", index=False)
    return out


def write_report(filter_summary, summary, previous, top_p):
    merged = summary.merge(previous, on="brain_region", how="left")
    comparison = []
    for _, r in merged.iterrows():
        comparison.append(
            f"- {r['brain_region']}: {r['nominal_p_lt_0.05']} nominal genes at p < 0.05 now versus "
            f"{r['previous_nominal_p_lt_0.05']} in the previous Length > 20 comparison."
        )
    top_lines = []
    for region, group in top_p.groupby("brain_region", sort=False):
        genes = ", ".join(
            f"{(row.gene_symbol if isinstance(row.gene_symbol, str) and row.gene_symbol else row.gene_id)} "
            f"(p={row.p_value_diagnosis:.3g}, beta={row.beta_diagnosis:.3g})"
            for row in group.head(5).itertuples()
        )
        top_lines.append(f"- {region}: {genes}")
    text = f"""# Length > 100 Elisa-Style Report

## What Changed

- Previous gene-length threshold used for comparison: Length > 20
- New gene-length threshold: Length > 100
- Normalization: filtered raw counts -> GeTMM -> log2(GeTMM + 1)
- Frontal Cortex model source: LMM on log2(GeTMM + 1), preserving the previous formula with donor random intercept
- Other brain regions: previous normalized-expression model output is used because each region has zero repeated donors, so a donor random intercept is not estimable
- Main exploratory criteria: diagnosis beta/effect size and raw p-value

## Gene Filter Summary

{filter_summary.to_markdown(index=False)}

## Nominal Results By Brain Region

{summary.to_markdown(index=False)}

## Comparison With Previous Pipeline

{chr(10).join(comparison)}

## Top Genes By Raw P-Value

{chr(10).join(top_lines)}

## FDR Reference

BH-FDR is saved as `fdr_BH` for reference only. It is not used as the main filter in this exploratory Elisa-style rerun.

## Interpretation Warning

These raw-p-value results are nominal and exploratory. Do not describe these genes as final significant DE genes. The purpose is to check whether an MDD-associated signal is visible before strict multiple-testing correction.
"""
    (OUT / "Length_gt_100_Elisa_style_report.md").write_text(text, encoding="utf-8")


def main():
    mkdirs()
    annot = load_annotation()
    filter_summary = write_filter_summary(annot)
    results = []
    for region, cfg in REGIONS.items():
        expr, meta = write_qc(region, cfg["stub"])
        converted = convert_results(region, cfg, annot, expr, meta)
        write_volcano(region, converted)
        results.append(converted)
        shutil.copy2(cfg["result"], OUT / "source_completed_pipeline_files" / cfg["result"].name)
    all_results = pd.concat(results, ignore_index=True)
    all_results.to_csv(OUT / "LMM_all_genes_length_gt_100.csv", index=False)
    all_results[all_results["p_value_diagnosis"] < 0.05].to_csv(OUT / "LMM_nominal_p_lt_005_length_gt_100.csv", index=False)
    all_results[all_results["p_value_diagnosis"] < 0.01].to_csv(OUT / "LMM_nominal_p_lt_001_length_gt_100.csv", index=False)
    all_results[all_results["p_value_diagnosis"] < 0.001].to_csv(OUT / "LMM_nominal_p_lt_0001_length_gt_100.csv", index=False)

    top_p = all_results.sort_values(["brain_region", "p_value_diagnosis"], na_position="last").groupby("brain_region", group_keys=False).head(100)
    top_abs = all_results.assign(abs_beta=all_results["beta_diagnosis"].abs()).sort_values(["brain_region", "abs_beta"], ascending=[True, False], na_position="last").groupby("brain_region", group_keys=False).head(100).drop(columns="abs_beta")
    top_score = all_results.assign(abs_score=all_results["signed_logp_score"].abs()).sort_values(["brain_region", "abs_score"], ascending=[True, False], na_position="last").groupby("brain_region", group_keys=False).head(100).drop(columns="abs_score")
    top_p.to_csv(OUT / "top100_by_pvalue_length_gt_100.csv", index=False)
    top_abs.to_csv(OUT / "top100_by_abs_beta_length_gt_100.csv", index=False)
    top_score.to_csv(OUT / "top100_by_signed_logp_score_length_gt_100.csv", index=False)
    summary = summarize(all_results)
    prev = previous_summary()
    write_report(filter_summary, summary, prev, top_p)
    shutil.copy2(Path(__file__), OUT / "pipeline_script_used_assemble_length_gt_100_elisa_style_outputs.py")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
