"""Test MDD candidate-gene enrichment in published AHBA canonical modules."""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path
from textwrap import wrap

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CLEAN_DIR = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
SOURCE_LENGTH100_DIR = PROJECT_ROOT / "outputs" / "length_gt_100_elisa_style"
CANONICAL_ML_RANKING = CLEAN_DIR / "10_ML_prioritized_Tier2_corrected" / "ML_prioritized_Tier2_genes_all_ranked.csv"
TIER2_DIR = CLEAN_DIR / "07_Tier2_candidate_gene_sets"
AHBA_DIR = CLEAN_DIR / "11_AHBA_modules_enrichment"
PUBLISHED_DIR = AHBA_DIR / "published_modules"
OUT_DIR = PUBLISHED_DIR / "enrichment"
FIG_DIR = PUBLISHED_DIR / "figures"
MODULE_TABLE = PUBLISHED_DIR / "ahba_published_gene_to_module_assignments.csv"
MIN_DIRECTION_SET_GENES = 5


HIGH_CLASS = "High-confidence ML-prioritized Tier 2"
MODERATE_CLASS = "Moderate-confidence ML-prioritized Tier 2"


def norm_symbol(x: object) -> str:
    return "" if pd.isna(x) else str(x).strip().upper()


def bh_fdr(pvals: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvals, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    n = len(pv)
    if n == 0:
        return pd.Series(out, index=pvals.index)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    restored = np.empty(n)
    restored[order] = np.minimum(adj, 1.0)
    out[valid] = restored
    return pd.Series(out, index=pvals.index)


def unique_nonempty_symbols(values: pd.Series) -> set[str]:
    return {x for x in values.map(norm_symbol) if x}


def direction_column(df: pd.DataFrame) -> str | None:
    for col in ("statistical_direction", "direction"):
        if col in df.columns:
            return col
    return None


def make_candidate_set(candidate_set: str, scope: str, frame: pd.DataFrame, reason: str) -> dict:
    return {
        "candidate_set": candidate_set,
        "scope": scope,
        "original_gene_region_pair_count": int(len(frame)),
        "genes": unique_nonempty_symbols(frame["gene_symbol"]),
        "reason_included": reason,
    }


def load_candidates() -> list[dict]:
    ranked_path = CANONICAL_ML_RANKING
    all_tier2_path = TIER2_DIR / "Tier2_candidate_genes_all_regions.csv"
    if not ranked_path.exists():
        raise FileNotFoundError(f"Missing ML-prioritized Tier 2 table: {ranked_path}")
    if not all_tier2_path.exists():
        raise FileNotFoundError(f"Missing all-Tier2 candidate table: {all_tier2_path}")

    ranked = pd.read_csv(ranked_path, dtype={"gene_id": str})
    all_tier2 = pd.read_csv(all_tier2_path, dtype={"gene_id": str})
    if "gene_symbol" not in ranked.columns or "gene_symbol" not in all_tier2.columns:
        raise ValueError("Candidate tables must include gene_symbol.")

    high = ranked[ranked["ML_priority_class"].eq(HIGH_CLASS)]
    high_moderate = ranked[ranked["ML_priority_class"].isin([HIGH_CLASS, MODERATE_CLASS])]

    sets: list[dict] = [
        make_candidate_set("high_confidence", "all_regions", high, "all high-confidence ML-prioritized Tier 2 gene-region rows"),
        make_candidate_set("high_plus_moderate", "all_regions", high_moderate, "all high- or moderate-confidence ML-prioritized Tier 2 gene-region rows"),
        make_candidate_set("all_Tier2", "all_regions", all_tier2, "all Tier 2 candidate gene-region rows"),
    ]

    for region, sub in high_moderate.groupby("region_or_subregion"):
        sets.append(
            make_candidate_set(
                "high_plus_moderate",
                f"region:{region}",
                sub,
                f"high- or moderate-confidence ML-prioritized Tier 2 rows in {region}",
            )
        )

    if "region_or_subregion" in all_tier2.columns:
        for region, sub in all_tier2.groupby("region_or_subregion"):
            sets.append(make_candidate_set("all_Tier2", f"region:{region}", sub, f"all Tier 2 rows in {region}"))

    for source_name, frame, set_name in [
        ("high_plus_moderate", high_moderate, "high_plus_moderate"),
        ("all_Tier2", all_tier2, "all_Tier2"),
    ]:
        dcol = direction_column(frame)
        if dcol is None or "region_or_subregion" not in frame.columns:
            continue
        for (region, direction), sub in frame.groupby(["region_or_subregion", dcol]):
            genes = unique_nonempty_symbols(sub["gene_symbol"])
            if len(genes) >= MIN_DIRECTION_SET_GENES:
                sets.append(
                    make_candidate_set(
                        set_name,
                        f"region:{region}|direction:{direction}",
                        sub,
                        f"{source_name} rows in {region} with direction {direction}; unique genes >= {MIN_DIRECTION_SET_GENES}",
                    )
                )
    return sets


def load_mdd_tested_universe() -> tuple[set[str], str]:
    candidates = [
        CLEAN_DIR / "01_inputs" / "tested_gene_universe_gene_symbols.csv",
        CLEAN_DIR / "03_broad_region_LMM" / "LMM_all_genes_length_gt_100.csv",
        CLEAN_DIR / "04_standardized_effect_size" / "LMM_all_genes_with_standardized_effect_size_4groups.csv",
        SOURCE_LENGTH100_DIR / "LMM_all_genes_length_gt_100.csv",
        SOURCE_LENGTH100_DIR / "standardized_effect_size_4groups" / "LMM_all_genes_with_standardized_effect_size_4groups.csv",
    ]
    for path in candidates:
        if not path.exists():
            continue
        df = pd.read_csv(path, usecols=lambda c: c in {"gene_symbol", "gene_id"})
        if "gene_symbol" in df.columns:
            genes = unique_nonempty_symbols(df["gene_symbol"])
            if genes:
                return genes, str(path)
    raise FileNotFoundError("Could not find a clean MDD tested-gene table with gene_symbol.")


def load_published_modules() -> pd.DataFrame:
    if not MODULE_TABLE.exists():
        raise FileNotFoundError(
            f"Published AHBA module table not found: {MODULE_TABLE}. "
            "Run code/ahba_modules/04_fetch_or_load_published_ahba_modules.py after placing Supplementary Table 6 in data/ahba/published_modules/."
        )
    modules = pd.read_csv(MODULE_TABLE)
    required = {"gene_symbol", "published_module"}
    missing = sorted(required - set(modules.columns))
    if missing:
        raise ValueError(f"Published module table is missing required columns: {missing}")
    modules["gene_symbol_norm"] = modules["gene_symbol"].map(norm_symbol)
    modules["published_module"] = modules["published_module"].astype(str).str.strip()
    modules = modules[modules["gene_symbol_norm"].ne("") & modules["published_module"].ne("")]
    modules = modules.drop_duplicates(["gene_symbol_norm", "published_module"])
    if modules.empty:
        raise ValueError("Published AHBA module table has no usable gene/module rows after cleaning.")
    return modules


def safe_label(text: str, width: int = 52) -> str:
    return "\n".join(wrap(text, width=width))


def write_plots(res: pd.DataFrame) -> list[str]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    created: list[str] = []
    suggestive = res[(res["fdr_BH"] < 0.10) | (res["p_value_hypergeometric"] < 0.05)].copy()
    if suggestive.empty:
        (FIG_DIR / "NO_SUGGESTIVE_PUBLISHED_MODULE_ENRICHMENT.txt").write_text(
            "No published AHBA module enrichment rows passed FDR < 0.10 or nominal p < 0.05, so no summary plot was created.\n",
            encoding="utf-8",
        )
        return created

    top = suggestive.sort_values(["fdr_BH", "p_value_hypergeometric", "overlap_genes"], ascending=[True, True, False]).head(30).copy()
    top["minus_log10_fdr"] = -np.log10(top["fdr_BH"].clip(lower=np.nextafter(0, 1)))
    labels = [safe_label(f"{r.candidate_set} | {r.scope} | {r.published_module}", 62) for r in top.itertuples()]

    fig, ax = plt.subplots(figsize=(12, max(6, 0.45 * len(top))))
    points = ax.scatter(
        top["minus_log10_fdr"],
        range(len(top)),
        s=35 + 30 * top["overlap_genes"],
        c=top["odds_ratio_approx"],
        cmap="viridis",
    )
    fig.colorbar(points, ax=ax, label="Approx. odds ratio")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("-log10(BH-FDR)")
    ax.set_title("Published AHBA module enrichment")
    fig.tight_layout()
    dotplot = FIG_DIR / "published_ahba_module_enrichment_dotplot.png"
    fig.savefig(dotplot, dpi=240)
    plt.close(fig)
    created.append(str(dotplot))

    fig, ax = plt.subplots(figsize=(12, max(6, 0.45 * len(top))))
    y = np.arange(len(top))
    ax.barh(y, top["overlap_genes"], color="#4C78A8")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Candidate genes overlapping module")
    ax.set_title("Published AHBA module candidate-gene overlaps")
    fig.tight_layout()
    barplot = FIG_DIR / "published_ahba_module_enrichment_barplot.png"
    fig.savefig(barplot, dpi=240)
    plt.close(fig)
    created.append(str(barplot))
    return created


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        modules = load_published_modules()
    except FileNotFoundError as exc:
        msg = str(exc)
        (OUT_DIR / "PUBLISHED_AHBA_ENRICHMENT_NOT_RUN.txt").write_text(msg + "\n", encoding="utf-8")
        print(msg)
        return
    candidate_sets = load_candidates()
    mdd_universe, mdd_source = load_mdd_tested_universe()

    module_universe = set(modules["gene_symbol_norm"])
    background = sorted(mdd_universe & module_universe)
    bg_set = set(background)
    module_groups = {
        module: set(group["gene_symbol_norm"]) & bg_set
        for module, group in modules.groupby("published_module")
    }
    module_groups = {module: genes for module, genes in module_groups.items() if genes}
    n_modules_tested = len(module_groups)
    M = len(bg_set)

    module_sizes = pd.DataFrame(
        [{"published_module": module, "background_gene_count": len(genes)} for module, genes in module_groups.items()]
    ).sort_values("background_gene_count", ascending=False)
    module_sizes.to_csv(OUT_DIR / "published_ahba_module_size_table.csv", index=False)
    pd.DataFrame({"background_gene_symbol": background}).to_csv(OUT_DIR / "published_ahba_mdd_background_genes.csv", index=False)

    audit_rows = []
    for cset in candidate_sets:
        unique_genes = set(cset["genes"])
        overlap = unique_genes & bg_set
        audit_rows.append({
            "candidate_set": cset["candidate_set"],
            "scope": cset["scope"],
            "original_gene_region_pair_count": cset["original_gene_region_pair_count"],
            "unique_gene_count": len(unique_genes),
            "unique_genes_overlapping_published_AHBA_background": len(overlap),
            "modules_tested": n_modules_tested,
            "reason_included": cset["reason_included"],
        })
    audit = pd.DataFrame(audit_rows).sort_values(["candidate_set", "scope"])
    audit_path = OUT_DIR / "published_ahba_candidate_set_audit.csv"
    audit.to_csv(audit_path, index=False)
    print("Candidate-set audit before published AHBA module enrichment:")
    print(audit.to_string(index=False))

    rows = []
    for cset in candidate_sets:
        genes = set(cset["genes"]) & bg_set
        N = len(genes)
        if M == 0 or N == 0:
            continue
        for module, module_genes in module_groups.items():
            K = len(module_genes)
            x = len(genes & module_genes)
            p_value = hypergeom.sf(x - 1, M, K, N) if x > 0 else 1.0
            rows.append({
                "candidate_set": cset["candidate_set"],
                "scope": cset["scope"],
                "original_gene_region_pair_count": cset["original_gene_region_pair_count"],
                "unique_candidate_genes": len(cset["genes"]),
                "published_module": module,
                "background_genes": M,
                "candidate_genes_in_background": N,
                "module_genes_in_background": K,
                "overlap_genes": x,
                "overlap_gene_symbols": ";".join(sorted(genes & module_genes)),
                "p_value_hypergeometric": p_value,
                "odds_ratio_approx": ((x + 0.5) / (N - x + 0.5)) / ((K - x + 0.5) / (M - N - K + x + 0.5)),
            })

    res = pd.DataFrame(rows)
    if res.empty:
        msg = "No enrichment rows generated after MDD tested-gene and published AHBA module background intersection."
        (OUT_DIR / "PUBLISHED_AHBA_ENRICHMENT_NOT_RUN.txt").write_text(msg + "\n", encoding="utf-8")
        print(msg)
        return

    res["fdr_BH"] = res.groupby(["candidate_set", "scope"], group_keys=False)["p_value_hypergeometric"].apply(bh_fdr)
    res = res.sort_values(["fdr_BH", "p_value_hypergeometric", "candidate_set", "scope"])
    results_path = OUT_DIR / "published_ahba_module_enrichment_results.csv"
    res.to_csv(results_path, index=False)
    plot_paths = write_plots(res)

    report = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "published_module_assignment_table": str(MODULE_TABLE),
        "mdd_tested_universe_source": mdd_source,
        "background_definition": "unique genes present in both the MDD clean-pipeline tested-gene universe and the published AHBA module assignment table",
        "background_gene_count": M,
        "published_modules_tested": n_modules_tested,
        "candidate_set_audit_file": str(audit_path),
        "results_file": str(results_path),
        "module_size_table": str(OUT_DIR / "published_ahba_module_size_table.csv"),
        "plots_created": plot_paths,
        "interpretation": "Published AHBA canonical modules provide healthy-brain coexpression context for ML-prioritized MDD candidate genes. They do not validate these genes as disease biomarkers.",
    }
    (OUT_DIR / "published_ahba_module_enrichment_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
