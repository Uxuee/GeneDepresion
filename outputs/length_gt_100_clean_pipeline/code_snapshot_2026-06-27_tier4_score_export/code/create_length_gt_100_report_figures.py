"""Organize and generate report figures for the Length > 100 MDD analysis.

This script only uses completed output tables and existing plots. It does not
rerun the LMM, enrichment, ML pilot, or ML-prioritized gene extraction.
"""

from __future__ import annotations

import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(r"D:\genes\GeneDepresion")
ANALYSIS = ROOT / "outputs" / "length_gt_100_elisa_style"
IMAGES = ANALYSIS / "images"
TIER2_DIR = IMAGES / "tier2"
MLV1_DIR = IMAGES / "ml_v1"
MLP_DIR = IMAGES / "ml_prioritized"
CONTACT_DIR = IMAGES / "contact_sheets"
INVENTORY_DIR = IMAGES / "inventory"

PLOT_EXTS = {".png", ".jpg", ".jpeg", ".svg", ".pdf"}

AREA_ORDER = [
    "Frontal Cortex",
    "Insula",
    "Ventral Striatum",
    "Hippocampal Formation",
    "Cingulate gyrus 25",
    "dlPFC BA8/9",
    "OFC BA11",
]

SLUG_LABELS = {
    "Frontal_Cortex": "Frontal Cortex",
    "Insula": "Insula",
    "Ventral_Striatum": "Ventral Striatum",
    "Hippocampal_Formation": "Hippocampal Formation",
    "Cingulate_gyrus_25_Cg25": "Cingulate gyrus 25",
    "dlPFC_BA8_9": "dlPFC BA8/9",
    "OFC_BA11": "OFC BA11",
}

records: list[dict[str, object]] = []
copied_count = 0
generated_count = 0


def ensure_dirs() -> None:
    for folder in [TIER2_DIR, MLV1_DIR, MLP_DIR, CONTACT_DIR, INVENTORY_DIR]:
        folder.mkdir(parents=True, exist_ok=True)


def rel(path: Path) -> str:
    return str(path.relative_to(ANALYSIS)).replace("\\", "/")


def record(
    category: str,
    dest: Path,
    source: Path | None,
    existed_already: bool,
    action: str,
    description: str,
    recommended: bool,
    notes: str = "",
) -> None:
    records.append(
        {
            "category": category,
            "filename": dest.name,
            "relative_path": rel(dest),
            "source_path": str(source) if source else "",
            "existed_already": existed_already,
            "copied_or_generated": action,
            "plot_description": description,
            "recommended_for_main_report": recommended,
            "notes": notes,
        }
    )


def copy_if_missing(
    source: Path,
    dest: Path,
    category: str,
    description: str,
    recommended: bool,
    notes: str = "",
) -> None:
    global copied_count
    existed = dest.exists()
    if source.exists() and not existed:
        shutil.copy2(source, dest)
        copied_count += 1
        action = "copied"
    elif source.exists():
        action = "already_present"
    else:
        action = "missing_source"
        notes = notes or "Missing; underlying data or plotting source not found."
    record(category, dest, source, existed, action, description, recommended, notes)


def savefig_if_missing(fig: plt.Figure, dest: Path, category: str, description: str, recommended: bool) -> None:
    global generated_count
    existed = dest.exists()
    if not existed:
        fig.savefig(dest, dpi=180, bbox_inches="tight")
        generated_count += 1
        action = "generated"
    else:
        action = "already_present"
    plt.close(fig)
    record(category, dest, None, existed, action, description, recommended)


def find_all_existing_plots() -> list[Path]:
    return [
        p
        for p in ANALYSIS.rglob("*")
        if p.is_file() and p.suffix.lower() in PLOT_EXTS
    ]


def load_tier2_summary() -> pd.DataFrame:
    summary = ANALYSIS / "volcano_plots_Tier2_all_areas" / "Tier2_volcano_all_areas_summary.csv"
    if summary.exists():
        df = pd.read_csv(summary)
        df["area_name"] = df["area_name"].replace({"Cingulate gyrus 25 (Cg25)": "Cingulate gyrus 25"})
        return df

    candidate_dir = ANALYSIS / "DE_candidate_gene_sets" / "Tier2_working_exploratory"
    rows = []
    for csv in candidate_dir.glob("*__Tier2_all.csv"):
        name = csv.name.replace("__Tier2_all.csv", "").replace("broad_", "").replace("_", " ")
        df = pd.read_csv(csv)
        beta_col = "standardized_effect_size_diagnosis"
        rows.append(
            {
                "area_name": name,
                "genes_tested": np.nan,
                "tier2_genes": len(df),
                "tier2_up_in_MDD": int((df[beta_col] > 0).sum()) if beta_col in df else np.nan,
                "tier2_down_in_MDD": int((df[beta_col] < 0).sum()) if beta_col in df else np.nan,
            }
        )
    return pd.DataFrame(rows)


def ordered(df: pd.DataFrame, col: str = "area_name") -> pd.DataFrame:
    order = {name: i for i, name in enumerate(AREA_ORDER)}
    return df.assign(_order=df[col].map(order).fillna(999)).sort_values("_order").drop(columns="_order")


def generate_tier2_counts() -> None:
    df = ordered(load_tier2_summary())
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    colors = ["#8b8f96" if x == "Frontal Cortex" else "#36688d" for x in df["area_name"]]
    ax.barh(df["area_name"], df["tier2_genes"], color=colors)
    ax.invert_yaxis()
    ax.set_xlabel("Tier 2 genes")
    ax.set_title("Tier 2 exploratory genes by region/subregion")
    for y, value in enumerate(df["tier2_genes"]):
        ax.text(value + max(df["tier2_genes"]) * 0.015, y, f"{int(value):,}", va="center", fontsize=9)
    ax.grid(axis="x", color="#e5e8eb")
    savefig_if_missing(fig, TIER2_DIR / "Tier2_counts_by_region.png", "tier2", "Tier 2 counts by region/subregion.", True)


def generate_tier2_up_down() -> None:
    df = ordered(load_tier2_summary())
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.barh(df["area_name"], df["tier2_down_in_MDD"], color="#2c7fb8", label="Down in MDD")
    ax.barh(df["area_name"], df["tier2_up_in_MDD"], left=df["tier2_down_in_MDD"], color="#d95f02", label="Up in MDD")
    ax.invert_yaxis()
    ax.set_xlabel("Tier 2 genes")
    ax.set_title("Direction of Tier 2 exploratory genes")
    ax.legend(loc="lower right")
    ax.grid(axis="x", color="#e5e8eb")
    savefig_if_missing(fig, TIER2_DIR / "Tier2_up_down_by_region.png", "tier2", "Tier 2 up/down counts by region/subregion.", True)


def copy_tier2_volcanoes() -> None:
    source_dir = ANALYSIS / "volcano_plots_Tier2_all_areas"
    mapping = {
        "Frontal_Cortex": "Tier2_volcano_FrontalCortex.png",
        "Insula": "Tier2_volcano_Insula.png",
        "Ventral_Striatum": "Tier2_volcano_VentralStriatum.png",
        "Hippocampal_Formation": "Tier2_volcano_HippocampalFormation.png",
        "Cingulate_gyrus_25_Cg25": "Tier2_volcano_Cg25.png",
        "dlPFC_BA8_9": "Tier2_volcano_dlPFC_BA8_9.png",
        "OFC_BA11": "Tier2_volcano_OFC_BA11.png",
    }
    for slug, filename in mapping.items():
        source = source_dir / f"{slug}_Tier2_volcano_standardized_effect.png"
        copy_if_missing(
            source,
            TIER2_DIR / filename,
            "tier2",
            f"Tier 2 volcano plot for {SLUG_LABELS[slug]}.",
            filename in {"Tier2_volcano_Cg25.png", "Tier2_volcano_VentralStriatum.png"},
        )


def generate_tier2_overlap() -> None:
    candidate_dir = ANALYSIS / "DE_candidate_gene_sets" / "Tier2_working_exploratory"
    sets = {}
    for path in candidate_dir.glob("*__Tier2_all.csv"):
        label = path.name.replace("__Tier2_all.csv", "").replace("broad_", "")
        label = label.replace("Cingulate_gyrus_25", "Cingulate gyrus 25").replace("dlPFC_BA8_9", "dlPFC BA8/9")
        label = label.replace("OFC_BA11", "OFC BA11").replace("_", " ")
        df = pd.read_csv(path, usecols=lambda c: c in {"gene_id", "gene_symbol"})
        sets[label] = set(df["gene_id"].astype(str))
    if not sets:
        record("tier2", TIER2_DIR / "Tier2_overlap_upset.png", None, False, "skipped", "Tier 2 overlap plot.", False, "Missing; underlying data or plotting source not found.")
        return
    labels = [x for x in AREA_ORDER if x in sets]
    intersections = []
    all_genes = sorted(set().union(*[sets[x] for x in labels]))
    for gene in all_genes:
        mask = tuple(gene in sets[label] for label in labels)
        intersections.append(mask)
    combo_counts = pd.Series(intersections).value_counts().head(20)
    fig = plt.figure(figsize=(10, 6.4))
    grid = fig.add_gridspec(2, 1, height_ratios=[2.2, 1.4], hspace=0.05)
    ax_bar = fig.add_subplot(grid[0])
    ax_dot = fig.add_subplot(grid[1], sharex=ax_bar)
    x = np.arange(len(combo_counts))
    ax_bar.bar(x, combo_counts.values, color="#36688d")
    ax_bar.set_ylabel("Shared Tier 2 genes")
    ax_bar.set_title("Largest Tier 2 gene overlaps across regions/subregions")
    ax_bar.grid(axis="y", color="#e5e8eb")
    for i, mask in enumerate(combo_counts.index):
        for y, present in enumerate(mask):
            ax_dot.scatter(i, y, s=34, color="#1f2528" if present else "#d1d6da")
        present_y = [y for y, present in enumerate(mask) if present]
        if len(present_y) > 1:
            ax_dot.plot([i, i], [min(present_y), max(present_y)], color="#1f2528", linewidth=1)
    ax_dot.set_yticks(range(len(labels)), labels)
    ax_dot.set_xlabel("Intersection combination")
    ax_dot.set_xticks([])
    savefig_if_missing(fig, TIER2_DIR / "Tier2_overlap_upset.png", "tier2", "UpSet-style Tier 2 overlap plot.", False)


def generate_mlv1_performance() -> None:
    perf = pd.read_csv(ANALYSIS / "ML_pilot_Tier2" / "best_model_by_region.csv")
    perf["label"] = perf["region_or_subregion"].str.replace("Hippocampal Formation", "Hippocampal\nFormation")
    perf = perf.sort_values("roc_auc_mean", ascending=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True)
    axes[0].barh(perf["label"], perf["balanced_accuracy_mean"], color="#4878a8")
    axes[0].set_title("Balanced accuracy")
    axes[0].set_xlim(0.45, 1.02)
    axes[1].barh(perf["label"], perf["roc_auc_mean"], color="#d95f02")
    axes[1].set_title("ROC-AUC")
    axes[1].set_xlim(0.45, 1.02)
    for ax in axes:
        ax.axvline(0.5, color="#555", linestyle="--", linewidth=1)
        ax.grid(axis="x", color="#e5e8eb")
    fig.suptitle("ML v1 exploratory best model performance by area")
    savefig_if_missing(fig, MLV1_DIR / "MLv1_performance_summary.png", "ml_v1", "ML v1 best-model balanced accuracy and ROC-AUC.", True)


def make_image_grid(sources: list[Path], dest: Path, title: str, category: str, description: str, recommended: bool) -> None:
    global generated_count
    existed = dest.exists()
    valid = [p for p in sources if p.exists()]
    if not valid:
        record(category, dest, None, existed, "skipped", description, recommended, "Missing; underlying data or plotting source not found.")
        return
    if not existed:
        thumbs = []
        for p in valid:
            img = Image.open(p).convert("RGB")
            img.thumbnail((520, 390))
            thumbs.append((p, img.copy()))
            img.close()
        cols = 2
        rows = math.ceil(len(thumbs) / cols)
        width = cols * 560
        height = rows * 440 + 60
        sheet = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(sheet)
        draw.text((20, 18), title, fill=(25, 31, 35))
        for i, (path, img) in enumerate(thumbs):
            x = (i % cols) * 560 + 20
            y = (i // cols) * 440 + 60
            sheet.paste(img, (x, y))
            draw.text((x, y + img.height + 4), path.name[:68], fill=(60, 68, 74))
        sheet.save(dest)
        generated_count += 1
        action = "generated"
    else:
        action = "already_present"
    record(category, dest, None, existed, action, description, recommended)


def copy_and_grid_mlv1() -> None:
    plot_dir = ANALYSIS / "ML_pilot_Tier2" / "plots"
    copy_if_missing(plot_dir / "top_lasso_elastic_net_feature_stability.png", MLV1_DIR / "MLv1_feature_stability.png", "ml_v1", "ML v1 feature stability/importance.", True)
    copy_if_missing(plot_dir / "heatmap_roc_auc_feature_sets_models.png", MLV1_DIR / "MLv1_performance_heatmap.png", "ml_v1", "ML v1 ROC-AUC heatmap across feature sets/models.", False)
    roc_sources = sorted(plot_dir.glob("*__ROC_curve.png"))
    pr_sources = sorted(plot_dir.glob("*__PR_curve.png"))
    make_image_grid(roc_sources, MLV1_DIR / "MLv1_best_ROC_curve.png", "Exploratory ML v1 ROC curves", "ml_v1", "Best ROC curves by region/subregion.", False)
    make_image_grid(pr_sources, MLV1_DIR / "MLv1_best_PR_curve.png", "Exploratory ML v1 PR curves", "ml_v1", "Best precision-recall curves by region/subregion.", False)
    if not (MLV1_DIR / "MLv1_workflow_summary.png").exists():
        fig, ax = plt.subplots(figsize=(9, 2.8))
        ax.axis("off")
        steps = ["Tier 2\ncandidate genes", "Classical\nML models", "Fold-level\nfeature importance", "Stable\nML-prioritized genes"]
        xs = np.linspace(0.1, 0.9, len(steps))
        for i, (x, step) in enumerate(zip(xs, steps)):
            ax.text(x, 0.5, step, ha="center", va="center", bbox={"boxstyle": "round,pad=0.45", "facecolor": "#eef3f5", "edgecolor": "#8091a0"})
            if i < len(steps) - 1:
                ax.annotate("", xy=(xs[i + 1] - 0.08, 0.5), xytext=(x + 0.08, 0.5), arrowprops={"arrowstyle": "->", "color": "#46525c", "lw": 1.5})
        ax.set_title("ML v1 exploratory gene-prioritization workflow")
        savefig_if_missing(fig, MLV1_DIR / "MLv1_workflow_summary.png", "ml_v1", "Simple ML v1 workflow summary.", False)
    else:
        record("ml_v1", MLV1_DIR / "MLv1_workflow_summary.png", None, True, "already_present", "Simple ML v1 workflow summary.", False)


def generate_ml_prioritized() -> None:
    ranked = pd.read_csv(ANALYSIS / "ML_prioritized_Tier2_genes" / "ML_prioritized_Tier2_genes_all_ranked.csv")
    summary = pd.read_csv(ANALYSIS / "ML_prioritized_Tier2_genes" / "ML_prioritized_Tier2_gene_summary_by_region.csv")
    summary = ordered(summary.rename(columns={"region_or_subregion": "area_name"}))

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.barh(summary["area_name"], summary["n_moderate_confidence"], color="#9bc2d1", label="Moderate confidence")
    ax.barh(summary["area_name"], summary["n_high_confidence"], color="#2b6f93", label="High confidence")
    ax.invert_yaxis()
    ax.set_xlabel("ML-prioritized Tier 2 genes")
    ax.set_title("ML-prioritized Tier 2 genes by area")
    ax.legend(loc="lower right")
    ax.grid(axis="x", color="#e5e8eb")
    savefig_if_missing(fig, MLP_DIR / "MLprioritized_counts_by_region.png", "ml_prioritized", "High- and moderate-confidence ML-prioritized Tier 2 counts.", True)

    top = ranked.sort_values(["region_or_subregion", "ML_priority_score"], ascending=[True, False]).groupby("region_or_subregion").head(8)
    top["label"] = top["region_or_subregion"] + " | " + top["gene_symbol"].fillna(top["gene_id"].astype(str))
    top = top.sort_values("ML_priority_score", ascending=True)
    fig, ax = plt.subplots(figsize=(9.5, 9.0))
    colors = top["ML_priority_class"].map(lambda x: "#1b6e8f" if str(x).startswith("High") else "#8ab6c6")
    ax.barh(top["label"], top["ML_priority_score"], color=colors)
    ax.set_xlabel("ML priority score")
    ax.set_title("Top ML-prioritized Tier 2 genes by area")
    ax.grid(axis="x", color="#e5e8eb")
    savefig_if_missing(fig, MLP_DIR / "MLprioritized_top_genes_by_region.png", "ml_prioritized", "Top ML-prioritized genes per area.", True)

    class_colors = {"High-confidence ML-prioritized Tier 2": "#1b6e8f", "Moderate-confidence ML-prioritized Tier 2": "#d95f02"}
    for xcol, xlabel, filename, desc in [
        ("standardized_effect_size_diagnosis", "Standardized effect size", "MLpriority_vs_standardized_effect_size.png", "ML priority score vs standardized effect size."),
        ("p_value_diagnosis", "-log10(p-value)", "MLpriority_vs_neglog10p.png", "ML priority score vs -log10 diagnosis p-value."),
    ]:
        fig, ax = plt.subplots(figsize=(7.8, 5.4))
        x = -np.log10(ranked[xcol].clip(lower=np.nextafter(0, 1))) if xcol == "p_value_diagnosis" else ranked[xcol]
        for cls, sub in ranked.assign(_x=x).groupby("ML_priority_class"):
            ax.scatter(sub["_x"], sub["ML_priority_score"], s=34, alpha=0.75, color=class_colors.get(cls, "#7d8790"), label=cls.replace(" ML-prioritized Tier 2", ""))
        ax.set_xlabel(xlabel)
        ax.set_ylabel("ML priority score")
        ax.set_title(desc)
        ax.grid(color="#e5e8eb")
        ax.legend(fontsize=8)
        savefig_if_missing(fig, MLP_DIR / filename, "ml_prioritized", desc, filename != "MLpriority_vs_neglog10p.png")

    support_cols = ["sparse_model_support", "linear_model_support", "rf_support", "v2_support_flag", "lmm_ml_sign_agreement"]
    heat = ranked.sort_values("ML_priority_score", ascending=False).head(30).copy()
    matrix = heat[support_cols].astype(float).to_numpy()
    fig, ax = plt.subplots(figsize=(8.2, 8.8))
    ax.imshow(matrix, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_yticks(range(len(heat)), heat["gene_symbol"].fillna(heat["gene_id"].astype(str)))
    ax.set_xticks(range(len(support_cols)), ["Sparse\nmodels", "Linear\nmodels", "Random\nforest", "Strict v2\nsupport", "LMM/ML\nsign agree"])
    ax.set_title("Model support for top ML-prioritized Tier 2 genes")
    savefig_if_missing(fig, MLP_DIR / "MLpriority_model_support_heatmap.png", "ml_prioritized", "Support heatmap for top ML-prioritized genes.", True)

    agreement = ranked.groupby(["region_or_subregion", "lmm_ml_sign_agreement"]).size().unstack(fill_value=0)
    agreement = agreement.rename(columns={0.0: "Disagree", 1.0: "Agree", False: "Disagree", True: "Agree"}).reset_index()
    agreement = ordered(agreement.rename(columns={"region_or_subregion": "area_name"}))
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    agree = agreement.get("Agree", pd.Series([0] * len(agreement)))
    disagree = agreement.get("Disagree", pd.Series([0] * len(agreement)))
    ax.barh(agreement["area_name"], disagree, color="#c55a4d", label="Disagree")
    ax.barh(agreement["area_name"], agree, left=disagree, color="#5c8f55", label="Agree")
    ax.invert_yaxis()
    ax.set_xlabel("ML-prioritized genes")
    ax.set_title("LMM and ML direction agreement")
    ax.legend(loc="lower right")
    savefig_if_missing(fig, MLP_DIR / "MLpriority_sign_agreement.png", "ml_prioritized", "LMM/ML sign agreement summary.", False)


def make_contact_sheets() -> None:
    make_image_grid(sorted(TIER2_DIR.glob("*.png")), CONTACT_DIR / "tier2_contact_sheet.png", "Tier 2 report figures", "contact_sheets", "Contact sheet for Tier 2 figures.", False)
    make_image_grid(sorted(MLV1_DIR.glob("*.png")), CONTACT_DIR / "ml_v1_contact_sheet.png", "ML v1 report figures", "contact_sheets", "Contact sheet for ML v1 figures.", False)
    make_image_grid(sorted(MLP_DIR.glob("*.png")), CONTACT_DIR / "ml_prioritized_contact_sheet.png", "ML-prioritized report figures", "contact_sheets", "Contact sheet for ML-prioritized figures.", False)


def write_inventory(found_count: int) -> None:
    inventory = pd.DataFrame(records).drop_duplicates(subset=["relative_path"], keep="last")
    inventory_path = INVENTORY_DIR / "image_inventory.csv"
    inventory.to_csv(inventory_path, index=False)

    main_figures = [
        "images/tier2/Tier2_counts_by_region.png",
        "images/tier2/Tier2_up_down_by_region.png",
        "images/tier2/Tier2_volcano_Cg25.png",
        "images/tier2/Tier2_volcano_VentralStriatum.png",
        "images/ml_v1/MLv1_performance_summary.png",
        "images/ml_v1/MLv1_feature_stability.png",
        "images/ml_prioritized/MLprioritized_counts_by_region.png",
        "images/ml_prioritized/MLprioritized_top_genes_by_region.png",
        "images/ml_prioritized/MLpriority_vs_standardized_effect_size.png",
        "images/ml_prioritized/MLpriority_model_support_heatmap.png",
    ]
    supplementary = [
        "images/ml_v1/MLv1_best_ROC_curve.png",
        "images/ml_v1/MLv1_best_PR_curve.png",
        "images/tier2/Tier2_overlap_upset.png",
        "images/ml_prioritized/MLpriority_vs_neglog10p.png",
        "images/ml_prioritized/MLpriority_sign_agreement.png",
    ]
    skipped = inventory[inventory["copied_or_generated"].isin(["skipped", "missing_source"])]

    readme = INVENTORY_DIR / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# Length > 100 Report Figure Inventory",
                "",
                f"Images folder: `{IMAGES}`",
                f"Total plots found recursively before organization: {found_count}",
                f"Total copied into central folder in this run: {copied_count}",
                f"Total newly generated in this run: {generated_count}",
                "",
                "No LMM, enrichment, ML v1, or stable ML-prioritized extraction was rerun. Figures were copied or generated from completed outputs only.",
                "",
                "## Tier 2 Plots",
                "Includes Tier 2 counts, direction counts, all-area standardized-effect volcano plots, and an UpSet-style overlap summary.",
                "",
                "## ML v1 Plots",
                "ML v1 figures summarize exploratory Tier2-guided ML performance and feature stability. These are not validated diagnostic classifier results.",
                "",
                "## ML-Prioritized Plots",
                "ML-prioritized figures summarize high/moderate confidence counts, top prioritized genes, score relationships, model support, and sign agreement.",
                "",
                "## Preferred Main Report Figure Set",
                *[f"- `{x}`" for x in main_figures],
                "",
                "## Supplementary Figure Candidates",
                *[f"- `{x}`" for x in supplementary],
                "",
                "## Missing or Skipped Figures",
                *(
                    [f"- `{row.relative_path}`: {row.notes}" for row in skipped.itertuples()]
                    if len(skipped)
                    else ["- None. Requested plots were present, copied, or generated from existing outputs."]
                ),
                "",
                "## Inventory",
                f"See `{inventory_path}` for source paths and copy/generation status.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    record("inventory", inventory_path, None, inventory_path.exists(), "generated", "CSV inventory of central report figures.", False)
    record("inventory", readme, None, readme.exists(), "generated", "Markdown guide to central report figures.", False)


def main() -> None:
    ensure_dirs()
    found = find_all_existing_plots()
    generate_tier2_counts()
    generate_tier2_up_down()
    copy_tier2_volcanoes()
    generate_tier2_overlap()
    generate_mlv1_performance()
    copy_and_grid_mlv1()
    generate_ml_prioritized()
    make_contact_sheets()
    write_inventory(len(found))

    print(f"plots_found={len(found)}")
    print(f"copied={copied_count}")
    print(f"generated={generated_count}")
    print(f"images_folder={IMAGES}")
    print(f"inventory_csv={INVENTORY_DIR / 'image_inventory.csv'}")
    print(f"inventory_readme={INVENTORY_DIR / 'README.md'}")


if __name__ == "__main__":
    main()
