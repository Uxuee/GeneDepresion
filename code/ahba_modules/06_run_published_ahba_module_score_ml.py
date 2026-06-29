"""Run exploratory ML using published AHBA module scores as predefined features."""

from __future__ import annotations

import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedGroupKFold, StratifiedKFold
from sklearn.svm import LinearSVC
import warnings

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", message="y_pred contains classes not in y_true")

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CLEAN_DIR = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
EXPR_DIR = CLEAN_DIR / "02_normalized_expression"
META_PATH = CLEAN_DIR / "01_inputs" / "sample_metadata_used.csv"
LMM_PATH = CLEAN_DIR / "03_broad_region_LMM" / "LMM_all_genes_length_gt_100.csv"
PUB_DIR = CLEAN_DIR / "ahba_modules" / "published_modules"
MODULE_PATH = PUB_DIR / "ahba_published_gene_to_module_assignments.csv"
ENR_DIR = PUB_DIR / "enrichment"
OUT_DIR = CLEAN_DIR / "ahba_modules" / "module_score_ml"
MATRIX_DIR = OUT_DIR / "module_score_matrices"
PLOT_DIR = OUT_DIR / "plots"
IMPORTANCE_DIR = OUT_DIR / "feature_importance"
N_REPEATS = 20
RANDOM_SEEDS = list(range(1001, 1001 + N_REPEATS))

REGION_FILES = {
    "Frontal_Cortex": EXPR_DIR / "normalized_GeTMM_log2_frontal_cortex.csv",
    "Hippocampal_Formation": EXPR_DIR / "normalized_GeTMM_log2_hippocampal_formation.csv",
    "Insula": EXPR_DIR / "normalized_GeTMM_log2_insula.csv",
    "Ventral_Striatum": EXPR_DIR / "normalized_GeTMM_log2_ventral_striatum.csv",
}


def norm_symbol(x: object) -> str:
    return "" if pd.isna(x) else str(x).strip().upper()


def load_modules() -> dict[str, set[str]]:
    modules = pd.read_csv(MODULE_PATH)
    modules["gene_symbol_norm"] = modules["gene_symbol"].map(norm_symbol)
    return {
        module: set(sub["gene_symbol_norm"]) - {""}
        for module, sub in modules.groupby("published_module")
    }


def load_gene_map() -> dict[str, str]:
    lmm = pd.read_csv(LMM_PATH, usecols=["gene_id", "gene_symbol"])
    lmm["gene_id"] = lmm["gene_id"].astype(str)
    lmm["gene_symbol_norm"] = lmm["gene_symbol"].map(norm_symbol)
    lmm = lmm[lmm["gene_symbol_norm"].ne("")].drop_duplicates("gene_id")
    return dict(zip(lmm["gene_id"], lmm["gene_symbol_norm"]))


def load_metadata() -> pd.DataFrame:
    meta = pd.read_csv(META_PATH)
    meta["sample_id"] = meta["geo_accession"].astype(str)
    meta["diagnosis_binary"] = meta["Diagnosis"].map({"MDD": 1, "CTRL": 0, "Control": 0, "CONTROL": 0})
    return meta[meta["diagnosis_binary"].notna()].copy()


def read_expression(path: Path, gene_map: dict[str, str]) -> pd.DataFrame:
    expr = pd.read_csv(path)
    gene_col = expr.columns[0]
    expr[gene_col] = expr[gene_col].astype(str)
    expr["gene_symbol_norm"] = expr[gene_col].map(gene_map).fillna(expr[gene_col].map(norm_symbol))
    expr = expr[expr["gene_symbol_norm"].ne("")]
    sample_cols = [c for c in expr.columns if c not in {gene_col, "gene_symbol_norm"}]
    expr[sample_cols] = expr[sample_cols].apply(pd.to_numeric, errors="coerce")
    expr = expr.groupby("gene_symbol_norm", as_index=True)[sample_cols].mean()
    return expr


def module_scores_from_scaled(gene_by_sample: pd.DataFrame, modules: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = {}
    coverage = []
    available = set(gene_by_sample.index)
    for module, genes in sorted(modules.items(), key=lambda kv: int(kv[0][1:]) if kv[0].startswith("M") and kv[0][1:].isdigit() else kv[0]):
        overlap = sorted(genes & available)
        coverage.append({
            "published_module": module,
            "module_gene_count": len(genes),
            "genes_available_in_expression": len(overlap),
            "fraction_available": len(overlap) / len(genes) if genes else np.nan,
        })
        if overlap:
            rows[module] = gene_by_sample.loc[overlap].mean(axis=0)
        else:
            rows[module] = pd.Series(np.nan, index=gene_by_sample.columns)
    scores = pd.DataFrame(rows, index=gene_by_sample.columns)
    return scores, pd.DataFrame(coverage)


def descriptive_module_scores(expr: pd.DataFrame, modules: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    means = expr.mean(axis=1)
    sds = expr.std(axis=1, ddof=0).replace(0, np.nan)
    scaled = expr.sub(means, axis=0).div(sds, axis=0).fillna(0.0)
    return module_scores_from_scaled(scaled, modules)


def get_feature_sets(all_modules: list[str]) -> dict[str, list[str]]:
    corrected = pd.read_csv(ENR_DIR / "published_ahba_enrichment_corrected_results_table.csv")
    fdr05 = sorted(set(corrected.loc[corrected["BH_FDR"] < 0.05, "AHBA_module_ID"].astype(str)))
    fdr10 = sorted(set(corrected.loc[corrected["BH_FDR"] < 0.10, "AHBA_module_ID"].astype(str)))
    sets = {"all32_predefined_modules": all_modules}
    if fdr05:
        sets["post_enrichment_FDR005_modules_sensitivity"] = fdr05
    if fdr10:
        sets["post_enrichment_FDR010_modules_sensitivity"] = fdr10
    if "M26" in all_modules:
        sets["post_enrichment_M26_only_sensitivity"] = ["M26"]
    return sets


def make_cv(y: np.ndarray, groups: np.ndarray | None, seed: int):
    min_class = int(pd.Series(y).value_counts().min())
    n_splits = 5 if min_class >= 5 else 3
    if groups is not None and pd.Series(groups).duplicated().any():
        return StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed), "StratifiedGroupKFold", n_splits
    if min_class >= 5:
        return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed), "StratifiedKFold_repeated_seeds", n_splits
    return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed), "StratifiedKFold_small_n_repeated_seeds", n_splits


def model_specs():
    return {
        "majority_baseline": DummyClassifier(strategy="most_frequent"),
        "logistic_L2": LogisticRegression(penalty="l2", solver="liblinear", max_iter=2000),
        "logistic_L1": LogisticRegression(penalty="l1", solver="liblinear", max_iter=2000),
        "elastic_net": LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, max_iter=5000),
        "linear_SVM": LinearSVC(C=1.0, max_iter=5000),
        "random_forest_exploratory": RandomForestClassifier(n_estimators=120, max_depth=None, random_state=0, class_weight="balanced"),
    }


def decision_scores(model, X):
    if hasattr(model, "predict_proba"):
        return model.predict_proba(X)[:, 1]
    if hasattr(model, "decision_function"):
        return model.decision_function(X)
    return model.predict(X)


def fold_module_scores(expr: pd.DataFrame, train_samples: list[str], test_samples: list[str], modules: dict[str, set[str]], feature_modules: list[str]):
    train_expr = expr[train_samples]
    test_expr = expr[test_samples]
    means = train_expr.mean(axis=1)
    sds = train_expr.std(axis=1, ddof=0).replace(0, np.nan)
    train_scaled = train_expr.sub(means, axis=0).div(sds, axis=0).fillna(0.0)
    test_scaled = test_expr.sub(means, axis=0).div(sds, axis=0).fillna(0.0)
    train_scores, _ = module_scores_from_scaled(train_scaled, {m: modules[m] for m in feature_modules})
    test_scores, _ = module_scores_from_scaled(test_scaled, {m: modules[m] for m in feature_modules})
    return train_scores[feature_modules].fillna(0.0), test_scores[feature_modules].fillna(0.0)


def summarize_metrics(rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    group_cols = ["region", "feature_set", "n_features", "model", "cv_type", "n_samples", "n_MDD", "n_Control", "leakage_status"]
    metric_cols = ["accuracy", "balanced_accuracy", "roc_auc", "pr_auc", "f1", "sensitivity", "specificity"]
    agg = df.groupby(group_cols)[metric_cols].agg(["mean", "std"]).reset_index()
    agg.columns = ["_".join(c).strip("_") for c in agg.columns.to_flat_index()]
    return agg


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MATRIX_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    IMPORTANCE_DIR.mkdir(parents=True, exist_ok=True)
    modules = load_modules()
    all_modules = sorted(modules, key=lambda x: int(x[1:]) if x.startswith("M") and x[1:].isdigit() else x)
    feature_sets = get_feature_sets(all_modules)
    gene_map = load_gene_map()
    meta = load_metadata()

    performance_rows = []
    coef_rows = []
    coverage_all = []
    best_curves = []

    for region, path in REGION_FILES.items():
        print(f"Loading {region}: {path}")
        expr = read_expression(path, gene_map)
        samples = [s for s in expr.columns if s in set(meta["sample_id"])]
        expr = expr[samples]
        meta_region = meta.set_index("sample_id").loc[samples].reset_index()
        y = meta_region["diagnosis_binary"].astype(int).to_numpy()
        groups = meta_region["BrainID"].astype(str).to_numpy() if "BrainID" in meta_region.columns else None
        desc_scores, coverage = descriptive_module_scores(expr, modules)
        coverage.insert(0, "region", region)
        coverage_all.append(coverage)
        desc = desc_scores.reset_index().rename(columns={"index": "sample_id"})
        desc = desc.merge(meta_region[["sample_id", "Diagnosis", "BrainID", "BrainArea", "AnalysisGroup"]], on="sample_id", how="left")
        desc.to_csv(MATRIX_DIR / f"{region}_published_AHBA_module_scores_descriptive.csv", index=False)
        meta_region.to_csv(MATRIX_DIR / f"{region}_module_score_sample_metadata.csv", index=False)

        for seed in RANDOM_SEEDS:
            cv, cv_type, n_splits = make_cv(y, groups, seed)
            split_iter = cv.split(np.zeros(len(y)), y, groups) if "Group" in cv_type else cv.split(np.zeros(len(y)), y)
            for fold_idx, (train_idx, test_idx) in enumerate(split_iter, start=1):
                train_samples = meta_region.iloc[train_idx]["sample_id"].tolist()
                test_samples = meta_region.iloc[test_idx]["sample_id"].tolist()
                # Compute all 32 predefined module scores once per fold, then subset for sensitivity analyses.
                all_X_train, all_X_test = fold_module_scores(expr, train_samples, test_samples, modules, all_modules)
                y_train, y_test = y[train_idx], y[test_idx]
                for feature_set, feature_modules in feature_sets.items():
                    leakage = "predefined_external_all32_modules" if feature_set == "all32_predefined_modules" else "post_enrichment_selected_sensitivity_exploratory"
                    X_train = all_X_train[feature_modules]
                    X_test = all_X_test[feature_modules]
                    for model_name, model in model_specs().items():
                        if model_name == "random_forest_exploratory":
                            model.set_params(random_state=seed)
                        model.fit(X_train, y_train)
                        pred = model.predict(X_test)
                        score = decision_scores(model, X_test)
                        tn, fp, fn, tp = confusion_matrix(y_test, pred, labels=[0, 1]).ravel()
                        row = {
                            "region": region,
                            "feature_set": feature_set,
                            "n_features": len(feature_modules),
                            "model": model_name,
                            "cv_type": cv_type,
                            "seed": seed,
                            "fold": fold_idx,
                            "n_samples": len(y),
                            "n_MDD": int(y.sum()),
                            "n_Control": int((y == 0).sum()),
                            "accuracy": accuracy_score(y_test, pred),
                            "balanced_accuracy": balanced_accuracy_score(y_test, pred),
                            "roc_auc": roc_auc_score(y_test, score) if len(set(y_test)) == 2 else np.nan,
                            "pr_auc": average_precision_score(y_test, score) if len(set(y_test)) == 2 else np.nan,
                            "f1": f1_score(y_test, pred, zero_division=0),
                            "sensitivity": tp / (tp + fn) if (tp + fn) else np.nan,
                            "specificity": tn / (tn + fp) if (tn + fp) else np.nan,
                            "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
                            "leakage_status": leakage,
                        }
                        performance_rows.append(row)
                        if hasattr(model, "coef_"):
                            for module, coef in zip(feature_modules, model.coef_.ravel()):
                                coef_rows.append({
                                    "region": region, "feature_set": feature_set, "model": model_name,
                                    "seed": seed, "fold": fold_idx, "published_module": module,
                                    "coefficient": float(coef), "abs_coefficient": abs(float(coef)),
                                })
                        elif hasattr(model, "feature_importances_"):
                            for module, imp in zip(feature_modules, model.feature_importances_):
                                coef_rows.append({
                                    "region": region, "feature_set": feature_set, "model": model_name,
                                    "seed": seed, "fold": fold_idx, "published_module": module,
                                    "coefficient": float(imp), "abs_coefficient": abs(float(imp)),
                                })
                        if feature_set == "all32_predefined_modules" and model_name == "logistic_L2" and seed == RANDOM_SEEDS[0] and fold_idx == 1 and len(set(y_test)) == 2:
                            fpr, tpr, _ = roc_curve(y_test, score)
                            precision, recall, _ = precision_recall_curve(y_test, score)
                            for a, b in zip(fpr, tpr):
                                best_curves.append({"region": region, "curve": "ROC", "x": a, "y": b})
                            for a, b in zip(recall, precision):
                                best_curves.append({"region": region, "curve": "PR", "x": a, "y": b})

    perf = pd.DataFrame(performance_rows)
    perf.to_csv(OUT_DIR / "module_score_ml_fold_metrics.csv", index=False)
    summary = summarize_metrics(performance_rows)
    summary.to_csv(OUT_DIR / "module_score_ml_performance_summary.csv", index=False)
    pd.concat(coverage_all, ignore_index=True).to_csv(OUT_DIR / "published_module_gene_coverage_by_region.csv", index=False)
    coef = pd.DataFrame(coef_rows)
    coef.to_csv(IMPORTANCE_DIR / "module_score_ml_feature_importance_by_fold.csv", index=False)
    if not coef.empty:
        imp = coef.groupby(["region", "feature_set", "model", "published_module"]).agg(
            mean_importance=("abs_coefficient", "mean"),
            mean_signed_coefficient=("coefficient", "mean"),
            sign_consistency=("coefficient", lambda s: max((s > 0).mean(), (s < 0).mean())),
        ).reset_index().sort_values(["region", "feature_set", "mean_importance"], ascending=[True, True, False])
        imp.to_csv(IMPORTANCE_DIR / "module_score_ml_feature_importance_summary.csv", index=False)

    best = summary.sort_values(["region", "feature_set", "balanced_accuracy_mean", "roc_auc_mean"], ascending=[True, True, False, False]).groupby(["region", "feature_set"]).head(1)
    best.to_csv(OUT_DIR / "module_score_ml_best_model_by_region_feature_set.csv", index=False)

    main = summary[summary["feature_set"].eq("all32_predefined_modules")].copy()
    fig, ax = plt.subplots(figsize=(10, 5))
    plot_df = main.sort_values("balanced_accuracy_mean", ascending=False)
    ax.barh(plot_df["region"] + " | " + plot_df["model"], plot_df["balanced_accuracy_mean"], xerr=plot_df["balanced_accuracy_std"].fillna(0), color="#4C78A8")
    ax.axvline(0.5, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Balanced accuracy mean across CV folds")
    ax.set_title("Published AHBA module-score ML: all 32 predefined modules")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "barplot_balanced_accuracy_all32_modules.png", dpi=220)
    plt.close(fig)

    pivot = main.pivot_table(index="region", columns="model", values="roc_auc_mean")
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    fig.colorbar(im, ax=ax, label="ROC-AUC mean")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=35, ha="right")
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index)
    ax.set_title("ROC-AUC across all 32 AHBA module-score models")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "heatmap_roc_auc_all32_modules.png", dpi=220)
    plt.close(fig)

    curves = pd.DataFrame(best_curves)
    if not curves.empty:
        curves.to_csv(PLOT_DIR / "representative_logistic_L2_curves_all32_modules.csv", index=False)
        for curve_type in ["ROC", "PR"]:
            fig, ax = plt.subplots(figsize=(6, 5))
            for region, sub in curves[curves["curve"].eq(curve_type)].groupby("region"):
                ax.plot(sub["x"], sub["y"], label=region)
            if curve_type == "ROC":
                ax.plot([0, 1], [0, 1], color="black", linestyle="--", linewidth=1)
                ax.set_xlabel("False positive rate")
                ax.set_ylabel("True positive rate")
            else:
                ax.set_xlabel("Recall")
                ax.set_ylabel("Precision")
            ax.set_title(f"Representative {curve_type} curves: logistic L2 all32 modules")
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(PLOT_DIR / f"representative_{curve_type.lower()}_curves_logistic_L2_all32_modules.png", dpi=220)
            plt.close(fig)

    report_lines = [
        "# Published AHBA Module-Score ML Report",
        "",
        "Published AHBA modules were used as external predefined biological feature sets. Module-score ML was exploratory and was used to assess whether AHBA-informed module-level expression patterns carried MDD/Control signal, not to establish a validated diagnostic classifier.",
        "",
        "## Methods",
        "",
        f"- Main feature set: all {len(all_modules)} published AHBA modules.",
        "- Module scores: within each CV fold, each gene was z-scored using training samples only; test samples were transformed using training means and SDs; module score was the mean standardized expression of available genes in the module.",
        "- Models: majority baseline, logistic regression L2/L1/elastic net, linear SVM, and exploratory random forest.",
        f"- CV repeats/seeds: {N_REPEATS}; grouped stratified folds were used when BrainID repeats occurred within a region.",
        "",
        "## Main all32-module performance",
        "",
        main.sort_values(["region", "balanced_accuracy_mean"], ascending=[True, False])[["region", "model", "cv_type", "n_samples", "n_MDD", "n_Control", "balanced_accuracy_mean", "balanced_accuracy_std", "roc_auc_mean", "roc_auc_std", "leakage_status"]].to_markdown(index=False),
        "",
        "## Best model per region/feature set",
        "",
        best[["region", "feature_set", "model", "balanced_accuracy_mean", "roc_auc_mean", "leakage_status"]].to_markdown(index=False),
        "",
        "## Limitations",
        "",
        "This is exploratory module-level ML. All 32-module models use external predefined AHBA modules and avoid feature-selection leakage, but they are still evaluated only within this dataset. Post-enrichment selected module subsets are sensitivity analyses and are more optimistic because the modules were selected after enrichment.",
    ]
    (OUT_DIR / "Module_Score_ML_Report.md").write_text("\n".join(report_lines), encoding="utf-8")
    manifest = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "n_repeats": N_REPEATS,
        "regions": list(REGION_FILES),
        "feature_sets": feature_sets,
        "main_feature_set": "all32_predefined_modules",
        "output_dir": str(OUT_DIR),
        "interpretation": "Exploratory module-level ML, not validated diagnostic prediction.",
    }
    (OUT_DIR / "module_score_ml_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
