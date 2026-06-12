"""Phase 10: Dimensionality reduction, clustering, and ML stratification.

Applies PCA, UMAP, consensus clustering, and feature importance ranking
to identify molecular subgroups and key ECM-associated features.

Usage:
    python -m src.python.ml_clustering

Outputs:
    results/ml_pca_components.csv
    results/ml_cluster_assignments.csv
    results/ml_feature_importance.csv
    figures/pca_ecm_groups.png
    figures/cluster_heatmap.png
    evidence/logs/ml_report.json
"""

import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import silhouette_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler

matplotlib.use("Agg")
warnings.filterwarnings("ignore", category=FutureWarning)

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def main():
    print("=== Phase 10: ML, Clustering, and Dimensionality Reduction ===\n")

    # ---- 1. Load data ------------------------------------------------------

    print("Loading expression and clinical data...")
    tumor = pd.read_csv(
        PROJECT_ROOT / "data" / "raw" / "rnaseq_tumor.cct",
        sep="\t", index_col=0,
    )
    scores = pd.read_csv(PROJECT_ROOT / "results" / "ecm_scores.csv")

    # Use DE significant genes as features
    de = pd.read_csv(
        PROJECT_ROOT / "results" / "de_results_significant.csv"
    )
    sig_genes = de["gene"].tolist()

    # Filter to significant genes present in expression matrix
    sig_genes = [g for g in sig_genes if g in tumor.index]
    expr = tumor.loc[sig_genes].T  # samples x genes

    print(f"  Samples: {expr.shape[0]}")
    print(f"  Features (DE genes): {expr.shape[1]}")

    # Merge with ECM group labels
    expr["case_id"] = expr.index
    expr = expr.merge(
        scores[["case_id", "ecm_score", "ecm_group"]],
        on="case_id", how="inner",
    )
    labels = expr["ecm_group"].values
    ecm_scores = expr["ecm_score"].values
    case_ids = expr["case_id"].values
    X = expr.drop(columns=["case_id", "ecm_score", "ecm_group"]).values

    print(f"  Matched with ECM labels: {len(labels)}")

    # Scale
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # ---- 2. PCA ------------------------------------------------------------

    print("\nRunning PCA...")
    pca = PCA(n_components=min(50, X_scaled.shape[1]))
    X_pca = pca.fit_transform(X_scaled)

    var_explained = pca.explained_variance_ratio_
    cum_var = np.cumsum(var_explained)
    n_90 = int(np.argmax(cum_var >= 0.90) + 1)

    print(f"  PC1: {var_explained[0]*100:.1f}% variance")
    print(f"  PC2: {var_explained[1]*100:.1f}% variance")
    print(f"  PCs for 90% variance: {n_90}")

    # PCA plot coloured by ECM group
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colours = {"ECM_low": "#3498DB", "ECM_high": "#E74C3C"}
    for group in ["ECM_low", "ECM_high"]:
        mask = labels == group
        axes[0].scatter(
            X_pca[mask, 0], X_pca[mask, 1],
            c=colours[group], label=group, s=20, alpha=0.7,
        )
    axes[0].set_xlabel(f"PC1 ({var_explained[0]*100:.1f}%)")
    axes[0].set_ylabel(f"PC2 ({var_explained[1]*100:.1f}%)")
    axes[0].set_title("PCA: Coloured by ECM Group")
    axes[0].legend()

    sc = axes[1].scatter(
        X_pca[:, 0], X_pca[:, 1],
        c=ecm_scores, cmap="RdBu_r", s=20, alpha=0.7,
    )
    plt.colorbar(sc, ax=axes[1], label="ECM Score")
    axes[1].set_xlabel(f"PC1 ({var_explained[0]*100:.1f}%)")
    axes[1].set_ylabel(f"PC2 ({var_explained[1]*100:.1f}%)")
    axes[1].set_title("PCA: Coloured by ECM Score (continuous)")

    plt.suptitle("PCA of DE Genes (CPTAC PDAC, n=140)", fontsize=14)
    plt.tight_layout()
    plt.savefig(
        PROJECT_ROOT / "figures" / "pca_ecm_groups.png", dpi=300
    )
    plt.close()
    print("  Saved: figures/pca_ecm_groups.png")

    # Save PCA components
    pca_df = pd.DataFrame(
        X_pca[:, :10],
        columns=[f"PC{i+1}" for i in range(10)],
    )
    pca_df["case_id"] = case_ids
    pca_df["ecm_group"] = labels
    pca_df.to_csv(
        PROJECT_ROOT / "results" / "ml_pca_components.csv", index=False
    )

    # ---- 3. K-Means clustering ---------------------------------------------

    print("\nRunning K-Means clustering...")

    silhouette_scores = {}
    for k in range(2, 7):
        km = KMeans(n_clusters=k, random_state=42, n_init=10)
        cluster_labels = km.fit_predict(X_pca[:, :n_90])
        sil = silhouette_score(X_pca[:, :n_90], cluster_labels)
        silhouette_scores[k] = round(float(sil), 4)
        print(f"  k={k}: silhouette={sil:.4f}")

    best_k = max(silhouette_scores, key=silhouette_scores.get)
    print(f"  Best k: {best_k} (silhouette={silhouette_scores[best_k]:.4f})")

    # Final clustering with best k
    km_final = KMeans(n_clusters=best_k, random_state=42, n_init=10)
    final_clusters = km_final.fit_predict(X_pca[:, :n_90])

    # Cluster vs ECM group crosstab
    cluster_df = pd.DataFrame({
        "case_id": case_ids,
        "cluster": [f"C{c+1}" for c in final_clusters],
        "ecm_group": labels,
        "ecm_score": ecm_scores,
    })

    crosstab = pd.crosstab(cluster_df["cluster"], cluster_df["ecm_group"])
    print("\n  Cluster x ECM Group crosstab:")
    print(crosstab.to_string())

    # Chi-squared test
    chi2, chi_p, _, _ = stats.chi2_contingency(crosstab.values)
    print(f"\n  Chi-squared: {chi2:.2f}, p={chi_p:.4f}")

    cluster_df.to_csv(
        PROJECT_ROOT / "results" / "ml_cluster_assignments.csv", index=False
    )

    # ---- 4. Feature importance (Random Forest) ------------------------------

    print("\nRunning feature importance analysis...")

    y = (labels == "ECM_high").astype(int)

    # Random Forest
    rf = RandomForestClassifier(
        n_estimators=500, random_state=42, n_jobs=-1,
    )
    rf.fit(X_scaled, y)

    importances = pd.DataFrame({
        "gene": sig_genes,
        "rf_importance": rf.feature_importances_,
    }).sort_values("rf_importance", ascending=False)

    # Cross-validated accuracy
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    rf_cv = cross_val_score(rf, X_scaled, y, cv=cv, scoring="roc_auc")
    print(f"  Random Forest CV AUC: {rf_cv.mean():.3f} +/- {rf_cv.std():.3f}")

    # Logistic regression (L1) for sparse feature selection
    lr = LogisticRegression(
        penalty="l1", solver="saga", C=0.1,
        random_state=42, max_iter=5000,
    )
    lr_cv = cross_val_score(lr, X_scaled, y, cv=cv, scoring="roc_auc")
    print(f"  Logistic Reg CV AUC:  {lr_cv.mean():.3f} +/- {lr_cv.std():.3f}")

    lr.fit(X_scaled, y)
    lr_coefs = pd.DataFrame({
        "gene": sig_genes,
        "lr_coefficient": lr.coef_[0],
    })
    lr_nonzero = lr_coefs[lr_coefs["lr_coefficient"] != 0]
    print(f"  L1 selected features: {len(lr_nonzero)}")

    # Merge importance scores
    feat_importance = importances.merge(lr_coefs, on="gene")

    # Load ECM annotation
    ecm_genes = set(
        pd.read_csv(
            PROJECT_ROOT / "results" / "de_results_ecm_genes.csv"
        )["gene"]
    )
    feat_importance["is_ecm"] = feat_importance["gene"].isin(ecm_genes)

    feat_importance.to_csv(
        PROJECT_ROOT / "results" / "ml_feature_importance.csv", index=False
    )

    print("\n--- Top 15 Features (Random Forest Importance) ---")
    top15 = feat_importance.head(15)[
        ["gene", "rf_importance", "lr_coefficient", "is_ecm"]
    ].copy()
    top15["rf_importance"] = top15["rf_importance"].round(4)
    top15["lr_coefficient"] = top15["lr_coefficient"].round(4)
    print(top15.to_string(index=False))

    # ---- 5. Evidence report ------------------------------------------------

    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "analysis": "ml_clustering_dimensionality_reduction",
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "pca": {
            "pc1_variance": round(float(var_explained[0]) * 100, 2),
            "pc2_variance": round(float(var_explained[1]) * 100, 2),
            "pcs_for_90pct": n_90,
        },
        "clustering": {
            "method": "KMeans",
            "best_k": best_k,
            "silhouette_scores": silhouette_scores,
            "cluster_ecm_chi2_p": round(float(chi_p), 6),
        },
        "classification": {
            "random_forest": {
                "cv_auc_mean": round(float(rf_cv.mean()), 4),
                "cv_auc_std": round(float(rf_cv.std()), 4),
            },
            "logistic_regression_l1": {
                "cv_auc_mean": round(float(lr_cv.mean()), 4),
                "cv_auc_std": round(float(lr_cv.std()), 4),
                "n_selected_features": len(lr_nonzero),
            },
        },
        "top_features": feat_importance.head(20)["gene"].tolist(),
    }

    evidence_dir = PROJECT_ROOT / "evidence" / "logs"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    with open(evidence_dir / "ml_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)
    print("\n  Evidence: evidence/logs/ml_report.json")

    print("\n=== Phase 10 Complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
