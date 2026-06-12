"""Phase 11: Candidate ECM target prioritisation.

Integrates evidence from all previous phases into a transparent,
multi-layer scoring framework to rank candidate ECM-associated targets.

Usage:
    python -m src.python.target_prioritisation

Outputs:
    results/target_prioritisation_table.csv
    results/target_shortlist_top30.csv
    figures/target_prioritisation_heatmap.png
    evidence/logs/prioritisation_report.json
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use("Agg")

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def score_column(series: pd.Series, ascending: bool = True) -> pd.Series:
    """Rank-based score normalised to 0-1 (1 = strongest evidence)."""
    ranked = series.rank(ascending=ascending, method="average", na_option="bottom")
    return (ranked - ranked.min()) / (ranked.max() - ranked.min())


def main():
    print("=== Phase 11: ECM Target Prioritisation ===\n")

    # ---- 1. Load all evidence layers ---------------------------------------

    print("Loading evidence layers...")

    # RNA DE
    rna_de = pd.read_csv(PROJECT_ROOT / "results" / "de_results_all.csv")
    rna_de = rna_de[rna_de["is_ecm"]].copy()
    rna_de = rna_de.rename(columns={
        "log2fc": "rna_log2fc",
        "padj": "rna_padj",
        "significant": "rna_significant",
    })
    print(f"  RNA DE (ECM genes): {len(rna_de)}")

    # Proteomics DA
    prot_da = pd.read_csv(PROJECT_ROOT / "results" / "proteomics_da_all.csv")
    prot_ecm = prot_da[prot_da["is_ecm"]].copy()
    prot_ecm = prot_ecm.rename(columns={
        "log2fc": "prot_log2fc",
        "padj": "prot_padj",
        "significant": "prot_significant",
    })
    print(f"  Proteomics DA (ECM proteins): {len(prot_ecm)}")

    # Multi-omics concordance
    multiomics = pd.read_csv(
        PROJECT_ROOT / "results" / "multiomics_integrated.csv"
    )
    concordant = multiomics[multiomics["concordant"]].copy()
    concordant_genes = set(concordant["gene"])
    print(f"  Multi-omics concordant: {len(concordant_genes)}")

    # Feature importance
    feat_imp = pd.read_csv(
        PROJECT_ROOT / "results" / "ml_feature_importance.csv"
    )
    print(f"  ML features: {len(feat_imp)}")

    # ---- 2. Build integrated evidence table --------------------------------

    print("\nBuilding integrated evidence table...")

    # Start with RNA DE ECM genes
    targets = rna_de[["gene", "rna_log2fc", "rna_padj",
                       "rna_significant"]].copy()

    # Merge proteomics
    targets = targets.merge(
        prot_ecm[["protein", "prot_log2fc", "prot_padj", "prot_significant"]],
        left_on="gene", right_on="protein", how="left",
    ).drop(columns=["protein"], errors="ignore")

    # Multi-omics concordance flag
    targets["multiomics_concordant"] = targets["gene"].isin(concordant_genes)

    # ML feature importance
    targets = targets.merge(
        feat_imp[["gene", "rf_importance", "lr_coefficient"]],
        on="gene", how="left",
    )

    print(f"  Integrated table: {len(targets)} ECM genes")

    # ---- 3. Multi-evidence scoring -----------------------------------------

    print("\nComputing prioritisation scores...")

    # Score 1: RNA significance (lower padj = better)
    targets["score_rna_significance"] = score_column(
        targets["rna_padj"], ascending=True
    )

    # Score 2: RNA effect size (higher |log2FC| = better)
    targets["score_rna_effect"] = score_column(
        targets["rna_log2fc"].abs(), ascending=False
    )

    # Score 3: Protein significance (lower padj = better)
    targets["score_prot_significance"] = score_column(
        targets["prot_padj"], ascending=True
    )

    # Score 4: Protein effect size
    targets["score_prot_effect"] = score_column(
        targets["prot_log2fc"].abs(), ascending=False
    )

    # Score 5: Multi-omics concordance (binary: 1 if concordant)
    targets["score_concordance"] = targets["multiomics_concordant"].astype(float)

    # Score 6: ML importance
    targets["score_ml_importance"] = score_column(
        targets["rf_importance"].fillna(0), ascending=False
    )

    # Composite score: weighted average
    weights = {
        "score_rna_significance": 0.20,
        "score_rna_effect": 0.10,
        "score_prot_significance": 0.20,
        "score_prot_effect": 0.10,
        "score_concordance": 0.25,
        "score_ml_importance": 0.15,
    }

    targets["composite_score"] = sum(
        targets[col] * weight for col, weight in weights.items()
    )

    # Rank
    targets = targets.sort_values("composite_score", ascending=False)
    targets["rank"] = range(1, len(targets) + 1)

    # Evidence tier
    targets["evidence_tier"] = "Tier 3: RNA only"
    targets.loc[
        targets["rna_significant"] & targets["prot_significant"].fillna(False),
        "evidence_tier"
    ] = "Tier 2: RNA + Protein"
    targets.loc[
        targets["multiomics_concordant"],
        "evidence_tier"
    ] = "Tier 1: Multi-omics concordant"

    tier_counts = targets["evidence_tier"].value_counts()
    print("\n  Evidence tiers:")
    for tier, count in tier_counts.items():
        print(f"    {tier}: {count}")

    # ---- 4. Top 30 shortlist -----------------------------------------------

    top30 = targets.head(30).copy()

    print("\n--- Top 30 Candidate ECM Targets ---")
    display = top30[[
        "rank", "gene", "rna_log2fc", "prot_log2fc",
        "composite_score", "evidence_tier",
    ]].copy()
    display["rna_log2fc"] = display["rna_log2fc"].round(2)
    display["prot_log2fc"] = display["prot_log2fc"].round(2)
    display["composite_score"] = display["composite_score"].round(3)
    print(display.to_string(index=False))

    # ---- 5. Figures ---------------------------------------------------------

    print("\nGenerating figures...")

    # Prioritisation score heatmap for top 30
    score_cols = [c for c in targets.columns if c.startswith("score_")]
    heatmap_data = top30[score_cols].copy()
    heatmap_data.index = top30["gene"]
    heatmap_data.columns = [c.replace("score_", "") for c in score_cols]

    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(heatmap_data.values, cmap="YlOrRd", aspect="auto")

    ax.set_xticks(range(len(heatmap_data.columns)))
    ax.set_xticklabels(heatmap_data.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(heatmap_data.index)))
    ax.set_yticklabels(heatmap_data.index, fontsize=8)

    plt.colorbar(im, ax=ax, label="Evidence Score (0-1)")
    ax.set_title("Top 30 ECM Target Candidates: Multi-Evidence Scores",
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(
        PROJECT_ROOT / "figures" / "target_prioritisation_heatmap.png",
        dpi=300,
    )
    plt.close()
    print("  Saved: figures/target_prioritisation_heatmap.png")

    # ---- 6. Save results ---------------------------------------------------

    print("\nSaving results...")

    targets.to_csv(
        PROJECT_ROOT / "results" / "target_prioritisation_table.csv",
        index=False,
    )
    print(f"  Saved: results/target_prioritisation_table.csv"
          f" ({len(targets)} rows)")

    top30.to_csv(
        PROJECT_ROOT / "results" / "target_shortlist_top30.csv",
        index=False,
    )
    print(f"  Saved: results/target_shortlist_top30.csv ({len(top30)} rows)")

    # ---- 7. Evidence report ------------------------------------------------

    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "analysis": "target_prioritisation",
        "method": "multi_evidence_weighted_scoring",
        "weights": weights,
        "total_ecm_genes_scored": len(targets),
        "evidence_tiers": {
            str(k): int(v) for k, v in tier_counts.items()
        },
        "top_10_targets": top30.head(10)["gene"].tolist(),
        "top_10_scores": top30.head(10)["composite_score"].round(3).tolist(),
    }

    evidence_dir = PROJECT_ROOT / "evidence" / "logs"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    with open(evidence_dir / "prioritisation_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)
    print("  Evidence: evidence/logs/prioritisation_report.json")

    print("\n=== Phase 11 Complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
