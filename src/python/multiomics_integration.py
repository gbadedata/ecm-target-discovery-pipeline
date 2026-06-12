"""Phase 8: Multi-omics integration (RNA + Protein).

Correlates differential expression (RNA) with differential abundance
(protein) to identify concordant ECM-associated targets supported
by both data layers.

Usage:
    python -m src.python.multiomics_integration

Outputs:
    results/multiomics_integrated.csv
    results/multiomics_concordant_ecm.csv
    figures/rna_protein_correlation.png
    figures/rna_protein_ecm_concordance.png
    evidence/logs/multiomics_report.json
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def main():
    print("=== Phase 8: Multi-Omics Integration ===\n")

    # ---- 1. Load DE and DA results -----------------------------------------

    print("Loading RNA DE results...")
    rna = pd.read_csv(PROJECT_ROOT / "results" / "de_results_all.csv")
    print(f"  RNA genes: {len(rna)}")

    print("Loading proteomics DA results...")
    prot = pd.read_csv(PROJECT_ROOT / "results" / "proteomics_da_all.csv")
    print(f"  Proteins: {len(prot)}")

    # ---- 2. Merge on gene/protein symbol -----------------------------------

    print("\nMerging RNA and protein data...")

    merged = pd.merge(
        rna[["gene", "log2fc", "padj", "significant", "direction", "is_ecm"]],
        prot[["protein", "log2fc", "padj", "significant", "direction"]],
        left_on="gene",
        right_on="protein",
        suffixes=("_rna", "_prot"),
        how="inner",
    )

    print(f"  Matched gene-protein pairs: {len(merged)}")

    # ---- 3. Concordance analysis -------------------------------------------

    print("\nAnalysing concordance...")

    # Both significant, same direction
    merged["concordant"] = (
        merged["significant_rna"]
        & merged["significant_prot"]
        & (merged["direction_rna"] == merged["direction_prot"])
    )

    # Both significant, opposite direction
    merged["discordant"] = (
        merged["significant_rna"]
        & merged["significant_prot"]
        & (merged["direction_rna"] != merged["direction_prot"])
    )

    # RNA only or protein only
    merged["rna_only"] = merged["significant_rna"] & ~merged["significant_prot"]
    merged["prot_only"] = ~merged["significant_rna"] & merged["significant_prot"]

    n_concordant = int(merged["concordant"].sum())
    n_discordant = int(merged["discordant"].sum())
    n_rna_only = int(merged["rna_only"].sum())
    n_prot_only = int(merged["prot_only"].sum())

    print(f"  Concordant (both sig, same direction): {n_concordant}")
    print(f"  Discordant (both sig, opposite dir):   {n_discordant}")
    print(f"  RNA-only significant:                  {n_rna_only}")
    print(f"  Protein-only significant:              {n_prot_only}")

    # Overall RNA-protein log2FC correlation
    valid = merged.dropna(subset=["log2fc_rna", "log2fc_prot"])
    rho, p_val = stats.spearmanr(valid["log2fc_rna"], valid["log2fc_prot"])
    print(f"\n  Spearman rho (all matched): {rho:.4f}")
    print(f"  P-value: {p_val:.2e}")

    # ---- 4. ECM concordance ------------------------------------------------

    print("\nECM-specific concordance...")

    ecm_merged = merged[merged["is_ecm"]].copy()
    n_ecm_matched = len(ecm_merged)
    n_ecm_concordant = int(ecm_merged["concordant"].sum())
    n_ecm_up = int(
        ecm_merged["concordant"]
        .multiply(ecm_merged["direction_rna"] == "Up in Tumour")
        .sum()
    )
    n_ecm_down = int(
        ecm_merged["concordant"]
        .multiply(ecm_merged["direction_rna"] == "Down in Tumour")
        .sum()
    )

    print(f"  ECM genes matched in both: {n_ecm_matched}")
    print(f"  ECM concordant: {n_ecm_concordant} ({n_ecm_up} up, {n_ecm_down} down)")

    # ECM correlation
    ecm_valid = ecm_merged.dropna(subset=["log2fc_rna", "log2fc_prot"])
    if len(ecm_valid) >= 10:
        ecm_rho, ecm_p = stats.spearmanr(
            ecm_valid["log2fc_rna"], ecm_valid["log2fc_prot"]
        )
        print(f"  ECM Spearman rho: {ecm_rho:.4f}, p={ecm_p:.2e}")
    else:
        ecm_rho, ecm_p = float("nan"), float("nan")

    # Top concordant ECM genes
    concordant_ecm = (
        ecm_merged[ecm_merged["concordant"]]
        .sort_values("padj_rna")
        .head(20)
    )

    print("\n--- Top 20 Concordant ECM Targets ---")
    display_cols = ["gene", "log2fc_rna", "log2fc_prot",
                    "padj_rna", "padj_prot", "direction_rna"]
    if len(concordant_ecm) > 0:
        display = concordant_ecm[display_cols].copy()
        display["padj_rna"] = display["padj_rna"].map("{:.2e}".format)
        display["padj_prot"] = display["padj_prot"].map("{:.2e}".format)
        display["log2fc_rna"] = display["log2fc_rna"].round(3)
        display["log2fc_prot"] = display["log2fc_prot"].round(3)
        print(display.to_string(index=False))

    # ---- 5. Scatter plot: RNA vs Protein log2FC ----------------------------

    print("\nGenerating figures...")

    fig, ax = plt.subplots(figsize=(8, 7))

    # Background: non-ECM
    non_ecm = merged[~merged["is_ecm"]]
    ax.scatter(
        non_ecm["log2fc_rna"], non_ecm["log2fc_prot"],
        c="grey", s=5, alpha=0.3, label="Non-ECM"
    )

    # ECM non-concordant
    ecm_nc = ecm_merged[~ecm_merged["concordant"]]
    ax.scatter(
        ecm_nc["log2fc_rna"], ecm_nc["log2fc_prot"],
        c="grey", s=15, alpha=0.5, edgecolors="black",
        linewidths=0.3, label="ECM (not concordant)"
    )

    # ECM concordant up
    ecm_up_df = ecm_merged[
        ecm_merged["concordant"]
        & (ecm_merged["direction_rna"] == "Up in Tumour")
    ]
    ax.scatter(
        ecm_up_df["log2fc_rna"], ecm_up_df["log2fc_prot"],
        c="#E74C3C", s=30, alpha=0.8, edgecolors="black",
        linewidths=0.3, label=f"ECM concordant up ({len(ecm_up_df)})"
    )

    # ECM concordant down
    ecm_down_df = ecm_merged[
        ecm_merged["concordant"]
        & (ecm_merged["direction_rna"] == "Down in Tumour")
    ]
    ax.scatter(
        ecm_down_df["log2fc_rna"], ecm_down_df["log2fc_prot"],
        c="#3498DB", s=30, alpha=0.8, edgecolors="black",
        linewidths=0.3, label=f"ECM concordant down ({len(ecm_down_df)})"
    )

    # Label top concordant ECM
    for _, row in concordant_ecm.head(10).iterrows():
        ax.annotate(
            row["gene"],
            (row["log2fc_rna"], row["log2fc_prot"]),
            fontsize=7, alpha=0.8,
            xytext=(5, 5), textcoords="offset points",
        )

    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")

    # Diagonal reference
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, "k--", alpha=0.2, linewidth=0.5)

    ax.set_xlabel("RNA log2(Fold Change)", fontsize=12)
    ax.set_ylabel("Protein log2(Fold Change)", fontsize=12)
    ax.set_title(
        f"RNA-Protein Concordance in PDAC\n"
        f"Spearman rho={rho:.3f}, p={p_val:.2e} (n={len(valid)})",
        fontsize=13,
    )
    ax.legend(fontsize=8, loc="upper left")
    plt.tight_layout()
    plt.savefig(
        PROJECT_ROOT / "figures" / "rna_protein_correlation.png", dpi=300
    )
    plt.close()
    print("  Saved: figures/rna_protein_correlation.png")

    # ---- 6. Save results ---------------------------------------------------

    print("\nSaving results...")

    results_dir = PROJECT_ROOT / "results"
    merged.to_csv(results_dir / "multiomics_integrated.csv", index=False)
    print(f"  Saved: results/multiomics_integrated.csv ({len(merged)} rows)")

    concordant_ecm_full = ecm_merged[ecm_merged["concordant"]].sort_values(
        "padj_rna"
    )
    concordant_ecm_full.to_csv(
        results_dir / "multiomics_concordant_ecm.csv", index=False
    )
    print(
        f"  Saved: results/multiomics_concordant_ecm.csv "
        f"({len(concordant_ecm_full)} rows)"
    )

    # ---- 7. Evidence report ------------------------------------------------

    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "analysis": "multiomics_integration",
        "method": "RNA_protein_log2FC_concordance",
        "matched_genes": len(merged),
        "correlation": {
            "spearman_rho": round(float(rho), 4),
            "p_value": float(p_val),
            "n_pairs": len(valid),
        },
        "concordance": {
            "concordant": n_concordant,
            "discordant": n_discordant,
            "rna_only": n_rna_only,
            "protein_only": n_prot_only,
        },
        "ecm_specific": {
            "ecm_matched": n_ecm_matched,
            "ecm_concordant": n_ecm_concordant,
            "ecm_concordant_up": n_ecm_up,
            "ecm_concordant_down": n_ecm_down,
            "ecm_spearman_rho": (
                round(float(ecm_rho), 4)
                if not np.isnan(ecm_rho) else None
            ),
            "top_concordant_ecm": concordant_ecm["gene"].tolist(),
        },
    }

    evidence_dir = PROJECT_ROOT / "evidence" / "logs"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    with open(evidence_dir / "multiomics_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)
    print("  Evidence: evidence/logs/multiomics_report.json")

    print("\n=== Phase 8 Complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
