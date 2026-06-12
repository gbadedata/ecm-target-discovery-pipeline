"""Clinical metadata cleaning and dictionary creation for CPTAC PDAC.

Reads raw clinical.tsv, standardises column names, creates derived
variables for survival and ECM analysis, and outputs:
  - data/processed/clinical_clean.csv
  - metadata/clinical_dictionary.csv
  - evidence/logs/clinical_cleaning_report.json

Usage:
    python -m src.python.clean_clinical
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
METADATA_DIR = PROJECT_ROOT / "metadata"
EVIDENCE_DIR = PROJECT_ROOT / "evidence" / "logs"

# ---------------------------------------------------------------------------
# Column dictionary: raw name -> (clean name, description, category)
# ---------------------------------------------------------------------------

COLUMN_DICTIONARY = {
    "case_id": (
        "case_id", "CPTAC patient identifier", "identifier"
    ),
    "tumor_included_for_the_study": (
        "tumour_in_study", "Whether tumour sample was included", "study"
    ),
    "normal_included_for_the_study": (
        "normal_in_study", "Whether normal sample was included", "study"
    ),
    "histology_diagnosis": (
        "histology", "Histological diagnosis", "pathology"
    ),
    "age": (
        "age", "Age at diagnosis (years)", "demographics"
    ),
    "sex": (
        "sex", "Biological sex (Male/Female)", "demographics"
    ),
    "race": (
        "race", "Self-reported race (33/140 complete)", "demographics"
    ),
    "participant_country": (
        "country", "Country of participant", "demographics"
    ),
    "tumor_site": (
        "tumour_site", "Anatomical site of tumour", "pathology"
    ),
    "tumor_focality": (
        "tumour_focality", "Tumour focality", "pathology"
    ),
    "tumor_size_cm": (
        "tumour_size_cm", "Tumour size in centimetres", "pathology"
    ),
    "tumor_necrosis": (
        "tumour_necrosis", "Presence of tumour necrosis", "pathology"
    ),
    "lymph_vascular_invasion": (
        "lymphovascular_invasion", "Lymphovascular invasion status", "pathology"
    ),
    "perineural_invasion": (
        "perineural_invasion", "Perineural invasion status", "pathology"
    ),
    "number_of_lymph_nodes_examined": (
        "lymph_nodes_examined", "Number of lymph nodes examined", "pathology"
    ),
    "number_of_lymph_nodes_positive_for_tumor": (
        "lymph_nodes_positive", "Number of positive lymph nodes", "pathology"
    ),
    "pathologic_staging_regional_lymph_nodes_pn": (
        "stage_pn", "Pathological N stage", "staging"
    ),
    "pathologic_staging_primary_tumor_pt": (
        "stage_pt", "Pathological T stage", "staging"
    ),
    "pathologic_staging_distant_metastasis_pm": (
        "stage_pm", "Pathological M stage", "staging"
    ),
    "clinical_staging_distant_metastasis_cm": (
        "stage_cm", "Clinical M stage", "staging"
    ),
    "residual_tumor": (
        "residual_tumour", "Residual tumour classification", "pathology"
    ),
    "tumor_stage_pathological": (
        "stage_overall", "Overall pathological stage", "staging"
    ),
    "additional_pathologic_findings": (
        "additional_pathology", "Additional pathological findings", "pathology"
    ),
    "bmi": (
        "bmi", "Body mass index", "demographics"
    ),
    "alcohol_consumption": (
        "alcohol", "Alcohol consumption history", "demographics"
    ),
    "tobacco_smoking_history": (
        "smoking", "Tobacco smoking history", "demographics"
    ),
    "medical_condition": (
        "medical_conditions", "Reported medical conditions", "demographics"
    ),
    "Neoplastic_cellularity": (
        "cellularity_neoplastic", "Pathologist-estimated neoplastic fraction",
        "tissue_composition"
    ),
    "Acinar_fraction": (
        "fraction_acinar", "Pathologist-estimated acinar fraction",
        "tissue_composition"
    ),
    "Islet_fraction": (
        "fraction_islet", "Pathologist-estimated islet fraction",
        "tissue_composition"
    ),
    "Stromal_fraction": (
        "fraction_stromal", "Pathologist-estimated stromal fraction",
        "tissue_composition"
    ),
    "Non_neoplastic_duct": (
        "fraction_duct_nonneoplastic",
        "Pathologist-estimated non-neoplastic duct fraction",
        "tissue_composition"
    ),
    "Fat_fraction": (
        "fraction_fat", "Pathologist-estimated fat fraction",
        "tissue_composition"
    ),
    "Inflammation_fraction": (
        "fraction_inflammation",
        "Pathologist-estimated inflammation fraction",
        "tissue_composition"
    ),
    "Muscle_fraction": (
        "fraction_muscle", "Pathologist-estimated muscle fraction",
        "tissue_composition"
    ),
    "follow_up_days": (
        "survival_days", "Overall survival time in days", "survival"
    ),
    "vital_status": (
        "vital_status", "Vital status at last follow-up", "survival"
    ),
    "is_this_patient_lost_to_follow_up": (
        "lost_to_followup", "Whether patient was lost to follow-up",
        "survival"
    ),
    "cause_of_death": (
        "cause_of_death", "Cause of death if deceased", "survival"
    ),
}


# ---------------------------------------------------------------------------
# Cleaning functions
# ---------------------------------------------------------------------------


def load_raw_clinical() -> pd.DataFrame:
    """Load the raw clinical TSV file."""
    filepath = RAW_DIR / "clinical.tsv"
    return pd.read_csv(filepath, sep="\t")


def rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to standardised names."""
    rename_map = {
        raw: clean for raw, (clean, _, _) in COLUMN_DICTIONARY.items()
    }
    return df.rename(columns=rename_map)


def clean_survival(df: pd.DataFrame) -> pd.DataFrame:
    """Create binary survival event and numeric survival time."""
    df = df.copy()

    # Binary event: 1 = dead, 0 = alive/censored
    df["survival_event"] = (
        df["vital_status"]
        .str.strip()
        .str.lower()
        .map({"deceased": 1, "dead": 1, "living": 0, "alive": 0})
    )

    # Survival months
    df["survival_months"] = pd.to_numeric(
        df["survival_days"], errors="coerce"
    ) / 30.44

    return df


def _parse_semicolon_values(value) -> float:
    """Parse semicolon-separated pathologist estimates and average them.

    CPTAC PDAC tissue composition fields contain multiple reviewer
    estimates separated by semicolons (e.g., '75;55;53'). This function
    parses them and returns the mean.
    """
    if pd.isna(value):
        return np.nan
    s = str(value).strip()
    if not s:
        return np.nan
    parts = s.split(";")
    numeric_parts = []
    for p in parts:
        p = p.strip()
        if p:
            try:
                numeric_parts.append(float(p))
            except ValueError:
                continue
    if not numeric_parts:
        return np.nan
    return np.mean(numeric_parts)


def clean_tissue_composition(df: pd.DataFrame) -> pd.DataFrame:
    """Parse and average multi-reviewer tissue composition estimates.

    CPTAC PDAC tissue composition fields may contain semicolon-separated
    values from multiple pathologist reviews. Each value is parsed,
    averaged, and converted to a fraction (0-1 scale).
    """
    df = df.copy()
    comp_cols = [c for c in df.columns if c.startswith("fraction_")]
    comp_cols.append("cellularity_neoplastic")

    for col in comp_cols:
        df[col] = df[col].apply(_parse_semicolon_values)
        # Convert percentage to fraction if values are > 1
        if df[col].dropna().max() > 1.0:
            df[col] = df[col] / 100.0

    return df


def create_stromal_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Create stromal-high/low groups based on median split."""
    df = df.copy()
    median_stromal = df["fraction_stromal"].median()

    df["stromal_group"] = np.where(
        df["fraction_stromal"] >= median_stromal,
        "stromal_high",
        "stromal_low",
    )
    df["stromal_group"] = df["stromal_group"].where(
        df["fraction_stromal"].notna(), other=np.nan
    )

    return df


def create_clinical_dictionary() -> pd.DataFrame:
    """Create a clinical data dictionary CSV."""
    rows = []
    for raw_name, (clean_name, description, category) in (
        COLUMN_DICTIONARY.items()
    ):
        rows.append({
            "raw_column": raw_name,
            "clean_column": clean_name,
            "description": description,
            "category": category,
        })

    # Add derived columns
    derived = [
        ("survival_event", "Binary survival event (1=dead, 0=alive)",
         "survival", "Derived from vital_status"),
        ("survival_months", "Overall survival in months",
         "survival", "Derived from follow_up_days / 30.44"),
        ("stromal_group", "Stromal-high or stromal-low (median split)",
         "tissue_composition", "Derived from fraction_stromal"),
    ]
    for clean_name, description, category, note in derived:
        rows.append({
            "raw_column": f"[derived: {note}]",
            "clean_column": clean_name,
            "description": description,
            "category": category,
        })

    return pd.DataFrame(rows)


def generate_report(
    df_raw: pd.DataFrame, df_clean: pd.DataFrame
) -> dict:
    """Generate a structured cleaning report."""
    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "raw_shape": {"rows": len(df_raw), "columns": len(df_raw.columns)},
        "clean_shape": {
            "rows": len(df_clean),
            "columns": len(df_clean.columns),
        },
        "columns_renamed": len(COLUMN_DICTIONARY),
        "derived_columns_added": 3,
    }

    # Survival summary
    surv = df_clean["survival_event"].value_counts().to_dict()
    report["survival"] = {
        "events_dead": int(surv.get(1, 0)),
        "events_alive": int(surv.get(0, 0)),
        "missing": int(df_clean["survival_event"].isna().sum()),
        "median_survival_days": (
            float(df_clean["survival_days"].median())
            if df_clean["survival_days"].notna().any() else None
        ),
    }

    # Tissue composition summary
    comp_cols = [c for c in df_clean.columns if c.startswith("fraction_")]
    comp_summary = {}
    for col in comp_cols:
        vals = df_clean[col].dropna()
        if len(vals) > 0:
            comp_summary[col] = {
                "mean": round(float(vals.mean()), 4),
                "median": round(float(vals.median()), 4),
                "min": round(float(vals.min()), 4),
                "max": round(float(vals.max()), 4),
            }
    report["tissue_composition"] = comp_summary

    # Stromal group counts
    sg = df_clean["stromal_group"].value_counts().to_dict()
    report["stromal_groups"] = {
        str(k): int(v) for k, v in sg.items()
    }

    # Demographics
    report["demographics"] = {
        "age_mean": round(float(df_clean["age"].mean()), 1),
        "age_range": [
            int(df_clean["age"].min()),
            int(df_clean["age"].max()),
        ],
        "sex": {
            str(k): int(v)
            for k, v in df_clean["sex"].value_counts().items()
        },
    }

    # Completeness
    completeness = {}
    for col in df_clean.columns:
        frac = round(float(df_clean[col].notna().mean()), 4)
        completeness[col] = frac
    report["completeness"] = completeness

    return report


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    """Run clinical metadata cleaning pipeline."""
    print("Cleaning clinical metadata...")

    # Load
    df_raw = load_raw_clinical()
    print(f"  Raw: {df_raw.shape[0]} patients, {df_raw.shape[1]} columns")

    # Clean
    df = rename_columns(df_raw)
    df = clean_survival(df)
    df = clean_tissue_composition(df)
    df = create_stromal_groups(df)
    print(f"  Clean: {df.shape[0]} patients, {df.shape[1]} columns")

    # Save cleaned clinical data
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    clean_path = PROCESSED_DIR / "clinical_clean.csv"
    df.to_csv(clean_path, index=False)
    print(f"  Saved: {clean_path}")

    # Save clinical dictionary
    dict_df = create_clinical_dictionary()
    dict_path = METADATA_DIR / "clinical_dictionary.csv"
    dict_df.to_csv(dict_path, index=False)
    print(f"  Dictionary: {dict_path}")

    # Generate and save report
    report = generate_report(df_raw, df)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    report_path = EVIDENCE_DIR / "clinical_cleaning_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"  Report: {report_path}")

    # Print summary
    print("\n--- Clinical Metadata Summary ---")
    print(f"Patients:       {report['clean_shape']['rows']}")
    print(f"Columns:        {report['clean_shape']['columns']}")
    s = report["survival"]
    print(f"Survival:       {s['events_dead']} dead, {s['events_alive']} alive"
          f", {s['missing']} missing")
    if s["median_survival_days"]:
        months = round(s["median_survival_days"] / 30.44, 1)
        print(f"Median OS:      {s['median_survival_days']:.0f} days"
              f" ({months} months)")
    print(f"Stromal groups: {report['stromal_groups']}")
    d = report["demographics"]
    print(f"Age:            mean {d['age_mean']}"
          f", range {d['age_range'][0]}-{d['age_range'][1]}")
    print(f"Sex:            {d['sex']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
