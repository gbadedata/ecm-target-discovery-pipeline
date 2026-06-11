"""Data validation for CPTAC PDAC datasets.

Validates downloaded files against checksums, checks data structure,
assesses completeness, and produces a structured validation report.

Usage:
    python src/python/validate_data.py

Output:
    evidence/logs/data_validation_report.json
    evidence/logs/data_validation_summary.txt
"""

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from src.python.config import load_params

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
EVIDENCE_DIR = PROJECT_ROOT / "evidence" / "logs"
METADATA_DIR = PROJECT_ROOT / "metadata"

EXPECTED_FILES = {
    "rnaseq_tumor.cct": {
        "description": "Bulk RNA-seq, RSEM-UQ log2, tumour samples",
        "expected_rows": 28057,
        "expected_samples_approx": 140,
        "data_type": "transcriptomics",
        "index_type": "gene_symbol",
    },
    "rnaseq_normal.cct": {
        "description": "Bulk RNA-seq, RSEM-UQ log2, normal adjacent tissue",
        "expected_rows": 28057,
        "expected_samples_approx": 21,
        "data_type": "transcriptomics",
        "index_type": "gene_symbol",
    },
    "proteomics_tumor.cct": {
        "description": "TMT proteomics, median-normalised intensity, tumour",
        "expected_rows": 11662,
        "expected_samples_approx": 140,
        "data_type": "proteomics",
        "index_type": "gene_symbol",
    },
    "proteomics_normal.cct": {
        "description": "TMT proteomics, median-normalised intensity, normal",
        "expected_rows": 11662,
        "expected_samples_approx": 75,
        "data_type": "proteomics",
        "index_type": "gene_symbol",
    },
    "clinical.tsv": {
        "description": "Clinical metadata for 140 PDAC patients",
        "expected_rows": 140,
        "expected_samples_approx": 140,
        "data_type": "clinical",
        "index_type": "case_id",
    },
}


# ---------------------------------------------------------------------------
# Validation functions
# ---------------------------------------------------------------------------


def compute_sha256(filepath: Path) -> str:
    """Compute SHA-256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def load_checksums(checksums_path: Path) -> dict[str, str]:
    """Load expected checksums from metadata/checksums.txt."""
    checksums = {}
    if not checksums_path.exists():
        return checksums
    with open(checksums_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) == 2:
                hash_val, filepath = parts
                filename = Path(filepath).name
                checksums[filename] = hash_val
    return checksums


def validate_file_presence(filename: str) -> dict:
    """Check if a file exists and return its size."""
    filepath = RAW_DIR / filename
    exists = filepath.exists()
    size_bytes = filepath.stat().st_size if exists else 0
    size_mb = round(size_bytes / (1024 * 1024), 2) if exists else 0
    return {
        "exists": exists,
        "size_bytes": size_bytes,
        "size_mb": size_mb,
    }


def validate_checksum(filename: str, expected_checksums: dict[str, str]) -> dict:
    """Validate file checksum against expected value."""
    filepath = RAW_DIR / filename
    if not filepath.exists():
        return {"valid": False, "reason": "file_not_found"}

    if filename not in expected_checksums:
        return {"valid": False, "reason": "no_expected_checksum"}

    actual = compute_sha256(filepath)
    expected = expected_checksums[filename]
    return {
        "valid": actual == expected,
        "expected": expected,
        "actual": actual,
    }


def validate_matrix_file(filename: str, spec: dict) -> dict:
    """Validate a .cct matrix file structure and content."""
    filepath = RAW_DIR / filename
    if not filepath.exists():
        return {"valid": False, "reason": "file_not_found"}

    try:
        df = pd.read_csv(filepath, sep="\t", index_col=0)

        n_features = len(df)
        n_samples = len(df.columns)
        missing_fraction = float(df.isna().sum().sum() / df.size)
        missing_per_feature = df.isna().mean(axis=1)
        missing_per_sample = df.isna().mean(axis=0)

        sample_ids = list(df.columns[:5])
        feature_ids = list(df.index[:5])

        # Check numeric
        numeric_cols = df.select_dtypes(include="number").columns
        all_numeric = len(numeric_cols) == n_samples

        # Value range
        if all_numeric and n_features > 0:
            val_min = float(df.min().min())
            val_max = float(df.max().max())
            val_mean = float(df.mean().mean())
        else:
            val_min = val_max = val_mean = None

        row_check = n_features == spec["expected_rows"]
        sample_check = abs(n_samples - spec["expected_samples_approx"]) <= 5

        return {
            "valid": row_check and sample_check and all_numeric,
            "n_features": n_features,
            "n_samples": n_samples,
            "expected_features": spec["expected_rows"],
            "expected_samples_approx": spec["expected_samples_approx"],
            "features_match": row_check,
            "samples_match": sample_check,
            "all_numeric": all_numeric,
            "missing_fraction": round(missing_fraction, 6),
            "features_with_no_missing": int((missing_per_feature == 0).sum()),
            "samples_with_no_missing": int((missing_per_sample == 0).sum()),
            "worst_feature_missing_pct": round(float(missing_per_feature.max()) * 100, 2),
            "worst_sample_missing_pct": round(float(missing_per_sample.max()) * 100, 2),
            "value_min": val_min,
            "value_max": val_max,
            "value_mean": round(val_mean, 4) if val_mean is not None else None,
            "sample_ids_preview": sample_ids,
            "feature_ids_preview": feature_ids,
        }

    except Exception as e:
        return {"valid": False, "reason": str(e)}


def validate_clinical_file(filename: str) -> dict:
    """Validate the clinical metadata file."""
    filepath = RAW_DIR / filename
    if not filepath.exists():
        return {"valid": False, "reason": "file_not_found"}

    try:
        df = pd.read_csv(filepath, sep="\t")
        n_patients = len(df)
        columns = list(df.columns)
        n_columns = len(columns)

        # Key clinical columns
        key_columns = [
            "case_id",
            "age",
            "gender",
            "histology_diagnosis",
            "tumor_stage_pathological",
            "vital_status",
            "overall_survival",
        ]
        present_key_columns = [c for c in key_columns if c in columns]
        missing_key_columns = [c for c in key_columns if c not in columns]

        # Completeness per column
        completeness = {}
        for col in columns:
            non_null = df[col].notna().sum()
            completeness[col] = round(non_null / n_patients, 4)

        # Summary stats for key fields
        field_summaries = {}
        for col in present_key_columns:
            if not pd.api.types.is_numeric_dtype(df[col]):
                field_summaries[col] = {
                    "unique_values": int(df[col].nunique()),
                    "top_values": {
                        str(k): int(v)
                        for k, v in df[col].value_counts().head(5).items()
                    },
                    "missing": int(df[col].isna().sum()),
                }
            else:
                field_summaries[col] = {
                    "min": float(df[col].min()) if df[col].notna().any() else None,
                    "max": float(df[col].max()) if df[col].notna().any() else None,
                    "mean": round(float(df[col].mean()), 2) if df[col].notna().any() else None,
                    "missing": int(df[col].isna().sum()),
                }

        return {
            "valid": n_patients == 140 and "case_id" in columns,
            "n_patients": n_patients,
            "n_columns": n_columns,
            "all_columns": columns,
            "key_columns_present": present_key_columns,
            "key_columns_missing": missing_key_columns,
            "completeness_summary": {
                "fully_complete": sum(1 for v in completeness.values() if v == 1.0),
                "above_90_pct": sum(1 for v in completeness.values() if v >= 0.9),
                "below_50_pct": sum(1 for v in completeness.values() if v < 0.5),
            },
            "field_summaries": field_summaries,
        }

    except Exception as e:
        return {"valid": False, "reason": str(e)}


# ---------------------------------------------------------------------------
# Main validation pipeline
# ---------------------------------------------------------------------------


def run_validation() -> dict:
    """Run full data validation and return structured report."""
    params = load_params()
    checksums = load_checksums(METADATA_DIR / "checksums.txt")

    report = {
        "validation_timestamp": datetime.now(timezone.utc).isoformat(),
        "project": params.get("project_name", "ecm-target-discovery-pipeline"),
        "dataset": params.get("dataset", {}),
        "python_version": sys.version,
        "pandas_version": pd.__version__,
        "files": {},
        "summary": {},
    }

    all_valid = True
    total_files = 0
    passed_files = 0

    for filename, spec in EXPECTED_FILES.items():
        total_files += 1
        file_report = {
            "description": spec["description"],
            "data_type": spec["data_type"],
            "presence": validate_file_presence(filename),
            "checksum": validate_checksum(filename, checksums),
        }

        if spec["data_type"] == "clinical":
            file_report["structure"] = validate_clinical_file(filename)
        else:
            file_report["structure"] = validate_matrix_file(filename, spec)

        # Overall file validity
        file_valid = (
            file_report["presence"]["exists"]
            and file_report["checksum"].get("valid", False)
            and file_report["structure"].get("valid", False)
        )
        file_report["valid"] = file_valid

        if file_valid:
            passed_files += 1
        else:
            all_valid = False

        report["files"][filename] = file_report

    report["summary"] = {
        "total_files": total_files,
        "passed": passed_files,
        "failed": total_files - passed_files,
        "all_valid": all_valid,
        "pass_rate": f"{passed_files}/{total_files}",
    }

    return report


def write_summary(report: dict, output_path: Path) -> None:
    """Write a human-readable validation summary."""
    lines = []
    lines.append("=" * 70)
    lines.append("ECM TARGET DISCOVERY PIPELINE - DATA VALIDATION REPORT")
    lines.append("=" * 70)
    lines.append(f"Timestamp: {report['validation_timestamp']}")
    source = report['dataset'].get('source', '')
    cohort = report['dataset'].get('cohort', '')
    lines.append(f"Dataset:   {source} {cohort}")
    lines.append(f"Result:    {report['summary']['pass_rate']} files passed validation")
    lines.append("")

    for filename, file_report in report["files"].items():
        status = "PASS" if file_report["valid"] else "FAIL"
        lines.append(f"[{status}] {filename}")
        lines.append(f"       {file_report['description']}")
        lines.append(f"       Size: {file_report['presence']['size_mb']} MB")
        cksum = 'valid' if file_report['checksum'].get('valid') else 'INVALID'
        lines.append(f"       Checksum: {cksum}")

        struct = file_report["structure"]
        if file_report["data_type"] == "clinical":
            if "n_patients" in struct:
                n_pat = struct['n_patients']
                n_col = struct['n_columns']
                lines.append(f"       Patients: {n_pat}, Columns: {n_col}")
                present = ', '.join(struct.get('key_columns_present', []))
                lines.append(f"       Key columns present: {present}")
                if struct.get("key_columns_missing"):
                    missing = ', '.join(struct['key_columns_missing'])
                    lines.append(f"       Key columns MISSING: {missing}")
        else:
            if "n_features" in struct:
                n_feat = struct['n_features']
                n_samp = struct['n_samples']
                lines.append(f"       Shape: {n_feat} features x {n_samp} samples")
                lines.append(f"       Missing: {struct['missing_fraction'] * 100:.2f}%")
                if struct.get("value_min") is not None:
                    vmin = struct['value_min']
                    vmax = struct['value_max']
                    lines.append(f"       Value range: [{vmin:.2f}, {vmax:.2f}]")
        lines.append("")

    lines.append("=" * 70)
    if report["summary"]["all_valid"]:
        lines.append("ALL FILES PASSED VALIDATION")
    else:
        lines.append(f"VALIDATION FAILURES: {report['summary']['failed']} file(s) failed")
    lines.append("=" * 70)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def main():
    """Run validation and write reports."""
    print("Running data validation...")
    report = run_validation()

    # Write JSON report (machine-readable evidence)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    json_path = EVIDENCE_DIR / "data_validation_report.json"
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"JSON report: {json_path}")

    # Write human-readable summary
    summary_path = EVIDENCE_DIR / "data_validation_summary.txt"
    write_summary(report, summary_path)
    print(f"Summary:     {summary_path}")

    # Print summary to terminal
    print()
    with open(summary_path) as f:
        print(f.read())

    return 0 if report["summary"]["all_valid"] else 1


if __name__ == "__main__":
    sys.exit(main())
