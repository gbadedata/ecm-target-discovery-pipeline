"""Tests for data validation module."""

import hashlib
from pathlib import Path

import pandas as pd
import pytest

import src.python.validate_data as vd
from src.python.validate_data import (
    compute_sha256,
    load_checksums,
)


@pytest.fixture
def tmp_raw_dir(tmp_path):
    """Create a temporary raw data directory with test files."""
    raw_dir = tmp_path / "data" / "raw"
    raw_dir.mkdir(parents=True)
    return raw_dir


@pytest.fixture
def sample_matrix(tmp_raw_dir):
    """Create a small sample matrix file."""
    df = pd.DataFrame(
        {"SAMPLE_01": [1.5, 2.3, 0.0], "SAMPLE_02": [3.1, 4.2, 1.1]},
        index=["GENE_A", "GENE_B", "GENE_C"],
    )
    filepath = tmp_raw_dir / "test_matrix.cct"
    df.to_csv(filepath, sep="\t")
    return filepath


@pytest.fixture
def sample_clinical(tmp_raw_dir):
    """Create a small sample clinical file."""
    df = pd.DataFrame(
        {
            "case_id": ["C3L-00001", "C3L-00002"],
            "age": [65, 72],
            "gender": ["Female", "Male"],
            "histology_diagnosis": ["PDAC", "PDAC"],
            "vital_status": ["Dead", "Alive"],
        }
    )
    filepath = tmp_raw_dir / "test_clinical.tsv"
    df.to_csv(filepath, sep="\t", index=False)
    return filepath


def test_compute_sha256(tmp_path):
    """Test SHA-256 computation."""
    test_file = tmp_path / "test.txt"
    test_file.write_text("hello world")
    result = compute_sha256(test_file)
    expected = hashlib.sha256(b"hello world").hexdigest()
    assert result == expected


def test_load_checksums(tmp_path):
    """Test checksum file loading."""
    checksums_file = tmp_path / "checksums.txt"
    checksums_file.write_text(
        "abc123  data/raw/file1.cct\ndef456  data/raw/file2.tsv\n"
    )
    result = load_checksums(checksums_file)
    assert result["file1.cct"] == "abc123"
    assert result["file2.tsv"] == "def456"


def test_load_checksums_missing_file(tmp_path):
    """Test that missing checksums file returns empty dict."""
    result = load_checksums(tmp_path / "nonexistent.txt")
    assert result == {}


def test_validate_matrix_file_structure(sample_matrix, monkeypatch):
    """Test matrix file validation with a small test file."""
    monkeypatch.setattr(vd, "RAW_DIR", sample_matrix.parent)

    spec = {
        "expected_rows": 3,
        "expected_samples_approx": 2,
        "data_type": "transcriptomics",
        "index_type": "gene_symbol",
    }
    result = vd.validate_matrix_file(sample_matrix.name, spec)
    assert result["valid"] is True
    assert result["n_features"] == 3
    assert result["n_samples"] == 2
    assert result["all_numeric"] is True
    assert result["missing_fraction"] == 0.0


def test_validate_matrix_missing_file(monkeypatch):
    """Test validation of a non-existent file."""
    monkeypatch.setattr(vd, "RAW_DIR", Path("/nonexistent"))
    spec = {"expected_rows": 10, "expected_samples_approx": 5}
    result = vd.validate_matrix_file("missing.cct", spec)
    assert result["valid"] is False


def test_validate_clinical_file_structure(sample_clinical, monkeypatch):
    """Test clinical file validation."""
    monkeypatch.setattr(vd, "RAW_DIR", sample_clinical.parent)

    result = vd.validate_clinical_file(sample_clinical.name)
    assert result["valid"] is False  # n_patients != 140
    assert result["n_patients"] == 2
    assert "case_id" in result["key_columns_present"]
    assert "age" in result["key_columns_present"]
