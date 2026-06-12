"""Tests for clinical metadata cleaning module."""

import pandas as pd
import pytest

import src.python.clean_clinical as cc


@pytest.fixture
def sample_raw_df():
    """Create a sample raw clinical DataFrame."""
    return pd.DataFrame({
        "case_id": ["C3L-001", "C3L-002", "C3L-003", "C3L-004"],
        "sex": ["Male", "Female", "Male", "Female"],
        "age": [65, 72, 58, 80],
        "histology_diagnosis": ["PDAC", "PDAC", "PDAC", "PDAC"],
        "tumor_stage_pathological": ["IIB", "IIA", "IIB", "III"],
        "Stromal_fraction": ["30", "70;60", "50", "10"],
        "Neoplastic_cellularity": ["60", "20;25", "40", "80"],
        "follow_up_days": [365, 180, 730, None],
        "vital_status": ["Deceased", "Living", "Deceased", None],
    })


def test_rename_columns(sample_raw_df):
    """Test that columns are renamed correctly."""
    result = cc.rename_columns(sample_raw_df)
    assert "sex" in result.columns
    assert "stage_overall" in result.columns
    assert "fraction_stromal" in result.columns
    assert "survival_days" in result.columns
    assert "tumor_stage_pathological" not in result.columns


def test_clean_survival(sample_raw_df):
    """Test survival variable derivation."""
    df = cc.rename_columns(sample_raw_df)
    result = cc.clean_survival(df)

    assert "survival_event" in result.columns
    assert "survival_months" in result.columns
    assert result.loc[0, "survival_event"] == 1  # Deceased
    assert result.loc[1, "survival_event"] == 0  # Living
    assert pd.isna(result.loc[3, "survival_event"])  # None

    expected_months = 365 / 30.44
    assert abs(result.loc[0, "survival_months"] - expected_months) < 0.1


def test_create_stromal_groups(sample_raw_df):
    """Test stromal high/low group creation."""
    df = cc.rename_columns(sample_raw_df)
    df = cc.clean_tissue_composition(df)
    result = cc.create_stromal_groups(df)

    assert "stromal_group" in result.columns
    # Raw values: 30, 70;60 (avg 65), 50, 10 -> as fractions: 0.3, 0.65, 0.5, 0.1
    # Median = 0.4, so >= 0.4: 0.65 and 0.5 are high; < 0.4: 0.3 and 0.1 are low
    groups = result["stromal_group"].value_counts()
    assert groups["stromal_high"] == 2
    assert groups["stromal_low"] == 2


def test_clean_tissue_composition_numeric(sample_raw_df):
    """Test tissue composition columns are numeric."""
    df = cc.rename_columns(sample_raw_df)
    result = cc.clean_tissue_composition(df)
    assert pd.api.types.is_numeric_dtype(result["fraction_stromal"])
    assert pd.api.types.is_numeric_dtype(result["cellularity_neoplastic"])


def test_create_clinical_dictionary():
    """Test dictionary creation."""
    dict_df = cc.create_clinical_dictionary()
    assert "raw_column" in dict_df.columns
    assert "clean_column" in dict_df.columns
    assert "description" in dict_df.columns
    assert "category" in dict_df.columns
    # Should include original + derived columns
    assert len(dict_df) >= len(cc.COLUMN_DICTIONARY)
    # Check derived columns exist
    clean_names = dict_df["clean_column"].tolist()
    assert "survival_event" in clean_names
    assert "survival_months" in clean_names
    assert "stromal_group" in clean_names


def test_parse_semicolon_values():
    """Test parsing of semicolon-separated pathologist estimates."""
    assert cc._parse_semicolon_values("75;55;53") == pytest.approx(61.0, abs=0.1)
    assert cc._parse_semicolon_values("30") == 30.0
    assert cc._parse_semicolon_values("50;50") == 50.0
    assert pd.isna(cc._parse_semicolon_values(None))
    assert pd.isna(cc._parse_semicolon_values(""))
