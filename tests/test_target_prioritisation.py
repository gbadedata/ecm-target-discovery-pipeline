"""Tests for target prioritisation module."""

import numpy as np
import pandas as pd

from src.python.target_prioritisation import score_column


def test_score_column_ascending_lower_is_better():
    """Test score_column with ascending=True (lower values score higher)."""
    series = pd.Series([0.001, 0.01, 0.05, 0.5, 1.0])
    scores = score_column(series, ascending=True)
    assert scores.iloc[0] == 1.0  # lowest padj = highest score
    assert scores.iloc[-1] == 0.0  # highest padj = lowest score


def test_score_column_descending_higher_is_better():
    """Test score_column with ascending=False (higher values score higher)."""
    series = pd.Series([0.1, 0.5, 1.0, 2.0, 5.0])
    scores = score_column(series, ascending=False)
    assert scores.iloc[-1] == 1.0  # highest value = highest score
    assert scores.iloc[0] == 0.0  # lowest value = lowest score


def test_score_column_handles_ties():
    """Test score_column handles tied values."""
    series = pd.Series([1.0, 1.0, 2.0, 3.0])
    scores = score_column(series, ascending=True)
    assert scores.iloc[0] == scores.iloc[1]  # ties get same score


def test_score_column_handles_na():
    """Test score_column handles NaN values."""
    series = pd.Series([0.01, 0.05, np.nan, 0.5])
    scores = score_column(series, ascending=True)
    assert not scores.iloc[0:2].isna().any()  # non-NA values scored
    assert scores.iloc[2] == 0.0  # NA gets lowest score (na_option=bottom)


def test_score_column_range():
    """Test score_column output is between 0 and 1."""
    series = pd.Series(np.random.rand(100))
    scores = score_column(series, ascending=True)
    assert scores.min() >= 0.0
    assert scores.max() <= 1.0
