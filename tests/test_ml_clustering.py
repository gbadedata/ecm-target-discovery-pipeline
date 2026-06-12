"""Tests for ML clustering module."""

import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def test_pca_reduces_dimensions():
    """Test PCA reduces a feature matrix correctly."""
    np.random.seed(42)
    X = np.random.randn(50, 100)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    pca = PCA(n_components=10)
    X_pca = pca.fit_transform(X_scaled)
    assert X_pca.shape == (50, 10)
    assert pca.explained_variance_ratio_.sum() > 0


def test_pca_variance_sums_to_one():
    """Test that all PCA components sum to ~1.0 total variance."""
    np.random.seed(42)
    X = np.random.randn(50, 20)
    pca = PCA(n_components=20)
    pca.fit(X)
    assert abs(pca.explained_variance_ratio_.sum() - 1.0) < 0.01
