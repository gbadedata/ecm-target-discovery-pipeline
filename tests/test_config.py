"""Tests for configuration loading."""

from src.python.config import get_param, load_params


def test_load_params():
    """Test that params.yaml loads successfully."""
    params = load_params()
    assert isinstance(params, dict)
    assert "project_name" in params
    assert params["project_name"] == "ecm-target-discovery-pipeline"


def test_load_params_dataset():
    """Test dataset parameters are present."""
    params = load_params()
    assert "dataset" in params
    assert params["dataset"]["source"] == "CPTAC"
    assert params["dataset"]["cohort"] == "PDAC"


def test_get_param_nested():
    """Test nested parameter access."""
    params = load_params()
    padj = get_param(params, "transcriptomics", "padj_threshold")
    assert padj == 0.05


def test_get_param_missing_returns_default():
    """Test that missing keys return the default value."""
    params = load_params()
    result = get_param(params, "nonexistent", "key", default="fallback")
    assert result == "fallback"


def test_get_param_ecm_scoring_method():
    """Test ECM signature scoring method parameter."""
    params = load_params()
    method = get_param(params, "ecm_signature", "scoring_method")
    assert method == "ssGSEA"


def test_reproducibility_seed():
    """Test random seed is set."""
    params = load_params()
    seed = get_param(params, "reproducibility", "random_seed")
    assert seed == 42
