"""Load pipeline parameters from params.yaml."""

from pathlib import Path

import yaml


def load_params(params_path: str | Path | None = None) -> dict:
    """Load parameters from params.yaml.

    Args:
        params_path: Path to params.yaml. Defaults to project root.

    Returns:
        Dictionary of pipeline parameters.
    """
    if params_path is None:
        params_path = Path(__file__).resolve().parents[2] / "params.yaml"
    else:
        params_path = Path(params_path)

    if not params_path.exists():
        raise FileNotFoundError(f"Parameters file not found: {params_path}")

    with open(params_path) as f:
        params = yaml.safe_load(f)

    return params


def get_param(params: dict, *keys: str, default=None):
    """Safely get a nested parameter value.

    Args:
        params: Parameters dictionary.
        keys: Sequence of keys for nested access.
        default: Default value if key path not found.

    Returns:
        Parameter value or default.

    Example:
        get_param(params, "transcriptomics", "padj_threshold")
    """
    current = params
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            return default
    return current
