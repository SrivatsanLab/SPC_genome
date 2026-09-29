"""Config loading with defaults, path resolution, and validation.

Configs are plain YAML. A single default-config file sits under
``phylo/utils/_defaults.yaml``; user configs override keys via a shallow merge.
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Any

import yaml


def _read_yaml(path: str | Path) -> dict[str, Any]:
    with open(path) as f:
        return yaml.safe_load(f) or {}


def _defaults() -> dict[str, Any]:
    with resources.files("phylo.utils").joinpath("_defaults.yaml").open() as f:
        return yaml.safe_load(f) or {}


def _deep_merge(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    """Recursive dict merge. Overlay wins for scalars and lists."""
    out = dict(base)
    for k, v in overlay.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a user YAML config and merge onto defaults.

    Parameters
    ----------
    path : str or Path
        Path to a user YAML file.

    Returns
    -------
    dict
        Resolved config dict.
    """
    user = _read_yaml(path)
    merged = _deep_merge(_defaults(), user)
    if "name" not in merged or merged["name"] is None:
        merged["name"] = Path(path).stem
    return resolve_config(merged, config_path=Path(path))


def resolve_config(cfg: dict[str, Any], config_path: Path | None = None) -> dict[str, Any]:
    """Expand relative paths against the config file's directory.

    Paths in ``inputs.*`` and ``output_root`` are resolved to absolute paths.
    """
    root = config_path.parent.resolve() if config_path is not None else Path.cwd()

    def _resolve(p: Any) -> Any:
        if not isinstance(p, str):
            return p
        pth = Path(p).expanduser()
        return str(pth if pth.is_absolute() else (root / pth).resolve())

    if "inputs" in cfg:
        cfg["inputs"] = {k: _resolve(v) for k, v in cfg["inputs"].items()}
    if "output_root" in cfg:
        cfg["output_root"] = _resolve(cfg["output_root"])
    # Additional path-holding sub-blocks: resolve any *_path keys.
    for section in ("clade_matching",):
        if section in cfg and isinstance(cfg[section], dict):
            for k, v in list(cfg[section].items()):
                if k.endswith("_path"):
                    cfg[section][k] = _resolve(v)
    return cfg
