"""Run directory management and per-stage completion markers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

STAGES = ("g0", "panels", "trees", "validate", "covariates", "bwm", "clade_matching")


@dataclass
class RunDir:
    """Handle to a run's output directory."""

    path: Path

    @classmethod
    def create(cls, output_root: str | Path, name: str, cfg: dict[str, Any]) -> RunDir:
        p = Path(output_root) / "runs" / name
        p.mkdir(parents=True, exist_ok=True)
        with open(p / "config.resolved.yaml", "w") as f:
            yaml.safe_dump(cfg, f, sort_keys=False)
        return cls(path=p)

    def subdir(self, stage: str) -> Path:
        d = self.path / stage
        d.mkdir(parents=True, exist_ok=True)
        return d


def stage_marker(run: RunDir, stage: str) -> Path:
    """Path to the completion sentinel for ``stage``."""
    return run.path / f".done_{stage}"


def stage_done(run: RunDir, stage: str) -> bool:
    """Whether ``stage`` has completed for this run."""
    return stage_marker(run, stage).exists()


def write_stage_marker(run: RunDir, stage: str) -> None:
    """Mark ``stage`` as complete."""
    stage_marker(run, stage).touch()


def clear_stage_marker(run: RunDir, stage: str) -> None:
    """Remove the completion marker for ``stage`` (used by --force)."""
    m = stage_marker(run, stage)
    if m.exists():
        m.unlink()
