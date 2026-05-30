"""Locate the conda interpreter that builds ligand force fields.

The MD dynamics run in the pip venv. One step — building a ligand's force
field (`md_relax --prepare-only`, `md_stage3`) — runs OpenFF parametrisation,
which lives in the `asyn-md` conda env (see `environment-md.yml`). This module
finds that env's python so callers don't have to be told to set anything.

Resolution order:
  1. $ASYN_MD_PYTHON, if set (explicit override).
  2. an env named `asyn-md` or `md` under a standard conda/mamba root.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

_ENV_NAMES = ("asyn-md", "md")


def _python_in(env_root: Path) -> Path:
    return env_root / ("python.exe" if os.name == "nt" else "bin/python")


def _conda_roots() -> list[Path]:
    roots: list[Path] = []
    for var in ("CONDA_PREFIX", "CONDA_ROOT", "MAMBA_ROOT_PREFIX"):
        v = os.environ.get(var)
        if v:
            p = Path(v)
            # CONDA_PREFIX may point at an *env*; its base is two levels up.
            roots += [p, p.parent.parent]
    home = Path.home()
    roots += [home / n for n in ("miniforge3", "mambaforge", "miniconda3", "anaconda3")]
    roots += [Path("/opt/conda")]
    # de-dup, keep order
    seen, out = set(), []
    for r in roots:
        if r not in seen:
            seen.add(r)
            out.append(r)
    return out


def md_python() -> Path:
    """Return the conda python that has OpenFF installed, or raise with the
    one command that creates it."""
    override = os.environ.get("ASYN_MD_PYTHON")
    if override:
        return Path(override)
    for root in _conda_roots():
        for name in _ENV_NAMES:
            cand = _python_in(root / "envs" / name)
            if cand.exists():
                return cand
    raise RuntimeError(
        "No MD conda env found. Create it once:\n"
        "    conda env create -f environment-md.yml\n"
        "(or set $ASYN_MD_PYTHON to a python with openff-toolkit installed)."
    )


if __name__ == "__main__":
    print(md_python())
    sys.exit(0)
