"""Weighted-feature activity score. Transparent and hand-tunable. Migrates
to a small NN once we have enough labelled data; not before.

Per the plan, contact_density is INVERSELY related to activity (fibrils are
stable and inert; oligomers are metastable). All other features are
positively related.

Each feature is z-scored across the input set before weighting, so the
weights express relative importance independent of feature units.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

WEIGHTS = {
    "exposed_hydrophobic_beta_sasa": 1.0,
    "membrane_insertion_propensity": 0.5,
    "nac_active_score": 1.0,
    "contact_density": -1.0,
    "disordered_hydrophobic_exposure": 0.5,
}


@dataclass
class Score:
    pdb_id: str
    label: str
    description: str
    raw: dict
    z: dict
    activity: float


def score_table(df: pd.DataFrame) -> pd.DataFrame:
    cols = list(WEIGHTS.keys())
    z = df[cols].apply(lambda c: (c - c.mean()) / c.std(ddof=0))
    activity = sum(z[c] * w for c, w in WEIGHTS.items())
    out = df.copy()
    for c in cols:
        out[f"z_{c}"] = z[c]
    out["activity"] = activity
    return out.sort_values("activity", ascending=False)
