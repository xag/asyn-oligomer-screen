"""Recompute Stage 3 dP(active) on already-docked poses.

Reuses the existing `<pair>_docked.pdbqt` (top pose unchanged because
Vina was run with seed=42, exhaustiveness=8) and the `<pair>_vina.log`
(for the affinity table) so no docking is repeated. Useful when a
Stage 2 feature changes, or when the multi-pose / weighting rule
changes, and we want the new dP without rerunning Vina across the
matrix.

Run after extending or modifying any feature in `features.py`, or after
changing the pose-aggregation rule in `stage3.score_all_poses`:

    python recompute_stage3.py

Overwrites each `<pair>_report.json` with the recomputed numbers and
prints a one-line summary per pair (top-pose vs Boltzmann-weighted
activity)."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from stage3 import (
    _parse_affinities,
    score_all_poses,
)

ROOT = Path(__file__).parent
STAGE3 = ROOT / "results" / "stage3"


def recompute(mol_id: str, pdb_id: str) -> dict:
    pair_tag = f"{mol_id}_{pdb_id.lower()}"
    report_path = STAGE3 / f"{pair_tag}_report.json"
    docked_pdbqt = STAGE3 / f"{pair_tag}_docked.pdbqt"
    vina_log = STAGE3 / f"{pair_tag}_vina.log"
    if not (report_path.exists() and docked_pdbqt.exists() and vina_log.exists()):
        raise FileNotFoundError(f"missing artifacts for {pair_tag}")

    old = json.loads(report_path.read_text(encoding="utf-8"))
    apo_chain = old["inner_chain_apo"]

    affinities = old.get("vina_all_affinities_kcal_per_mol")
    if not affinities:
        affinities = _parse_affinities(vina_log.read_text(encoding="utf-8"))

    anchor_df = pd.read_csv(ROOT / "results" / "anchor_features.csv")
    multi = score_all_poses(
        pdb_id, apo_chain, docked_pdbqt, affinities, anchor_df,
    )

    new = dict(old)
    new["vina_all_affinities_kcal_per_mol"] = [float(a) for a in affinities]
    new["n_poses"] = multi["n_poses"]
    new["apo_features"] = multi["apo_features"]
    new["complex_features_top"] = multi["complex_features_top"]
    new["delta_features_top"] = multi["delta_features_top"]
    new["activity_apo"] = multi["activity_apo"]
    new["activity_complex_top"] = multi["activity_complex_top"]
    new["delta_activity_top"] = multi["delta_activity_top"]
    new["activity_complex_weighted"] = multi["activity_complex_weighted"]
    new["delta_activity_weighted"] = multi["delta_activity_weighted"]
    new["affinity_gate_threshold_kcal_per_mol"] = multi["affinity_gate_threshold_kcal_per_mol"]
    new["affinity_gate"] = multi["affinity_gate"]
    new["delta_activity_gated"] = multi["delta_activity_gated"]
    new["poses"] = multi["poses"]
    # Strip the pre-multi-pose keys so the report shape is consistent.
    for legacy in ("complex_features", "delta_features", "activity_complex",
                   "delta_activity", "z_apo", "z_complex"):
        new.pop(legacy, None)
    report_path.write_text(json.dumps(new, indent=2))
    return new


PAIRS = [
    ("curcumin", "6PEO"),
    ("curcumin", "6UFR"),
    ("quercetin", "6PEO"),
    ("methylglyoxal", "6PEO"),
    ("curcumin", "2N0A"),
    ("methylglyoxal", "2N0A"),
]


def main() -> None:
    header = (
        f"{'molecule':14}{'anchor':7}"
        f"{'aff_top':>9}{'act_apo':>9}"
        f"{'dact_top':>10}{'dact_wtd':>10}"
        f"{'gate':>8}{'dact_gat':>10}"
        f"{'w_top':>8}"
    )
    print(header)
    print("-" * len(header))
    for mol, pdb in PAIRS:
        r = recompute(mol, pdb)
        w_top = r["poses"][0]["weight"]
        print(
            f"{mol:14}{pdb:7}"
            f"{r['vina_top_affinity_kcal_per_mol']:>+9.2f}"
            f"{r['activity_apo']:>+9.3f}"
            f"{r['delta_activity_top']:>+10.3f}"
            f"{r['delta_activity_weighted']:>+10.3f}"
            f"{r['affinity_gate']:>8.3f}"
            f"{r['delta_activity_gated']:>+10.3f}"
            f"{w_top:>8.3f}"
        )


if __name__ == "__main__":
    main()
