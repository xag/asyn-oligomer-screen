"""Scan vicinity-molecule list against a generated oligomer target.

Loads all vicinity molecules with non-null SMILES from lib/vicinity_molecules.js,
docks each against the specified oligomer PDB via stage3.perturb_oligomer(),
and writes a ranked summary CSV.

Usage:
    python sweep_oligomer.py
    python sweep_oligomer.py --oligomer results/oligomers/fusco_parallel_3mer_core70-88_relaxed.pdb
    python sweep_oligomer.py --exhaustiveness 16 --skip-existing

Outputs:
    results/sweep/<oligomer-stem>_sweep.csv  — full results, ranked by dact_gated
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import traceback
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))

from adduct_score import LIGAND_REACTIVITY, aspr_score
from stage3 import (
    RESULTS as STAGE3_RESULTS,
    _ENTRY_RE,
    _NULL_FIELD_RE,
    _STR_FIELD_RE,
    load_structure_from_file,
    perturb_oligomer,
)

VICINITY_JS = ROOT / "data" / "vicinity_molecules.js"
SWEEP_DIR = ROOT / "results" / "sweep"

DEFAULT_OLIGOMER = ROOT / "results" / "oligomers" / "fusco_parallel_3mer_core70-88_relaxed.pdb"


def load_all_vicinity_molecules() -> list[dict]:
    """Return all entries with non-null SMILES from vicinity_molecules.js."""
    text = VICINITY_JS.read_text(encoding="utf-8")
    molecules = []
    for m in _ENTRY_RE.finditer(text):
        mol_id = m.group(1)
        body = m.group(2)
        entry: dict = {"id": mol_id}
        for fm in _STR_FIELD_RE.finditer(body):
            entry.setdefault(fm.group(1), fm.group(2))
        for fm in _NULL_FIELD_RE.finditer(body):
            entry.setdefault(fm.group(1), None)
        if entry.get("smiles"):
            molecules.append(entry)
    return molecules


def sweep(
    oligo_pdb: Path,
    exhaustiveness: int = 8,
    n_poses: int = 5,
    seed: int = 42,
    skip_existing: bool = True,
) -> pd.DataFrame:
    oligo_pdb = Path(oligo_pdb)
    if not oligo_pdb.exists():
        raise FileNotFoundError(f"oligomer PDB not found: {oligo_pdb}")

    SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    molecules = load_all_vicinity_molecules()
    print(f"Sweeping {len(molecules)} molecules with SMILES against {oligo_pdb.name}")
    print(f"Exhaustiveness={exhaustiveness}, n_poses={n_poses}, seed={seed}")
    print(f"Skip-existing: {skip_existing}\n")

    # Pre-load the receptor once for the aspr (covalent / adduct) channel.
    # aspr depends only on the apo structure (per-residue SASA × ligand
    # reactivity table), not on the Vina pose, so backfilling cached
    # reports is cheap and adds no Vina calls.
    aspr_struct = load_structure_from_file(oligo_pdb)
    aspr_chain_ids = [
        ch.id for ch in next(iter(aspr_struct))
        if any(r.id[0] == " " for r in ch)
    ]

    stem = oligo_pdb.stem
    rows = []
    for i, mol in enumerate(molecules, 1):
        mol_id = mol["id"]
        pair_tag = f"{mol_id}_{stem}"
        report_path = STAGE3_RESULTS / f"{pair_tag}_report.json"

        if skip_existing and report_path.exists():
            print(f"[{i:3d}/{len(molecules)}] {mol_id:30s} (skip — exists)")
            report = json.loads(report_path.read_text(encoding="utf-8"))
        else:
            print(f"[{i:3d}/{len(molecules)}] {mol_id:30s} ...", end=" ", flush=True)
            try:
                report = perturb_oligomer(
                    mol_id, oligo_pdb,
                    exhaustiveness=exhaustiveness, n_poses=n_poses, seed=seed,
                )
                aff = report["vina_top_affinity_kcal_per_mol"]
                dg = report["delta_activity_gated"]
                print(f"aff={aff:+.2f}  dact_gated={dg:+.4f}")
            except Exception as exc:
                print(f"FAILED: {exc}")
                traceback.print_exc()
                rows.append({
                    "mol_id": mol_id,
                    "mol_name": mol.get("name", ""),
                    "status": "error",
                    "error": str(exc),
                })
                continue

        # aspr (covalent / adduct) channel — recomputed unconditionally so
        # the column is present even on cached reports written before the
        # channel landed. Cheap: one SASA pass shared across all ligands.
        aspr_val = float(aspr_score(aspr_struct, mol_id, chain_ids=aspr_chain_ids))
        is_reactive = mol_id in LIGAND_REACTIVITY
        if "aspr_score" not in report:
            report["aspr_score"] = aspr_val
            report["aspr_reactive"] = is_reactive
            report_path.write_text(json.dumps(report, indent=2))

        rows.append({
            "mol_id": mol_id,
            "mol_name": report.get("molecule_name", mol.get("name", "")),
            "smiles": report.get("smiles", mol.get("smiles", "")),
            "status": "ok",
            "apo_activity": report.get("activity_apo"),
            "vina_top_affinity_kcal_per_mol": report.get("vina_top_affinity_kcal_per_mol"),
            "n_poses": report.get("n_poses"),
            "delta_activity_top": report.get("delta_activity_top"),
            "delta_activity_weighted": report.get("delta_activity_weighted"),
            "affinity_gate": report.get("affinity_gate"),
            "delta_activity_gated": report.get("delta_activity_gated"),
            "aspr_reactive": is_reactive,
            "aspr_score": aspr_val,
        })

    df = pd.DataFrame(rows)
    ok = df[df["status"] == "ok"].copy()
    if not ok.empty:
        ok = ok.sort_values("delta_activity_gated")  # most negative = most protective first
    err = df[df["status"] != "ok"]
    out = pd.concat([ok, err], ignore_index=True)

    csv_path = SWEEP_DIR / f"{stem}_sweep.csv"
    out.to_csv(csv_path, index=False)
    print(f"\nSweep complete. {len(ok)} succeeded, {len(err)} failed.")
    print(f"Results: {csv_path}")

    if not ok.empty:
        print("\nTop 10 most protective (by delta_activity_gated):")
        top = ok.head(10)[["mol_id", "mol_name", "vina_top_affinity_kcal_per_mol",
                            "delta_activity_gated", "affinity_gate"]]
        print(top.to_string(index=False))

        reactive = ok[ok["aspr_reactive"] == True].sort_values("aspr_score", ascending=False)
        if not reactive.empty:
            print("\nReactive metabolites by covalent / adduct channel (aspr_score):")
            cols = ["mol_id", "vina_top_affinity_kcal_per_mol",
                    "delta_activity_gated", "aspr_score"]
            print(reactive[cols].to_string(index=False))

    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--oligomer", type=Path, default=DEFAULT_OLIGOMER,
        help="path to relaxed oligomer PDB (default: reference trimer)",
    )
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--num_modes", type=int, default=5)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--skip-existing", action="store_true", default=True,
                   help="skip molecules whose report JSON already exists")
    p.add_argument("--no-skip", dest="skip_existing", action="store_false",
                   help="re-dock all molecules even if results exist")
    args = p.parse_args()

    sweep(
        args.oligomer,
        exhaustiveness=args.exhaustiveness,
        n_poses=args.num_modes,
        seed=args.seed,
        skip_existing=args.skip_existing,
    )


if __name__ == "__main__":
    main()
