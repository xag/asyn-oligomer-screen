"""MD-relax the top-pose Stage 3 complexes and re-score Stage 2 features.

For each (mol, anchor) pair with an existing `<pair>_complex.pdb` in
`results/stage3/`, this script:

  1. Invokes `md_relax.py` in the `md` conda env (OpenMM + OpenFF live
     there; the pip venv that runs the rest of the pipeline doesn't).
     Inputs: the complex PDB + ligand SMILES from vicinity_molecules.js.
     Output: `<pair>_relaxed.pdb` (protein with hydrogens + ligand on
     chain Z residue LIG; waters / ions stripped).

  2. Loads the relaxed PDB with Biopython, runs Stage 2 features on the
     inner chain (which is the same chain id pdbfixer preserved from the
     complex), and computes `delta_activity_relaxed` against the Stage 2
     anchor z-axis.

  3. Updates `<pair>_report.json` with the relaxed feature block and
     `activity_complex_relaxed` / `delta_activity_relaxed` alongside the
     existing top-pose, Boltzmann-weighted, and gated columns.

The point of this step is to test whether MD relaxation produces
`delta_activity_relaxed > 0` for harm-side ligands. Under the current
static-pose feature weights, Δactivity is bounded above by zero
(any docked pose occludes SASA and adds Cα–LIG contacts, both of
which subtract from activity). The only path to a sign flip is
conformational rearrangement of the receptor around the bound pose,
which is exactly what NPT MD provides. See STATUS.md "Recommended
next moves" #1.

Usage:
    python md_stage3.py                    # all 6 PAIRS, defaults
    python md_stage3.py --equil-ps 200 --prod-ps 1000   # 1 ns production
    python md_stage3.py --only curcumin 6PEO            # one pair

The `md_relax.py` step is the slow one. On OpenCL/GPU the default
100 ps equil + 100 ps prod takes ~3 min/pair; 1 ns production is ~35
min/pair.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser

from stage3 import (
    compute_features_on,
    delta_activity,
    load_vicinity_molecule,
)

ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "results" / "stage3"
_MD_PY_ENV = os.environ.get("ASYN_MD_PYTHON")
if not _MD_PY_ENV:
    raise RuntimeError(
        "Set ASYN_MD_PYTHON to the python interpreter of a conda env that has "
        "OpenMM + openff-toolkit + openmmforcefields installed. "
        "See README.md (Reproduction / MD relaxation)."
    )
MD_PYTHON = Path(_MD_PY_ENV)
MD_RELAX = Path(__file__).parent / "md_relax.py"


def _run_md(cmd: list[str], out_pdb: Path, tag: str) -> None:
    print("  $ " + " ".join(cmd), flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.stdout:
        for line in res.stdout.splitlines():
            print(f"  md > {line}", flush=True)
    if res.returncode != 0 or not out_pdb.exists():
        if res.stderr:
            print(res.stderr, file=sys.stderr)
        raise RuntimeError(f"md_relax failed (rc={res.returncode}) for {tag}")


def run_md_complex(complex_pdb: Path, ligand_smiles: str, out_pdb: Path,
                   equil_ps: float, prod_ps: float) -> None:
    cmd = [
        str(MD_PYTHON), str(MD_RELAX),
        "--complex-pdb", str(complex_pdb),
        "--ligand-smiles", ligand_smiles,
        "--out-pdb", str(out_pdb),
        "--equil-ps", str(equil_ps),
        "--prod-ps", str(prod_ps),
    ]
    _run_md(cmd, out_pdb, complex_pdb.name)


def run_md_apo(apo_pdb: Path, out_pdb: Path,
               equil_ps: float, prod_ps: float) -> Path:
    """Cache by output path — apo MD is per-anchor, shared across pairs."""
    if out_pdb.exists():
        print(f"  apo cache hit: {out_pdb.name}", flush=True)
        return out_pdb
    cmd = [
        str(MD_PYTHON), str(MD_RELAX),
        "--apo-pdb", str(apo_pdb),
        "--out-pdb", str(out_pdb),
        "--equil-ps", str(equil_ps),
        "--prod-ps", str(prod_ps),
    ]
    _run_md(cmd, out_pdb, apo_pdb.name)
    return out_pdb


def relax_and_rescore(mol_id: str, pdb_id: str,
                      equil_ps: float, prod_ps: float,
                      anchor_df: pd.DataFrame) -> dict:
    pair_tag = f"{mol_id}_{pdb_id.lower()}"
    pdb_tag = pdb_id.lower()
    report_path = STAGE3 / f"{pair_tag}_report.json"
    complex_pdb = STAGE3 / f"{pair_tag}_complex.pdb"
    apo_pdb_input = STAGE3 / f"{pair_tag}_apo.pdb"
    complex_relaxed_pdb = STAGE3 / f"{pair_tag}_relaxed.pdb"
    apo_relaxed_pdb = STAGE3 / f"{pdb_tag}_apo_relaxed.pdb"

    if not report_path.exists() or not complex_pdb.exists() or not apo_pdb_input.exists():
        raise FileNotFoundError(f"missing artifacts for {pair_tag}")

    report = json.loads(report_path.read_text(encoding="utf-8"))
    dock_chain = report["inner_chain_dock"]
    smiles = report.get("smiles") or load_vicinity_molecule(mol_id)["smiles"]
    apo_feats_static = report["apo_features"]

    # MD on apo and complex, both with the same equil/prod budget so the
    # baseline-MD sampling drift cancels in delta_activity_relaxed =
    # activity(complex_relaxed) - activity(apo_relaxed).
    run_md_apo(apo_pdb_input, apo_relaxed_pdb, equil_ps=equil_ps, prod_ps=prod_ps)
    run_md_complex(complex_pdb, smiles, complex_relaxed_pdb,
                   equil_ps=equil_ps, prod_ps=prod_ps)

    parser = PDBParser(QUIET=True)
    apo_struct = parser.get_structure(f"{pair_tag}_apo", str(apo_relaxed_pdb))
    complex_struct = parser.get_structure(pair_tag, str(complex_relaxed_pdb))
    apo_feats_md = compute_features_on(apo_struct, chain_ids=[dock_chain])
    complex_feats_md = compute_features_on(complex_struct, chain_ids=[dock_chain])

    d = delta_activity(apo_feats_md, complex_feats_md, anchor_df)

    report["md_equil_ps"] = float(equil_ps)
    report["md_prod_ps"] = float(prod_ps)
    report["apo_relaxed_pdb_path"] = apo_relaxed_pdb.name
    report["complex_relaxed_pdb_path"] = complex_relaxed_pdb.name
    report["apo_features_relaxed"] = apo_feats_md
    report["complex_features_relaxed"] = complex_feats_md
    report["delta_features_relaxed"] = {k: complex_feats_md[k] - apo_feats_md[k] for k in complex_feats_md}
    report["activity_apo_relaxed"] = d["activity_apo"]
    report["activity_complex_relaxed"] = d["activity_complex"]
    report["delta_activity_relaxed"] = d["delta_activity"]
    # Also keep an apo-static baseline delta for diagnostic comparison.
    d_static_apo = delta_activity(apo_feats_static, complex_feats_md, anchor_df)
    report["delta_activity_relaxed_vs_static_apo"] = d_static_apo["delta_activity"]
    report_path.write_text(json.dumps(report, indent=2))
    return report


PAIRS = [
    ("curcumin", "6PEO"),
    ("curcumin", "6UFR"),
    ("quercetin", "6PEO"),
    ("methylglyoxal", "6PEO"),
    ("curcumin", "2N0A"),
    ("methylglyoxal", "2N0A"),
]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--equil-ps", type=float, default=100.0,
                   help="ps of NPT equilibration before production (default 100)")
    p.add_argument("--prod-ps", type=float, default=100.0,
                   help="ps of NPT production used for the relaxed snapshot (default 100)")
    p.add_argument("--only", nargs=2, metavar=("MOL", "PDB"),
                   help="run a single (molecule, anchor) pair instead of all")
    args = p.parse_args()

    pairs = [(args.only[0], args.only[1])] if args.only else PAIRS

    anchor_df = pd.read_csv(ROOT / "results" / "anchor_features.csv")

    rows: list[dict] = []
    for mol, pdb in pairs:
        print(f"\n=== {mol} × {pdb} ===", flush=True)
        r = relax_and_rescore(mol, pdb, equil_ps=args.equil_ps, prod_ps=args.prod_ps, anchor_df=anchor_df)
        rows.append({
            "mol": mol,
            "pdb": pdb,
            "aff_top": r["vina_top_affinity_kcal_per_mol"],
            "act_apo": r["activity_apo"],
            "act_apo_md": r["activity_apo_relaxed"],
            "dact_top": r["delta_activity_top"],
            "dact_wtd": r["delta_activity_weighted"],
            "dact_gated": r["delta_activity_gated"],
            "dact_relax": r["delta_activity_relaxed"],
        })

    print("\nMD-relaxation summary (equil={:.0f} ps, prod={:.0f} ps)".format(args.equil_ps, args.prod_ps))
    header = (
        f"{'molecule':14}{'anchor':7}"
        f"{'aff_top':>9}{'act_apo':>9}{'act_apoMD':>11}"
        f"{'dact_top':>10}{'dact_wtd':>10}{'dact_gat':>10}"
        f"{'dact_rlx':>10}"
    )
    print(header)
    print("-" * len(header))
    for r in rows:
        print(
            f"{r['mol']:14}{r['pdb']:7}"
            f"{r['aff_top']:>+9.2f}"
            f"{r['act_apo']:>+9.3f}"
            f"{r['act_apo_md']:>+11.3f}"
            f"{r['dact_top']:>+10.3f}"
            f"{r['dact_wtd']:>+10.3f}"
            f"{r['dact_gated']:>+10.3f}"
            f"{r['dact_relax']:>+10.3f}"
        )


if __name__ == "__main__":
    main()
