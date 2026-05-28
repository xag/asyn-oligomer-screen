"""Covalent / adduct channel for reactive metabolites.

Vina's flexible-ligand docking sees only reversible non-covalent binding.
Reactive electrophiles — methylglyoxal, 4-HNE, acrolein, malondialdehyde —
exert their α-syn effect through covalent adducts (Lys-CEL, His-Michael,
Lys-Schiff, …), which Vina cannot model. The absolute-affinity gate in
stage3 correctly collapses their Δactivity_gated to ~0; that is the
right default, but it leaves a real channel of biology unscored.

This module is the cheap first pass called out in STATUS.md
"Recommended next moves" item 1:

    aspr_score = (1/n_chains) · Σ_chains Σ_relevant_r
                     min(1, sasa(r) / SASA_REF) · rxty(ligand, restype(r))

It is purely receptor- and chemistry-derived; no Vina, no pose. Per-chain
mean keeps it on the same scale as the Stage 3 all-chain mean Δactivity
(both are averages over the chains of the assembly).

Sign convention: aspr_score is non-negative. Reactive adducts push α-syn
toward toxic conformers (crosslinking, charge neutralisation on Lys), so
a high aspr_score is *harm-leaning*. It is reported as an orthogonal
column, not folded into delta_activity_gated. Non-reactive ligands have
no entry in LIGAND_REACTIVITY and receive aspr_score = 0.

Usage:
    python adduct_score.py methylglyoxal results/oligomers/fusco_parallel_3mer_core70-88_relaxed.pdb
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio.PDB import PDBParser, SASA

ROOT = Path(__file__).parent

# -----------------------------------------------------------------------------
# Ligand reactivity table. Weights are intrinsic per-residue chemical rates,
# ordinalised on a [0, 1] scale. Primary target = 1.0, secondary = 0.3–0.7,
# minor = 0.1–0.3. Sources: the standard adduct-chemistry literature on each
# electrophile (Vicente-Miranda 2017 for MGO–Lys; Esterbauer 1991 / Bae 2013
# for 4-HNE; Uchida 1999 for acrolein; LoPachin 2014 for SBT reactivity rules).
#
# α-syn-specific consequence (WT, no mutation): 15 Lys, 0 Arg, 0 Cys, 1 His,
# 4 Tyr per chain. So MGO and glyoxal — Arg-primary chemically — fall back to
# their secondary Lys target (the empirically observed dominant adduct on
# α-syn, per Vicente-Miranda 2017). 4-HNE and acrolein — Cys-primary
# chemically — fall back to His/Lys with attenuated efficiency.
# -----------------------------------------------------------------------------

LIGAND_REACTIVITY: dict[str, dict[str, float]] = {
    # Reactive dicarbonyls (glycation). Primary: Arg (MGH adduct). Secondary:
    # Lys (CEL/CML). Minor: Cys (rapid but reversible thiohemiketal).
    "methylglyoxal": {"R": 1.0, "K": 0.4, "C": 0.3},
    "glyoxal":       {"R": 1.0, "K": 0.4, "C": 0.3},

    # α,β-unsaturated aldehydes (Michael acceptors). Primary: Cys (fastest
    # by 2–3 orders of magnitude). Secondary: His. Tertiary: Lys via
    # Schiff base + Michael combination.
    "4-hne":         {"C": 1.0, "H": 0.6, "K": 0.4},
    "acrolein":      {"C": 1.0, "H": 0.5, "K": 0.6},  # K bumped: acrolein-Lys is well-characterised
    "crotonaldehyde":{"C": 1.0, "H": 0.5, "K": 0.5},

    # Dialdehydes — Schiff base on Lys, with 1,3-crosslinking. Cys side
    # reactions exist but minor.
    "malondialdehyde": {"K": 1.0, "C": 0.2},
    "formaldehyde":    {"K": 0.7, "C": 0.5, "R": 0.3},

    # Reactive nitrogen / sulfur species — included for completeness but
    # tagged "low confidence": NO acts on Tyr (nitration) and Cys
    # (S-nitrosylation); H2S persulfidates Cys.
    "nitric-oxide":     {"Y": 0.8, "C": 1.0},
    "hydrogen-sulfide": {"C": 1.0},
}

# Reference per-residue accessible surface area in Å². Roughly the
# fully-exposed-in-Gly-X-Gly tripeptide SASA of a typical exposed side
# chain; used to normalise per-residue SASA to a fractional accessibility
# in [0, 1]. Single value across residue types keeps the scoring
# interpretable — the residue-specific differences live in the
# reactivity weights, not in the accessibility normaliser.
SASA_REF = 200.0


# Three-letter to one-letter map (subset relevant for adduct chemistry).
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _compute_sasa_per_residue(structure) -> None:
    """Attach Shrake-Rupley per-residue SASA to the first model in-place."""
    if getattr(structure, "_aspr_sasa_done", False):
        return
    model = next(iter(structure))
    SASA.ShrakeRupley().compute(model, level="R")
    structure._aspr_sasa_done = True


def _protein_chain_ids(structure) -> list[str]:
    out = []
    for ch in next(iter(structure)):
        if any(r.id[0] == " " and r.get_resname() in _THREE_TO_ONE for r in ch):
            out.append(ch.id)
    return out


def aspr_score(
    structure,
    mol_id: str,
    chain_ids: list[str] | None = None,
    sasa_ref: float = SASA_REF,
    breakdown: bool = False,
) -> float | dict:
    """Per-chain mean adduct surface propensity for ligand mol_id.

    If `breakdown=True`, returns a dict {'score', 'per_residue_type', 'rxty'}
    instead of the scalar, useful for inspection.

    Returns 0.0 (or {'score': 0.0, ...}) for any ligand not listed in
    LIGAND_REACTIVITY — including all the polyphenols, steroids,
    neurotransmitters, etc. that act through non-covalent channels.
    """
    rxty = LIGAND_REACTIVITY.get(mol_id)
    if rxty is None:
        return {"score": 0.0, "per_residue_type": {}, "rxty": {}} if breakdown else 0.0

    _compute_sasa_per_residue(structure)
    if chain_ids is None:
        chain_ids = _protein_chain_ids(structure)
    if not chain_ids:
        return {"score": 0.0, "per_residue_type": {}, "rxty": rxty} if breakdown else 0.0

    per_type = {one: 0.0 for one in rxty}  # accumulator across all chains
    for ch in next(iter(structure)):
        if ch.id not in chain_ids:
            continue
        for residue in ch:
            if residue.id[0] != " ":
                continue
            one = _THREE_TO_ONE.get(residue.get_resname())
            if one is None or one not in rxty:
                continue
            sasa = float(getattr(residue, "sasa", 0.0))
            frac = min(1.0, sasa / sasa_ref) if sasa_ref > 0 else 0.0
            per_type[one] += frac * rxty[one]

    score = sum(per_type.values()) / len(chain_ids)
    if breakdown:
        return {
            "score": float(score),
            "per_residue_type": {k: float(v / len(chain_ids)) for k, v in per_type.items()},
            "rxty": rxty,
        }
    return float(score)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def _load_structure(spec: str):
    p = Path(spec)
    if p.exists() and p.suffix.lower() == ".pdb":
        return PDBParser(QUIET=True).get_structure(p.stem, str(p))
    # 4-letter PDB id — load via assembly module to get BIOMT-expanded
    # assembly consistent with Stage 2 / Stage 3 fibril scoring.
    import sys
    sys.path.insert(0, str(ROOT))
    from assembly import load_assembly
    return load_assembly(spec)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("molecule", help="vicinity-molecule id, e.g. 'methylglyoxal'")
    p.add_argument("receptor", help="anchor PDB id (e.g. '6PEO') or path to oligomer PDB")
    p.add_argument("--json", action="store_true", help="emit breakdown as JSON")
    args = p.parse_args()

    structure = _load_structure(args.receptor)
    out = aspr_score(structure, args.molecule, breakdown=True)

    if args.json:
        print(json.dumps(out, indent=2))
        return

    rxty = out["rxty"]
    if not rxty:
        print(f"{args.molecule!r} not in LIGAND_REACTIVITY (aspr_score = 0)")
        print("Known reactive ligands:", ", ".join(sorted(LIGAND_REACTIVITY)))
        return
    print(f"=== aspr_score: {args.molecule} × {args.receptor} ===")
    print(f"  reactivity weights: {rxty}")
    print(f"  per-residue-type contributions (per-chain mean):")
    for k, v in out["per_residue_type"].items():
        print(f"    {k}: {v:+.4f}")
    print(f"  aspr_score = {out['score']:+.4f}")


if __name__ == "__main__":
    main()
