"""Shape metrics for the dwell-time channel (issue #14).

Two structure-vs-structure measurements that quantify *how far a frame has
moved away from a reference toxic shape*:

  1. ``beta_core_rmsd`` — Cα RMSD of the β-core (residues 70-88 by
     default, the Fusco strand) after optimal superposition. Measures
     whether the load-bearing β-sheet is still where the toxic reference
     put it. A ligand that loosens the sheet shows up as larger RMSD.

  2. ``contact_jaccard`` — Jaccard similarity of the *inter-chain* Cα
     contact map over the NAC region (60-100 by default). Measures
     whether the chains are still packed against each other the same
     way. A ligand that peels chains apart shows up as lower Jaccard.

Both are pure geometry on Biopython structures — no MD, no force field —
so they run in the pip venv and are testable on the existing
``results/oligomers/*_relaxed.pdb`` structures. ``dwell_time.py`` calls
them once per MD trajectory frame to build the dwell-time distribution
that issue #14 specifies.

Why these two: they are orthogonal. RMSD sees the sheet *deforming in
place*; Jaccard sees chains *dissociating* without the surviving sheet
necessarily deforming. A toxic shape needs both an intact β-core and
intact inter-chain packing, so a frame "dwells" in the toxic basin only
when it scores low RMSD *and* high Jaccard against the reference (see
``in_toxic_basin``).

Usage (standalone sanity check — compares two PDBs):
    python screen/shape_metrics.py REF.pdb MOBILE.pdb
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Structure import Structure

# Standard amino acids — guards against scoring HETATM / ligand residues.
_AA_3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

# Fusco β-core strand (residues 70-88) — the restrained sheet in the
# Track A build. RMSD is measured here because this is the part the
# topology prior pins down and the part a destabiliser has to disrupt.
BETA_CORE_RANGE = (70, 88)

# NAC region (60-100) — the structured, aggregation-driving stretch.
# Inter-chain contacts are counted here; the disordered tails (1-59,
# 101-140) contribute only sampling noise to a contact map.
NAC_RANGE = (60, 100)

# Default Cα–Cα cutoff for an inter-chain contact. 9 Å is the standard
# coarse Cα contact radius (one residue's Cα to another's, allowing for
# side chains reaching ~4-5 Å each). Matches the 8 Å heavy-atom cutoff
# used in features.contact_density without being so tight that thermal
# wobble flips contacts on and off frame-to-frame.
CONTACT_CUTOFF = 9.0

# Default toxic-basin thresholds. A frame is "in the toxic basin" when
# its β-core has not drifted past TOXIC_RMSD_MAX Å and it retains at
# least TOXIC_JACCARD_MIN of the reference inter-chain contacts. These
# are deliberately lenient starting points; dwell_time.py recalibrates
# the operating point from the apo distribution (issue #30's threshold
# work depends on the apo dwell spread, not on these constants).
TOXIC_RMSD_MAX = 3.0
TOXIC_JACCARD_MIN = 0.5


def load_pdb(path: Path | str) -> Structure:
    path = Path(path)
    return PDBParser(QUIET=True).get_structure(path.stem, str(path))


# Ligand convention shared with md_relax / stage3: docked ligand sits on
# chain Z as residue LIG. The occupancy check looks for it there.
LIGAND_CHAIN = "Z"
LIGAND_RESNAME = "LIG"

# Default Cα–ligand cutoff for "the ligand is still at the β-core site".
# A frame whose nearest β-core Cα is further than this from every ligand
# heavy atom is counted as unbound (the ligand diffused off).
OCCUPANCY_CUTOFF = 6.0


def _first_model(structure):
    """Return the model to iterate chains over. Accepts a Structure
    (returns model 0) or a Model directly (returns it) so callers can
    pass one trajectory frame (a Model) without rewrapping it."""
    if getattr(structure, "level", None) == "M":
        return structure
    return next(iter(structure))


def _ca_by_residue(
    structure: Structure,
    resnum_range: tuple[int, int] | None,
    chain_ids: set[str] | None,
) -> dict[tuple[str, int], np.ndarray]:
    """Map (chain_id, resnum) → Cα coordinate for standard-AA residues
    in ``resnum_range`` (inclusive) and, if given, ``chain_ids``."""
    lo, hi = resnum_range if resnum_range is not None else (-(10**9), 10**9)
    out: dict[tuple[str, int], np.ndarray] = {}
    for chain in _first_model(structure):
        if chain_ids is not None and chain.id not in chain_ids:
            continue
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if residue.get_resname() not in _AA_3:
                continue
            resnum = residue.id[1]
            if not (lo <= resnum <= hi):
                continue
            if "CA" not in {a.get_name() for a in residue}:
                continue
            out[(chain.id, resnum)] = residue["CA"].coord
    return out


# -----------------------------------------------------------------------------
# β-core Cα RMSD after optimal superposition.
# -----------------------------------------------------------------------------

def beta_core_rmsd(
    mobile: Structure,
    reference: Structure,
    core_range: tuple[int, int] = BETA_CORE_RANGE,
    chain_ids: set[str] | None = None,
) -> float:
    """Cα RMSD (Å) of the β-core between ``mobile`` and ``reference`` after
    optimal rigid-body superposition on those same β-core Cα atoms.

    Only residues present in *both* structures (same chain id + residue
    number) are used, so a frame that lost or gained a residue during MD
    is still comparable on the shared core. Raises if fewer than 3 common
    Cα atoms remain (superposition is undefined)."""
    ref_ca = _ca_by_residue(reference, core_range, chain_ids)
    mob_ca = _ca_by_residue(mobile, core_range, chain_ids)
    common = sorted(ref_ca.keys() & mob_ca.keys())
    if len(common) < 3:
        raise ValueError(
            f"β-core superposition needs ≥3 shared Cα; got {len(common)} "
            f"for range {core_range} (ref={len(ref_ca)}, mobile={len(mob_ca)})"
        )
    ref_coords = np.array([ref_ca[k] for k in common], dtype=float)
    mob_coords = np.array([mob_ca[k] for k in common], dtype=float)
    sup = Superimposer()
    # Superimposer works on Atom objects; feed lightweight shims carrying
    # coordinates so we don't depend on Atom identity across structures.
    sup.set_atoms(_coord_atoms(ref_coords), _coord_atoms(mob_coords))
    return float(sup.rms)


class _CoordAtom:
    """Minimal stand-in exposing the Atom interface Superimposer needs
    (``get_coord`` / ``coord`` / ``set_coord``)."""

    __slots__ = ("coord",)

    def __init__(self, coord: np.ndarray):
        self.coord = np.asarray(coord, dtype=float)

    def get_coord(self) -> np.ndarray:
        return self.coord

    def set_coord(self, c: np.ndarray) -> None:
        self.coord = np.asarray(c, dtype=float)


def _coord_atoms(coords: np.ndarray) -> list[_CoordAtom]:
    return [_CoordAtom(c) for c in coords]


# -----------------------------------------------------------------------------
# Inter-chain Cα contact-map Jaccard.
# -----------------------------------------------------------------------------

def interchain_contact_set(
    structure: Structure,
    cutoff: float = CONTACT_CUTOFF,
    resnum_range: tuple[int, int] | None = NAC_RANGE,
) -> set[tuple[str, int, str, int]]:
    """Set of inter-chain Cα contacts in ``structure``. Each contact is a
    canonical key ``(chain_a, resnum_a, chain_b, resnum_b)`` with the two
    (chain, resnum) endpoints sorted, so the relation is symmetric and the
    same physical contact yields one key regardless of iteration order.
    Only pairs on *different* chains within ``cutoff`` Å are kept."""
    ca = _ca_by_residue(structure, resnum_range, chain_ids=None)
    keys = list(ca.keys())
    coords = np.array([ca[k] for k in keys], dtype=float)
    contacts: set[tuple[str, int, str, int]] = set()
    if len(keys) < 2:
        return contacts
    # Pairwise distances; vectorised over the (typically few-hundred) Cα.
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    iu, ju = np.triu_indices(len(keys), k=1)
    close = dist[iu, ju] < cutoff
    for i, j in zip(iu[close], ju[close]):
        (ca_i, ri), (ca_j, rj) = keys[i], keys[j]
        if ca_i == ca_j:
            continue  # intra-chain — not an inter-chain packing contact
        a, b = sorted(((ca_i, ri), (ca_j, rj)))
        contacts.add((a[0], a[1], b[0], b[1]))
    return contacts


def contact_jaccard(
    a: Structure,
    b: Structure,
    cutoff: float = CONTACT_CUTOFF,
    resnum_range: tuple[int, int] | None = NAC_RANGE,
) -> float:
    """Jaccard similarity |A∩B|/|A∪B| of the two inter-chain contact maps.
    1.0 = identical packing; 0.0 = no shared inter-chain contact. Two
    contact-free structures return 1.0 (vacuously identical)."""
    sa = interchain_contact_set(a, cutoff, resnum_range)
    sb = interchain_contact_set(b, cutoff, resnum_range)
    union = sa | sb
    if not union:
        return 1.0
    return len(sa & sb) / len(union)


# -----------------------------------------------------------------------------
# Combined dwell-basin test.
# -----------------------------------------------------------------------------

def frame_metrics(
    frame: Structure,
    reference: Structure,
    core_range: tuple[int, int] = BETA_CORE_RANGE,
    contact_range: tuple[int, int] | None = NAC_RANGE,
    cutoff: float = CONTACT_CUTOFF,
) -> dict[str, float]:
    """Both shape metrics for one frame against the toxic reference."""
    return {
        "beta_core_rmsd": beta_core_rmsd(frame, reference, core_range),
        "contact_jaccard": contact_jaccard(frame, reference, cutoff, contact_range),
    }


def in_toxic_basin(
    metrics: dict[str, float],
    rmsd_max: float = TOXIC_RMSD_MAX,
    jaccard_min: float = TOXIC_JACCARD_MIN,
) -> bool:
    """A frame dwells in the toxic basin when the β-core has not drifted
    past ``rmsd_max`` Å *and* it retains at least ``jaccard_min`` of the
    reference inter-chain contacts. Both conditions are required: the
    toxic shape is an intact sheet *and* intact inter-chain packing."""
    return (
        metrics["beta_core_rmsd"] <= rmsd_max
        and metrics["contact_jaccard"] >= jaccard_min
    )


# -----------------------------------------------------------------------------
# Ligand-occupancy check (issue #14: don't score a diffused-off ligand as
# "no effect"). A frame counts as "bound" when some ligand heavy atom is
# within `cutoff` Å of some β-core Cα.
# -----------------------------------------------------------------------------

def _ligand_heavy_coords(structure) -> np.ndarray:
    coords: list[np.ndarray] = []
    for chain in _first_model(structure):
        for residue in chain:
            if residue.get_resname() != LIGAND_RESNAME:
                continue
            for atom in residue:
                el = (atom.element or atom.get_name()[0:1]).strip().upper()
                if el == "H":
                    continue
                coords.append(atom.coord)
    return np.asarray(coords, dtype=float)


def ligand_bound(
    frame,
    core_range: tuple[int, int] = BETA_CORE_RANGE,
    cutoff: float = OCCUPANCY_CUTOFF,
) -> bool | None:
    """Is the ligand still at the β-core in this frame? True/False for a
    frame that has a LIG residue; None when there is no ligand (an apo
    frame), so callers can distinguish "apo, occupancy not applicable"
    from "complex, ligand left"."""
    lig = _ligand_heavy_coords(frame)
    if len(lig) == 0:
        return None
    core = _ca_by_residue(frame, core_range, chain_ids=None)
    if not core:
        return False
    core_coords = np.array(list(core.values()), dtype=float)
    dmin = np.linalg.norm(lig[:, None, :] - core_coords[None, :, :], axis=-1).min()
    return bool(dmin <= cutoff)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("reference", type=Path, help="reference toxic-shape PDB")
    p.add_argument("mobile", type=Path, help="PDB to compare against the reference")
    p.add_argument("--core", default="70-88", help="β-core residue range (default 70-88)")
    p.add_argument("--nac", default="60-100", help="inter-chain contact range (default 60-100)")
    p.add_argument("--cutoff", type=float, default=CONTACT_CUTOFF)
    args = p.parse_args()

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    lo_c, hi_c = (int(x) for x in args.core.split("-"))
    lo_n, hi_n = (int(x) for x in args.nac.split("-"))
    ref = load_pdb(args.reference)
    mob = load_pdb(args.mobile)
    m = frame_metrics(mob, ref, (lo_c, hi_c), (lo_n, hi_n), args.cutoff)
    print(f"reference: {args.reference.name}")
    print(f"mobile:    {args.mobile.name}")
    print(f"  β-core Cα RMSD (res {args.core}):      {m['beta_core_rmsd']:.3f} Å")
    print(f"  inter-chain contact Jaccard (res {args.nac}): {m['contact_jaccard']:.3f}")
    print(f"  in toxic basin: {in_toxic_basin(m)}")


if __name__ == "__main__":
    main()
