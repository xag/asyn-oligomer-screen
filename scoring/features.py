"""Stage 2 features. Each function takes a Biopython Structure and returns a
scalar. Combined by a weighted score in classifier.py.

The features are the per-conformer quantities listed in IN_SILICO_PLAN.md:
exposed hydrophobic β-sheet area, membrane-insertion propensity (amphipathic
moment), NAC interface exposure × β-content, contact density (stability
proxy, inversely related to activity), and disordered hydrophobic exposure.
"""
from __future__ import annotations

import math

import numpy as np
import pydssp
from Bio.PDB import SASA
from Bio.PDB.Structure import Structure

HYDROPHOBIC = set("AVLIMFWYPC")

EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
    "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
    "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}

NAC_RESIDUES = set(range(61, 96))

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _first_model(structure: Structure):
    """Use model 0 only — NMR ensembles otherwise multiply every residue
    count by the number of models. Fibril chains within model 0 are kept
    (the contact pattern is the structure)."""
    return next(iter(structure))


def _iter_residues(structure: Structure, chain_ids=None, mask=None):
    for chain in _first_model(structure):
        if chain_ids is not None and chain.id not in chain_ids:
            continue
        for residue in chain:
            if residue.id[0] != " ":
                continue
            name = residue.get_resname()
            if name not in THREE_TO_ONE:
                continue
            if mask is not None and residue.get_full_id() not in mask:
                continue
            yield residue, THREE_TO_ONE[name]


def _count_residues(structure: Structure, chain_ids=None, mask=None) -> int:
    return sum(1 for _ in _iter_residues(structure, chain_ids, mask=mask))


def ordered_core_full_ids(
    structure: Structure,
    contact_radius: float = 8.0,
    min_contacts: int = 6,
    min_size: int = 20,
) -> set | None:
    """Set of residue full_ids forming the structurally ordered core of the
    full assembly. A residue qualifies when its Cα has at least
    `min_contacts` non-sequential Cα atoms (across any chain) within
    `contact_radius` Å.

    Purpose: equalise NMR full-length (with modelled disordered tails) and
    cryo-EM fibril cores (tails absent). Without the mask, NMR features
    are inflated by tail SASA / disorder that has no cryo-EM counterpart.

    Returns None when fewer than `min_size` residues meet the threshold —
    callers should then fall back to no mask (otherwise tiny assemblies
    like SDS-bound monomer 1XQ8 would have all features zeroed out).
    Cached on the structure to avoid recomputation across features."""
    cache_key = (contact_radius, min_contacts, min_size)
    cached = getattr(structure, "_core_mask_cache", None)
    if cached is not None and cached[0] == cache_key:
        return cached[1]

    coords: list = []
    full_ids: list = []
    chains_of: list = []
    for chain in _first_model(structure):
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if residue.get_resname() not in THREE_TO_ONE:
                continue
            if "CA" not in {atom.get_name() for atom in residue}:
                continue
            coords.append(residue["CA"].coord)
            full_ids.append(residue.get_full_id())
            chains_of.append(chain.id)

    if len(coords) < 3:
        structure._core_mask_cache = (cache_key, None)
        return None

    coords_np = np.asarray(coords)
    diffs = coords_np[:, None, :] - coords_np[None, :, :]
    dists = np.linalg.norm(diffs, axis=-1)
    contact_mask = (dists > 0) & (dists < contact_radius)
    # Exclude sequential neighbours, but only when on the same chain (so
    # the residue immediately before/after a chain break still counts).
    n = len(full_ids)
    for i in range(n):
        for j in (i - 1, i + 1):
            if 0 <= j < n and chains_of[i] == chains_of[j]:
                contact_mask[i, j] = False

    counts = contact_mask.sum(axis=1)
    core = {full_ids[i] for i in range(n) if counts[i] >= min_contacts}

    if len(core) < min_size:
        result = None
    else:
        result = core
    structure._core_mask_cache = (cache_key, result)
    return result


def _attach_sasa(structure: Structure) -> None:
    if getattr(structure, "_sasa_done", False):
        return
    model = _first_model(structure)
    SASA.ShrakeRupley().compute(model, level="R")
    structure._sasa_done = True


def dssp_ss_map(structure: Structure) -> dict:
    """Per-residue DSSP secondary-structure code ('H', 'E', '-') keyed by
    residue full_id. Replaces the previous φ/ψ box heuristic — much less
    noisy on NMR ensembles whose dihedrals are dispersed (STATUS.md
    problem 2). pydssp's strand detection relies on the H-bond pattern
    in a 3×3 neighbourhood window, which only resolves correctly when
    all chains of the biological assembly are passed together: a single
    fibril chain has no intra-chain H-bonds, so per-chain DSSP would
    classify every residue as loop. Cached on the structure."""
    cached = getattr(structure, "_dssp_ss_map", None)
    if cached is not None:
        return cached

    coords: list = []
    full_ids: list = []
    for chain in _first_model(structure):
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if residue.get_resname() not in THREE_TO_ONE:
                continue
            try:
                atoms = [
                    residue["N"].coord,
                    residue["CA"].coord,
                    residue["C"].coord,
                    residue["O"].coord,
                ]
            except KeyError:
                continue
            coords.append(atoms)
            full_ids.append(residue.get_full_id())

    out: dict = {}
    if len(coords) >= 3:
        codes = pydssp.assign(np.asarray(coords), out_type="c3")
        for fid, code in zip(full_ids, codes):
            out[fid] = str(code)

    structure._dssp_ss_map = out
    return out


def _is_beta_ss(ss_map: dict, full_id) -> bool:
    return ss_map.get(full_id) == "E"


def _is_helix_ss(ss_map: dict, full_id) -> bool:
    return ss_map.get(full_id) == "H"


def exposed_hydrophobic_beta_sasa(structure: Structure, chain_ids=None, core_mask=None) -> float:
    """Per-residue average SASA on hydrophobic residues in β-conformation.

    SASA is always computed on the full assembly so inter-chain burial is
    accounted for. Accumulation runs over `chain_ids` only (or all
    chains when None), restricted to `core_mask` residues when given so
    that disordered tails do not dilute the mean."""
    _attach_sasa(structure)
    ss_map = dssp_ss_map(structure)
    total = 0.0
    for residue, one in _iter_residues(structure, chain_ids, mask=core_mask):
        if one not in HYDROPHOBIC:
            continue
        if not _is_beta_ss(ss_map, residue.get_full_id()):
            continue
        total += getattr(residue, "sasa", 0.0)
    n = _count_residues(structure, chain_ids, mask=core_mask)
    return total / n if n else 0.0


def membrane_insertion_propensity(structure: Structure, chain_ids=None, core_mask=None, window: int = 11) -> float:
    """Maximum α-helix amphipathic moment over sliding windows. Eisenberg μ_H
    with 100° per residue. Sequence-only; intentionally ignores
    `core_mask` because membrane-binding propensity is a property of the
    full chain sequence (the N-terminal amphipathic region typically
    falls outside the ordered fibril core but still defines membrane
    affinity)."""
    seq = [one for _, one in _iter_residues(structure, chain_ids)]
    if len(seq) < window:
        return 0.0
    delta = math.radians(100.0)
    best = 0.0
    for i in range(len(seq) - window + 1):
        sx = sy = 0.0
        for j in range(window):
            h = EISENBERG.get(seq[i + j], 0.0)
            sx += h * math.cos(j * delta)
            sy += h * math.sin(j * delta)
        mu = math.hypot(sx, sy) / window
        if mu > best:
            best = mu
    return best


def nac_active_score(structure: Structure, chain_ids=None, core_mask=None) -> float:
    """Per-chain mean NAC SASA × β-fraction within NAC. High when residues
    61–95 are simultaneously surface-exposed and in extended/β conformation
    — the state that recruits further monomers. `core_mask` filters NAC
    residues to those that are also structurally ordered, so unstructured
    NAC residues in NMR ensembles don't drag down the β-fraction."""
    _attach_sasa(structure)
    ss_map = dssp_ss_map(structure)
    chain_scores: list[float] = []
    for chain in _first_model(structure):
        if chain_ids is not None and chain.id not in chain_ids:
            continue
        nac_total = 0.0
        nac_beta = 0.0
        nac_count = 0
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if residue.get_resname() not in THREE_TO_ONE:
                continue
            if residue.id[1] not in NAC_RESIDUES:
                continue
            if core_mask is not None and residue.get_full_id() not in core_mask:
                continue
            nac_count += 1
            s = getattr(residue, "sasa", 0.0)
            nac_total += s
            if _is_beta_ss(ss_map, residue.get_full_id()):
                nac_beta += s
        if nac_count == 0 or nac_total == 0:
            continue
        beta_fraction = nac_beta / nac_total
        chain_scores.append((nac_beta / nac_count) * beta_fraction)
    return sum(chain_scores) / len(chain_scores) if chain_scores else 0.0


def _stage3_ligand_heavy_coords(structure: Structure, resname: str = "LIG") -> np.ndarray | None:
    """Heavy-atom coordinates of a Stage-3-appended ligand (HETATM residue
    with resname `"LIG"`, added by `stage3.add_ligand_to_structure`).
    Returns None when no such residue is present — apo and anchor
    structures (which only carry crystallographic waters as HETATM, never
    a `LIG` residue) take the fast path, leaving contact_density
    behaviour unchanged."""
    coords: list = []
    for chain in _first_model(structure):
        for residue in chain:
            if residue.get_resname() != resname:
                continue
            for atom in residue:
                element = (atom.element or "").strip().upper()
                if element in ("H", "D"):
                    continue
                coords.append(atom.coord)
    return np.asarray(coords) if coords else None


def contact_density(structure: Structure, chain_ids=None, core_mask=None, cutoff: float = 8.0) -> float:
    """Mean non-sequential Cα contacts within `cutoff` Å for residues in
    `chain_ids` (and `core_mask` when given). Stability proxy — fibrils
    high, disordered low. Combined into activity with a negative weight.

    Coordinates for neighbour-search come from the full assembly so
    inter-chain contacts count. The mask restricts only the residues
    counted as "of interest" in the numerator and denominator; their
    neighbours can be anywhere in the assembly.

    When a Stage-3-appended ligand (resname `"LIG"`) is present in the
    structure, each ligand heavy atom within `cutoff` of an of-interest
    Cα adds one to that residue's contact count. This lets binders show
    up in the stability proxy proportional to how many ligand heavy
    atoms fall in the residue's first shell — larger / better-packed
    poses raise contact_density more, which (with the −1 classifier
    weight) lowers activity. Apo and anchor structures contain no `LIG`
    residue and are unaffected."""
    all_coords: list = []
    all_full_ids: list = []
    all_chains: list = []
    for chain in _first_model(structure):
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if residue.get_resname() not in THREE_TO_ONE:
                continue
            if "CA" not in {atom.get_name() for atom in residue}:
                continue
            all_coords.append(residue["CA"].coord)
            all_full_ids.append(residue.get_full_id())
            all_chains.append(chain.id)

    if len(all_coords) < 3:
        return 0.0
    all_coords_np = np.asarray(all_coords)

    of_interest: list[int] = []
    for i, fid in enumerate(all_full_ids):
        if chain_ids is not None and all_chains[i] not in chain_ids:
            continue
        if core_mask is not None and fid not in core_mask:
            continue
        of_interest.append(i)

    if not of_interest:
        return 0.0

    sub_coords = all_coords_np[of_interest]
    diffs = sub_coords[:, None, :] - all_coords_np[None, :, :]
    dists = np.linalg.norm(diffs, axis=-1)
    contacts = ((dists > 0) & (dists < cutoff)).sum(axis=1)

    # Subtract within-chain sequential neighbours.
    for k, i in enumerate(of_interest):
        for j in (i - 1, i + 1):
            if 0 <= j < len(all_full_ids) and all_chains[i] == all_chains[j]:
                if 0 < np.linalg.norm(all_coords_np[i] - all_coords_np[j]) < cutoff:
                    contacts[k] -= 1

    # Stage-3 ligand contribution. Each LIG heavy atom within cutoff of an
    # of-interest Cα adds one contact to that residue. Apo / anchor
    # structures have no LIG residue and skip this entirely.
    ligand_coords = _stage3_ligand_heavy_coords(structure)
    if ligand_coords is not None:
        of_interest_coords = all_coords_np[of_interest]
        lig_diffs = of_interest_coords[:, None, :] - ligand_coords[None, :, :]
        lig_dists = np.linalg.norm(lig_diffs, axis=-1)
        contacts = contacts + (lig_dists < cutoff).sum(axis=1).astype(contacts.dtype)

    return float(contacts.mean())


def disordered_hydrophobic_exposure(structure: Structure, chain_ids=None, core_mask=None) -> float:
    """Per-residue average exposed hydrophobic SASA on residues outside
    regular secondary structure. Captures sticky disorder — exposed greasy
    patches the structure cannot bury. With `core_mask`, only ordered-core
    residues that *still* fail DSSP's H/E test are counted (so 2N0A's
    disordered tails, which have no cryo-EM counterpart in modelling
    coverage, drop out)."""
    _attach_sasa(structure)
    ss_map = dssp_ss_map(structure)
    total = 0.0
    for residue, one in _iter_residues(structure, chain_ids, mask=core_mask):
        if one not in HYDROPHOBIC:
            continue
        fid = residue.get_full_id()
        if _is_beta_ss(ss_map, fid) or _is_helix_ss(ss_map, fid):
            continue
        total += getattr(residue, "sasa", 0.0)
    n = _count_residues(structure, chain_ids, mask=core_mask)
    return total / n if n else 0.0


FEATURES = {
    "exposed_hydrophobic_beta_sasa": exposed_hydrophobic_beta_sasa,
    "membrane_insertion_propensity": membrane_insertion_propensity,
    "nac_active_score": nac_active_score,
    "contact_density": contact_density,
    "disordered_hydrophobic_exposure": disordered_hydrophobic_exposure,
}
