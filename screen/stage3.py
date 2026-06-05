"""Stage 3 perturbation engine — single (molecule, conformer) pair.

For one vicinity molecule M from data/vicinity_molecules.js and one anchor
PDB id C, this script:

  1. Loads C as its REMARK 350 biological assembly.
  2. Picks the inner (most-buried) chain and computes its Stage 2 features
     (apo baseline).
  3. Prepares an AutoDock-ready PDBQT receptor (mk_prepare_receptor) and a
     PDBQT ligand (rdkit + meeko).
  4. Runs AutoDock Vina 1.2 over a box centred on the inner chain.
  5. Builds a top-pose complex PDB by appending the rank-1 docked pose
     to the assembly as a HETATM residue on its own chain (for
     inspection / back-compat).
  6. Recomputes the Stage 2 features on every Vina pose (each pose is
     added to a fresh copy of the assembly). SASA picks up ligand
     occlusion; contact_density picks up ligand heavy-atom contacts
     inside the 8 Å cutoff.
  7. Aggregates per-pose Δactivity values with Boltzmann weights
     exp(-aff_i / RT) / Z at T=300 K, so a clearly best pose dominates
     and a spread of similar affinities averages over poses.
  8. Writes a report with ΔP(active | M, C) anchored against the Stage 2
     anchor z-table — both the top-pose Δactivity and the
     Boltzmann-weighted Δactivity are emitted.

This is the prototype called for by IN_SILICO_PLAN.md. It deliberately skips MD
relaxation: a static docked complex catches the SASA-occlusion and
contact-density components of the perturbation but not conformational
rearrangement. MD comes later.

Usage:
    python stage3.py <molecule_id> <anchor_pdb> [--exhaustiveness N]
    python stage3.py curcumin 6PEO
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, Select
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from Bio.PDB import PDBParser as _PDBParser

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scoring"))

from adduct_score import LIGAND_REACTIVITY, aspr_score
from assembly import inner_chain_ids, load_assembly
from classifier import WEIGHTS
from features import FEATURES, ordered_core_full_ids

RESULTS = ROOT / "results" / "stage3"


def _tool_path(name: str, bundled: Path) -> str:
    """Resolve an external docking tool across platforms so the `dock` chunk runs
    on a volunteer's Linux/Mac box, not only the Windows dev machine. Order:
    1. explicit override env var (e.g. ``VINA_BIN`` / ``MK_PREPARE_RECEPTOR_BIN``);
    2. a console script next to the running Python (meeko installs
       ``mk_prepare_receptor`` into the venv's bin/Scripts dir);
    3. anything on ``PATH`` (``conda install vina`` / a downloaded binary);
    4. the bundled Windows binary as a last-resort fallback.
    """
    env = os.environ.get(name.upper().replace("-", "_") + "_BIN")
    if env:
        return env
    exe = name + (".exe" if os.name == "nt" else "")
    beside = Path(sys.executable).parent / exe
    if beside.exists():
        return str(beside)
    return shutil.which(name) or shutil.which(exe) or str(bundled)


def _vina_bin() -> str:
    return _tool_path("vina", ROOT / "bin" / "vina.exe")


def _mk_prepare_receptor() -> str:
    return _tool_path("mk_prepare_receptor", ROOT / ".venv" / "Scripts" / "mk_prepare_receptor.exe")
VICINITY_JS = ROOT / "data" / "vicinity_molecules.js"

# -----------------------------------------------------------------------------
# Vicinity molecule loader — parses data/vicinity_molecules.js with regex. The
# JS file is a hand-authored list of object literals; full JS parsing would
# need a Node round-trip. The regex is tolerant enough for the seeded fields.
# -----------------------------------------------------------------------------

_ENTRY_RE = re.compile(
    r"\{\s*id:\s*'([^']+)',(.*?)\n\s*\},",
    re.DOTALL,
)
_STR_FIELD_RE = re.compile(r"^\s*(\w+):\s*'([^']*)'", re.MULTILINE)
_NULL_FIELD_RE = re.compile(r"^\s*(\w+):\s*null", re.MULTILINE)


def load_vicinity_molecule(mol_id: str) -> dict:
    text = VICINITY_JS.read_text(encoding="utf-8")
    for m in _ENTRY_RE.finditer(text):
        if m.group(1) != mol_id:
            continue
        body = m.group(2)
        entry: dict = {"id": mol_id}
        for fm in _STR_FIELD_RE.finditer(body):
            entry.setdefault(fm.group(1), fm.group(2))
        for fm in _NULL_FIELD_RE.finditer(body):
            entry.setdefault(fm.group(1), None)
        return entry
    raise KeyError(f"vicinity molecule {mol_id!r} not found in {VICINITY_JS}")


# -----------------------------------------------------------------------------
# Ligand prep — SMILES → rdkit 3D conformer → meeko PDBQT.
# -----------------------------------------------------------------------------

def prepare_ligand(smiles: str, out_pdbqt: Path, random_seed: int = 42) -> Path:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"rdkit could not parse SMILES {smiles!r}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=random_seed) != 0:
        raise RuntimeError(f"3D embedding failed for SMILES {smiles!r}")
    AllChem.MMFFOptimizeMolecule(mol)
    setups = MoleculePreparation().prepare(mol)
    pdbqt, ok, errs = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(f"meeko ligand PDBQT failed: {errs}")
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    out_pdbqt.write_text(pdbqt, encoding="utf-8")
    return out_pdbqt


# -----------------------------------------------------------------------------
# Receptor prep — assembly PDB → single-letter chains, AA-only → PDBQT via the
# meeko mk_prepare_receptor CLI.
# -----------------------------------------------------------------------------

_AA_3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


class _AASelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " " and residue.get_resname() in _AA_3


def _renumber_chains(structure: Structure) -> None:
    """Rename chains to single ASCII letters. PDB chain id is one character;
    the assembly builder may produce IDs like 'A_2' which break PDB writers.
    Detach + re-add with new ids in order. Coordinates are not touched."""
    model = next(iter(structure))
    available = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    chains = list(model)
    for ch in chains:
        model.detach_child(ch.id)
    for ch, new_id in zip(chains, available):
        ch.id = new_id
        model.add(ch)


def write_apo_pdb(structure: Structure, out_pdb: Path) -> Path:
    _renumber_chains(structure)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_pdb), select=_AASelect())
    return out_pdb


def prepare_receptor(apo_pdb: Path, out_pdbqt: Path) -> Path:
    cmd = [_mk_prepare_receptor(), "--read_pdb", str(apo_pdb), "-p", str(out_pdbqt)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(
            f"mk_prepare_receptor failed (rc={res.returncode}):\n"
            f"STDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    return out_pdbqt


# -----------------------------------------------------------------------------
# Docking box — centred on the inner chain centroid, sized to the chain
# bounding box + padding but capped at Vina's recommended 30 Å per side.
# -----------------------------------------------------------------------------

def docking_box_for_chain(
    structure: Structure,
    chain_id: str,
    padding: float = 6.0,
    max_side: float = 30.0,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    coords = []
    for chain in next(iter(structure)):
        if chain.id != chain_id:
            continue
        for residue in chain:
            for atom in residue:
                coords.append(atom.coord)
    if not coords:
        raise ValueError(f"chain {chain_id!r} has no atoms")
    arr = np.asarray(coords)
    lo = arr.min(axis=0) - padding
    hi = arr.max(axis=0) + padding
    center = (lo + hi) / 2.0
    size = np.minimum(hi - lo, max_side)
    return tuple(float(x) for x in center), tuple(float(x) for x in size)


# -----------------------------------------------------------------------------
# Oligomer support: load from local file, per-chain mean features, multi-chain box
# -----------------------------------------------------------------------------

def load_structure_from_file(pdb_path: Path):
    """Load a PDB from a local file without BIOMT assembly expansion."""
    return _PDBParser(QUIET=True).get_structure(pdb_path.stem, str(pdb_path))


def compute_mean_features_on(
    structure, chain_ids: list[str], use_core_mask: bool = True
) -> dict[str, float]:
    """Score each chain individually, return mean feature vector. Matches the
    all-chain-mean mode used by score_oligomer.py."""
    per_chain = [compute_features_on(structure, [cid], use_core_mask) for cid in chain_ids]
    return {k: float(np.mean([f[k] for f in per_chain])) for k in per_chain[0]}


def docking_box_for_chains(
    structure,
    chain_ids: list[str],
    padding: float = 6.0,
    max_side: float = 30.0,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Docking box centred on the centroid of all listed chains."""
    ids = set(chain_ids)
    coords = []
    for chain in next(iter(structure)):
        if chain.id not in ids:
            continue
        for residue in chain:
            for atom in residue:
                coords.append(atom.coord)
    if not coords:
        raise ValueError(f"chains {chain_ids!r} have no atoms")
    arr = np.asarray(coords)
    lo = arr.min(axis=0) - padding
    hi = arr.max(axis=0) + padding
    center = (lo + hi) / 2.0
    size = np.minimum(hi - lo, max_side)
    return tuple(float(x) for x in center), tuple(float(x) for x in size)


# Residue-number range of the α-syn NAC region used to anchor the docking
# box for generated oligomers. The Fusco β-core (70-88) lives inside this
# range; using the broader 60-100 window captures any supported core range
# (65-83, 70-88, 73-91) and is robust to exact boundary choices.
_NAC_CORE_RESNUMS = set(range(60, 101))


def docking_box_for_nac_core(
    structure,
    padding: float = 6.0,
    max_side: float = 30.0,
    min_residues: int = 6,
) -> tuple[tuple[float, float, float], tuple[float, float, float]] | None:
    """Docking box centred on Cα of NAC-region residues (60-100) across all chains.

    For α-syn oligomers built from the Fusco topology, the β-core is always
    within residues 60-100 by construction. This is more reliable than
    DSSP-based selection which can pick up tail β-structure formed during
    OBC2 relaxation. Returns None if fewer than min_residues qualify."""
    coords = []
    for chain in next(iter(structure)):
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if residue.id[1] not in _NAC_CORE_RESNUMS:
                continue
            if "CA" in {a.get_name() for a in residue}:
                coords.append(residue["CA"].coord)
    if len(coords) < min_residues:
        return None
    arr = np.asarray(coords)
    lo = arr.min(axis=0) - padding
    hi = arr.max(axis=0) + padding
    center = (lo + hi) / 2.0
    size = np.minimum(hi - lo, max_side)
    return tuple(float(x) for x in center), tuple(float(x) for x in size)


# -----------------------------------------------------------------------------
# Vina docking via the standalone Windows binary.
# -----------------------------------------------------------------------------

def run_vina(
    receptor: Path,
    ligand: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    out_pdbqt: Path,
    exhaustiveness: int = 8,
    n_poses: int = 5,
    seed: int = 42,
) -> tuple[Path, list[float], str]:
    cmd = [
        _vina_bin(),
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--out", str(out_pdbqt),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(n_poses),
        "--seed", str(seed),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if res.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(
            f"vina failed (rc={res.returncode}):\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    return out_pdbqt, _parse_affinities(res.stdout), res.stdout


def _parse_affinities(log: str) -> list[float]:
    """All mode affinities from Vina's stdout mode table, in rank order.
    Format per data row: "   <rank>       <aff>     <rmsd_lb>    <rmsd_ub>"."""
    affs: list[float] = []
    for line in log.splitlines():
        toks = line.strip().split()
        if len(toks) >= 4 and toks[0].isdigit():
            try:
                affs.append(float(toks[1]))
            except ValueError:
                continue
    return affs


def _parse_top_affinity(log: str) -> float:
    affs = _parse_affinities(log)
    return affs[0] if affs else float("nan")


# -----------------------------------------------------------------------------
# Pose extraction — read MODEL 1 from the vina output PDBQT as heavy atoms.
# -----------------------------------------------------------------------------

_AD_ELEM = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "O": "O", "OA": "O",
    "S": "S", "SA": "S", "HD": "H", "H": "H", "F": "F", "Cl": "Cl",
    "Br": "Br", "I": "I", "P": "P",
}


def parse_pdbqt_all_poses(pdbqt_path: Path) -> list[list[tuple[str, np.ndarray, str]]]:
    """All MODEL blocks from a Vina output PDBQT, in rank order."""
    poses: list[list[tuple[str, np.ndarray, str]]] = []
    current: list[tuple[str, np.ndarray, str]] | None = None
    for line in pdbqt_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL"):
            current = []
            continue
        if line.startswith("ENDMDL"):
            if current is not None:
                poses.append(current)
                current = None
            continue
        if current is None:
            continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        name = line[12:16].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        ad_type = line[77:].strip() if len(line) > 77 else ""
        elem = _AD_ELEM.get(ad_type, name[0])
        current.append((name, np.array([x, y, z]), elem))
    if not poses:
        raise RuntimeError(f"no MODEL blocks in {pdbqt_path}")
    return poses


def parse_pdbqt_top_pose(pdbqt_path: Path) -> list[tuple[str, np.ndarray, str]]:
    return parse_pdbqt_all_poses(pdbqt_path)[0]


def add_ligand_to_structure(
    structure: Structure,
    ligand_atoms: list[tuple[str, np.ndarray, str]],
    ligand_resname: str = "LIG",
    ligand_chain_id: str = "Z",
) -> None:
    model = next(iter(structure))
    used = {ch.id for ch in model}
    chain_id = ligand_chain_id
    suffix = 0
    while chain_id in used:
        suffix += 1
        chain_id = f"{ligand_chain_id}{suffix}"
    chain = Chain(chain_id)
    res = Residue(("H_LIG", 1, " "), ligand_resname, "")
    # Biopython residues key atoms by name and reject duplicates. The Vina
    # output PDBQT keeps the original PDB-style names (mostly element
    # letters), so multiple "C" / "O" / "N" atoms collide. Re-number per
    # element to keep names unique while still parseable.
    elem_counters: dict[str, int] = {}
    for i, (_orig_name, coord, elem) in enumerate(ligand_atoms, start=1):
        elem_counters[elem] = elem_counters.get(elem, 0) + 1
        unique_name = f"{elem}{elem_counters[elem]}"
        atom = Atom(
            name=unique_name,
            coord=coord,
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=f"{unique_name:<4}",
            serial_number=i,
            element=elem,
        )
        res.add(atom)
    chain.add(res)
    model.add(chain)
    for atom in structure.get_atoms():
        atom.full_id = atom.get_full_id()


# -----------------------------------------------------------------------------
# Stage 2 feature recomputation on an arbitrary structure (apo or complex).
# -----------------------------------------------------------------------------

def compute_features_on(structure: Structure, chain_ids: list[str], use_core_mask: bool = True) -> dict[str, float]:
    core_mask = ordered_core_full_ids(structure) if use_core_mask else None
    feats: dict[str, float] = {}
    for name, fn in FEATURES.items():
        try:
            feats[name] = float(fn(structure, chain_ids=chain_ids, core_mask=core_mask))
        except Exception as exc:
            print(f"    feature {name} failed: {exc}")
            feats[name] = float("nan")
    return feats


# -----------------------------------------------------------------------------
# Boltzmann pose weighting at T=300 K. RT in kcal/mol.
# -----------------------------------------------------------------------------

RT_KCAL_300K = 1.987204e-3 * 300.0  # ~0.59616 kcal/mol

# Absolute-affinity gate. Ligands whose top-pose affinity is weaker
# (less negative) than this threshold get exponentially damped. Set
# near the "drug-like binder" boundary so curcumin / quercetin pass
# unchanged and metabolite-scale binders (methylglyoxal at ~-2.9) are
# strongly suppressed. Damping factor: min(1, exp((thr - aff_top)/RT)).
AFFINITY_GATE_THRESHOLD_KCAL = -6.0


def affinity_gate(aff_top: float, threshold: float = AFFINITY_GATE_THRESHOLD_KCAL,
                  temperature_kcal: float = RT_KCAL_300K) -> float:
    """Damping factor in [0, 1] for a ligand's pose ensemble.
    Returns 1 when aff_top is at or below `threshold` (real binder)
    and exp((threshold - aff_top) / RT) otherwise. Orthogonal to the
    intra-ligand softmax weights — this penalises ligands whose best
    pose is itself weak."""
    if not np.isfinite(aff_top):
        return 0.0
    if aff_top <= threshold:
        return 1.0
    return float(np.exp((threshold - aff_top) / temperature_kcal))


def softmax_weights(affinities: list[float], temperature_kcal: float = RT_KCAL_300K) -> list[float]:
    """Boltzmann weights from binding affinities (kcal/mol; lower is tighter).
    w_i = exp(-aff_i / RT) / Z. Shifted by the max logit for numerical
    stability so positive-affinity (clashing) poses collapse to ~0 weight."""
    if not affinities:
        return []
    arr = np.asarray(affinities, dtype=float)
    logits = -arr / temperature_kcal
    logits -= np.max(logits)
    expv = np.exp(logits)
    s = float(expv.sum())
    if s == 0.0:
        return [1.0 / len(affinities)] * len(affinities)
    return [float(x) for x in (expv / s)]


# -----------------------------------------------------------------------------
# ΔP(active) — activity scores under the Stage 2 anchor z-axis, so the apo
# baseline matches anchor_scores.csv to within numerical noise.
# -----------------------------------------------------------------------------

def delta_activity(apo_feats: dict, complex_feats: dict, anchor_df: pd.DataFrame) -> dict:
    z_apo: dict[str, float] = {}
    z_complex: dict[str, float] = {}
    for col, w in WEIGHTS.items():
        mu = anchor_df[col].mean()
        sd = anchor_df[col].std(ddof=0)
        if sd == 0:
            z_apo[col] = 0.0
            z_complex[col] = 0.0
        else:
            z_apo[col] = float((apo_feats[col] - mu) / sd)
            z_complex[col] = float((complex_feats[col] - mu) / sd)
    activity_apo = sum(z_apo[c] * w for c, w in WEIGHTS.items())
    activity_complex = sum(z_complex[c] * w for c, w in WEIGHTS.items())
    return {
        "activity_apo": float(activity_apo),
        "activity_complex": float(activity_complex),
        "delta_activity": float(activity_complex - activity_apo),
        "z_apo": z_apo,
        "z_complex": z_complex,
    }


def score_all_poses(
    pdb_id: str | None,
    apo_chain: str | None,
    docked_pdbqt: Path,
    affinities: list[float],
    anchor_df: pd.DataFrame,
    apo_feats: dict | None = None,
    *,
    load_struct_fn=None,
    mean_chain_ids: list[str] | None = None,
) -> dict:
    """Score every Vina pose under Stage 2 features and aggregate with
    Boltzmann weights. Each pose is built by adding its ligand atoms to a
    fresh copy of the assembly so the receptor is identical across poses
    (we do not relax) — only the ligand-aware features (SASA occlusion,
    contact_density) vary between poses.

    For oligomer targets, pass load_struct_fn (callable → Structure from file)
    and mean_chain_ids (all protein chains) to use per-chain-mean scoring
    instead of the single inner-chain mode used for deposited anchors."""
    _load = load_struct_fn if load_struct_fn is not None else (lambda: load_assembly(pdb_id))
    if mean_chain_ids is not None:
        def _compute_feats(struct):
            return compute_mean_features_on(struct, mean_chain_ids)
    else:
        def _compute_feats(struct):
            return compute_features_on(struct, chain_ids=[apo_chain])

    poses = parse_pdbqt_all_poses(docked_pdbqt)
    if len(poses) < len(affinities):
        # Vina occasionally drops trailing clashing poses (positive affinity)
        # from the output PDBQT while still listing them in the log table.
        # Truncate the affinity tail to match — the dropped poses are the
        # worst-ranked ones and would carry vanishing Boltzmann weight anyway.
        print(
            f"  note: log lists {len(affinities)} affinities but PDBQT has "
            f"{len(poses)} poses (dropped: {affinities[len(poses):]}); truncating."
        )
        affinities = affinities[: len(poses)]
    elif len(poses) > len(affinities):
        raise RuntimeError(
            f"more poses ({len(poses)}) than affinities ({len(affinities)}) — log parse failure"
        )

    if apo_feats is None:
        apo_feats = _compute_feats(_load())

    weights = softmax_weights(affinities)
    pose_records: list[dict] = []
    for rank, (ligand_atoms, aff, w) in enumerate(zip(poses, affinities, weights), start=1):
        complex_struct = _load()
        add_ligand_to_structure(complex_struct, ligand_atoms)
        complex_feats = _compute_feats(complex_struct)
        d = delta_activity(apo_feats, complex_feats, anchor_df)
        pose_records.append({
            "rank": rank,
            "affinity_kcal_per_mol": float(aff),
            "weight": float(w),
            "n_ligand_atoms": len(ligand_atoms),
            "complex_features": complex_feats,
            "activity_complex": d["activity_complex"],
            "delta_activity": d["delta_activity"],
        })

    activity_apo = pose_records[0]["activity_complex"] - pose_records[0]["delta_activity"]
    activity_complex_weighted = sum(p["weight"] * p["activity_complex"] for p in pose_records)
    delta_activity_weighted = sum(p["weight"] * p["delta_activity"] for p in pose_records)
    aff_top = float(pose_records[0]["affinity_kcal_per_mol"])
    gate = affinity_gate(aff_top)
    delta_activity_gated = float(gate * delta_activity_weighted)
    top = pose_records[0]
    return {
        "apo_features": apo_feats,
        "activity_apo": float(activity_apo),
        "n_poses": len(pose_records),
        "poses": pose_records,
        "complex_features_top": top["complex_features"],
        "delta_features_top": {k: top["complex_features"][k] - apo_feats[k] for k in apo_feats},
        "activity_complex_top": float(top["activity_complex"]),
        "delta_activity_top": float(top["delta_activity"]),
        "activity_complex_weighted": float(activity_complex_weighted),
        "delta_activity_weighted": float(delta_activity_weighted),
        "affinity_gate_threshold_kcal_per_mol": float(AFFINITY_GATE_THRESHOLD_KCAL),
        "affinity_gate": float(gate),
        "delta_activity_gated": delta_activity_gated,
    }


# -----------------------------------------------------------------------------
# Chain id matching across renumbered/original assembly copies.
# -----------------------------------------------------------------------------

def _match_inner_chain(dock_struct: Structure, apo_struct: Structure, apo_chain_id: str) -> str:
    apo_coords = [
        r["CA"].coord
        for ch in next(iter(apo_struct))
        if ch.id == apo_chain_id
        for r in ch
        if "CA" in {a.get_name() for a in r}
    ]
    if not apo_coords:
        raise ValueError(f"no Cα in chain {apo_chain_id!r} of apo")
    apo_center = np.mean(apo_coords, axis=0)
    best = None
    best_dist = float("inf")
    for ch in next(iter(dock_struct)):
        coords = [r["CA"].coord for r in ch if "CA" in {a.get_name() for a in r}]
        if not coords:
            continue
        d = float(np.linalg.norm(np.mean(coords, axis=0) - apo_center))
        if d < best_dist:
            best_dist = d
            best = ch.id
    if best is None:
        raise RuntimeError("could not match inner chain to dock structure")
    return best


# -----------------------------------------------------------------------------
# Top-level: dock one (molecule, anchor) pair and emit a report.
# -----------------------------------------------------------------------------

def perturb(
    mol_id: str,
    pdb_id: str,
    exhaustiveness: int = 8,
    n_poses: int = 5,
    padding: float = 6.0,
    seed: int = 42,
) -> dict:
    RESULTS.mkdir(parents=True, exist_ok=True)
    pair_tag = f"{mol_id}_{pdb_id.lower()}"
    out = RESULTS

    print(f"\n=== stage 3 perturb: {mol_id} × {pdb_id} ===")

    # 1) Apo features under Stage 2 conventions
    apo_struct = load_assembly(pdb_id)
    apo_chain = inner_chain_ids(apo_struct, top_k=1)[0]
    print(f"  apo inner chain: {apo_chain!r}")
    apo_feats = compute_features_on(apo_struct, chain_ids=[apo_chain])

    # 2) Renumbered assembly for docking I/O (PDB single-letter chains)
    dock_struct = load_assembly(pdb_id)
    apo_pdb = out / f"{pair_tag}_apo.pdb"
    write_apo_pdb(dock_struct, apo_pdb)
    dock_inner = _match_inner_chain(dock_struct, apo_struct, apo_chain)
    print(f"  dock inner chain: {dock_inner!r}")

    # 3) Receptor PDBQT
    receptor_pdbqt = out / f"{pair_tag}_receptor.pdbqt"
    prepare_receptor(apo_pdb, receptor_pdbqt)
    print(f"  receptor.pdbqt: {receptor_pdbqt.name}")

    # 4) Ligand PDBQT
    mol = load_vicinity_molecule(mol_id)
    smiles = mol.get("smiles")
    if not smiles:
        raise ValueError(
            f"vicinity molecule {mol_id!r} has smiles=null in vicinity_molecules.js. "
            f"Stage 3 needs a structure — fill in the SMILES before docking this entry."
        )
    ligand_pdbqt = out / f"{pair_tag}_ligand.pdbqt"
    prepare_ligand(smiles, ligand_pdbqt, random_seed=seed)
    print(f"  ligand SMILES: {smiles}")
    print(f"  ligand.pdbqt: {ligand_pdbqt.name}")

    # 5) Docking box
    center, size = docking_box_for_chain(dock_struct, dock_inner, padding=padding)
    print(f"  box center: {tuple(round(x, 2) for x in center)} size: {tuple(round(x, 2) for x in size)}")

    # 6) Dock
    docked_pdbqt = out / f"{pair_tag}_docked.pdbqt"
    _, affinities, vina_log = run_vina(
        receptor_pdbqt, ligand_pdbqt, center, size, docked_pdbqt,
        exhaustiveness=exhaustiveness, n_poses=n_poses, seed=seed,
    )
    (out / f"{pair_tag}_vina.log").write_text(vina_log, encoding="utf-8")
    print(f"  vina affinities (kcal/mol): {[f'{a:.2f}' for a in affinities]}")

    # 7) Top-pose complex artifact (for inspection / back-compat). Multi-pose
    #    scoring builds its own per-pose complexes in memory; this PDB stays
    #    the top-pose representative.
    complex_struct = load_assembly(pdb_id)
    top_ligand_atoms = parse_pdbqt_top_pose(docked_pdbqt)
    add_ligand_to_structure(complex_struct, top_ligand_atoms)
    complex_pdb = out / f"{pair_tag}_complex.pdb"
    io = PDBIO()
    io.set_structure(complex_struct)
    io.save(str(complex_pdb))

    # 8) Per-pose features + Boltzmann-weighted aggregate
    anchor_df = pd.read_csv(ROOT / "results" / "anchor_features.csv")
    multi = score_all_poses(
        pdb_id, apo_chain, docked_pdbqt, affinities, anchor_df, apo_feats=apo_feats,
    )

    # 9) Covalent / adduct channel. Independent of Vina; scored on the
    #    inner chain alone for consistency with the apo features above.
    aspr = aspr_score(apo_struct, mol_id, chain_ids=[apo_chain])

    report = {
        "molecule_id": mol_id,
        "molecule_name": mol.get("name"),
        "anchor_pdb_id": pdb_id,
        "inner_chain_apo": apo_chain,
        "inner_chain_dock": dock_inner,
        "smiles": smiles,
        "vina_top_affinity_kcal_per_mol": float(affinities[0]),
        "vina_all_affinities_kcal_per_mol": [float(a) for a in affinities],
        "n_ligand_atoms": len(top_ligand_atoms),
        "n_poses": multi["n_poses"],
        "apo_features": apo_feats,
        "complex_features_top": multi["complex_features_top"],
        "delta_features_top": multi["delta_features_top"],
        "activity_apo": multi["activity_apo"],
        "activity_complex_top": multi["activity_complex_top"],
        "delta_activity_top": multi["delta_activity_top"],
        "activity_complex_weighted": multi["activity_complex_weighted"],
        "delta_activity_weighted": multi["delta_activity_weighted"],
        "affinity_gate_threshold_kcal_per_mol": multi["affinity_gate_threshold_kcal_per_mol"],
        "affinity_gate": multi["affinity_gate"],
        "delta_activity_gated": multi["delta_activity_gated"],
        "aspr_score": float(aspr),
        "aspr_reactive": mol_id in LIGAND_REACTIVITY,
        "poses": multi["poses"],
    }
    (out / f"{pair_tag}_report.json").write_text(json.dumps(report, indent=2))
    print(f"  apo activity              {report['activity_apo']:+.3f}")
    print(f"  complex activity (top)    {report['activity_complex_top']:+.3f}")
    print(f"  complex activity (w'tied) {report['activity_complex_weighted']:+.3f}")
    direction_top = "protective" if report["delta_activity_top"] < 0 else "harm-leaning"
    direction_w = "protective" if report["delta_activity_weighted"] < 0 else "harm-leaning"
    print(f"  delta (top)               {report['delta_activity_top']:+.3f}  ({direction_top})")
    print(f"  delta (Boltzmann)         {report['delta_activity_weighted']:+.3f}  ({direction_w})")
    print(
        f"  delta (gated, thr={report['affinity_gate_threshold_kcal_per_mol']:.1f}) "
        f"{report['delta_activity_gated']:+.3f}  (gate={report['affinity_gate']:.3f})"
    )
    if report["aspr_reactive"]:
        print(f"  aspr_score (covalent)     {report['aspr_score']:+.3f}  (harm-leaning if > 0)")
    return report


# -----------------------------------------------------------------------------
# Oligomer-mode perturb: generated PDB as receptor, all-chain mean Δactivity.
# -----------------------------------------------------------------------------

def perturb_oligomer(
    mol_id: str,
    oligo_pdb: Path,
    exhaustiveness: int = 8,
    n_poses: int = 5,
    padding: float = 6.0,
    seed: int = 42,
) -> dict:
    """Dock mol_id against a generated oligomer PDB and score ΔP(active).

    Differs from perturb() in two ways:
      - receptor loaded from a local file (no RCSB fetch, no BIOMT)
      - apo/complex activity computed as per-chain mean over all protein chains
        (matching score_oligomer.py default), not single inner-chain
    """
    oligo_pdb = Path(oligo_pdb)
    RESULTS.mkdir(parents=True, exist_ok=True)
    pair_tag = f"{mol_id}_{oligo_pdb.stem}"
    out = RESULTS

    print(f"\n=== stage 3 perturb (oligomer): {mol_id} × {oligo_pdb.name} ===")

    # 1) Load structure, identify protein chains
    apo_struct = load_structure_from_file(oligo_pdb)
    model0 = next(iter(apo_struct))
    protein_chain_ids = [
        ch.id for ch in model0
        if any(r.id[0] == " " and r.get_resname() in _AA_3 for r in ch)
    ]
    if not protein_chain_ids:
        raise ValueError(f"no standard-AA chains found in {oligo_pdb}")
    print(f"  protein chains: {protein_chain_ids}")

    # 2) Apo features: per-chain mean (consistent with score_oligomer.py default)
    apo_feats = compute_mean_features_on(apo_struct, protein_chain_ids)
    anchor_df = pd.read_csv(ROOT / "results" / "anchor_features.csv")

    # 3) Write receptor PDB (renumbers chains to single-letter IDs for Vina)
    dock_struct = load_structure_from_file(oligo_pdb)
    apo_pdb = out / f"{pair_tag}_apo.pdb"
    write_apo_pdb(dock_struct, apo_pdb)
    # Grab chain IDs after renaming (for box centering)
    dock_chain_ids = [ch.id for ch in next(iter(dock_struct))]
    print(f"  dock chains (after renumber): {dock_chain_ids}")

    # 4) Receptor PDBQT
    receptor_pdbqt = out / f"{pair_tag}_receptor.pdbqt"
    prepare_receptor(apo_pdb, receptor_pdbqt)
    print(f"  receptor.pdbqt: {receptor_pdbqt.name}")

    # 5) Ligand PDBQT
    mol = load_vicinity_molecule(mol_id)
    smiles = mol.get("smiles")
    if not smiles:
        raise ValueError(
            f"vicinity molecule {mol_id!r} has smiles=null in vicinity_molecules.js. "
            f"Stage 3 needs a structure — fill in the SMILES before docking this entry."
        )
    ligand_pdbqt = out / f"{pair_tag}_ligand.pdbqt"
    prepare_ligand(smiles, ligand_pdbqt, random_seed=seed)
    print(f"  ligand SMILES: {smiles}")

    # 6) Docking box: NAC-region residues (60-100) are always the beta-core
    #    in Fusco-topology oligomers; fall back to all-chains when the
    #    residue numbering doesn't match (e.g. non-standard constructs).
    box = docking_box_for_nac_core(dock_struct, padding=padding)
    if box is not None:
        center, size = box
        print(f"  box source: NAC-core residues 60-100")
    else:
        center, size = docking_box_for_chains(dock_struct, dock_chain_ids, padding=padding)
        print(f"  box source: all chains (NAC range not found)")
    print(f"  box center: {tuple(round(x, 2) for x in center)} size: {tuple(round(x, 2) for x in size)}")

    # 7) Dock
    docked_pdbqt = out / f"{pair_tag}_docked.pdbqt"
    _, affinities, vina_log = run_vina(
        receptor_pdbqt, ligand_pdbqt, center, size, docked_pdbqt,
        exhaustiveness=exhaustiveness, n_poses=n_poses, seed=seed,
    )
    (out / f"{pair_tag}_vina.log").write_text(vina_log, encoding="utf-8")
    print(f"  vina affinities (kcal/mol): {[f'{a:.2f}' for a in affinities]}")

    # 8) Top-pose complex artifact for inspection
    complex_struct = load_structure_from_file(oligo_pdb)
    top_ligand_atoms = parse_pdbqt_top_pose(docked_pdbqt)
    add_ligand_to_structure(complex_struct, top_ligand_atoms)
    complex_pdb = out / f"{pair_tag}_complex.pdb"
    io = PDBIO()
    io.set_structure(complex_struct)
    io.save(str(complex_pdb))

    # 9) Per-pose scoring with all-chain mean Δactivity
    load_fn = lambda: load_structure_from_file(oligo_pdb)  # noqa: E731
    multi = score_all_poses(
        pdb_id=None, apo_chain=None,
        docked_pdbqt=docked_pdbqt, affinities=affinities, anchor_df=anchor_df,
        apo_feats=apo_feats,
        load_struct_fn=load_fn,
        mean_chain_ids=protein_chain_ids,
    )

    # 10) Covalent / adduct channel — receptor accessibility × ligand
    #     reactivity, per-chain mean across all protein chains. Independent
    #     of Vina pose.
    aspr = aspr_score(apo_struct, mol_id, chain_ids=protein_chain_ids)

    report = {
        "molecule_id": mol_id,
        "molecule_name": mol.get("name"),
        "receptor_type": "oligomer",
        "receptor_pdb": str(oligo_pdb),
        "receptor_stem": oligo_pdb.stem,
        "protein_chain_ids": protein_chain_ids,
        "smiles": smiles,
        "vina_top_affinity_kcal_per_mol": float(affinities[0]),
        "vina_all_affinities_kcal_per_mol": [float(a) for a in affinities],
        "n_ligand_atoms": len(top_ligand_atoms),
        "n_poses": multi["n_poses"],
        "apo_features": apo_feats,
        "complex_features_top": multi["complex_features_top"],
        "delta_features_top": multi["delta_features_top"],
        "activity_apo": multi["activity_apo"],
        "activity_complex_top": multi["activity_complex_top"],
        "delta_activity_top": multi["delta_activity_top"],
        "activity_complex_weighted": multi["activity_complex_weighted"],
        "delta_activity_weighted": multi["delta_activity_weighted"],
        "affinity_gate_threshold_kcal_per_mol": multi["affinity_gate_threshold_kcal_per_mol"],
        "affinity_gate": multi["affinity_gate"],
        "delta_activity_gated": multi["delta_activity_gated"],
        "aspr_score": float(aspr),
        "aspr_reactive": mol_id in LIGAND_REACTIVITY,
        "poses": multi["poses"],
    }
    (out / f"{pair_tag}_report.json").write_text(json.dumps(report, indent=2))

    print(f"  apo activity              {report['activity_apo']:+.3f}")
    print(f"  complex activity (top)    {report['activity_complex_top']:+.3f}")
    print(f"  complex activity (w'tied) {report['activity_complex_weighted']:+.3f}")
    direction_top = "protective" if report["delta_activity_top"] < 0 else "harm-leaning"
    direction_w = "protective" if report["delta_activity_weighted"] < 0 else "harm-leaning"
    print(f"  delta (top)               {report['delta_activity_top']:+.3f}  ({direction_top})")
    print(f"  delta (Boltzmann)         {report['delta_activity_weighted']:+.3f}  ({direction_w})")
    print(
        f"  delta (gated, thr={report['affinity_gate_threshold_kcal_per_mol']:.1f}) "
        f"{report['delta_activity_gated']:+.3f}  (gate={report['affinity_gate']:.3f})"
    )
    if report["aspr_reactive"]:
        print(f"  aspr_score (covalent)     {report['aspr_score']:+.3f}  (harm-leaning if > 0)")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("molecule", help="vicinity-molecule id, e.g. 'curcumin'")
    parser.add_argument(
        "anchor",
        help="anchor PDB id (e.g. '6PEO') OR path to a generated oligomer PDB file",
    )
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num_modes", type=int, default=5)
    parser.add_argument("--padding", type=float, default=6.0)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    receptor = Path(args.anchor)
    if receptor.exists() and receptor.suffix.lower() == ".pdb":
        perturb_oligomer(
            args.molecule, receptor,
            exhaustiveness=args.exhaustiveness, n_poses=args.num_modes,
            padding=args.padding, seed=args.seed,
        )
    else:
        perturb(
            args.molecule, args.anchor,
            exhaustiveness=args.exhaustiveness, n_poses=args.num_modes,
            padding=args.padding, seed=args.seed,
        )


if __name__ == "__main__":
    main()


# -----------------------------------------------------------------------------
# What this prototype detects, and what it doesn't:
#
#   Detects:
#     - SASA-occlusion from a docked ligand: hydrophobic β-residues hidden
#       under the pose lose surface area and the activity score drops.
#     - Ligand-aware contact_density: Cα–LIG heavy-atom contacts inside the
#       8 Å cutoff inflate the stability proxy in proportion to binding
#       footprint.
#     - Multi-pose aggregation: each of Vina's --num_modes poses is
#       individually scored on Stage 2 features and the per-pose
#       Δactivity values are combined with Boltzmann weights
#       exp(-aff_i / RT) / Z at T=300 K. A tightly-binding pose with a
#       single dominant affinity reproduces the top-pose Δactivity; a
#       weak-affinity ligand whose poses span a small affinity range
#       sees its score smeared across the rank-2..5 contributions.
#       Top-pose Δactivity is still reported for comparison.
#
#   Misses (deliberately, for now):
#     - Conformational rearrangement on binding. Without MD, the receptor
#       cannot move and the per-feature deltas are static-occlusion +
#       ligand-Cα-contact only.
#
#   Wired in but reported as an orthogonal column (not folded into the
#   primary delta_activity_weighted):
#     - Absolute-affinity gate. delta_activity_gated = gate * weighted
#       where gate = min(1, exp((thr - aff_top) / RT)) at thr=-6 kcal/mol.
#       Penalises ligands whose top-pose affinity is itself weak, on top
#       of the intra-ligand Boltzmann pose weighting.
# -----------------------------------------------------------------------------
