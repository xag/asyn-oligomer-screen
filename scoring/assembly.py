"""Build the biological assembly from REMARK 350 BIOMT operations.

The asymmetric unit deposited in a PDB entry is often a fraction of the
biologically relevant multimer — three chains of a nine-chain fibril
block (8A9L), for example. Stage 2 features must compare equivalent
assemblies across anchors, otherwise rankings reflect depositor
formatting rather than biology.

This module reconstructs the BIOMOLECULE 1 assembly from BIOMT records.
It also exposes `inner_chain_ids`, which picks the most buried chain(s)
of an assembly so that features can be averaged over chains whose
environment is dominated by inter-chain contacts rather than free
surface.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB import SASA
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure


@dataclass(frozen=True)
class Biomol:
    chains: tuple[str, ...]
    operators: tuple[tuple[np.ndarray, np.ndarray], ...]  # (R, t) pairs


_APPLY_RE = re.compile(r"APPLY THE FOLLOWING TO CHAINS:\s*(.*)", re.IGNORECASE)
_AND_RE = re.compile(r"AND CHAINS:\s*(.*)", re.IGNORECASE)
_BIOMT_RE = re.compile(
    r"^REMARK 350\s+BIOMT([123])\s+(\d+)\s+"
    r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
)


def parse_biomt(pdb_path: Path) -> Biomol | None:
    """Read REMARK 350 BIOMOLECULE 1 and return the chains and operators
    that build the author-determined biological unit. Returns None when
    no BIOMT records are present (asymmetric unit IS the assembly)."""
    chains: list[str] = []
    rows: dict[int, dict[int, list[float]]] = {}
    in_biomol_1 = False
    saw_biomol = False
    reading_chains = False

    for raw in pdb_path.read_text(errors="ignore").splitlines():
        if not raw.startswith("REMARK 350"):
            continue
        line = raw[10:].strip()

        if line.startswith("BIOMOLECULE:"):
            num = int(line.split(":")[1].strip())
            if saw_biomol and not in_biomol_1:
                break
            in_biomol_1 = num == 1
            saw_biomol = True
            reading_chains = False
            continue
        if saw_biomol and not in_biomol_1:
            continue

        m = _APPLY_RE.search(line)
        if m:
            chains = [s.strip() for s in m.group(1).split(",") if s.strip()]
            reading_chains = True
            continue
        m = _AND_RE.search(line)
        if m and reading_chains:
            chains.extend(s.strip() for s in m.group(1).split(",") if s.strip())
            continue

        m = _BIOMT_RE.match(raw)
        if m:
            reading_chains = False
            row_n = int(m.group(1))
            op_n = int(m.group(2))
            rows.setdefault(op_n, {})[row_n] = [float(x) for x in m.groups()[2:6]]

    if not rows or not chains:
        return None

    ops: list[tuple[np.ndarray, np.ndarray]] = []
    for op_n in sorted(rows):
        rec = rows[op_n]
        if not all(k in rec for k in (1, 2, 3)):
            continue
        mat = np.array([rec[1], rec[2], rec[3]])
        ops.append((mat[:, :3], mat[:, 3]))
    return Biomol(chains=tuple(chains), operators=tuple(ops))


def _copy_chain(chain: Chain, new_id: str, R: np.ndarray, t: np.ndarray) -> Chain:
    """Manual clone of a chain with coordinates transformed by R, t.
    Avoids deepcopy's parent-pointer traversal which is fragile in
    Biopython."""
    new_chain = Chain(new_id)
    for residue in chain:
        new_res = Residue(residue.id, residue.resname, residue.segid)
        for atom in residue:
            coord = R @ atom.coord + t
            new_atom = Atom(
                atom.get_name(),
                coord,
                atom.bfactor,
                atom.occupancy,
                atom.altloc,
                atom.fullname,
                atom.serial_number,
                element=atom.element,
            )
            new_res.add(new_atom)
        new_chain.add(new_res)
    return new_chain


def build_assembly(structure: Structure, biomol: Biomol | None) -> Structure:
    """Apply BIOMT operators to model 0 and return a new Structure whose
    single model holds one chain per (original-chain, operator) pair.

    Identity operator (#1) preserves original chain IDs. Subsequent
    operators get suffixed IDs (A → A_2, A_3, ...). Chains in model 0
    that aren't listed under APPLY THE FOLLOWING TO CHAINS are dropped.
    """
    src_model = next(iter(structure))
    new_struct = Structure(structure.id + "_asm")
    new_model = Model(0)
    new_struct.add(new_model)

    if biomol is None:
        for chain in src_model:
            R = np.eye(3)
            t = np.zeros(3)
            new_model.add(_copy_chain(chain, chain.id, R, t))
        for atom in new_struct.get_atoms():
            atom.full_id = atom.get_full_id()
        return new_struct

    src_by_id = {c.id: c for c in src_model}
    for op_idx, (R, t) in enumerate(biomol.operators, start=1):
        for chain_id in biomol.chains:
            if chain_id not in src_by_id:
                continue
            new_id = chain_id if op_idx == 1 else f"{chain_id}_{op_idx}"
            new_model.add(_copy_chain(src_by_id[chain_id], new_id, R, t))

    # Biopython caches atom.full_id when set_parent() runs (i.e. when the
    # atom is first attached to its residue). At that moment the residue
    # has no chain parent yet, so every atom's cached full_id collapses
    # to (residue_id, (name, altloc)) — identical across chains. This
    # corrupts SASA aggregation (atomdict in Bio.PDB.SASA keys on
    # full_id and silently overwrites duplicates). Refresh full_ids
    # now that the full hierarchy exists.
    for atom in new_struct.get_atoms():
        atom.full_id = atom.get_full_id()
    return new_struct


def load_assembly(pdb_id: str) -> Structure:
    """Fetch (if needed) and return the biological assembly Structure."""
    from anchors import fetch

    path = fetch(pdb_id)
    structure = PDBParser(QUIET=True).get_structure(pdb_id, path)
    biomol = parse_biomt(path)
    return build_assembly(structure, biomol)


_STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def inner_chain_ids(structure: Structure, top_k: int = 1, min_residues: int = 20) -> list[str]:
    """Return the chain IDs with the smallest per-residue SASA in the
    assembly — i.e. the most buried chains. These chains see a
    chain-stack environment dominated by inter-chain contacts and best
    approximate the bulk fibril/oligomer interior.

    Chains with fewer than `min_residues` recognized amino acids are
    excluded; many fibril depositions include short UNK fragments
    (modelled ligand-like density) which would otherwise game the
    ranking. Falls back to the longest available chain if no chain
    meets the threshold."""
    model = next(iter(structure))
    if not getattr(structure, "_sasa_done", False):
        SASA.ShrakeRupley().compute(model, level="R")
        structure._sasa_done = True

    per_chain: list[tuple[str, int, float]] = []
    for chain in model:
        residues = [r for r in chain if r.id[0] == " " and r.get_resname() in _STANDARD_AA]
        if not residues:
            continue
        total = sum(getattr(r, "sasa", 0.0) for r in residues)
        per_chain.append((chain.id, len(residues), total / len(residues)))

    eligible = [c for c in per_chain if c[1] >= min_residues]
    if not eligible:
        eligible = sorted(per_chain, key=lambda x: x[1], reverse=True)[: max(top_k, 1)]
    eligible.sort(key=lambda x: x[2])
    return [cid for cid, _, _ in eligible[: max(top_k, 1)]]
