"""Load anchor structures from the PDB. See ANCHORS.md for the curated list
and confidence levels. Active anchors are deferred until verified — initial
Stage 2 validation runs on inert anchors only, with mutant fibrils added as
graded-active controls once their PDB IDs are confirmed."""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import PDBList, PDBParser

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "anchors"


@dataclass(frozen=True)
class Anchor:
    pdb_id: str
    label: str
    description: str
    # Literature-curated biological protofilament count. 0 = not applicable
    # (monomer 1XQ8). Note: this may differ from the count of chains
    # actually present in the deposited biological assembly — several
    # cryo-EM α-syn fibrils deposit only one protofilament and leave the
    # second as a symmetry mate (see protofilaments.count_protofilaments
    # for the deposited count).
    n_protofilaments: int


ANCHORS: list[Anchor] = [
    # inert: physiological state + disease endpoints
    Anchor("1XQ8", "inert", "micelle-bound helical monomer (Ulmer & Bax 2005)", 0),
    Anchor("2N0A", "inert", "Greek-key ssNMR fibril (Tuttle 2016)", 1),
    Anchor("6CU7", "inert", "rod cryo-EM polymorph (Li 2018)", 2),
    Anchor("6CU8", "inert", "twister cryo-EM polymorph (Li 2018)", 2),
    Anchor("6H6B", "inert", "recombinant cytotoxic fibril (Guerrero-Ferreira 2018)", 2),
    Anchor("6XYO", "inert", "MSA Type I brain-derived fibril (Schweighauser 2020)", 2),
    Anchor("6XYP", "inert", "MSA Type II-1 brain-derived fibril (Schweighauser 2020)", 2),
    Anchor("8A9L", "inert", "Lewy fold from PD/PDD/DLB brain (Yang 2022)", 1),
    Anchor("8A4L", "inert", "lipidic fibril polymorph L2A (Antonschmidt 2022)", 2),
    # graded-active: familial mutant fibrils; should rank ≥ WT fibrils ordinally
    Anchor("6LRQ", "graded-active", "A53T fibril (Sun 2020)", 2),
    Anchor("7WO0", "graded-active", "A53T fibril induced by Ca2+ (vicinity-molecule context)", 2),
    Anchor("6UFR", "graded-active", "E46K pathogenic fibril (Boyer 2020)", 2),
    Anchor("6PEO", "graded-active", "H50Q narrow fibril (Boyer 2019)", 1),
    Anchor("6PES", "graded-active", "H50Q wide fibril (Boyer 2019)", 2),
    # direct active: none yet — see ANCHORS.md for the gap
]


def fetch(pdb_id: str) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    target = DATA_DIR / f"{pdb_id.lower()}.pdb"
    if target.exists():
        return target
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        PDBList(verbose=False).retrieve_pdb_file(
            pdb_id, pdir=str(DATA_DIR), file_format="pdb"
        )
    ent = DATA_DIR / f"pdb{pdb_id.lower()}.ent"
    if ent.exists():
        ent.rename(target)
    if not target.exists():
        raise FileNotFoundError(f"PDB fetch failed for {pdb_id}")
    return target


def load(pdb_id: str):
    path = fetch(pdb_id)
    return PDBParser(QUIET=True).get_structure(pdb_id, path)
