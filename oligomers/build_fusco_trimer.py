"""Build a Fusco 2017 Type B* α-syn oligomer starting structure.

Topology prior (from Fusco et al., Science 2017, "Type B*" toxic oligomer):

  - 2-4 α-syn monomers (full length, residues 1-140)
  - rigid β-core in the C-terminal NAC region, one extended β-strand
    per monomer (default residues 70-88)
  - strands form a small β-sheet (parallel by default; antiparallel
    available)
  - N-terminal (1-69) and C-terminal (89-140) tails as extended/PPII
    coil — will collapse to disordered conformations during MD
    refinement; here we set the initial coords far enough apart that
    they don't clash

This script writes ONLY the starting structure. The intended next
step is `md_relax.py --apo-pdb <out> --restrain-residues 70-88` to
relax side chains and tails while preserving the β-core topology.
That relaxed ensemble is the candidate active-side anchor the
deposited PDB does not provide; the original Stage 2 framework
calibrates with mutant fibrils only because no toxic-oligomer
structure exists.

Many topologies can be tried by varying CLI args:
  --n-mers {2,3,4}              dimer / trimer / tetramer
  --core-start, --core-end      β-core residue range (default 70-88)
  --arrangement {parallel, antiparallel}
  --spacing 4.7                 inter-strand Cα-Cα distance (Å)
  --tag NAME                    output filename stem

Output: results/oligomers/<tag>.pdb (one PDB per topology).
"""
from __future__ import annotations

import argparse
import copy
import random
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
from PeptideBuilder import Geometry, PeptideBuilder

# Docstrings and help strings contain non-ASCII (α, β, …) which would
# otherwise crash --help on Windows consoles defaulting to cp1252.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

ALPHA_SYN_SEQ = (
    "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYV"   # 1-40
    "GSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQK"   # 41-80
    "TVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP"   # 81-120
    "DNEAYEMPSEEGYQDYEPEA"                       # 121-140
)
assert len(ALPHA_SYN_SEQ) == 140, f"α-syn must be 140 aa, got {len(ALPHA_SYN_SEQ)}"

# β-strand (rigid core): φ=-130°, ψ=+130°.
BETA_PHI = -130.0
BETA_PSI = +130.0

# Disordered-coil φ/ψ mixture for the unrestrained tails.
#
# Per-residue independent draws from two extended basins (β and PPII).
# A 50/50 mix gives a moderately random walk in 3D — kinks where the
# residue type changes — but both basins are EXTENDED, so the chain
# does not loop back through itself.
#
# We deliberately exclude the α-helix basin: α produces tight 3.6-residue
# turns and a 70-residue tail sampling α with reasonable weight has high
# probability of self-intersecting. Empirically (built and inspected),
# including α at 30% weight produced >30 atom pairs within 1.0 Å in the
# built oligomer (chain folding back on its own backbone), which OpenMM
# minimization could not resolve and crashed at the first integration
# step with NaN coordinates.
#
# Each entry: (φ_mean°, ψ_mean°, sigma°, weight).
COIL_BASINS: list[tuple[float, float, float, float]] = [
    (-130.0, +130.0, 20.0, 0.40),   # β-region
    (-75.0,  +145.0, 20.0, 0.60),   # PPII
]


def sample_coil(rng: random.Random) -> tuple[float, float]:
    r = rng.random()
    cum = 0.0
    for phi, psi, sigma, w in COIL_BASINS:
        cum += w
        if r <= cum:
            return rng.gauss(phi, sigma), rng.gauss(psi, sigma)
    phi, psi, _, _ = COIL_BASINS[-1]
    return phi, psi


def build_monomer(core_start: int, core_end: int, seed: int = 42) -> "Structure":
    """Build one α-syn chain. Residues [core_start, core_end] take a uniform
    β-strand backbone (preserved by Cα restraints during MD); all other
    residues draw φ/ψ independently from a coil mixture so the chain
    random-walks into a compact disordered shape."""
    rng = random.Random(seed)

    # Pre-sample every residue's (φ, ψ).
    phis: list[float] = []
    psis: list[float] = []
    for r in range(1, len(ALPHA_SYN_SEQ) + 1):
        if core_start <= r <= core_end:
            phis.append(BETA_PHI)
            psis.append(BETA_PSI)
        else:
            phi, psi = sample_coil(rng)
            phis.append(phi)
            psis.append(psi)

    aa0 = ALPHA_SYN_SEQ[0]
    geo = Geometry.geometry(aa0)
    geo.phi = phis[0]
    geo.psi_im1 = psis[0]
    structure = PeptideBuilder.initialize_res(geo)

    for idx, aa in enumerate(ALPHA_SYN_SEQ[1:], start=2):
        g = Geometry.geometry(aa)
        g.phi = phis[idx - 1]
        g.psi_im1 = psis[idx - 2]
        structure = PeptideBuilder.add_residue(structure, g)
    return structure


def ca_coords(structure, resnum_range: tuple[int, int] | None = None) -> np.ndarray:
    chain = next(iter(next(iter(structure))))
    out: list[np.ndarray] = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        if resnum_range is not None and not (resnum_range[0] <= residue.id[1] <= resnum_range[1]):
            continue
        if "CA" in {a.get_name() for a in residue}:
            out.append(residue["CA"].coord)
    return np.asarray(out)


def orient_to_strand_frame(structure, core_start: int, core_end: int) -> None:
    """Rigid-body transform so the β-core's long axis aligns with +x,
    sheet normal aligns with +z, and core centroid sits at the origin.

    When core_start > core_end (impossible range — used for all-coil controls),
    the core Cα set is empty; fall back to using all-Cα for orientation so
    the strands still point roughly along +x before assembly translation."""
    core = ca_coords(structure, (core_start, core_end))
    if len(core) == 0:
        # all-coil control: orient using the full chain
        core = ca_coords(structure)
    center = core.mean(axis=0)
    centered = core - center
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    R = vh.copy()  # rows: axes in decreasing variance
    if np.linalg.det(R) < 0:
        R[2] = -R[2]  # enforce right-handed coordinate system
    for atom in structure.get_atoms():
        atom.coord = (R @ (atom.coord - center)).astype(np.float32)


def clone(structure, new_chain_id: str):
    s2 = copy.deepcopy(structure)
    chain = next(iter(next(iter(s2))))
    chain.id = new_chain_id
    return s2


def translate(structure, dx: float, dy: float, dz: float) -> None:
    delta = np.array([dx, dy, dz], dtype=np.float32)
    for atom in structure.get_atoms():
        atom.coord = atom.coord + delta


def rotate_z180(structure) -> None:
    """180° around z-axis: flips a strand for antiparallel pairing while
    keeping the sheet normal pointing +z."""
    for atom in structure.get_atoms():
        c = atom.coord
        atom.coord = np.array([-c[0], -c[1], c[2]], dtype=np.float32)


def rotate_x_deg(structure, angle_deg: float) -> None:
    """Rotate around the x-axis (β-strand axis). Keeps the β-core Cα
    atoms approximately on the x-axis but spins the rest of the chain
    (tails) into a different azimuthal direction. Used to separate the
    monomer tails in an oligomer so they don't all emerge into the same
    volume of space."""
    angle = np.radians(angle_deg)
    c, s = np.cos(angle), np.sin(angle)
    R = np.array([[1.0, 0.0, 0.0],
                  [0.0,   c,   -s],
                  [0.0,   s,    c]], dtype=np.float32)
    for atom in structure.get_atoms():
        atom.coord = (R @ atom.coord).astype(np.float32)


def assemble(monomers: list, chain_ids: list[str]):
    sb = StructureBuilder()
    sb.init_structure("oligomer")
    sb.init_model(0)
    out = sb.get_structure()
    model = next(iter(out))
    for mono, cid in zip(monomers, chain_ids):
        ch = next(iter(next(iter(mono))))
        ch.detach_parent()
        ch.id = cid
        model.add(ch)
    # Refresh full_ids so atom keys reflect the new hierarchy. SASA / DSSP
    # in features.py cache by full_id and break otherwise.
    for atom in out.get_atoms():
        atom.full_id = atom.get_full_id()
    return out


def build(
    n_mers: int,
    core_start: int,
    core_end: int,
    arrangement: str,
    spacing: float,
    seed: int = 42,
):
    chain_ids = list("ABCDEFGH"[:n_mers])
    monomers = []
    for i, cid in enumerate(chain_ids):
        # Independent random walk per chain so the three tails don't
        # superimpose on top of each other in identical conformations.
        # NOTE: this still gives all tails roughly the same expected
        # direction; inter-chain tail clashes are expected and are
        # resolved by vacuum minimization before solvation (see
        # md_relax.py --vacuum-min-iter).
        m = build_monomer(core_start, core_end, seed=seed + i)
        orient_to_strand_frame(m, core_start, core_end)
        chain = next(iter(next(iter(m))))
        chain.id = cid
        if arrangement == "antiparallel" and i % 2 == 1:
            rotate_z180(m)
        translate(m, 0.0, i * spacing, 0.0)
        monomers.append(m)
    return assemble(monomers, chain_ids)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-mers", type=int, default=3)
    p.add_argument("--core-start", type=int, default=70)
    p.add_argument("--core-end", type=int, default=88)
    p.add_argument("--arrangement", choices=["parallel", "antiparallel"], default="parallel")
    p.add_argument("--spacing", type=float, default=4.7,
                   help="inter-strand Cα-Cα distance in Å (β-sheet H-bond spacing)")
    p.add_argument("--seed", type=int, default=42,
                   help="random seed for the coil φ/ψ draws (vary across "
                        "ensemble members to get different tail conformations)")
    p.add_argument("--tag", type=str, default=None)
    p.add_argument("--out-dir", type=Path,
                   default=Path(__file__).resolve().parents[1] / "results" / "oligomers")
    args = p.parse_args()

    tag = args.tag or (
        f"fusco_{args.arrangement}_{args.n_mers}mer_"
        f"core{args.core_start}-{args.core_end}"
    )
    out_pdb = args.out_dir / f"{tag}.pdb"
    out_pdb.parent.mkdir(parents=True, exist_ok=True)

    print(f"Building {args.n_mers}-mer, beta-core {args.core_start}-{args.core_end}, "
          f"{args.arrangement}, spacing {args.spacing} A")
    s = build(
        n_mers=args.n_mers,
        core_start=args.core_start,
        core_end=args.core_end,
        arrangement=args.arrangement,
        spacing=args.spacing,
        seed=args.seed,
    )
    io = PDBIO()
    io.set_structure(s)
    io.save(str(out_pdb))
    n_atoms = sum(1 for _ in s.get_atoms())
    print(f"Wrote {out_pdb}  ({n_atoms} atoms, chains {[ch.id for ch in next(iter(s))]})")


if __name__ == "__main__":
    main()
