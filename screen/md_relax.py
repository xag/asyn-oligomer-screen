"""MD relaxation of a Stage 3 docked complex (OpenMM + OpenFF).

The `--prepare-only` build path runs in the conda MD env (`environment-md.yml`,
located by `md_env.py`); the OpenFF packages it imports aren't in the pip venv.
The `--system-xml` dynamics path runs in the pip venv. `stage3.py` and
`md_stage3.py` invoke this script via subprocess.

Inputs:
  --complex-pdb   apo + docked LIG on chain Z (what stage3 emits as
                  `<pair>_complex.pdb`)
  --ligand-smiles  SMILES used to dock the ligand; needed to assign
                   bond orders to the docked PDB block
  --out-pdb       where to write the relaxed protein + ligand PDB

Defaults to 100 ps total equilibration (NVT + NPT) + 100 ps production
NPT at 300 K, 1 atm, TIP3P, 0.15 M NaCl, 1 nm padding, 2 fs steps with
HBonds constraints. Bumpable via flags.

The point of this step: MD relaxation around the docked pose so the receptor side-chains
and backbone can rearrange, which is the only mechanism through which
Δactivity can flip positive under the Stage 2 feature weights
(static-pose features only ever subtract from activity).

Distributed / crowdsourced split (issue #34). The OpenFF parametrisation the
docked-complex path needs is GPU-free and one-time, so it is separable from
the per-replica GPU dynamics. Two extra modes implement that split so a
docked-complex dwell chunk can run on a volunteer's basic GPU with no conda:

  --prepare-only PREFIX   build + solvate + parametrise, serialise the OpenMM
      System to PREFIX_system.xml and the solvated topology/positions to
      PREFIX_solvated.pdb, then exit without running dynamics. This is the
      *only* step that uses OpenFF (runs in the conda MD env).
  --system-xml/--solvated-pdb  run the dynamics from a serialised System
      instead of building one. Needs *only* pip-installed openmm — no OpenFF,
      no conda — because every force-field term (including the ligand's
      SMIRNOFF parameters) is already baked into system.xml. This is the
      chunk a contributor's GPU runs; identical to the apo chunk's runtime.
"""
from __future__ import annotations

import argparse
import io
import sys
import time
from pathlib import Path

import numpy as np

import openmm
import openmm.app as app
import openmm.unit as u
from pdbfixer import PDBFixer

# OpenFF / openmmforcefields are imported lazily inside the ligand-parameterising
# functions (offmol_from_rdkit, build_system). They are only needed for the
# *complex* path; the *apo* baseline (off_ligand=None) uses amber14 + TIP3P
# only, so it runs in a pip-only env (openmm + pdbfixer) without the OpenFF
# stack — which, unlike OpenMM, is not currently installable from PyPI.

from rdkit import Chem
from rdkit.Chem import AllChem


NAGL_MODEL = "openff-gnn-am1bcc-0.1.0-rc.3.pt"
OFF_FF = "openff-2.2.0.offxml"

LIGAND_CHAIN = "Z"
LIGAND_RESNAME = "LIG"


# -----------------------------------------------------------------------------
# Split a complex PDB into protein-only and ligand-only PDB blocks.
# -----------------------------------------------------------------------------

def split_complex_pdb(complex_pdb: Path) -> tuple[str, str]:
    """Return (protein_pdb_text, ligand_pdb_text). Splits on chain Z LIG."""
    prot_lines: list[str] = []
    lig_lines: list[str] = []
    for raw in complex_pdb.read_text(encoding="utf-8").splitlines():
        if raw.startswith(("ATOM", "HETATM", "TER")):
            chain = raw[21:22]
            if chain == LIGAND_CHAIN:
                if raw.startswith("HETATM") or raw.startswith("ATOM"):
                    lig_lines.append(raw)
                # drop TER on Z (we wrap it ourselves)
            else:
                prot_lines.append(raw)
        elif raw.startswith(("HEADER", "CRYST1", "REMARK")):
            prot_lines.append(raw)
        # END / others: drop, we wrap our own
    prot_text = "\n".join(prot_lines) + "\nEND\n"
    lig_text = "\n".join(lig_lines) + "\nEND\n"
    return prot_text, lig_text


def load_apo_pdb(apo_pdb: Path) -> str:
    """Return protein-only PDB text from an apo PDB file (no LIG chain).
    Used for the apo-MD baseline."""
    lines: list[str] = []
    for raw in apo_pdb.read_text(encoding="utf-8").splitlines():
        if raw.startswith(("ATOM", "HETATM", "TER", "HEADER", "CRYST1", "REMARK")):
            lines.append(raw)
    return "\n".join(lines) + "\nEND\n"


# -----------------------------------------------------------------------------
# Build a properly bonded ligand RDKit Mol from (PDB block, reference SMILES).
# The PDB block has only heavy atoms + a couple of polar Hs (Vina output) and
# no CONECT records. We rely on RDKit's distance-based bond perception, then
# transplant bond orders + formal charges from the SMILES template, then add
# the missing hydrogens back with sensible 3D positions.
# -----------------------------------------------------------------------------

def ligand_from_pdb_and_smiles(ligand_pdb_text: str, smiles: str) -> Chem.Mol:
    # Strip the PDB-block Hs — AssignBondOrdersFromTemplate matches heavy atoms.
    pdb_mol_with_h = Chem.MolFromPDBBlock(ligand_pdb_text, removeHs=False, sanitize=False)
    if pdb_mol_with_h is None:
        raise ValueError("RDKit could not parse the ligand PDB block")
    pdb_mol = Chem.RemoveHs(pdb_mol_with_h, sanitize=False)

    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError(f"RDKit could not parse SMILES {smiles!r}")

    fixed = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
    fixed = Chem.AddHs(fixed, addCoords=True)
    # Light MMFF cleanup of H positions only — heavy atoms are frozen at the
    # docked coords so the pose stays put while Hs settle into a sane geometry.
    AllChem.MMFFOptimizeMolecule(fixed, maxIters=200)
    # Drop per-atom PDB monomer info. AddHs leaves the new H atoms with no
    # MonomerInfo, which makes openff.toolkit.Molecule.from_rdkit split the
    # molecule into multiple residues (one per metadata cluster) on the way
    # to OpenMM Topology — and SMIRNOFFTemplateGenerator then can't match
    # the partial residue. Clearing the info forces a single-residue topology.
    for atom in fixed.GetAtoms():
        atom.SetMonomerInfo(None)
    return fixed


# -----------------------------------------------------------------------------
# Build an OpenFF Molecule with NAGL charges from an RDKit Mol.
# -----------------------------------------------------------------------------

def offmol_from_rdkit(rdkit_mol: Chem.Mol) -> "Molecule":
    from openff.toolkit import Molecule
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

    off = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
    off.name = "LIG"
    off.assign_partial_charges(NAGL_MODEL, toolkit_registry=NAGLToolkitWrapper())
    return off


# -----------------------------------------------------------------------------
# Build the system: amber14 protein + TIP3P water + NaCl + SMIRNOFF ligand.
# -----------------------------------------------------------------------------

def build_system(
    protein_pdb_text: str,
    off_ligand: Molecule | None,
    padding_nm: float,
    salt_mol: float,
    rectangular_box: bool = False,
):
    # Use PDBFixer to add missing Hs to the protein. PDB chains here are
    # rebuilt by stage3 and already have all heavy atoms.
    fixer = PDBFixer(pdbfile=io.StringIO(protein_pdb_text))
    fixer.findMissingResidues()
    fixer.missingResidues.clear()  # don't model in loops we don't have
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    protein_topology = fixer.topology
    protein_positions = fixer.positions

    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    modeller = app.Modeller(protein_topology, protein_positions)

    if off_ligand is not None:
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
        smirnoff_gen = SMIRNOFFTemplateGenerator(molecules=[off_ligand], forcefield=OFF_FF.replace(".offxml", ""))
        forcefield.registerTemplateGenerator(smirnoff_gen.generator)
        lig_top = off_ligand.to_topology().to_openmm()
        lig_pos = off_ligand.conformers[0].to_openmm()
        modeller.add(lig_top, lig_pos)

    # Box + solvent + ions.
    #
    # padding= makes addSolvent build a *cubic* box sized to the largest
    # solute dimension. For an elongated construct (an extended β-strand is
    # ~14 nm long but only ~4 nm wide) that cube is mostly water: the
    # core58-102 chunk solvates to ~390k atoms cubic vs ~55k rectangular.
    # rectangular_box=True instead fits a tight box to the solute's actual
    # bounding box + padding on each axis — same 1 nm minimum-image clearance,
    # ~7× fewer atoms, which is what keeps one dwell chunk on a basic GPU.
    if rectangular_box:
        pos_nm = np.array(modeller.positions.value_in_unit(u.nanometer))
        extent = pos_nm.max(axis=0) - pos_nm.min(axis=0)
        box = extent + 2.0 * padding_nm
        modeller.addSolvent(
            forcefield,
            model="tip3p",
            boxSize=openmm.Vec3(*box) * u.nanometers,
            ionicStrength=salt_mol * u.molar,
            positiveIon="Na+",
            negativeIon="Cl-",
        )
    else:
        modeller.addSolvent(
            forcefield,
            model="tip3p",
            padding=padding_nm * u.nanometers,
            ionicStrength=salt_mol * u.molar,
            positiveIon="Na+",
            negativeIon="Cl-",
        )

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * u.nanometers,
        constraints=app.HBonds,
        rigidWater=True,
    )

    n_atoms = modeller.topology.getNumAtoms()
    n_prot = protein_topology.getNumAtoms()
    n_lig = off_ligand.to_topology().to_openmm().getNumAtoms() if off_ligand is not None else 0
    print(f"  system: {n_atoms} atoms (protein {n_prot}, ligand {n_lig}, solvent+ions {n_atoms - n_prot - n_lig})", flush=True)
    return modeller, system


# -----------------------------------------------------------------------------
# Position restraints + implicit-solvent collapse (for oligomer hand-builds).
# Extended-coil tails on a from-sequence-built oligomer would otherwise force
# a huge explicit-solvent box; OBC2 GBSA collapses tails in a few ns at
# minimal cost. Restraints on the β-core Cα atoms preserve the prescribed
# topology while everything else relaxes.
# -----------------------------------------------------------------------------


def parse_residue_range(spec: str | None) -> tuple[int, int] | None:
    if not spec:
        return None
    if "-" not in spec:
        n = int(spec)
        return (n, n)
    lo, hi = spec.split("-", 1)
    return (int(lo), int(hi))


def parse_chain_set(spec: str | None) -> set[str] | None:
    if not spec:
        return None
    return {c.strip() for c in spec.split(",") if c.strip()}


def add_position_restraints(
    system,
    topology,
    positions,
    restrain_range: tuple[int, int] | None,
    restrain_chains: set[str] | None,
    k_kj_per_mol_nm2: float,
) -> int:
    """Pull Cα atoms in `restrain_range` (and `restrain_chains` if given)
    back to their current positions with a harmonic spring of strength
    `k_kj_per_mol_nm2` kJ/mol/nm². Returns the number of restrained
    particles. No-op when `restrain_range is None`."""
    if restrain_range is None:
        return 0
    lo, hi = restrain_range
    force = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    k = k_kj_per_mol_nm2 * u.kilojoule_per_mole / u.nanometer**2
    k_val = k.value_in_unit_system(u.md_unit_system)
    n = 0
    for atom in topology.atoms():
        if atom.name != "CA":
            continue
        res = atom.residue
        if restrain_chains is not None and res.chain.id not in restrain_chains:
            continue
        try:
            resid = int(res.id)
        except (ValueError, TypeError):
            continue
        if not (lo <= resid <= hi):
            continue
        pos = positions[atom.index].value_in_unit(u.nanometer)
        force.addParticle(atom.index, [k_val, pos[0], pos[1], pos[2]])
        n += 1
    if n > 0:
        system.addForce(force)
    return n


def write_heavy_only_pdb(pdb_text: str, out_pdb: Path) -> None:
    """Strip H atoms from a PDB text block and write to disk. Stage 2
    features were calibrated on RCSB heavy-atom-only structures, so
    relaxed-PDB outputs from md_relax need to match that convention
    for downstream scoring."""
    kept: list[str] = []
    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            kept.append(line)
            continue
        elem = line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[0]
        if elem == "H":
            continue
        kept.append(line)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(kept) + "\n", encoding="utf-8")


def vacuum_minimize_apo(
    prot_pdb_text: str,
    max_iterations: int,
    restrain_range: tuple[int, int] | None,
    restrain_chains: set[str] | None,
    restrain_k: float,
) -> str:
    """Energy-minimize the protein in vacuum (no waters) before solvation.
    A hand-built oligomer with random-coil tails typically has both
    intra-chain and inter-chain clashes that minimization-after-solvation
    cannot resolve (the surrounding waters block atom rearrangement).
    Vacuum minimization gives the tails free space to slide out of each
    other's way. The β-core is held by Cα restraints so the prescribed
    topology survives. Returns a PDB text block in the relaxed frame."""
    fixer = PDBFixer(pdbfile=io.StringIO(prot_pdb_text))
    fixer.findMissingResidues()
    fixer.missingResidues.clear()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    forcefield = app.ForceField("amber14-all.xml")
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )
    n_restr = add_position_restraints(
        system, fixer.topology, fixer.positions,
        restrain_range, restrain_chains, restrain_k,
    )
    print(f"  vacuum-minimize: {fixer.topology.getNumAtoms()} atoms, "
          f"{n_restr} Cα restraints, up to {max_iterations} iter", flush=True)

    integrator = openmm.VerletIntegrator(1.0 * u.femtosecond)  # unused
    platform, props = pick_platform()
    sim = app.Simulation(fixer.topology, system, integrator, platform, props)
    sim.context.setPositions(fixer.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    print(f"  vacuum pre-min energy: {e0:.1f} kJ/mol", flush=True)
    sim.minimizeEnergy(maxIterations=max_iterations)
    e1 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    print(f"  vacuum post-min energy: {e1:.1f} kJ/mol (delta {e1-e0:+.1f})", flush=True)

    final = sim.context.getState(getPositions=True)
    buf = io.StringIO()
    app.PDBFile.writeFile(fixer.topology, final.getPositions(), buf)
    return buf.getvalue()


def implicit_collapse_apo(
    prot_pdb_text: str,
    collapse_ps: float,
    temperature_k: float,
    timestep_fs: float,
    restrain_range: tuple[int, int] | None,
    restrain_chains: set[str] | None,
    restrain_k: float,
) -> str:
    """OBC2 implicit-solvent MD on protein only — collapses extended coil
    regions before explicit solvation. Returns a PDB text block in the
    collapsed coordinate frame, ready for the standard explicit-solvent
    `build_system` + production path."""
    fixer = PDBFixer(pdbfile=io.StringIO(prot_pdb_text))
    fixer.findMissingResidues()
    fixer.missingResidues.clear()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=2.0 * u.nanometers,
        constraints=app.HBonds,
    )
    n_restr = add_position_restraints(
        system, fixer.topology, fixer.positions,
        restrain_range, restrain_chains, restrain_k,
    )
    print(f"  implicit-collapse: {fixer.topology.getNumAtoms()} atoms, "
          f"{n_restr} Cα restraints, {collapse_ps:.0f} ps OBC2 MD", flush=True)

    integrator = openmm.LangevinMiddleIntegrator(
        temperature_k * u.kelvin, 1.0 / u.picosecond, timestep_fs * u.femtosecond,
    )
    platform, props = pick_platform()
    sim = app.Simulation(fixer.topology, system, integrator, platform, props)
    sim.context.setPositions(fixer.positions)
    sim.minimizeEnergy(maxIterations=2000)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    print(f"  implicit post-min energy: {e0:.1f} kJ/mol", flush=True)
    sim.context.setVelocitiesToTemperature(temperature_k * u.kelvin)
    steps = int(round(collapse_ps * 1000.0 / timestep_fs))
    sim.reporters.append(app.StateDataReporter(
        sys.stdout, max(1, steps // 10),
        step=True, potentialEnergy=True, temperature=True, speed=True,
    ))
    sim.step(steps)

    final = sim.context.getState(getPositions=True)
    buf = io.StringIO()
    app.PDBFile.writeFile(fixer.topology, final.getPositions(), buf)
    return buf.getvalue()


# -----------------------------------------------------------------------------
# Pick the fastest OpenMM platform actually available.
# -----------------------------------------------------------------------------

def pick_platform() -> tuple[openmm.Platform, dict]:
    priority = ["CUDA", "OpenCL", "CPU", "Reference"]
    available = {openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())}
    for name in priority:
        if name in available:
            plat = openmm.Platform.getPlatformByName(name)
            props = {}
            if name in ("CUDA", "OpenCL"):
                props["Precision"] = "mixed"
            return plat, props
    raise RuntimeError("no OpenMM platforms available")


# -----------------------------------------------------------------------------
# Serialise a prepared (solvated, parametrised) system so the GPU integration
# can run on a machine with only pip-installed OpenMM — no OpenFF, no conda.
# This is the split (issue #34) that lets docked-complex dwell chunks be
# crowdsourced: the one-time, GPU-free SMIRNOFF parametrisation (build_system)
# happens centrally; the per-replica dynamics (run_dynamics) load the
# serialised System and run anywhere.
# -----------------------------------------------------------------------------

def serialize_prepared_system(system, topology, positions, prefix: Path) -> tuple[Path, Path]:
    """Write {prefix}_system.xml (the fully parametrised OpenMM System — ligand
    SMIRNOFF terms baked in) + {prefix}_solvated.pdb (topology + positions,
    including water/ions) so load_prepared_system can rebuild a runnable
    Simulation with pip-only OpenMM. Returns (system_xml, solvated_pdb)."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    sys_xml = prefix.with_name(prefix.name + "_system.xml")
    solv_pdb = prefix.with_name(prefix.name + "_solvated.pdb")
    sys_xml.write_text(openmm.XmlSerializer.serialize(system), encoding="utf-8")
    with solv_pdb.open("w") as f:
        app.PDBFile.writeFile(topology, positions, f, keepIds=True)
    return sys_xml, solv_pdb


def load_prepared_system(system_xml: Path, solvated_pdb: Path):
    """Inverse of serialize_prepared_system: load the serialised System plus the
    solvated topology/positions. Needs only pip-installed openmm (no OpenFF) —
    every force-field term, including the ligand's SMIRNOFF parameters, is
    already inside system_xml. Returns (topology, positions, system)."""
    pdb = app.PDBFile(str(solvated_pdb))
    system = openmm.XmlSerializer.deserialize(Path(system_xml).read_text(encoding="utf-8"))
    if system.getNumParticles() != pdb.topology.getNumAtoms():
        raise ValueError(
            f"system/topology mismatch: {system.getNumParticles()} particles vs "
            f"{pdb.topology.getNumAtoms()} atoms — {Path(system_xml).name} and "
            f"{Path(solvated_pdb).name} are not a matched pair"
        )
    return pdb.topology, pdb.positions, system


# -----------------------------------------------------------------------------
# Dynamics: minimise → NVT warm-up → NPT equilibration → NPT production.
# Shared by the build path (build_system) and the prepared-system path
# (load_prepared_system), so a docked-complex chunk runs identically whether
# OpenFF parametrised it in-process or it arrived as a serialised System on a
# volunteer's GPU.
# -----------------------------------------------------------------------------

def run_dynamics(
    out_pdb: Path,
    topology,
    positions,
    system,
    *,
    equil_ps: float,
    prod_ps: float,
    temperature_k: float,
    pressure_atm: float,
    timestep_fs: float,
    report_interval_ps: float,
    restrain_range: tuple[int, int] | None,
    restrain_chains: set[str] | None,
    restrain_k: float,
    seed: int | None,
    traj_out: Path | None,
    traj_interval_ps: float,
    t0: float | None = None,
) -> None:
    if t0 is None:
        t0 = time.time()

    # Re-apply position restraints on the explicit-solvent system so the
    # β-core topology survives equilibration + production at full T.
    n_restr = add_position_restraints(
        system, topology, positions,
        restrain_range, restrain_chains, restrain_k,
    )
    if n_restr:
        print(f"  explicit-solvent: {n_restr} Cα restraints carried over", flush=True)

    # Build integrator + simulation WITHOUT a barostat first — we minimize and
    # NVT-thermalize in a constant-volume box, then add the barostat for the
    # actual NPT phases. Doing NPT immediately on a freshly-solvated, possibly
    # strained protein (especially a hand-built oligomer with random-coil
    # tails) makes the barostat react to a non-equilibrium pressure spike
    # while the integrator tries to ramp T, and the combination crashes
    # ("Particle coordinate is NaN").
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_k * u.kelvin,
        1.0 / u.picosecond,
        timestep_fs * u.femtosecond,
    )
    # Seed the integrator's random stream so velocity-seeded replicas
    # (dwell_time.py) diverge deterministically: same seed → same
    # trajectory, different seed → an independent thermal realisation.
    if seed is not None:
        integrator.setRandomNumberSeed(int(seed))
    platform, props = pick_platform()
    print(f"  platform: {platform.getName()} (props={props})", flush=True)
    sim = app.Simulation(topology, system, integrator, platform, props)
    sim.context.setPositions(positions)

    print("  minimizing...", flush=True)
    sim.minimizeEnergy(maxIterations=10000)
    state = sim.context.getState(getEnergy=True)
    print(f"  post-min energy: {state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole):.1f} kJ/mol", flush=True)

    # NVT thermalization ramp at small dt: 50 K → 150 K → 300 K, 1 ps each
    # at 0.5 fs. Catches strained bonds before they NaN. Cheap (~2 sec).
    warmup_dt_fs = 0.5
    warmup_steps_per_ps = int(round(1000.0 / warmup_dt_fs))
    integrator.setStepSize(warmup_dt_fs * u.femtosecond)
    for warmup_t in (50.0, 150.0, temperature_k):
        integrator.setTemperature(warmup_t * u.kelvin)
        sim.context.setVelocitiesToTemperature(warmup_t * u.kelvin)
        sim.step(warmup_steps_per_ps)  # 1 ps each
    print(f"  warmup: 3 ps NVT ramp 50→150→{temperature_k:.0f} K at 0.5 fs", flush=True)

    # Add barostat for NPT and restore the production timestep.
    barostat = openmm.MonteCarloBarostat(pressure_atm * u.atmosphere, temperature_k * u.kelvin, 25)
    if seed is not None:
        barostat.setRandomNumberSeed(int(seed))
    system.addForce(barostat)
    sim.context.reinitialize(preserveState=True)
    integrator.setStepSize(timestep_fs * u.femtosecond)

    # Fresh production velocities — seeded per replica so each replica is an
    # independent draw from the Maxwell-Boltzmann distribution.
    if seed is not None:
        sim.context.setVelocitiesToTemperature(temperature_k * u.kelvin, int(seed))
    else:
        sim.context.setVelocitiesToTemperature(temperature_k * u.kelvin)

    steps_per_ps = int(round(1000.0 / timestep_fs))
    report_every = max(1, int(round(report_interval_ps * steps_per_ps)))

    equil_steps = int(round(equil_ps * steps_per_ps))
    prod_steps = int(round(prod_ps * steps_per_ps))

    sim.reporters.append(
        app.StateDataReporter(
            sys.stdout,
            report_every,
            step=True,
            potentialEnergy=True,
            temperature=True,
            volume=True,
            speed=True,
            elapsedTime=True,
        )
    )

    print(f"  equilibration NPT ({equil_ps:.0f} ps, {equil_steps} steps)", flush=True)
    sim.step(equil_steps)

    top = topology
    if traj_out is not None and traj_interval_ps > 0:
        # Production with periodic frame dumps — the dwell-time channel
        # scores the shape at each Δt, not just the endpoint. Heavy-atom
        # protein+ligand frames only (waters/ions stripped per frame).
        frame_steps = max(1, int(round(traj_interval_ps * steps_per_ps)))
        n_frames = max(1, prod_steps // frame_steps)
        writer = HeavyTrajectoryWriter(traj_out, top)
        print(
            f"  production NPT ({prod_ps:.0f} ps, {prod_steps} steps) → "
            f"{n_frames} frames every {traj_interval_ps:.0f} ps to {traj_out.name}",
            flush=True,
        )
        for fi in range(n_frames):
            sim.step(frame_steps)
            fr = sim.context.getState(getPositions=True, enforcePeriodicBox=False)
            writer.write_frame(fr.getPositions(asNumpy=True))
        writer.close()
    else:
        print(f"  production NPT ({prod_ps:.0f} ps, {prod_steps} steps)", flush=True)
        sim.step(prod_steps)

    # Strip waters/ions/barostat, dump the final relaxed protein + ligand.
    final = sim.context.getState(getPositions=True, enforcePeriodicBox=False)
    positions_out = final.getPositions(asNumpy=True)

    write_protein_ligand_pdb(top, positions_out, out_pdb)
    print(f"  wrote {out_pdb} ({(time.time() - t0)/60:.1f} min total)", flush=True)


# -----------------------------------------------------------------------------
# Resumable chunk dynamics (portable State; distributable chunks — issue #34).
#
# A replica's production is split so it can be advanced in small, independent
# steps that any machine can run with pip-only OpenMM:
#   build (--prepare-only)  → system.xml + solvated.pdb        (once per shape/ligand)
#   equilibrate_chunk       → state_0.xml                       (per replica seed)
#   segment_chunk × K       → state_{i+1}.xml + seg_i frames    (chain within a replica)
# Continuation state travels as a serialised OpenMM State (positions/velocities/
# box) — portable, unlike a machine-specific binary Checkpoint. The System
# (forces) travels separately as system.xml. The barostat is re-added identically
# in both steps so the physics matches; the Langevin noise is memoryless, so a
# fresh per-step --seed is a valid continuation, not a discontinuity.
# -----------------------------------------------------------------------------

def serialize_state(sim, state_xml: Path) -> Path:
    """Serialise positions + velocities + box to portable State XML so another
    machine can continue the run with pip-only OpenMM."""
    state = sim.context.getState(getPositions=True, getVelocities=True)
    Path(state_xml).parent.mkdir(parents=True, exist_ok=True)
    Path(state_xml).write_text(openmm.XmlSerializer.serialize(state), encoding="utf-8")
    return Path(state_xml)


def restore_state(sim, state_xml: Path) -> None:
    state = openmm.XmlSerializer.deserialize(Path(state_xml).read_text(encoding="utf-8"))
    sim.context.setState(state)


def add_production_barostat(system, *, pressure_atm: float, temperature_k: float, seed: int | None):
    """Add the same MonteCarloBarostat run_dynamics uses, so equilibrate/segment
    chunks share identical NPT physics."""
    barostat = openmm.MonteCarloBarostat(pressure_atm * u.atmosphere, temperature_k * u.kelvin, 25)
    if seed is not None:
        barostat.setRandomNumberSeed(int(seed))
    system.addForce(barostat)
    return barostat


def equilibrate_chunk(
    state_out: Path,
    system_xml: Path,
    solvated_pdb: Path,
    *,
    equil_ps: float,
    temperature_k: float,
    pressure_atm: float,
    timestep_fs: float,
    report_interval_ps: float,
    seed: int | None,
    t0: float | None = None,
) -> None:
    """Chunk step 'equilibrate': from a built system (system_xml + solvated_pdb),
    minimise → NVT warm-up → add barostat → NPT equilibrate, then serialise the
    end State to state_out. Velocity-seeded per replica. Pip-only OpenMM."""
    if t0 is None:
        t0 = time.time()
    print(f"=== md_relax (equilibrate chunk): {Path(system_xml).name} ===", flush=True)
    topology, positions, system = load_prepared_system(Path(system_xml), Path(solvated_pdb))
    print(f"  loaded {system.getNumParticles()} particles", flush=True)

    integrator = openmm.LangevinMiddleIntegrator(
        temperature_k * u.kelvin, 1.0 / u.picosecond, timestep_fs * u.femtosecond)
    if seed is not None:
        integrator.setRandomNumberSeed(int(seed))
    platform, props = pick_platform()
    print(f"  platform: {platform.getName()} (props={props})", flush=True)
    sim = app.Simulation(topology, system, integrator, platform, props)
    sim.context.setPositions(positions)

    print("  minimizing...", flush=True)
    sim.minimizeEnergy(maxIterations=10000)

    # NVT thermalization ramp at small dt (mirrors run_dynamics), before NPT.
    warmup_dt_fs = 0.5
    warmup_steps_per_ps = int(round(1000.0 / warmup_dt_fs))
    integrator.setStepSize(warmup_dt_fs * u.femtosecond)
    for warmup_t in (50.0, 150.0, temperature_k):
        integrator.setTemperature(warmup_t * u.kelvin)
        sim.context.setVelocitiesToTemperature(warmup_t * u.kelvin)
        sim.step(warmup_steps_per_ps)
    print(f"  warmup: 3 ps NVT ramp 50→150→{temperature_k:.0f} K at 0.5 fs", flush=True)

    add_production_barostat(system, pressure_atm=pressure_atm, temperature_k=temperature_k, seed=seed)
    sim.context.reinitialize(preserveState=True)
    integrator.setStepSize(timestep_fs * u.femtosecond)
    if seed is not None:
        sim.context.setVelocitiesToTemperature(temperature_k * u.kelvin, int(seed))
    else:
        sim.context.setVelocitiesToTemperature(temperature_k * u.kelvin)

    steps_per_ps = int(round(1000.0 / timestep_fs))
    report_every = max(1, int(round(report_interval_ps * steps_per_ps)))
    equil_steps = int(round(equil_ps * steps_per_ps))
    sim.reporters.append(app.StateDataReporter(
        sys.stdout, report_every, step=True, potentialEnergy=True,
        temperature=True, volume=True, speed=True, elapsedTime=True))
    print(f"  equilibration NPT ({equil_ps:.0f} ps, {equil_steps} steps)", flush=True)
    sim.step(equil_steps)

    serialize_state(sim, Path(state_out))
    print(f"  wrote {Path(state_out).name} ({(time.time() - t0)/60:.1f} min total)", flush=True)


def segment_chunk(
    state_out: Path,
    seg_out: Path,
    system_xml: Path,
    solvated_pdb: Path,
    state_in: Path,
    *,
    segment_ps: float,
    traj_interval_ps: float,
    temperature_k: float,
    pressure_atm: float,
    timestep_fs: float,
    report_interval_ps: float,
    seed: int | None,
    t0: float | None = None,
) -> None:
    """Chunk step 'segment': continue a replica from state_in for segment_ps,
    dumping heavy-atom frames to seg_out, then serialise state_out. No minimise/
    warm-up/equilibrate — velocities come from state_in. Pip-only OpenMM."""
    if t0 is None:
        t0 = time.time()
    print(f"=== md_relax (segment chunk): {Path(state_in).name} → {Path(state_out).name} ===", flush=True)
    topology, positions, system = load_prepared_system(Path(system_xml), Path(solvated_pdb))
    add_production_barostat(system, pressure_atm=pressure_atm, temperature_k=temperature_k, seed=seed)

    integrator = openmm.LangevinMiddleIntegrator(
        temperature_k * u.kelvin, 1.0 / u.picosecond, timestep_fs * u.femtosecond)
    if seed is not None:
        integrator.setRandomNumberSeed(int(seed))
    platform, props = pick_platform()
    print(f"  platform: {platform.getName()} (props={props})", flush=True)
    sim = app.Simulation(topology, system, integrator, platform, props)
    sim.context.setPositions(positions)
    restore_state(sim, Path(state_in))

    steps_per_ps = int(round(1000.0 / timestep_fs))
    seg_steps = int(round(segment_ps * steps_per_ps))
    frame_steps = max(1, int(round(traj_interval_ps * steps_per_ps)))
    n_frames = max(1, seg_steps // frame_steps)
    report_every = max(1, int(round(report_interval_ps * steps_per_ps)))
    sim.reporters.append(app.StateDataReporter(
        sys.stdout, report_every, step=True, potentialEnergy=True,
        temperature=True, volume=True, speed=True, elapsedTime=True))

    writer = HeavyTrajectoryWriter(Path(seg_out), topology)
    print(f"  segment NPT ({segment_ps:.0f} ps, {seg_steps} steps) → "
          f"{n_frames} frames every {traj_interval_ps:.0f} ps to {Path(seg_out).name}", flush=True)
    for _ in range(n_frames):
        sim.step(frame_steps)
        fr = sim.context.getState(getPositions=True, enforcePeriodicBox=False)
        writer.write_frame(fr.getPositions(asNumpy=True))
    writer.close()

    serialize_state(sim, Path(state_out))
    print(f"  wrote {Path(seg_out).name} + {Path(state_out).name} "
          f"({(time.time() - t0)/60:.1f} min total)", flush=True)


# -----------------------------------------------------------------------------
# Main relaxation pipeline.
# -----------------------------------------------------------------------------

def relax(
    out_pdb: Path,
    complex_pdb: Path | None = None,
    ligand_smiles: str | None = None,
    apo_pdb: Path | None = None,
    equil_ps: float = 100.0,
    prod_ps: float = 100.0,
    temperature_k: float = 300.0,
    pressure_atm: float = 1.0,
    salt_mol: float = 0.15,
    padding_nm: float = 1.0,
    timestep_fs: float = 2.0,
    report_interval_ps: float = 10.0,
    collapse_ps: float = 0.0,
    vacuum_min_iter: int = 0,
    no_explicit: bool = False,
    restrain_range: tuple[int, int] | None = None,
    restrain_chains: set[str] | None = None,
    restrain_k: float = 1000.0,
    seed: int | None = None,
    traj_out: Path | None = None,
    traj_interval_ps: float = 0.0,
    rectangular_box: bool = False,
    system_xml: Path | None = None,
    solvated_pdb: Path | None = None,
    prepare_prefix: Path | None = None,
):
    t0 = time.time()

    # Edge path: run a pre-parametrised, serialised system with pip-only
    # OpenMM — no OpenFF, no conda, even for a docked complex (issue #34).
    if system_xml is not None:
        print(f"=== md_relax (prepared system): {Path(system_xml).name} ===", flush=True)
        topology, positions, system = load_prepared_system(Path(system_xml), Path(solvated_pdb))
        print(f"  loaded {system.getNumParticles()} particles from serialised "
              f"System (no OpenFF)", flush=True)
        run_dynamics(
            out_pdb, topology, positions, system,
            equil_ps=equil_ps, prod_ps=prod_ps, temperature_k=temperature_k,
            pressure_atm=pressure_atm, timestep_fs=timestep_fs,
            report_interval_ps=report_interval_ps, restrain_range=restrain_range,
            restrain_chains=restrain_chains, restrain_k=restrain_k, seed=seed,
            traj_out=traj_out, traj_interval_ps=traj_interval_ps, t0=t0,
        )
        return

    if complex_pdb is not None:
        if collapse_ps > 0:
            raise ValueError("--collapse-ps is apo-only (ligand has no OBC2 GB parameters)")
        print(f"=== md_relax (complex): {complex_pdb.name} ===", flush=True)
        prot_pdb_text, lig_pdb_text = split_complex_pdb(complex_pdb)
        print("  parsed complex PDB", flush=True)
        rdkit_lig = ligand_from_pdb_and_smiles(lig_pdb_text, ligand_smiles)
        print(f"  ligand: {rdkit_lig.GetNumAtoms()} atoms (incl Hs); SMILES = {ligand_smiles}", flush=True)
        off_lig = offmol_from_rdkit(rdkit_lig)
        qsum = float(sum(c.m for c in off_lig.partial_charges))
        print(f"  NAGL charges assigned, sum = {qsum:+.3f} e", flush=True)
    elif apo_pdb is not None:
        print(f"=== md_relax (apo): {apo_pdb.name} ===", flush=True)
        prot_pdb_text = load_apo_pdb(apo_pdb)
        print("  parsed apo PDB (no ligand)", flush=True)
        off_lig = None
        if vacuum_min_iter > 0:
            prot_pdb_text = vacuum_minimize_apo(
                prot_pdb_text,
                max_iterations=vacuum_min_iter,
                restrain_range=restrain_range,
                restrain_chains=restrain_chains,
                restrain_k=restrain_k,
            )
        if collapse_ps > 0:
            prot_pdb_text = implicit_collapse_apo(
                prot_pdb_text,
                collapse_ps=collapse_ps,
                temperature_k=temperature_k,
                timestep_fs=timestep_fs,
                restrain_range=restrain_range,
                restrain_chains=restrain_chains,
                restrain_k=restrain_k,
            )
        if no_explicit:
            print("  --no-explicit: writing post-collapse structure (heavy atoms only)", flush=True)
            write_heavy_only_pdb(prot_pdb_text, out_pdb)
            print(f"  wrote {out_pdb} ({(time.time() - t0)/60:.1f} min total)", flush=True)
            return
    else:
        raise ValueError("md_relax requires --complex-pdb (+ --ligand-smiles) or --apo-pdb")

    modeller, system = build_system(prot_pdb_text, off_lig, padding_nm=padding_nm, salt_mol=salt_mol,
                                    rectangular_box=rectangular_box)

    # Prepare-only: serialise the parametrised, solvated system and stop. The
    # GPU-free, OpenFF-using half of a docked-complex chunk (issue #34). The
    # serialised System then runs anywhere with pip-only OpenMM via
    # --system-xml/--solvated-pdb.
    if prepare_prefix is not None:
        sys_xml, solv_pdb = serialize_prepared_system(
            system, modeller.topology, modeller.positions, Path(prepare_prefix),
        )
        print(f"  prepared (no dynamics): wrote {sys_xml.name} + {solv_pdb.name}; "
              f"run anywhere with --system-xml {sys_xml.name} "
              f"--solvated-pdb {solv_pdb.name}", flush=True)
        return

    run_dynamics(
        out_pdb, modeller.topology, modeller.positions, system,
        equil_ps=equil_ps, prod_ps=prod_ps, temperature_k=temperature_k,
        pressure_atm=pressure_atm, timestep_fs=timestep_fs,
        report_interval_ps=report_interval_ps, restrain_range=restrain_range,
        restrain_chains=restrain_chains, restrain_k=restrain_k, seed=seed,
        traj_out=traj_out, traj_interval_ps=traj_interval_ps, t0=t0,
    )


SOLVENT_RESNAMES = {"HOH", "WAT", "TIP3", "TIP", "NA", "CL", "Na+", "Cl-", "K+", "MG", "ZN", "CA"}


def _is_ligand_residue(residue) -> bool:
    """Detect the LIG residue regardless of whatever name openff/openmm gave it.
    Anything that isn't a standard amino acid (recognised by having a Cα) and
    isn't solvent/ion is treated as the ligand."""
    if residue.name in SOLVENT_RESNAMES:
        return False
    atom_names = {a.name for a in residue.atoms()}
    return "CA" not in atom_names


def _heavy(atom) -> bool:
    # Stage 2 anchor features were calibrated on RCSB PDBs that have no
    # explicit hydrogens (cryo-EM α-syn entries are heavy-atom only). Strip
    # Hs from the relaxed output so SASA / contact computations remain on
    # the same axis. MD ran with explicit Hs (correct physics); the saved
    # PDB just drops them for downstream feature recompute.
    return atom.element is not None and atom.element.symbol != "H"


def build_heavy_subtopology(topology: app.Topology) -> tuple[app.Topology, list[int]]:
    """Build the heavy-atom protein + ligand sub-topology and the ordered
    list of original atom indices that map into it. Protein chains keep
    their ids; the ligand (any non-AA, non-solvent residue) is collected
    onto chain Z residue LIG. Returns (sub_top, add_order) so a single
    frame or a whole trajectory can be written against one fixed topology
    (atom order is identical frame-to-frame).

    The Stage 2 / Stage 3 feature recompute keys off chain Z residue 'LIG'
    to detect appended ligand atoms (see features.py
    _stage3_ligand_heavy_coords), so we restore that convention here."""
    ligand_res_atoms: list = []
    for residue in topology.residues():
        if residue.name in SOLVENT_RESNAMES:
            continue
        if _is_ligand_residue(residue):
            ligand_res_atoms.extend(residue.atoms())

    sub_top = app.Topology()
    atom_map: dict[int, app.topology.Atom] = {}
    add_order: list[int] = []
    used_chain_ids: set[str] = set()
    for chain in topology.chains():
        chain_residues = [r for r in chain.residues() if not _is_ligand_residue(r) and r.name not in SOLVENT_RESNAMES]
        if not chain_residues:
            continue
        new_chain = sub_top.addChain(chain.id)
        used_chain_ids.add(chain.id)
        for residue in chain_residues:
            new_res = sub_top.addResidue(residue.name, new_chain, residue.id, residue.insertionCode)
            for a in residue.atoms():
                if not _heavy(a):
                    continue
                atom_map[a.index] = sub_top.addAtom(a.name, a.element, new_res)
                add_order.append(a.index)

    if ligand_res_atoms:
        lig_chain_id = LIGAND_CHAIN
        suffix = 0
        while lig_chain_id in used_chain_ids:
            suffix += 1
            lig_chain_id = f"{LIGAND_CHAIN}{suffix}"
        lig_chain = sub_top.addChain(lig_chain_id)
        lig_res = sub_top.addResidue(LIGAND_RESNAME, lig_chain, "1", " ")
        for a in ligand_res_atoms:
            if not _heavy(a):
                continue
            atom_map[a.index] = sub_top.addAtom(a.name, a.element, lig_res)
            add_order.append(a.index)

    for bond in topology.bonds():
        a, b = bond[0], bond[1]
        if a.index in atom_map and b.index in atom_map:
            sub_top.addBond(atom_map[a.index], atom_map[b.index])

    return sub_top, add_order


def _subset_positions(positions, add_order: list[int]):
    return np.array([positions[i].value_in_unit(u.nanometer) for i in add_order]) * u.nanometer


def write_protein_ligand_pdb(topology: app.Topology, positions, out_pdb: Path) -> None:
    """Write only protein chains + ligand renamed to chain Z residue LIG; drop
    water/ions."""
    sub_top, add_order = build_heavy_subtopology(topology)
    sub_pos = _subset_positions(positions, add_order)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    with out_pdb.open("w") as f:
        app.PDBFile.writeFile(sub_top, sub_pos, f, keepIds=True)


class HeavyTrajectoryWriter:
    """Append heavy-atom protein+ligand frames to a multi-MODEL PDB during
    production MD. The sub-topology is built once from the full system
    topology, so every MODEL shares the same atom order — dwell_time.py
    reads each MODEL back as one frame and scores it with shape_metrics.

    Frames carry only heavy atoms on the same chain/residue convention as
    write_protein_ligand_pdb (chain Z LIG ligand), so per-frame Stage 2 /
    shape scoring sees exactly what the single-frame relaxed PDB shows."""

    def __init__(self, out_pdb: Path, topology: app.Topology):
        self.sub_top, self.add_order = build_heavy_subtopology(topology)
        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        self._fh = out_pdb.open("w")
        app.PDBFile.writeHeader(self.sub_top, self._fh)
        self._n = 0

    def write_frame(self, positions) -> None:
        self._n += 1
        sub_pos = _subset_positions(positions, self.add_order)
        app.PDBFile.writeModel(self.sub_top, sub_pos, self._fh, modelIndex=self._n, keepIds=True)

    def close(self) -> None:
        app.PDBFile.writeFooter(self.sub_top, self._fh)
        self._fh.close()


def main() -> None:
    # Progress lines use non-ASCII (β, Å, →); force UTF-8 so a redirected or
    # subprocess-captured stdout on a Windows cp1252 console doesn't crash the
    # whole MD run on a print. dwell_time.py captures this stdout, so the same
    # fix is needed there for the pilot.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if hasattr(sys.stderr, "reconfigure"):
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--complex-pdb", type=Path, default=None,
                   help="apo + docked LIG on chain Z (mutually exclusive with --apo-pdb)")
    p.add_argument("--ligand-smiles", type=str, default=None,
                   help="SMILES of the ligand; required with --complex-pdb")
    p.add_argument("--apo-pdb", type=Path, default=None,
                   help="protein-only PDB for apo MD baseline (mutually exclusive with --complex-pdb)")
    p.add_argument("--system-xml", type=Path, default=None,
                   help="run dynamics from a serialised OpenMM System (from "
                        "--prepare-only) instead of building one. Needs only "
                        "pip-installed openmm — no OpenFF/conda, even for a "
                        "docked complex. Requires --solvated-pdb.")
    p.add_argument("--solvated-pdb", type=Path, default=None,
                   help="solvated topology+positions paired with --system-xml")
    p.add_argument("--prepare-only", type=Path, default=None, metavar="PREFIX",
                   help="build + solvate + parametrise, then serialise to "
                        "PREFIX_system.xml + PREFIX_solvated.pdb and exit "
                        "(no dynamics). The OpenFF/conda half of a chunk; the "
                        "serialised System runs with pip-only OpenMM.")
    p.add_argument("--out-pdb", type=Path, default=None,
                   help="output relaxed PDB (required unless --prepare-only)")
    p.add_argument("--equil-ps", type=float, default=100.0)
    p.add_argument("--prod-ps", type=float, default=100.0)
    p.add_argument("--temperature-k", type=float, default=300.0)
    p.add_argument("--pressure-atm", type=float, default=1.0)
    p.add_argument("--salt-mol", type=float, default=0.15)
    p.add_argument("--padding-nm", type=float, default=1.0)
    p.add_argument("--timestep-fs", type=float, default=2.0)
    p.add_argument("--report-interval-ps", type=float, default=10.0)
    p.add_argument("--vacuum-min-iter", type=int, default=0,
                   help="iterations of vacuum (no solvent) energy minimization "
                        "before solvation; resolves clashes that solvent-trapped "
                        "minimization can't fix. Apo-mode only.")
    p.add_argument("--no-explicit", action="store_true",
                   help="skip explicit-solvent solvation + production. Outputs "
                        "the structure after vacuum-min + implicit-collapse "
                        "stages. Use for oligomer scoring where Stage 2 features "
                        "only need a relaxed structure, not equilibrium dynamics.")
    p.add_argument("--collapse-ps", type=float, default=0.0,
                   help="ps of OBC2 implicit-solvent MD before explicit "
                        "solvation; apo-mode only (no GB params for ligands)")
    p.add_argument("--restrain-residues", type=str, default=None,
                   help="restrain Cα of these residues (e.g. '70-88') to "
                        "their initial positions during collapse and explicit "
                        "MD; preserves a hand-built β-core")
    p.add_argument("--restrain-chains", type=str, default=None,
                   help="comma-separated chain ids the restraint applies to "
                        "(default: all chains)")
    p.add_argument("--restrain-k", type=float, default=1000.0,
                   help="harmonic restraint constant in kJ/mol/nm²")
    p.add_argument("--seed", type=int, default=None,
                   help="random seed for the Langevin integrator, barostat, and "
                        "production velocities. Same seed → identical trajectory; "
                        "different seed → an independent thermal replica. Used by "
                        "the dwell-time channel (dwell_time.py) for velocity-seeded "
                        "replicas. Default: OpenMM picks its own (single-shot mode).")
    p.add_argument("--traj-out", type=Path, default=None,
                   help="write a multi-MODEL PDB of heavy-atom protein+ligand "
                        "frames sampled every --traj-interval-ps during production. "
                        "Each MODEL is one dwell-time frame.")
    p.add_argument("--traj-interval-ps", type=float, default=0.0,
                   help="ps between trajectory frames (requires --traj-out)")
    p.add_argument("--rect-box", action="store_true",
                   help="solvate in a tight rectangular box (solute bounding box "
                        "+ padding per axis) instead of a cube sized to the largest "
                        "dimension. ~7× fewer atoms for an elongated β-strand "
                        "construct — keeps one dwell chunk on a basic GPU.")
    # Resumable chunk modes (distributable per-replica steps; issue #34). Both
    # consume a built system (--system-xml + --solvated-pdb from --prepare-only).
    p.add_argument("--equilibrate", type=Path, default=None, metavar="STATE_OUT",
                   help="chunk step: minimise + NVT warm-up + NPT equilibrate the "
                        "built system, serialise the end State to STATE_OUT. "
                        "Requires --system-xml/--solvated-pdb; velocity-seeded via --seed.")
    p.add_argument("--segment", action="store_true",
                   help="chunk step: continue a replica from --state-in for "
                        "--segment-ps, dump frames to --seg-out, serialise --state-out. "
                        "Requires --system-xml/--solvated-pdb/--state-in/--state-out/--seg-out.")
    p.add_argument("--state-in", type=Path, default=None,
                   help="input serialised State for --segment (from equilibrate or a prior segment)")
    p.add_argument("--state-out", type=Path, default=None,
                   help="output serialised State for --segment")
    p.add_argument("--seg-out", type=Path, default=None,
                   help="output multi-MODEL frame PDB for --segment")
    p.add_argument("--segment-ps", type=float, default=100.0,
                   help="ps of production per --segment step")
    args = p.parse_args()

    # --- Resumable chunk dispatch (pip-only; needs a built --system-xml). ---
    if args.equilibrate is not None or args.segment:
        if args.system_xml is None or args.solvated_pdb is None:
            p.error("--equilibrate/--segment require --system-xml and --solvated-pdb")
        if args.equilibrate is not None:
            equilibrate_chunk(
                args.equilibrate, args.system_xml, args.solvated_pdb,
                equil_ps=args.equil_ps, temperature_k=args.temperature_k,
                pressure_atm=args.pressure_atm, timestep_fs=args.timestep_fs,
                report_interval_ps=args.report_interval_ps, seed=args.seed,
            )
            return
        if not all((args.state_in, args.state_out, args.seg_out)):
            p.error("--segment requires --state-in, --state-out and --seg-out")
        if args.traj_interval_ps <= 0:
            p.error("--segment requires a positive --traj-interval-ps")
        segment_chunk(
            args.state_out, args.seg_out, args.system_xml, args.solvated_pdb, args.state_in,
            segment_ps=args.segment_ps, traj_interval_ps=args.traj_interval_ps,
            temperature_k=args.temperature_k, pressure_atm=args.pressure_atm,
            timestep_fs=args.timestep_fs, report_interval_ps=args.report_interval_ps,
            seed=args.seed,
        )
        return

    if args.traj_out is not None and args.traj_interval_ps <= 0:
        p.error("--traj-out requires a positive --traj-interval-ps")

    n_modes = sum(x is not None for x in (args.complex_pdb, args.apo_pdb, args.system_xml))
    if n_modes != 1:
        p.error("exactly one of --complex-pdb, --apo-pdb, or --system-xml must be provided")
    if args.complex_pdb is not None and not args.ligand_smiles:
        p.error("--ligand-smiles is required with --complex-pdb")
    if args.system_xml is not None and args.solvated_pdb is None:
        p.error("--system-xml requires --solvated-pdb")
    if args.prepare_only is not None and args.system_xml is not None:
        p.error("--prepare-only builds a system; it cannot combine with --system-xml")
    if args.prepare_only is None and args.out_pdb is None:
        p.error("--out-pdb is required (unless --prepare-only)")

    relax(
        complex_pdb=args.complex_pdb,
        ligand_smiles=args.ligand_smiles,
        apo_pdb=args.apo_pdb,
        out_pdb=args.out_pdb,
        system_xml=args.system_xml,
        solvated_pdb=args.solvated_pdb,
        prepare_prefix=args.prepare_only,
        equil_ps=args.equil_ps,
        prod_ps=args.prod_ps,
        temperature_k=args.temperature_k,
        pressure_atm=args.pressure_atm,
        salt_mol=args.salt_mol,
        padding_nm=args.padding_nm,
        timestep_fs=args.timestep_fs,
        report_interval_ps=args.report_interval_ps,
        collapse_ps=args.collapse_ps,
        vacuum_min_iter=args.vacuum_min_iter,
        no_explicit=args.no_explicit,
        restrain_range=parse_residue_range(args.restrain_residues),
        restrain_chains=parse_chain_set(args.restrain_chains),
        restrain_k=args.restrain_k,
        seed=args.seed,
        traj_out=args.traj_out,
        traj_interval_ps=args.traj_interval_ps,
        rectangular_box=args.rect_box,
    )


if __name__ == "__main__":
    main()
