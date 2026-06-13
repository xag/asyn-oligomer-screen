"""Serialised-System round-trip for md_relax.

The crowdsourcing claim for docked-complex dwell chunks rests on one fact:
OpenFF is needed only to *build* the force field (one-time, GPU-free), not to
*run* dynamics. So the central side serialises a fully-parametrised OpenMM
``System`` and any pip-only-OpenMM machine integrates it. These tests pin that
contract:

  * serialisation is lossless — the deserialised System gives identical
    energy/forces at identical coordinates (no parameter is dropped);
  * the deserialised System is actually runnable through the same
    ``run_dynamics`` path the build side uses, producing the expected output +
    trajectory frames;
  * a system.xml / solvated.pdb pair that don't match is rejected, not run.

They build a tiny toy System directly with OpenMM, so they need neither OpenFF
nor pdbfixer/amber — only ``openmm`` (the ``md`` dependency group). Importing
md_relax pulls in pdbfixer at module load, so that is import-skipped too; the
tests run under ``uv run --group md`` and skip cleanly in the bare pip venv.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

openmm = pytest.importorskip("openmm")
pytest.importorskip("pdbfixer")  # md_relax imports it at module load
import openmm.app as app  # noqa: E402
import openmm.unit as u  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "screen"))
import md_relax  # noqa: E402


def _toy_system_and_topology(n: int = 8, box_nm: float = 2.5):
    """A small periodic LJ system + matching single-chain Topology. Particles
    sit on a loose grid (no steep contacts) so PDB coordinate truncation does
    not swing the LJ energy. Carbon atoms named C1.. so none reads as a Cα."""
    system = openmm.System()
    top = app.Topology()
    chain = top.addChain("A")
    carbon = app.element.carbon

    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(0.9 * u.nanometer)

    # grid positions, spacing 0.5 nm, inside the box
    coords = []
    side = int(np.ceil(n ** (1 / 3)))
    i = 0
    for a in range(side):
        for b in range(side):
            for c in range(side):
                if i >= n:
                    break
                coords.append((0.4 + 0.5 * a, 0.4 + 0.5 * b, 0.4 + 0.5 * c))
                i += 1
    for idx in range(n):
        system.addParticle(12.011 * u.amu)
        res = top.addResidue("LIG", chain)
        top.addAtom(f"C{idx}", carbon, res)
        nb.addParticle(0.0 * u.elementary_charge, 0.34 * u.nanometer, 0.4 * u.kilojoule_per_mole)
    system.addForce(nb)

    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box_nm, 0, 0) * u.nanometer,
        openmm.Vec3(0, box_nm, 0) * u.nanometer,
        openmm.Vec3(0, 0, box_nm) * u.nanometer,
    )
    # Topology wants a single Quantity wrapping three Vec3 (unit on the tuple).
    top.setPeriodicBoxVectors(
        (openmm.Vec3(box_nm, 0, 0), openmm.Vec3(0, box_nm, 0), openmm.Vec3(0, 0, box_nm)) * u.nanometer
    )
    positions = np.array(coords) * u.nanometer
    return system, top, positions


def _potential_kj(system, positions) -> float:
    integ = openmm.VerletIntegrator(1.0 * u.femtosecond)
    ctx = openmm.Context(system, integ, openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(positions)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    del ctx, integ
    return e


def test_serialized_system_is_lossless(tmp_path):
    system, top, pos = _toy_system_and_topology()
    sx, sp = md_relax.serialize_prepared_system(system, top, pos, tmp_path / "chunk")
    assert sx.exists() and sp.exists()
    assert sx.name == "chunk_system.xml" and sp.name == "chunk_solvated.pdb"

    top2, pos2, system2 = md_relax.load_prepared_system(sx, sp)
    assert system2.getNumParticles() == system.getNumParticles()
    assert top2.getNumAtoms() == top.getNumAtoms()

    # Same parameters → identical energy at identical coordinates. Evaluate both
    # systems at pos2 to isolate serialisation from PDB coordinate truncation:
    # if any force-field term were lost, this would differ.
    assert _potential_kj(system, pos2) == pytest.approx(_potential_kj(system2, pos2), rel=1e-9, abs=1e-7)

    # Box vectors survive the round-trip on both the System and the Topology.
    a0 = system.getDefaultPeriodicBoxVectors()[0][0].value_in_unit(u.nanometer)
    a2 = system2.getDefaultPeriodicBoxVectors()[0][0].value_in_unit(u.nanometer)
    assert a2 == pytest.approx(a0, abs=1e-3)


def test_prepared_system_runs_through_run_dynamics(tmp_path):
    system, top, pos = _toy_system_and_topology(n=8)
    sx, sp = md_relax.serialize_prepared_system(system, top, pos, tmp_path / "chunk")
    top2, pos2, system2 = md_relax.load_prepared_system(sx, sp)

    out = tmp_path / "relaxed.pdb"
    traj = tmp_path / "traj.pdb"
    # Tiny budget so the test is fast: 0.04 ps prod / 0.02 ps interval = 2 frames.
    md_relax.run_dynamics(
        out, top2, pos2, system2,
        equil_ps=0.02, prod_ps=0.04, temperature_k=120.0, pressure_atm=1.0,
        timestep_fs=1.0, report_interval_ps=0.02,
        restrain_range=None, restrain_chains=None, restrain_k=1000.0,
        seed=7, traj_out=traj, traj_interval_ps=0.02,
    )
    assert out.exists()
    # All 8 heavy atoms make it into the relaxed output (chain Z LIG).
    n_atoms = sum(1 for ln in out.read_text().splitlines() if ln.startswith(("ATOM", "HETATM")))
    assert n_atoms == 8
    assert traj.read_text().count("MODEL ") == 2


def test_mismatched_pair_is_rejected(tmp_path):
    sys_a, top_a, pos_a = _toy_system_and_topology(n=6)
    sys_b, top_b, pos_b = _toy_system_and_topology(n=9)
    sx_a, _ = md_relax.serialize_prepared_system(sys_a, top_a, pos_a, tmp_path / "a")
    _, sp_b = md_relax.serialize_prepared_system(sys_b, top_b, pos_b, tmp_path / "b")
    with pytest.raises(ValueError, match="not a matched pair"):
        md_relax.load_prepared_system(sx_a, sp_b)  # 6-particle system vs 9-atom pdb
