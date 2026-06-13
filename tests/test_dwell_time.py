"""Tests for the dwell-time channel: shape metrics + bootstrap.

No pytest dependency required — run directly with the pipeline venv:

    .venv/bin/python tests/test_dwell_time.py

The functions are also named ``test_*`` so pytest discovers them if it is
ever added. They use the committed ``results/oligomers/*_relaxed.pdb``
structures as fixtures, so they are deterministic and need no MD.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure as BioStructure

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))

import shape_metrics as sm  # noqa: E402
from dwell_time import bootstrap_dwell_shift, score_trajectory, summarise_pair  # noqa: E402

OLIG = ROOT / "results" / "oligomers"
PARALLEL = OLIG / "fusco_parallel_3mer_core70-88_relaxed.pdb"
PARALLEL_S123 = OLIG / "fusco_parallel_3mer_core70-88_s123_relaxed.pdb"
ANTIPARALLEL = OLIG / "fusco_antiparallel_3mer_core70-88_relaxed.pdb"


# ---- shape metrics ---------------------------------------------------------

def test_identity_rmsd_zero_jaccard_one():
    s = sm.load_pdb(PARALLEL)
    assert sm.beta_core_rmsd(s, s) < 1e-6
    assert sm.contact_jaccard(s, s) == 1.0
    assert sm.in_toxic_basin(sm.frame_metrics(s, s)) is True


def test_same_topology_small_drift_stays_in_basin():
    ref = sm.load_pdb(PARALLEL)
    mob = sm.load_pdb(PARALLEL_S123)
    m = sm.frame_metrics(mob, ref)
    # Same build, different MD seed: real but modest drift.
    assert 0.0 < m["beta_core_rmsd"] < sm.TOXIC_RMSD_MAX
    assert m["contact_jaccard"] > sm.TOXIC_JACCARD_MIN
    assert sm.in_toxic_basin(m) is True


def test_cross_topology_leaves_basin():
    ref = sm.load_pdb(PARALLEL)
    mob = sm.load_pdb(ANTIPARALLEL)
    m = sm.frame_metrics(mob, ref)
    # Genuinely different shape: far in both metrics, out of the basin.
    assert m["beta_core_rmsd"] > 10.0
    assert m["contact_jaccard"] < 0.2
    assert sm.in_toxic_basin(m) is False


def test_contact_set_symmetric_and_interchain_only():
    s = sm.load_pdb(PARALLEL)
    contacts = sm.interchain_contact_set(s)
    assert len(contacts) > 0
    for ca, ra, cb, rb in contacts:
        assert ca != cb, "contacts must be inter-chain"
        assert (ca, ra) <= (cb, rb), "endpoints must be canonically ordered"


def test_apo_structure_has_no_ligand():
    s = sm.load_pdb(PARALLEL)
    assert sm.ligand_bound(next(iter(s))) is None  # no LIG residue → not applicable


# ---- trajectory scoring ----------------------------------------------------
#
# No MD has run yet (the pilot needs a GPU host), so there are no real
# multi-MODEL trajectory PDBs to score. We synthesise them from the
# committed single-frame relaxed structures: each source PDB becomes one
# MODEL of a trajectory, and a ligand residue (chain Z / resname LIG) is
# optionally grafted so the occupancy branch of score_trajectory is also
# exercised. This covers the exact code path the GPU output will flow
# through — multi-MODEL parse, per-frame dwell, bound-frame masking.


def _graft_ligand(model, mode: str) -> None:
    """Add a chain-Z LIG residue to ``model``. ``mode='bound'`` places a
    heavy atom on one of the model's own β-core Cα (so it reads as bound);
    ``mode='left'`` places it 100 Å away (diffused off the site)."""
    core = sm._ca_by_residue(model, sm.BETA_CORE_RANGE, chain_ids=None)
    anchor = np.array(next(iter(core.values())), dtype=float)
    coord = anchor if mode == "bound" else anchor + 100.0
    chain = Chain("Z")
    res = Residue((" ", 1, " "), "LIG", "")
    res.add(Atom("C1", coord, 0.0, 1.0, " ", "C1", 1, element="C"))
    chain.add(res)
    model.add(chain)


def _make_traj(tmp: Path, name: str, sources: list[Path], ligand: str = "none") -> Path:
    """Write a multi-MODEL trajectory PDB: one MODEL per source structure,
    optionally with a grafted ligand on every frame."""
    out = BioStructure("traj")
    for i, src in enumerate(sources):
        model = next(iter(sm.load_pdb(src))).copy()
        model.id = i
        model.serial_num = i + 1
        if ligand != "none":
            _graft_ligand(model, ligand)
        out.add(model)
    path = tmp / name
    io = PDBIO()
    io.set_structure(out)
    io.save(str(path))
    return path


def test_score_apo_trajectory_all_in_basin():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        traj = _make_traj(tmp, "apo.pdb", [PARALLEL, PARALLEL_S123, PARALLEL])
        ref = sm.load_pdb(PARALLEL)
        s = score_trajectory(traj, ref)
        assert s["n_frames"] == 3
        assert s["dwell_fraction"] == 1.0          # every frame stays in the basin
        assert np.isnan(s["occupancy"])            # apo: occupancy not applicable
        assert all(b is None for b in s["per_frame"]["bound"])


def test_score_mixed_trajectory_partial_dwell():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        traj = _make_traj(tmp, "mixed.pdb", [PARALLEL, ANTIPARALLEL, PARALLEL])
        ref = sm.load_pdb(PARALLEL)
        s = score_trajectory(traj, ref)
        # parallel frames dwell, the antiparallel frame does not → 2/3.
        assert s["per_frame"]["in_basin"] == [True, False, True]
        assert abs(s["dwell_fraction"] - 2 / 3) < 1e-9


def test_score_complex_bound_dwell_over_bound_frames():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        # Ligand stays bound but holds the oligomer in the (non-toxic)
        # antiparallel shape → bound everywhere, dwell over bound frames = 0.
        traj = _make_traj(tmp, "cpx.pdb", [ANTIPARALLEL, ANTIPARALLEL], ligand="bound")
        ref = sm.load_pdb(PARALLEL)
        s = score_trajectory(traj, ref)
        assert s["occupancy"] == 1.0
        assert s["n_bound"] == 2
        assert s["dwell_fraction"] == 0.0
        assert all(b is True for b in s["per_frame"]["bound"])


def test_score_complex_ligand_left_is_nan_dwell():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        traj = _make_traj(tmp, "left.pdb", [PARALLEL, PARALLEL], ligand="left")
        ref = sm.load_pdb(PARALLEL)
        s = score_trajectory(traj, ref)
        assert s["occupancy"] == 0.0               # ligand off-site in every frame
        assert np.isnan(s["dwell_fraction"])       # no bound frame → undefined


def test_summarise_pair_end_to_end_destabiliser():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        apo = [
            _make_traj(tmp, "apo0.pdb", [PARALLEL, PARALLEL_S123]),
            _make_traj(tmp, "apo1.pdb", [PARALLEL_S123, PARALLEL]),
        ]
        cpx = [
            _make_traj(tmp, "cpx0.pdb", [ANTIPARALLEL, ANTIPARALLEL], ligand="bound"),
            _make_traj(tmp, "cpx1.pdb", [ANTIPARALLEL, ANTIPARALLEL], ligand="bound"),
        ]
        pair = summarise_pair(PARALLEL, apo, cpx)
        b = pair["bootstrap"]
        # apo dwells ~1.0, bound complex dwells 0.0 → clean destabiliser, full occupancy.
        assert pair["apo_dwell"] == [1.0, 1.0]
        assert pair["complex_dwell"] == [0.0, 0.0]
        assert b["shift"] == -1.0
        assert b["classification"] == "destabiliser"
        assert b["occupancy_ok"] is True
        assert b["mean_complex_occupancy"] == 1.0


def test_summarise_pair_complex_without_ligand_flags_not_diffusion():
    # A complex trajectory with no LIG residue must report "no ligand
    # present", not be mistaken for a ligand that diffused off (occupancy
    # stays nan, ligand_present is False, and no empty-slice warning).
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        apo = [_make_traj(tmp, "apo.pdb", [PARALLEL, PARALLEL])]
        cpx = [_make_traj(tmp, "cpx.pdb", [PARALLEL, PARALLEL])]   # no ligand grafted
        b = summarise_pair(PARALLEL, apo, cpx)["bootstrap"]
        assert b["ligand_present"] is False
        assert np.isnan(b["mean_complex_occupancy"])
        assert b["occupancy_ok"] is False


# ---- bootstrap -------------------------------------------------------------

def test_bootstrap_destabiliser():
    apo = [0.95, 0.9, 0.92, 0.88, 0.97]
    cpx = [0.2, 0.3, 0.25, 0.15, 0.35]
    b = bootstrap_dwell_shift(apo, cpx, seed=1)
    assert b["shift"] < 0
    assert b["ci_high"] < 0
    assert b["classification"] == "destabiliser"
    assert b["prob_destabiliser"] > 0.95


def test_bootstrap_stabiliser_symmetric():
    # A positive shift is detectable too.
    apo = [0.45, 0.5, 0.55, 0.48, 0.52]
    cpx = [0.9, 0.95, 0.88, 0.92, 0.97]
    b = bootstrap_dwell_shift(apo, cpx, seed=1)
    assert b["shift"] > 0
    assert b["ci_low"] > 0
    assert b["classification"] == "stabiliser"
    assert b["prob_stabiliser"] > 0.95


def test_bootstrap_inconclusive_for_neutral():
    apo = [0.5, 0.55, 0.45, 0.52, 0.48]
    cpx = [0.49, 0.53, 0.47, 0.5, 0.51]
    b = bootstrap_dwell_shift(apo, cpx, seed=1)
    assert b["classification"] == "inconclusive"
    assert b["ci_low"] < 0 < b["ci_high"]


def test_bootstrap_drops_nan_replicas():
    # NaN dwell (ligand never stayed bound) must be excluded, not crash.
    apo = [0.9, 0.92, 0.88]
    cpx = [float("nan"), 0.3, 0.25]
    b = bootstrap_dwell_shift(apo, cpx, seed=1)
    assert b["n_complex"] == 2
    assert b["n_dropped_complex"] == 1
    assert b["classification"] == "destabiliser"


# ---- runner ----------------------------------------------------------------

def _run_all() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failures = 0
    for t in tests:
        try:
            t()
            print(f"  ok   {t.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"  FAIL {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            failures += 1
            print(f"  ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failures}/{len(tests)} passed")
    return 1 if failures else 0


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.exit(_run_all())
