"""Tests for the dwell-time channel (#14): shape metrics + bootstrap.

No pytest dependency required — run directly with the pipeline venv:

    .venv/bin/python tests/test_dwell_time.py

The functions are also named ``test_*`` so pytest discovers them if it is
ever added. They use the committed ``results/oligomers/*_relaxed.pdb``
structures as fixtures, so they are deterministic and need no MD.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))

import shape_metrics as sm  # noqa: E402
from dwell_time import bootstrap_dwell_shift, score_trajectory  # noqa: E402

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
    # The whole point of #14: a positive shift is detectable too (#30).
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
