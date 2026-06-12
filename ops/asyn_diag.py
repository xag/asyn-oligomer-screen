#!/usr/bin/env python3
"""CPU-only diagnostic on existing dwell-time trajectories — why did the gate
come back null? Pure geometry on the *_rep*.pdb files already on disk (no GPU,
no MD). Reports, per label:
  * occupancy  — did the ligand stay bound to the beta-core, or diffuse off?
  * dwell distribution + bootstrapped shift vs apo (reproduces the gate),
  * mean beta-core RMSD / contact Jaccard,
then how the basin thresholds sit relative to the actual apo frame distribution,
and a frame-resolved apo-vs-control time-course. Together these separate
"ligand never bound" from "basin/metric mis-set" from "underpowered".

Env: ASYN_REPO (required), SHAPE, LIGANDS (default "silibinin dhea caffeine").
"""
from __future__ import annotations

import glob
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(os.environ["ASYN_REPO"]).resolve()
sys.path.insert(0, str(REPO / "screen"))
os.chdir(REPO)

import dwell_time as dt  # noqa: E402
from dwell_time import score_trajectory, bootstrap_dwell_shift  # noqa: E402
from shape_metrics import (  # noqa: E402
    load_pdb, TOXIC_RMSD_MAX, TOXIC_JACCARD_MIN, OCCUPANCY_CUTOFF,
    BETA_CORE_RANGE, NAC_RANGE,
)

SHAPE = os.environ.get("SHAPE", "fusco_parallel_3mer_core70-88_relaxed")
LIGANDS = [x for x in os.environ.get("LIGANDS", "silibinin dhea caffeine").split() if x]
LO, HI = dt.CHUNK_RANGE
SD = dt.DWELL_DIR / SHAPE
REF = load_pdb(SD / f"apo_core{LO}-{HI}.pdb")
TRAJ_INTERVAL_PS = float(os.environ.get("TRAJ_INTERVAL_PS", "20"))

_cache: dict[str, list[dict]] = {}


def scores(label: str) -> list[dict]:
    if label not in _cache:
        out = []
        for p in sorted(glob.glob(str(SD / f"{label}_rep*.pdb"))):
            if p.endswith("_final.pdb"):
                continue  # single-frame endpoint, not a trajectory
            try:
                out.append(score_trajectory(Path(p), REF))
            except Exception as e:  # noqa: BLE001
                print(f"  ! {Path(p).name}: {type(e).__name__}: {e}")
        _cache[label] = out
    return _cache[label]


def _fin(xs):
    return [x for x in xs if np.isfinite(x)]


def summarise(label: str) -> list[float]:
    sc = scores(label)
    if not sc:
        print(f"[{label}] no trajectories found")
        return []
    dwell = [s["dwell_fraction"] for s in sc]
    occ = _fin([s["occupancy"] for s in sc])
    rmsd = _fin([s["mean_rmsd"] for s in sc])
    jac = _fin([s["mean_jaccard"] for s in sc])
    fd = _fin(dwell)
    print(f"[{label}] n={len(sc)} replicas")
    print(f"   dwell fraction: mean={np.mean(fd):.3f}  per-rep={[round(x, 2) for x in dwell]}")
    if occ:
        lo = min(occ)
        flag = "  <-- LOW: ligand left the core" if np.mean(occ) < 0.5 else "  (stays bound)"
        print(f"   occupancy (bound frac): mean={np.mean(occ):.2f} min={lo:.2f}{flag}")
    print(f"   mean beta-core RMSD={np.mean(rmsd):.2f} A   mean contact Jaccard={np.mean(jac):.2f}")
    return dwell


def mean_curve(sc: list[dict], key: str):
    arrs = [s["per_frame"][key] for s in sc if s["per_frame"][key]]
    if not arrs:
        return np.array([])
    n = min(len(a) for a in arrs)
    return np.array([a[:n] for a in arrs]).mean(axis=0)


def _slope_per_ns(values: list[float], dt_ns: float) -> float:
    """Least-squares slope of a per-frame series vs time (units/ns)."""
    v = np.asarray(values, dtype=float)
    if v.size < 3:
        return float("nan")
    t = np.arange(v.size) * dt_ns
    return float(np.polyfit(t, v, 1)[0])


def _verdict(lo: float, hi: float, higher_is_destab: bool) -> str:
    """Direction-aware label from a bootstrap CI of (complex - apo)."""
    if lo > 0:
        return "destabiliser" if higher_is_destab else "stabiliser"
    if hi < 0:
        return "stabiliser" if higher_is_destab else "destabiliser"
    return "inconclusive"


def rescore() -> None:
    """Ligand-blind rate observables on the existing trajectories — no basin
    threshold, no tuning. A destabiliser should make the beta-core RMSD rise
    FASTER (more positive slope) and the contact Jaccard decay FASTER (more
    negative slope) than apo. Reports the bootstrapped complex-minus-apo
    difference per axis with a direction-aware verdict."""
    dt_ns = TRAJ_INTERVAL_PS / 1000.0

    def rates(label: str, key: str) -> list[float]:
        return _fin([_slope_per_ns(s["per_frame"][key], dt_ns) for s in scores(label)])

    apo_r = rates("apo", "beta_core_rmsd")
    apo_j = rates("apo", "contact_jaccard")
    print("\n=== RESCORE: ligand-blind rate observables (no basin threshold) ===")
    print(f"   sampling dt={TRAJ_INTERVAL_PS:.0f} ps/frame; "
          f"apo RMSD-rise {np.mean(apo_r):+.3f} A/ns, apo Jaccard-change {np.mean(apo_j):+.4f}/ns")
    print("   destabiliser = RMSD rises faster (shift>0) OR Jaccard decays faster (shift<0)\n")
    print(f"   {'ligand':12} | {'dRMSD/ns shift [CI95]':28} {'verdict':12} | "
          f"{'dJacc/ns shift [CI95]':30} {'verdict':12}")
    for lig in LIGANDS:
        lr = rates(lig, "beta_core_rmsd")
        lj = rates(lig, "contact_jaccard")
        if not lr:
            print(f"   {lig:12} | (no trajectories)")
            continue
        br = bootstrap_dwell_shift(apo_r, lr)
        bj = bootstrap_dwell_shift(apo_j, lj)
        vr = _verdict(br["ci_low"], br["ci_high"], higher_is_destab=True)
        vj = _verdict(bj["ci_low"], bj["ci_high"], higher_is_destab=False)
        rcol = f"{br['shift']:+.3f} [{br['ci_low']:+.3f},{br['ci_high']:+.3f}]"
        jcol = f"{bj['shift']:+.4f} [{bj['ci_low']:+.4f},{bj['ci_high']:+.4f}]"
        print(f"   {lig:12} | {rcol:28} {vr:12} | {jcol:30} {vj:12}")
    print("\n   For the channel to be trustworthy: silibinin (control) should read "
          "'destabiliser' on at least one axis,\n   and caffeine (decoy) 'inconclusive'.")


def main() -> None:
    print(f"shape = {SHAPE}")
    print(f"thresholds: basin requires RMSD < {TOXIC_RMSD_MAX} A AND Jaccard > {TOXIC_JACCARD_MIN}; "
          f"ligand-bound distance cutoff {OCCUPANCY_CUTOFF} A")
    print(f"beta-core range {BETA_CORE_RANGE}, contact range {NAC_RANGE}\n")

    print("=== per-label summary ===")
    apo_dwell = summarise("apo")
    for lig in LIGANDS:
        d = summarise(lig)
        if d and apo_dwell:
            b = bootstrap_dwell_shift(apo_dwell, d)
            print(f"   shift vs apo = {b['shift']:+.3f}  "
                  f"CI[{b['ci_low']:+.3f},{b['ci_high']:+.3f}]  {b['classification']}")
        print()

    # Where do apo frames actually sit relative to the basin thresholds?
    apo_sc = scores("apo")
    allr = np.array([v for s in apo_sc for v in s["per_frame"]["beta_core_rmsd"]])
    allj = np.array([v for s in apo_sc for v in s["per_frame"]["contact_jaccard"]])
    if allr.size:
        print("=== apo per-frame distribution vs basin thresholds ===")
        print(f"   RMSD(A):  min={allr.min():.2f}  median={np.median(allr):.2f}  max={allr.max():.2f}   "
              f"(in-basin needs < {TOXIC_RMSD_MAX})")
        print(f"   Jaccard:  min={allj.min():.2f}  median={np.median(allj):.2f}  max={allj.max():.2f}   "
              f"(in-basin needs > {TOXIC_JACCARD_MIN})")
        in_basin = np.mean((allr < TOXIC_RMSD_MAX) & (allj > TOXIC_JACCARD_MIN))
        print(f"   apo frames in basin: {in_basin:.2f}\n")

    # Frame-resolved mean for apo vs the control (first ligand).
    ctrl = LIGANDS[0] if LIGANDS else None
    if ctrl:
        ar, aj = mean_curve(apo_sc, "beta_core_rmsd"), mean_curve(apo_sc, "contact_jaccard")
        cs = scores(ctrl)
        cr, cj = mean_curve(cs, "beta_core_rmsd"), mean_curve(cs, "contact_jaccard")
        n = min(len(ar), len(cr)) if len(cr) else len(ar)
        if n:
            print(f"=== frame-resolved mean: apo vs {ctrl} (does the control drift differently?) ===")
            print(f"   {'frame':>5}  {'apo_RMSD':>8}  {ctrl+'_RMSD':>14}  {'apo_Jacc':>8}  {ctrl+'_Jacc':>14}")
            step = max(1, n // 12)
            for i in range(0, n, step):
                print(f"   {i:5d}  {ar[i]:8.2f}  {cr[i]:14.2f}  {aj[i]:8.2f}  {cj[i]:14.2f}")


if __name__ == "__main__":
    main()
    if "--rescore" in sys.argv:
        rescore()
