#!/usr/bin/env python3
"""Packed, GPU-pinned dwell-time replica launcher for a shared multi-GPU node.

The stock `dwell_time.py pilot` runs replicas serially and unpinned — on a box
where other users hold GPUs 0-3, that lands on a busy card. This controller
instead:
  * prepares each (apo / docked-complex) system ONCE via md_relax --prepare-only
    (CPU, conda env) -> system.xml + solvated.pdb,
  * fans the velocity-seeded replicas across a pool of idle GPUs, pinning each
    md_relax process to one card with CUDA_VISIBLE_DEVICES and packing PER_GPU
    per card,
  * scores the trajectories with the dwell-time channel and prints the gate.

Runs under the PIP VENV python (needs numpy/biopython, and rdkit/meeko/Vina for
docking). It spawns md_relax via the CONDA env python ($ASYN_MD_PYTHON), which
carries CUDA OpenMM + OpenFF.

Config via env (all optional except ASYN_REPO + ASYN_MD_PYTHON):
  ASYN_REPO, ASYN_MD_PYTHON, SHAPE, GPUS="4 5 6 7", REPLICAS=10, PER_GPU=2,
  PROD_NS=2.0, EQUIL_PS=100, TRAJ_INTERVAL_PS=20, LIGANDS="" (space-sep; empty
  => apo only, no Vina/OpenFF needed).
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(os.environ["ASYN_REPO"]).resolve()
sys.path.insert(0, str(REPO / "screen"))
os.chdir(REPO)

MDPY = os.environ["ASYN_MD_PYTHON"]
MD_RELAX = str(REPO / "screen" / "md_relax.py")
SHAPE = os.environ.get("SHAPE", "fusco_parallel_3mer_core70-88_relaxed")
GPUS = os.environ.get("GPUS", "4 5 6 7").split()
REPLICAS = int(os.environ.get("REPLICAS", "10"))
PER_GPU = int(os.environ.get("PER_GPU", "2"))
PROD_NS = float(os.environ.get("PROD_NS", "2.0"))
EQUIL_PS = float(os.environ.get("EQUIL_PS", "100"))
TRAJ_INT = float(os.environ.get("TRAJ_INTERVAL_PS", "20"))
LIGANDS = [x for x in os.environ.get("LIGANDS", "").split() if x]

import dwell_time as dt  # noqa: E402
from dwell_time import _truncate_chunk, score_trajectory, bootstrap_dwell_shift  # noqa: E402
from shape_metrics import load_pdb  # noqa: E402

LO, HI = dt.CHUNK_RANGE
SHAPE_DIR = dt.DWELL_DIR / SHAPE
SHAPE_PDB = REPO / "results" / "oligomers" / f"{SHAPE}.pdb"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def prepare(prefix: Path, build_args: list[str]) -> tuple[str, str]:
    """Run md_relax --prepare-only once; return (system_xml, solvated_pdb)."""
    sysxml = prefix.with_name(prefix.name + "_system.xml")
    solv = prefix.with_name(prefix.name + "_solvated.pdb")
    if sysxml.exists() and solv.exists():
        log(f"prep cache hit: {prefix.name}")
        return str(sysxml), str(solv)
    cmd = [MDPY, MD_RELAX, *build_args, "--rect-box", "--prepare-only", str(prefix)]
    log("prepare: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    if not (sysxml.exists() and solv.exists()):
        raise RuntimeError(f"prepare produced no system for {prefix.name}")
    return str(sysxml), str(solv)


def build_jobs() -> list[dict]:
    """Prepare apo (+ each docked complex) and expand into per-replica jobs."""
    SHAPE_DIR.mkdir(parents=True, exist_ok=True)
    jobs: list[dict] = []

    # --- apo ---
    apo_chunk = _truncate_chunk(SHAPE_PDB, SHAPE_DIR / f"apo_core{LO}-{HI}.pdb", LO, HI)
    apo_sys, apo_solv = prepare(SHAPE_DIR / f"apo_core{LO}-{HI}", ["--apo-pdb", str(apo_chunk)])
    for i in range(REPLICAS):
        jobs.append(_job("apo", apo_sys, apo_solv, 1000 + i, i))

    # --- complexes (need Vina + OpenFF) ---
    if LIGANDS:
        import stage3  # noqa: F401  (pulls rdkit/meeko/Vina; only needed when docking)
        from dwell_time import _ensure_complex
        skipped: list[str] = []
        for lig in LIGANDS:
            # One ligand's dock/parametrise failure must not sink the whole gate —
            # log it and carry on with the rest.
            try:
                complex_pdb, smiles = _ensure_complex(lig, SHAPE)
                complex_chunk = _truncate_chunk(
                    complex_pdb, SHAPE_DIR / f"{lig}_core{LO}-{HI}_complex.pdb", LO, HI)
                sysxml, solv = prepare(
                    SHAPE_DIR / f"{lig}_core{LO}-{HI}",
                    ["--complex-pdb", str(complex_chunk), "--ligand-smiles", smiles])
            except Exception as e:
                log(f"SKIP ligand {lig}: prep failed ({type(e).__name__}: {e})")
                skipped.append(lig)
                continue
            for i in range(REPLICAS):
                jobs.append(_job(lig, sysxml, solv, 2000 + i, i))
        if skipped:
            log(f"NOTE: {len(skipped)} ligand(s) skipped at prep: {skipped}")

    return jobs


def _job(label: str, sysxml: str, solv: str, seed: int, i: int) -> dict:
    return {
        "label": label, "rep": i, "sysxml": sysxml, "solv": solv, "seed": seed,
        "traj": str(SHAPE_DIR / f"{label}_rep{i:02d}.pdb"),
        "final": str(SHAPE_DIR / f"{label}_rep{i:02d}_final.pdb"),
        "log": str(SHAPE_DIR / f"{label}_rep{i:02d}.mdlog"),
    }


def run_packed(jobs: list[dict]) -> None:
    """Launch jobs across the GPU pool, PER_GPU concurrent per card."""
    load = {g: 0 for g in GPUS}
    running: list[tuple] = []
    todo = list(jobs)
    n_done = n_skip = 0

    def pick_gpu():
        free = [g for g in GPUS if load[g] < PER_GPU]
        return min(free, key=lambda g: load[g]) if free else None

    while todo or running:
        while todo:
            if Path(todo[0]["final"]).exists():
                j = todo.pop(0); n_skip += 1
                log(f"skip cached {j['label']} rep{j['rep']:02d}")
                continue
            g = pick_gpu()
            if g is None:
                break
            j = todo.pop(0)
            env = dict(os.environ); env["CUDA_VISIBLE_DEVICES"] = g
            cmd = [MDPY, MD_RELAX, "--system-xml", j["sysxml"], "--solvated-pdb", j["solv"],
                   "--out-pdb", j["final"], "--seed", str(j["seed"]),
                   "--equil-ps", str(EQUIL_PS), "--prod-ps", str(PROD_NS * 1000.0),
                   "--traj-out", j["traj"], "--traj-interval-ps", str(TRAJ_INT)]
            fh = open(j["log"], "w")
            p = subprocess.Popen(cmd, env=env, stdout=fh, stderr=subprocess.STDOUT)
            running.append((p, g, j, fh)); load[g] += 1
            log(f"launch {j['label']} rep{j['rep']:02d} -> GPU{g} (pid {p.pid}); "
                f"{len(running)} running, {len(todo)} queued")
        time.sleep(5)
        for tup in running[:]:
            p, g, j, fh = tup
            if p.poll() is not None:
                fh.close(); load[g] -= 1; running.remove(tup); n_done += 1
                ok = (p.returncode == 0 and Path(j["final"]).exists())
                log(f"{'done' if ok else 'FAILED'} {j['label']} rep{j['rep']:02d} "
                    f"rc={p.returncode} ({n_done} done, {len(todo)} queued)")
    log(f"all replicas finished: {n_done} ran, {n_skip} cached")


def score(jobs: list[dict]) -> None:
    ref = load_pdb(SHAPE_DIR / f"apo_core{LO}-{HI}.pdb")
    by_label: dict[str, list[str]] = {}
    for j in jobs:
        by_label.setdefault(j["label"], []).append(j["traj"])

    def dwell(trajs):
        out = []
        for t in trajs:
            if Path(t).exists():
                out.append(score_trajectory(Path(t), ref)["dwell_fraction"])
        return out

    apo = dwell(by_label.get("apo", []))
    log(f"apo dwell fractions: {[round(x, 3) for x in apo]}")
    if not LIGANDS:
        import numpy as np
        finite = [x for x in apo if x == x]
        log(f"apo mean dwell (toxic-basin occupancy of the bare oligomer): "
            f"{np.mean(finite):.3f}" if finite else "apo mean dwell: n/a")
        log("apo-only run: no ligand gate to score. Re-run with LIGANDS set for the gate.")
        return

    print("\n=== GATE: dwell shift vs apo (shift<0 destabiliser, >0 stabiliser) ===", flush=True)
    print(f"{'ligand':14} {'shift':>8}  {'CI95':>20}  class", flush=True)
    for lig in LIGANDS:
        d = dwell(by_label.get(lig, []))
        if not d:
            print(f"{lig:14} {'--':>8}  {'(skipped at prep)':>20}", flush=True)
            continue
        b = bootstrap_dwell_shift(apo, d)
        print(f"{lig:14} {b['shift']:+8.3f}  "
              f"[{b['ci_low']:+.3f},{b['ci_high']:+.3f}]  {b['classification']}", flush=True)


def main() -> None:
    log(f"repo={REPO} shape={SHAPE} gpus={GPUS} replicas={REPLICAS} "
        f"per_gpu={PER_GPU} prod_ns={PROD_NS} ligands={LIGANDS or '(apo only)'}")
    jobs = build_jobs()
    log(f"{len(jobs)} replica jobs across {len(GPUS)} GPUs (cap {len(GPUS) * PER_GPU} concurrent)")
    run_packed(jobs)
    score(jobs)
    log("ALL DONE")


if __name__ == "__main__":
    main()
