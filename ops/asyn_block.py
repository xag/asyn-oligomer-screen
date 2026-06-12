#!/usr/bin/env python3
"""Blocked, matched-seed dwell-time runner for the low-cost validation.

One invocation = one BLOCK = one starting conformer. Within the block, apo /
test / decoy are run from MATCHED velocity seeds, so the per-block contrast
cancels the (large) structural + thermal variance instead of letting it drown
the signal. Blocks accumulate across short opportunistic sessions: each night
pick a conformer (and the launcher skips any replica already on disk), so the
store under results/blocks/<conformer>/ grows until asyn_pool.py calls it.

Runs under the PIP VENV python; spawns md_relax via the CONDA env ($ASYN_MD_PYTHON,
CUDA OpenMM + OpenFF). GPU-pinned + packed exactly like the pilot launcher.

Env: ASYN_REPO, ASYN_MD_PYTHON (required); CONFORMER (shape stem), LIGANDS
("silibinin caffeine"), REPLICAS (5), SEED_BASE (7000), GPUS ("4 5 6 7"),
PER_GPU (2), PROD_NS (2.0), TRAJ_INTERVAL_PS (20). Apo is shared per conformer;
test/decoy share apo's seeds for the matched contrast.
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
CONFORMER = os.environ.get("CONFORMER", "fusco_parallel_3mer_core70-88_relaxed")
LIGANDS = [x for x in os.environ.get("LIGANDS", "silibinin caffeine").split() if x]
REPLICAS = int(os.environ.get("REPLICAS", "5"))
SEED_BASE = int(os.environ.get("SEED_BASE", "7000"))
GPUS = os.environ.get("GPUS", "4 5 6 7").split()
PER_GPU = int(os.environ.get("PER_GPU", "2"))
PROD_NS = float(os.environ.get("PROD_NS", "2.0"))
EQUIL_PS = float(os.environ.get("EQUIL_PS", "100"))
TRAJ_INT = float(os.environ.get("TRAJ_INTERVAL_PS", "20"))

import dwell_time as dt  # noqa: E402
from dwell_time import _truncate_chunk  # noqa: E402

LO, HI = dt.CHUNK_RANGE
BLOCK_DIR = REPO / "results" / "blocks" / CONFORMER
SHAPE_PDB = REPO / "results" / "oligomers" / f"{CONFORMER}.pdb"
SEEDS = [SEED_BASE + i for i in range(REPLICAS)]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def prepare(prefix: Path, build_args: list[str]) -> tuple[str, str]:
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


def _job(label: str, sysxml: str, solv: str, seed: int) -> dict:
    return {
        "label": label, "seed": seed, "sysxml": sysxml, "solv": solv,
        "traj": str(BLOCK_DIR / f"{label}_s{seed}.pdb"),
        "final": str(BLOCK_DIR / f"{label}_s{seed}_final.pdb"),
        "log": str(BLOCK_DIR / f"{label}_s{seed}.mdlog"),
    }


def build_jobs() -> list[dict]:
    if not SHAPE_PDB.exists():
        raise FileNotFoundError(f"conformer PDB not found: {SHAPE_PDB}")
    BLOCK_DIR.mkdir(parents=True, exist_ok=True)
    jobs: list[dict] = []

    apo_chunk = _truncate_chunk(SHAPE_PDB, BLOCK_DIR / f"apo_core{LO}-{HI}.pdb", LO, HI)
    apo_sys, apo_solv = prepare(BLOCK_DIR / "apo", ["--apo-pdb", str(apo_chunk)])
    for s in SEEDS:
        jobs.append(_job("apo", apo_sys, apo_solv, s))

    import stage3  # noqa: F401
    from dwell_time import _ensure_complex
    skipped: list[str] = []
    for lig in LIGANDS:
        try:
            complex_pdb, smiles = _ensure_complex(lig, CONFORMER)
            chunk = _truncate_chunk(complex_pdb, BLOCK_DIR / f"{lig}_complex_core{LO}-{HI}.pdb", LO, HI)
            sysxml, solv = prepare(BLOCK_DIR / lig, ["--complex-pdb", str(chunk), "--ligand-smiles", smiles])
        except Exception as e:  # noqa: BLE001
            log(f"SKIP ligand {lig}: prep failed ({type(e).__name__}: {e})")
            skipped.append(lig)
            continue
        for s in SEEDS:
            jobs.append(_job(lig, sysxml, solv, s))
    if skipped:
        log(f"NOTE: {len(skipped)} ligand(s) skipped at prep: {skipped}")
    return jobs


def run_packed(jobs: list[dict]) -> None:
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
                log(f"skip cached {j['label']} s{j['seed']}")
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
            log(f"launch {j['label']} s{j['seed']} -> GPU{g} (pid {p.pid}); "
                f"{len(running)} running, {len(todo)} queued")
        time.sleep(5)
        for tup in running[:]:
            p, g, j, fh = tup
            if p.poll() is not None:
                fh.close(); load[g] -= 1; running.remove(tup); n_done += 1
                ok = (p.returncode == 0 and Path(j["final"]).exists())
                log(f"{'done' if ok else 'FAILED'} {j['label']} s{j['seed']} "
                    f"rc={p.returncode} ({n_done} done, {len(todo)} queued)")
    log(f"block finished: {n_done} ran, {n_skip} cached")


def main() -> None:
    log(f"BLOCK conformer={CONFORMER} ligands={LIGANDS} seeds={SEEDS} "
        f"gpus={GPUS} per_gpu={PER_GPU} prod_ns={PROD_NS}")
    jobs = build_jobs()
    log(f"{len(jobs)} matched replicas (cap {len(GPUS) * PER_GPU} concurrent); store={BLOCK_DIR}")
    run_packed(jobs)
    log("BLOCK DONE")


if __name__ == "__main__":
    main()
