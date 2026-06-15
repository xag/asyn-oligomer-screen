"""E8 node runner: a swarm of short, velocity-seeded MD replicas for MSM occupancy.

Manifest-driven (manifest.json in this dir). Unlike the E7 runner (a few long
single quenches), this fans MANY SHORT replicas from a library of distinct seed
conformations across this node's idle GPUs, to feed a Markov state model of the
oligomer's occupancy (issue #57). One process per replica, packed PER_GPU per
card via a slot scheduler, so #replicas may exceed #GPUs.

Each replica: velocity-seeded equilibrate -> N short segments (unbiased, no
restraints), frames every traj_interval_ps. Resumable: a finished replica drops
a DONE marker and is skipped on re-run; per-segment state files let an
interrupted replica continue from its last completed segment.

All temp / JIT / HOME writes are kept under the work tree. Runs under the conda
`md` python (CUDA OpenMM); shells out to md_relax.py for every MD step.
"""
import json
import os
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
M = json.loads((HERE / "manifest.json").read_text())
WORK = Path(M["work"])
PY = M["py"]
MD = M.get("md_relax", str(WORK / "md_relax.py"))
GPUS = [str(g) for g in M["gpus"]]
PER_GPU = int(M.get("per_gpu", 1))
EQUIL = float(M.get("equil_ps", 200.0))
SEG = float(M.get("segment_ps", 5000.0))
NSEG = int(M.get("n_segments", 3))
TI = float(M.get("traj_interval_ps", 50.0))
PREP_WORKERS = int(M.get("prep_workers", 4))
JOBS = M["jobs"]
BASE = str(WORK.parent)

TMP = WORK / "tmp"
NVC = WORK / "nvcache"
for d in ("prep", "traj", "logs"):
    (WORK / d).mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)
NVC.mkdir(parents=True, exist_ok=True)

_print_lock = threading.Lock()


def log(msg):
    with _print_lock:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def base_env(gpu=None):
    e = dict(os.environ)
    e.update(HOME=BASE, TMPDIR=str(TMP), CUDA_CACHE_PATH=str(NVC))
    if gpu is not None:
        e["CUDA_VISIBLE_DEVICES"] = str(gpu)
    return e


def sh(cmd, logpath, env):
    with open(logpath, "a") as f:
        f.write("\n$ " + " ".join(cmd) + "\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env).returncode


def prep_prefix(seed):
    return WORK / "prep" / seed


def prepare(seed, pdb):
    """Build + solvate one seed conformation once (CPU). Returns True on success."""
    prefix = prep_prefix(seed)
    sysxml = Path(str(prefix) + "_system.xml")
    solv = Path(str(prefix) + "_solvated.pdb")
    if sysxml.exists() and solv.exists():
        return True
    chunk = str(WORK / pdb) if not os.path.isabs(pdb) else pdb
    logpath = WORK / "logs" / f"prepare_{seed}.log"
    rc = sh([PY, MD, "--apo-pdb", chunk, "--rect-box", "--prepare-only", str(prefix)],
            logpath, base_env())
    ok = rc == 0 and sysxml.exists() and solv.exists()
    log(f"prepare {seed}: {'ok' if ok else 'FAIL rc=' + str(rc)}")
    return ok


def replica(job, gpu):
    """Velocity-seeded equilibrate + N short segments for one replica, pinned to gpu."""
    seed, rep = job["seed"], job["rep"]
    name = f"{seed}_rep{rep:02d}"
    prefix = prep_prefix(seed)
    sysxml = str(Path(str(prefix) + "_system.xml"))
    solv = str(Path(str(prefix) + "_solvated.pdb"))
    env = base_env(gpu)
    logpath = WORK / "logs" / f"{name}.log"
    repdir = WORK / "traj" / name
    repdir.mkdir(parents=True, exist_ok=True)
    if (repdir / "DONE").exists():
        return True

    s_e = repdir / "state_equil.xml"
    if not s_e.exists():
        rc = sh([PY, MD, "--equilibrate", str(s_e), "--system-xml", sysxml,
                 "--solvated-pdb", solv, "--equil-ps", str(EQUIL),
                 "--seed", str(1000 + rep)], logpath, env)
        if rc or not s_e.exists():
            log(f"[gpu{gpu}] FAIL equil {name} rc={rc}")
            return False

    prev = s_e
    for i in range(NSEG):
        so = repdir / f"state_{i}.xml"
        seg = repdir / f"seg_{i}.pdb"
        if seg.exists() and so.exists():
            prev = so
            continue
        rc = sh([PY, MD, "--segment", "--system-xml", sysxml, "--solvated-pdb", solv,
                 "--state-in", str(prev), "--state-out", str(so), "--seg-out", str(seg),
                 "--segment-ps", str(SEG), "--traj-interval-ps", str(TI),
                 "--checkpoint-s", "600", "--seed", str(10000 + 100 * rep + i)],
                logpath, env)
        if rc or not seg.exists():
            log(f"[gpu{gpu}] FAIL seg{i} {name} rc={rc}")
            return False
        prev = so
    (repdir / "DONE").write_text("ok\n")
    log(f"[gpu{gpu}] DONE {name}")
    return True


def run_packed(jobs):
    """Schedule replicas across GPUS, PER_GPU concurrent per card."""
    load = {g: 0 for g in GPUS}
    lock = threading.Lock()
    todo = list(jobs)
    running = []
    n_done = n_fail = 0

    def pick_gpu():
        free = [g for g in GPUS if load[g] < PER_GPU]
        return min(free, key=lambda g: load[g]) if free else None

    def worker(job, gpu):
        nonlocal n_done, n_fail
        ok = False
        try:
            ok = replica(job, gpu)
        finally:
            with lock:
                load[gpu] -= 1
                if ok:
                    n_done += 1
                else:
                    n_fail += 1

    while todo or running:
        with lock:
            while todo:
                g = pick_gpu()
                if g is None:
                    break
                job = todo.pop(0)
                load[g] += 1
                t = threading.Thread(target=worker, args=(job, g), daemon=True)
                t.start()
                running.append(t)
                log(f"launch {job['seed']} rep{job['rep']:02d} -> GPU{g} "
                    f"({len(todo)} queued)")
        running = [t for t in running if t.is_alive()]
        time.sleep(5)
    log(f"swarm finished: {n_done} done, {n_fail} failed")


def main():
    log(f"work={WORK} gpus={GPUS} per_gpu={PER_GPU} jobs={len(JOBS)} "
        f"equil={EQUIL}ps seg={SEG}ps x{NSEG} ti={TI}ps")
    seeds = {}
    for j in JOBS:
        seeds.setdefault(j["seed"], j["pdb"])
    log(f"preparing {len(seeds)} unique seed conformations ({PREP_WORKERS}-wide)")
    with ThreadPoolExecutor(max_workers=PREP_WORKERS) as ex:
        results = dict(zip(seeds, ex.map(lambda kv: prepare(*kv), seeds.items())))
    ready = {s for s, ok in results.items() if ok}
    jobs = [j for j in JOBS if j["seed"] in ready]
    dropped = len(JOBS) - len(jobs)
    if dropped:
        log(f"NOTE: dropping {dropped} jobs whose seed prep failed: "
            f"{sorted(set(results) - ready)}")
    log(f"{len(jobs)} replica jobs across {len(GPUS)} GPUs "
        f"(cap {len(GPUS) * PER_GPU} concurrent)")
    run_packed(jobs)
    log("ALL DONE")


if __name__ == "__main__":
    main()
