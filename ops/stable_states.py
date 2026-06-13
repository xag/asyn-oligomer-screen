"""Stable-states probe (E5): which conformations does the toxic oligomer
actually occupy, and does any stable basin remain toxic-looking?

Every screen channel scores against one hand-built oligomer shape that E2
showed is not metastable. This experiment releases that shape under *unbiased*
MD (no beta-core restraints), runs until the beta-core RMSD plateaus, clusters
the sampled shapes, and scores each basin's Stage-2 activity.

Two subcommands:

  run      Unbiased seeded MD from a relaxed oligomer. Reuses md_relax.py's
           resumable segment path (which applies NO position restraints, so the
           dynamics are unbiased) on the NAC-core chunk. Per shape: build the
           apo system once, then per velocity seed equilibrate + N x segment.
           MD runs in the conda MD env (located by screen/md_env.py).

  analyze  Score every frame with shape_metrics (beta-core Ca RMSD + inter-chain
           contact Jaccard vs the shape's own relaxed reference), report the
           per-replica plateau (last-20 ns RMSD slope), pool post-equilibration
           frames, cluster in standardized (RMSD, Jaccard) space (Ward; k by
           silhouette), and score each basin medoid with score_oligomer.py.
           Baselines: each relaxed reference chunk + the coil control chunk.

Usage:
    python ops/stable_states.py run                 # both default shapes, 3 seeds x 50 ns
    python ops/stable_states.py analyze             # cluster + score what run produced
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))

from dwell_time import _truncate_chunk, CHUNK_RANGE          # noqa: E402
from shape_metrics import (                                   # noqa: E402
    beta_core_rmsd, contact_jaccard, in_toxic_basin, load_pdb,
    TOXIC_RMSD_MAX, TOXIC_JACCARD_MIN,
)

OLIGO = ROOT / "results" / "oligomers"
OUT = ROOT / "results" / "stable_states"
SCORE = ROOT / "oligomers" / "score_oligomer.py"
VENV_PY = sys.executable

# Reference shape (what the screen docks against) + the structurally distinct
# alternative register. Coil trimer is the non-toxic floor for the activity scale.
SHAPES = ["fusco_parallel_3mer_core70-88", "fusco_antiparallel_3mer_core70-88"]
COIL_CTRL = "ctrl_coil_3mer"

EQUIL_PS = 200.0
SEGMENT_PS = 10000.0
N_SEGMENTS = 5            # 50 ns production / replica
TRAJ_INTERVAL_PS = 50.0
N_SEEDS = 3
FRAME_PS = TRAJ_INTERVAL_PS


# ---------------------------------------------------------------------------
# run: unbiased seeded MD
# ---------------------------------------------------------------------------

def _md_python() -> str:
    from md_env import md_python
    return str(md_python())


def _chunk(shape: str) -> Path:
    lo, hi = CHUNK_RANGE
    out = OUT / "inputs" / f"{shape}_core{lo}-{hi}.pdb"
    if not out.exists():
        _truncate_chunk(OLIGO / f"{shape}_relaxed.pdb", out, lo, hi)
    return out


def run(shapes, n_seeds=N_SEEDS, n_segments=N_SEGMENTS):
    md = _md_python()
    md_relax = str(ROOT / "screen" / "md_relax.py")
    for shape in shapes:
        chunk = _chunk(shape)
        prefix = OUT / "prep" / shape
        prefix.parent.mkdir(parents=True, exist_ok=True)
        sysxml = Path(str(prefix) + "_system.xml")
        if not sysxml.exists():
            subprocess.run([md, md_relax, "--apo-pdb", str(chunk), "--rect-box",
                            "--prepare-only", str(prefix)], check=True)
        solv = str(Path(str(prefix) + "_solvated.pdb"))
        for seed in range(n_seeds):
            rep = OUT / f"{shape}_rep{seed}"
            rep.mkdir(parents=True, exist_ok=True)
            s_e = rep / "state_equil.xml"
            if not s_e.exists():
                subprocess.run([md, md_relax, "--equilibrate", str(s_e),
                                "--system-xml", str(sysxml), "--solvated-pdb", solv,
                                "--equil-ps", str(EQUIL_PS), "--seed", str(1000 + seed)], check=True)
            prev = s_e
            for i in range(n_segments):
                state_out, seg_out = rep / f"state_{i}.xml", rep / f"seg_{i}.pdb"
                if seg_out.exists() and state_out.exists():
                    prev = state_out
                    continue
                subprocess.run([md, md_relax, "--segment",
                                "--system-xml", str(sysxml), "--solvated-pdb", solv,
                                "--state-in", str(prev), "--state-out", str(state_out),
                                "--seg-out", str(seg_out), "--segment-ps", str(SEGMENT_PS),
                                "--traj-interval-ps", str(TRAJ_INTERVAL_PS),
                                "--seed", str(10000 + 100 * seed + i)], check=True)
                prev = state_out
            (rep / "DONE").write_text("ok\n")


# ---------------------------------------------------------------------------
# analyze: plateau -> cluster -> score
# ---------------------------------------------------------------------------

def _seg_files(rep_dir: Path):
    return sorted(rep_dir.glob("seg_*.pdb"),
                  key=lambda p: int(re.search(r"seg_(\d+)", p.name).group(1)))


def _frame_series(rep_dir: Path, reference):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    rows, g = [], 0
    for seg in _seg_files(rep_dir):
        for li, model in enumerate(parser.get_structure(seg.stem, str(seg))):
            rows.append((str(seg), li, (g + 1) * FRAME_PS,
                         beta_core_rmsd(model, reference), contact_jaccard(model, reference)))
            g += 1
    return rows


def _slope_per_ns(t_ps, y):
    if len(t_ps) < 3:
        return float("nan")
    return float(np.polyfit(np.asarray(t_ps) / 1000.0, np.asarray(y), 1)[0])


def _silhouette(X, labels):
    if len(set(labels)) < 2:
        return -1.0
    D = squareform(pdist(X))
    uniq, sl = np.unique(labels), []
    for i in range(len(X)):
        same = labels == labels[i]; same[i] = False
        a = D[i, same].mean() if same.any() else 0.0
        b = min((D[i, labels == u].mean() for u in uniq if u != labels[i]), default=np.inf)
        sl.append(0.0 if max(a, b) == 0 else (b - a) / max(a, b))
    return float(np.mean(sl))


def _cluster(pts):
    if len(pts) < 10:
        return np.ones(len(pts), dtype=int), 1
    Z = linkage(pts, method="ward")
    idx = np.arange(len(pts))
    if len(pts) > 1500:
        idx = np.random.default_rng(0).choice(len(pts), 1500, replace=False)
    best = (1, -1.0, np.ones(len(pts), dtype=int))
    for k in (2, 3):
        labels = fcluster(Z, k, criterion="maxclust")
        s = _silhouette(pts[idx], labels[idx])
        if s > best[1]:
            best = (k, s, labels)
    return (best[2], best[0]) if best[1] >= 0.5 else (np.ones(len(pts), dtype=int), 1)


def _write_medoid(seg_path, local_idx, out_pdb):
    from Bio.PDB import PDBParser, PDBIO
    model = list(PDBParser(QUIET=True).get_structure("m", seg_path))[local_idx]
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO(); io.set_structure(model); io.save(str(out_pdb))


def _activity(pdb: Path):
    r = subprocess.run([VENV_PY, str(SCORE), str(pdb)], capture_output=True, text=True)
    m = re.search(r"activity score:\s*([+-]?\d+\.\d+)", r.stdout)
    return float(m.group(1)) if m else None


def analyze(shapes):
    ana = OUT / "analysis"; ana.mkdir(parents=True, exist_ok=True)
    summary = {"frame_ps": FRAME_PS, "toxic_rmsd_max": TOXIC_RMSD_MAX,
               "toxic_jaccard_min": TOXIC_JACCARD_MIN, "baselines": {}, "shapes": {}}

    print("== baselines ==")
    for shape in shapes:
        a = _activity(_chunk(shape))
        summary["baselines"][shape + "_ref"] = a
        print(f"  {shape}_ref activity={a}")
    a = _activity(_chunk(COIL_CTRL)); summary["baselines"]["coil_ctrl"] = a
    print(f"  coil_ctrl activity={a}")

    for shape in shapes:
        reps = [d for d in sorted(OUT.glob(f"{shape}_rep*")) if _seg_files(d)]
        if not reps:
            print(f"[{shape}] no trajectories; run first"); continue
        reference = load_pdb(_chunk(shape))
        print(f"\n== {shape}: {len(reps)} replicas ==")
        per_rep, pool, plateau = {}, [], []
        for d in reps:
            rows = _frame_series(d, reference); per_rep[d.name] = rows
            t = [r[2] for r in rows]; rr = [r[3] for r in rows]; tmax = t[-1]
            tail = [(tt, y) for tt, y in zip(t, rr) if tt >= tmax - 20000]
            plateau.append({"rep": d.name, "t_ns": tmax / 1000.0,
                            "slope_last20ns": _slope_per_ns([x for x, _ in tail], [y for _, y in tail]),
                            "mean_rmsd_last20ns": float(np.mean([y for _, y in tail]))})
            print(f"  {d.name}: {tmax/1000:.0f} ns, slope_last20ns={plateau[-1]['slope_last20ns']:+.4f} A/ns, "
                  f"mean_rmsd_tail={plateau[-1]['mean_rmsd_last20ns']:.2f} A")
            with open(ana / f"{d.name}_timeseries.csv", "w") as f:
                f.write("t_ps,beta_core_rmsd,contact_jaccard\n")
                for r in rows:
                    f.write(f"{r[2]:.0f},{r[3]:.4f},{r[4]:.4f}\n")
            half = rows[-1][2] / 2.0
            pool += [(d.name,) + r for r in rows if r[2] >= half]

        rmsd = np.array([p[4] for p in pool]); jacc = np.array([p[5] for p in pool])
        X = np.column_stack([rmsd, jacc]); Xs = (X - X.mean(0)) / (X.std(0) + 1e-9)
        labels, k = _cluster(Xs)
        print(f"  post-equil pool: {len(pool)} frames -> {k} basin(s)")
        basins = []
        for cl in sorted(set(labels)):
            m = labels == cl; cen = Xs[m].mean(0)
            medoid = np.where(m)[0][int(np.argmin(((Xs[m] - cen) ** 2).sum(1)))]
            rep, seg, li, t, rr, jj = pool[medoid]
            med_pdb = ana / f"{shape}_basin{cl}_medoid.pdb"; _write_medoid(seg, li, med_pdb)
            frac = float(np.mean([in_toxic_basin({"beta_core_rmsd": pool[i][4], "contact_jaccard": pool[i][5]})
                                  for i in np.where(m)[0]]))
            act = _activity(med_pdb)
            basins.append({"basin": int(cl), "pop_frac": float(m.mean()),
                           "mean_rmsd": float(rmsd[m].mean()), "mean_jaccard": float(jacc[m].mean()),
                           "frac_in_toxic_basin": frac, "activity": act,
                           "medoid": {"rep": rep, "seg": Path(seg).name, "model": li, "t_ns": t / 1000.0}})
            print(f"    basin {cl}: pop={m.mean()*100:.0f}% mean_rmsd={rmsd[m].mean():.2f}A "
                  f"mean_jacc={jacc[m].mean():.2f} toxic_frac={frac:.2f} activity={act}")
        summary["shapes"][shape] = {"plateau": plateau, "n_basins": int(k), "basins": basins}

    (ana / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nwrote {ana / 'summary.json'}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("run"); pr.add_argument("--shapes", nargs="+", default=SHAPES)
    pr.add_argument("--seeds", type=int, default=N_SEEDS)
    pr.add_argument("--segments", type=int, default=N_SEGMENTS)
    pa = sub.add_parser("analyze"); pa.add_argument("--shapes", nargs="+", default=SHAPES)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    args = p.parse_args()
    if args.cmd == "run":
        run(args.shapes, args.seeds, args.segments)
    else:
        analyze(args.shapes)


if __name__ == "__main__":
    main()
