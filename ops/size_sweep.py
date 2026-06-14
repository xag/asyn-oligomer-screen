"""Size sweep (E7): does the register-split, history-dependent landscape (E6)
hold, sharpen, or break with oligomer order?

  run      Unbiased 50 ns MD (no restraints, md_relax.py segment path) from
           2-mer and 4-mer cores, both registers, N seeds each. MD env via md_env.
  analyze  Per (size,register) cell: plateau slope, drift from own start,
           end-basin activity (score_oligomer; per-chain mean => size-comparable),
           and within-cell basin count. Cross-size comparison is via activity,
           since whole-assembly RMSD is not comparable across chain counts.

Usage:
    python ops/size_sweep.py run [--seeds 3]
    python ops/size_sweep.py analyze
"""
from __future__ import annotations
import argparse, json, re, subprocess, sys
from pathlib import Path
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))
from dwell_time import _truncate_chunk, CHUNK_RANGE
from shape_metrics import beta_core_rmsd, contact_jaccard, load_pdb

OLIGO = ROOT / "results" / "oligomers"
OUT = ROOT / "results" / "size_sweep_states"
SCORE = ROOT / "oligomers" / "score_oligomer.py"
VENV_PY = sys.executable
FRAME_PS = 50.0
EQUIL_PS, SEGMENT_PS, N_SEGMENTS, TRAJ_PS = 200.0, 10000.0, 5, 50.0

CELLS = {  # (size, register): shape stem ; seed count per cell for run()
 ("2mer", "parallel"): "fusco_parallel_2mer_core70-88",
 ("2mer", "antiparallel"): "fusco_antiparallel_2mer_core70-88",
 ("4mer", "parallel"): "fusco_parallel_4mer_core70-88",
 ("4mer", "antiparallel"): "fusco_antiparallel_4mer_core70-88",
}


def _chunk(shape):
    lo, hi = CHUNK_RANGE
    out = OUT / "inputs" / f"{shape}_core{lo}-{hi}.pdb"
    if not out.exists():
        _truncate_chunk(OLIGO / f"{shape}_relaxed.pdb", out, lo, hi)
    return out


def run(seeds=3):
    from md_env import md_python
    md = str(md_python()); mdrelax = str(ROOT / "screen" / "md_relax.py")
    for shape in CELLS.values():
        chunk = str(_chunk(shape)); prefix = OUT / "prep" / shape
        prefix.parent.mkdir(parents=True, exist_ok=True)
        sysxml = Path(str(prefix) + "_system.xml")
        if not sysxml.exists():
            subprocess.run([md, mdrelax, "--apo-pdb", chunk, "--rect-box", "--prepare-only", str(prefix)], check=True)
        solv = str(Path(str(prefix) + "_solvated.pdb"))
        for seed in range(seeds):
            rep = OUT / f"{shape}_rep{seed}"; rep.mkdir(parents=True, exist_ok=True)
            s_e = rep / "state_equil.xml"
            if not s_e.exists():
                subprocess.run([md, mdrelax, "--equilibrate", str(s_e), "--system-xml", str(sysxml),
                                "--solvated-pdb", solv, "--equil-ps", str(EQUIL_PS), "--seed", str(1000 + seed)], check=True)
            prev = s_e
            for i in range(N_SEGMENTS):
                so, seg = rep / f"state_{i}.xml", rep / f"seg_{i}.pdb"
                if seg.exists() and so.exists(): prev = so; continue
                subprocess.run([md, mdrelax, "--segment", "--system-xml", str(sysxml), "--solvated-pdb", solv,
                                "--state-in", str(prev), "--state-out", str(so), "--seg-out", str(seg),
                                "--segment-ps", str(SEGMENT_PS), "--traj-interval-ps", str(TRAJ_PS),
                                "--seed", str(10000 + 100 * seed + i)], check=True)
                prev = so
            (rep / "DONE").write_text("ok\n")


def _segs(d): return sorted(d.glob("seg_*.pdb"), key=lambda p: int(re.search(r"seg_(\d+)", p.name).group(1)))
def _activity(pdb):
    r = subprocess.run([VENV_PY, str(SCORE), str(pdb)], capture_output=True, text=True)
    m = re.search(r"activity score:\s*([+-]?\d+\.\d+)", r.stdout); return float(m.group(1)) if m else None
def _silh(X, lab):
    if len(set(lab)) < 2: return -1.0
    D = squareform(pdist(X)); s = []
    for i in range(len(X)):
        same = lab == lab[i]; same[i] = False
        a = D[i, same].mean() if same.any() else 0.0
        b = min((D[i, lab == u].mean() for u in np.unique(lab) if u != lab[i]), default=np.inf)
        s.append(0.0 if max(a, b) == 0 else (b - a) / max(a, b))
    return float(np.mean(s))
def _nclust(pts):
    if len(pts) < 10: return 1
    Z = linkage(pts, "ward"); idx = np.arange(len(pts)); best = (1, -1.0)
    if len(pts) > 1500: idx = np.random.default_rng(0).choice(len(pts), 1500, replace=False)
    for k in (2, 3):
        lab = fcluster(Z, k, criterion="maxclust"); s = _silh(pts[idx], lab[idx])
        if s > best[1]: best = (k, s)
    return best[0] if best[1] >= 0.5 else 1


def analyze():
    from Bio.PDB import PDBParser, PDBIO
    ana = OUT / "analysis"; ana.mkdir(parents=True, exist_ok=True)
    parser = PDBParser(QUIET=True); summary = {"cells": []}
    print(f"{'cell':22}{'seeds':6}{'start':7}{'end(mean[min,max])':22}{'drift':8}{'slope':9}basins")
    for (size, reg), shape in CELLS.items():
        ref = load_pdb(_chunk(shape)); start = _activity(_chunk(shape))
        reps = [d for d in sorted(OUT.glob(f"{shape}_rep*")) if _segs(d)]
        ends, drifts, slopes, pool = [], [], [], []
        for d in reps:
            rows, g = [], 0
            for s in _segs(d):
                for li, m in enumerate(parser.get_structure(s.stem, str(s))):
                    rows.append((str(s), li, (g + 1) * FRAME_PS, beta_core_rmsd(m, ref), contact_jaccard(m, ref))); g += 1
            t = [r[2] for r in rows]; rr = [r[3] for r in rows]; tmax = t[-1]
            tail = [(x, y) for x, y in zip(t, rr) if x >= tmax - 20000]
            slopes.append(float(np.polyfit(np.array([x for x, _ in tail]) / 1000., np.array([y for _, y in tail]), 1)[0]))
            drifts.append(float(np.mean([y for _, y in tail])))
            post = [r for r in rows if r[2] >= tmax / 2]
            X = np.array([[r[3], r[4]] for r in post]); med = post[int(np.argmin(((X - X.mean(0)) ** 2).sum(1)))]
            mp = ana / f"{d.name}_medoid.pdb"; mo = list(parser.get_structure("m", med[0]))[med[1]]
            io = PDBIO(); io.set_structure(mo); io.save(str(mp)); ends.append(_activity(mp))
            pool += [(r[3], r[4]) for r in post]
            with open(ana / f"{d.name}_timeseries.csv", "w") as f:
                f.write("t_ps,beta_core_rmsd_vs_own_start,contact_jaccard_vs_own_start\n")
                for r in rows: f.write(f"{r[2]:.0f},{r[3]:.4f},{r[4]:.4f}\n")
        P = np.array(pool); nb = _nclust((P - P.mean(0)) / (P.std(0) + 1e-9)); ea = np.array([a for a in ends if a is not None])
        cell = {"size": size, "register": reg, "shape": shape, "n_seeds": len(reps), "start_activity": start,
                "end_activity_mean": float(ea.mean()), "end_activity_min": float(ea.min()), "end_activity_max": float(ea.max()),
                "end_activity_per_seed": ends, "drift_mean": float(np.mean(drifts)), "slope_mean": float(np.mean(slopes)),
                "within_cell_basins": int(nb)}
        summary["cells"].append(cell)
        print(f"{size+'/'+reg:22}{len(reps):<6}{start:+6.1f} {ea.mean():+6.1f}[{ea.min():+.1f},{ea.max():+.1f}]{'':5}"
              f"{np.mean(drifts):7.2f}{np.mean(slopes):+9.4f} {nb}")
    (ana / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\ncoil floor ~3.7. wrote {ana/'summary.json'}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0]); sub = p.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("run"); pr.add_argument("--seeds", type=int, default=3); sub.add_parser("analyze")
    if hasattr(sys.stdout, "reconfigure"): sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    a = p.parse_args()
    run(a.seeds) if a.cmd == "run" else analyze()


if __name__ == "__main__":
    main()
