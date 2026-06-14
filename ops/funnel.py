"""Funnel test (E6): do diverse independent starting conformations relax to a
common ensemble, or to distinct, history-dependent basins?

E5 showed the one hand-built shape is not metastable but relaxes to a still-
toxic-scoring basin. That used a single starting region. This test launches
many *independent* starts — different coil draws, both registers, two shifted
beta-core windows, and an all-coil (disordered) start — runs each 50 ns
unbiased, and asks whether they converge.

  run      Unbiased seeded MD (no restraints, md_relax.py segment path) from
           each start's NAC-core chunk; one replica each. MD env via md_env.py.
  analyze  Per start: plateau slope (last 20 ns), drift from its own start,
           end-basin activity (score_oligomer). Then the cross-start question:
           pairwise beta-core RMSD among all relaxed medoids, hierarchical
           cluster counts at several cuts, and a classical-MDS map. Tight
           single cluster => a common attractor; register-split macro-basins
           => rugged, history-dependent landscape.

Usage:
    python ops/funnel.py run
    python ops/funnel.py analyze
"""
from __future__ import annotations
import argparse, json, re, subprocess, sys
from pathlib import Path
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))
from dwell_time import _truncate_chunk, CHUNK_RANGE
from shape_metrics import beta_core_rmsd, contact_jaccard, load_pdb

OLIGO = ROOT / "results" / "oligomers"
OUT = ROOT / "results" / "funnel_states"
SCORE = ROOT / "oligomers" / "score_oligomer.py"
VENV_PY = sys.executable
FRAME_PS = 50.0

STARTS = [
 "fusco_parallel_3mer_core70-88",
 "fusco_parallel_3mer_core70-88_s123",
 "fusco_parallel_3mer_core70-88_s777",
 "fusco_antiparallel_3mer_core70-88",
 "fusco_antiparallel_3mer_core70-88_s123",
 "fusco_parallel_3mer_core65-83",
 "fusco_parallel_3mer_core73-91",
 "ctrl_coil_3mer",
]
EQUIL_PS, SEGMENT_PS, N_SEGMENTS, TRAJ_PS = 200.0, 10000.0, 5, 50.0


def _chunk(shape):
    lo, hi = CHUNK_RANGE
    out = OUT / "inputs" / f"{shape}_core{lo}-{hi}.pdb"
    if not out.exists():
        _truncate_chunk(OLIGO / f"{shape}_relaxed.pdb", out, lo, hi)
    return out


def run(starts):
    from md_env import md_python
    md = str(md_python()); mdrelax = str(ROOT / "screen" / "md_relax.py")
    for shape in starts:
        chunk = str(_chunk(shape)); prefix = OUT / "prep" / shape
        prefix.parent.mkdir(parents=True, exist_ok=True)
        sysxml = Path(str(prefix) + "_system.xml")
        if not sysxml.exists():
            subprocess.run([md, mdrelax, "--apo-pdb", chunk, "--rect-box", "--prepare-only", str(prefix)], check=True)
        solv = str(Path(str(prefix) + "_solvated.pdb"))
        rep = OUT / f"{shape}_rep0"; rep.mkdir(parents=True, exist_ok=True)
        s_e = rep / "state_equil.xml"
        if not s_e.exists():
            subprocess.run([md, mdrelax, "--equilibrate", str(s_e), "--system-xml", str(sysxml),
                            "--solvated-pdb", solv, "--equil-ps", str(EQUIL_PS), "--seed", "1000"], check=True)
        prev = s_e
        for i in range(N_SEGMENTS):
            so, seg = rep / f"state_{i}.xml", rep / f"seg_{i}.pdb"
            if seg.exists() and so.exists(): prev = so; continue
            subprocess.run([md, mdrelax, "--segment", "--system-xml", str(sysxml), "--solvated-pdb", solv,
                            "--state-in", str(prev), "--state-out", str(so), "--seg-out", str(seg),
                            "--segment-ps", str(SEGMENT_PS), "--traj-interval-ps", str(TRAJ_PS),
                            "--seed", str(10000 + i)], check=True)
            prev = so
        (rep / "DONE").write_text("ok\n")


def _segs(d):
    return sorted(d.glob("seg_*.pdb"), key=lambda p: int(re.search(r"seg_(\d+)", p.name).group(1)))


def _activity(pdb):
    r = subprocess.run([VENV_PY, str(SCORE), str(pdb)], capture_output=True, text=True)
    m = re.search(r"activity score:\s*([+-]?\d+\.\d+)", r.stdout)
    return float(m.group(1)) if m else None


def analyze(starts):
    from Bio.PDB import PDBParser, PDBIO
    ana = OUT / "analysis"; ana.mkdir(parents=True, exist_ok=True)
    parser = PDBParser(QUIET=True)
    per, medoids = {}, {}
    for shape in starts:
        d = OUT / f"{shape}_rep0"
        if not _segs(d):
            print(f"[{shape}] no data; skip"); continue
        ref = load_pdb(_chunk(shape))
        rows, g = [], 0
        for seg in _segs(d):
            for li, m in enumerate(parser.get_structure(seg.stem, str(seg))):
                rows.append((str(seg), li, (g + 1) * FRAME_PS, beta_core_rmsd(m, ref), contact_jaccard(m, ref))); g += 1
        t = [r[2] for r in rows]; rr = [r[3] for r in rows]; tmax = t[-1]
        tail = [(x, y) for x, y in zip(t, rr) if x >= tmax - 20000]
        sl = float(np.polyfit(np.array([x for x, _ in tail]) / 1000., np.array([y for _, y in tail]), 1)[0])
        post = [r for r in rows if r[2] >= tmax / 2]
        X = np.array([[r[3], r[4]] for r in post]); cen = X.mean(0)
        med = post[int(np.argmin(((X - cen) ** 2).sum(1)))]
        mp = ana / f"{shape}_medoid.pdb"
        mo = list(parser.get_structure("m", med[0]))[med[1]]
        io = PDBIO(); io.set_structure(mo); io.save(str(mp)); medoids[shape] = mp
        per[shape] = {"t_ns": tmax / 1000., "slope_last20ns": sl,
                      "own_rmsd_mean_tail": float(np.mean([y for _, y in tail])), "activity": _activity(mp)}
        with open(ana / f"{shape}_rep0_timeseries.csv", "w") as f:
            f.write("t_ps,beta_core_rmsd_vs_own_start,contact_jaccard_vs_own_start\n")
            for r in rows: f.write(f"{r[2]:.0f},{r[3]:.4f},{r[4]:.4f}\n")
        print(f"  {shape:40} drift={per[shape]['own_rmsd_mean_tail']:.2f}A slope={sl:+.4f} act={per[shape]['activity']}")

    labels = list(medoids) + ["HANDBUILT_parallel_target"]
    st = {k: load_pdb(v) for k, v in medoids.items()}
    st["HANDBUILT_parallel_target"] = load_pdb(_chunk("fusco_parallel_3mer_core70-88"))
    n = len(labels); R = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            R[i, j] = R[j, i] = beta_core_rmsd(st[labels[i]], st[labels[j]])
    md_idx = [i for i, l in enumerate(labels) if l in medoids]
    sub = R[np.ix_(md_idx, md_idx)]
    Z = linkage(squareform(sub, checks=False), "average")
    clusters = {str(c): int(len(set(fcluster(Z, c, criterion="distance")))) for c in (4, 8, 12)}
    off = sub[np.triu_indices(len(md_idx), 1)]
    # classical MDS for the map
    J = np.eye(n) - np.ones((n, n)) / n; B = -0.5 * J @ (R ** 2) @ J
    w, V = np.linalg.eigh(B); idx = np.argsort(w)[::-1]
    mds = (V[:, idx[:2]] * np.sqrt(np.clip(w[idx[:2]], 0, None))).tolist()
    print(f"\nmean pairwise RMSD among {len(md_idx)} medoids: {off.mean():.2f} A "
          f"(min {off.min():.2f} max {off.max():.2f}); clusters {clusters}")
    summary = {"per_start": per, "labels": labels, "pairwise_rmsd": R.tolist(),
               "mds": mds, "clusters_by_cut": clusters, "mean_pairwise_rmsd": float(off.mean())}
    (ana / "funnel_summary.json").write_text(json.dumps(summary, indent=2))
    _map(labels, np.array(mds), ana / "funnel_map.png")
    print(f"wrote {ana/'funnel_summary.json'} and funnel_map.png")


def _map(labels, X, out):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    def short(l): return (l.replace("fusco_", "").replace("3mer_", "").replace("_3mer", "")
                          .replace("HANDBUILT_parallel_target", "hand-built target"))
    def col(l): return "#c1660a" if "antiparallel" in l else ("#2e7d32" if "coil" in l
                       else ("#000" if "HANDBUILT" in l else "#1f5fa8"))
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, l in enumerate(labels):
        ax.scatter(X[i, 0], X[i, 1], c=col(l), marker="*" if "HANDBUILT" in l else "o",
                   s=260 if "HANDBUILT" in l else 90, zorder=3, edgecolor="white", linewidth=.6)
        ax.annotate(short(l), (X[i, 0], X[i, 1]), fontsize=8.5, xytext=(6, 4), textcoords="offset points")
    ax.set_title("Funnel test — relaxed end-states (classical MDS of β-core RMSD)", fontsize=11)
    ax.set_xlabel("MDS-1 (Å)"); ax.set_ylabel("MDS-2 (Å)"); ax.grid(alpha=.15)
    ax.legend(handles=[Line2D([0], [0], marker='o', color='w', markerfacecolor=c, label=t, markersize=9)
              for c, t in [("#1f5fa8", "parallel"), ("#c1660a", "antiparallel"),
                           ("#2e7d32", "coil"), ("#000", "hand-built target")]], frameon=False, fontsize=8.5)
    fig.tight_layout(); fig.savefig(out, dpi=130)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("run").add_argument("--starts", nargs="+", default=STARTS)
    sub.add_parser("analyze").add_argument("--starts", nargs="+", default=STARTS)
    if hasattr(sys.stdout, "reconfigure"): sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    args = p.parse_args()
    (run if args.cmd == "run" else analyze)(args.starts)


if __name__ == "__main__":
    main()
