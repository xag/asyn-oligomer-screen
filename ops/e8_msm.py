#!/usr/bin/env python3
"""E8 analysis: occupancy MSM over the oligomer shape coordinates (issue #57).

From the swarm of short, velocity-seeded replicas (ops/e8_launch.ps1, collected
by ops/e8_collect.ps1), build a per-register Markov state model and read off the
*occupancy* of each shape basin — the populations single-quench MD (E5-E7) could
not give. Pipeline, per register (parallel / antiparallel):

  1. featurize every frame of every replica against that register's relaxed
     reference, reusing screen/shape_metrics (β-core Cα RMSD, inter-chain contact
     Jaccard, inter-chain contact count, β-core radius of gyration);
  2. cluster frames into microstates (k-means), count transitions at lag τ
     *within* each replica, build a row-stochastic T (Laplace-smoothed for
     ergodicity), solve for the stationary distribution π = microstate occupancy;
  3. lump microstates into macrostates (Ward on standardized centroids); a
     macrostate's occupancy is Σπ over its microstates;
  4. for each macrostate, extract the medoid frame as a PDB and score its toxic
     activity with oligomers/score_oligomer.py.

Output: results/msm_states/analysis/<register>_macrostates.csv + medoid PDBs, and
a combined occupancy_summary.json — the occupancy-weighted basins on which
channels 1 & 2 should be re-grounded.

Run: uv run python ops/e8_msm.py [--lag 4] [--micro 80] [--macro 5]
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.cluster.vq import kmeans2, whiten  # noqa: F401  (whiten unused; explicit std below)

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "screen"))
from shape_metrics import (  # noqa: E402
    BETA_CORE_RANGE, NAC_RANGE, CONTACT_CUTOFF,
    beta_core_rmsd, contact_jaccard, interchain_contact_set, _ca_by_residue, load_pdb,
)

TRAJ_ROOT = REPO / "results" / "msm_states" / "traj"
OUT = REPO / "results" / "msm_states" / "analysis"
REF = {
    "parallel": REPO / "results" / "oligomers" / "fusco_parallel_3mer_core70-88_relaxed.pdb",
    "antiparallel": REPO / "results" / "oligomers" / "fusco_antiparallel_3mer_core70-88_relaxed.pdb",
}
SCORE = REPO / "oligomers" / "score_oligomer.py"


def split_models(seg_pdb: Path) -> list[str]:
    """Return the text of each MODEL...ENDMDL block in a multi-model PDB."""
    blocks, cur, inmodel = [], [], False
    head = []
    for line in seg_pdb.read_text().splitlines(keepends=True):
        if line.startswith("MODEL"):
            inmodel, cur = True, [line]
        elif line.startswith("ENDMDL"):
            cur.append(line)
            blocks.append("".join(head) + "".join(cur))
            inmodel = False
        elif inmodel:
            cur.append(line)
        elif line.startswith(("CRYST1", "REMARK", "HEADER")):
            head.append(line)
    if not blocks:  # single-model PDB
        blocks = [seg_pdb.read_text()]
    return blocks


def frame_features(model, ref) -> np.ndarray:
    """[β-core RMSD, contact Jaccard, #inter-chain contacts, β-core Rg]."""
    rmsd = beta_core_rmsd(model, ref, BETA_CORE_RANGE)
    jac = contact_jaccard(model, ref, CONTACT_CUTOFF, NAC_RANGE)
    nc = len(interchain_contact_set(model, CONTACT_CUTOFF, NAC_RANGE))
    ca = _ca_by_residue(model, BETA_CORE_RANGE, chain_ids=None)
    coords = np.array(list(ca.values()), dtype=float)
    rg = float(np.sqrt(((coords - coords.mean(0)) ** 2).sum(1).mean())) if len(coords) else 0.0
    return np.array([rmsd, jac, nc, rg], dtype=float)


def replica_dirs() -> list[Path]:
    return sorted({p.parent for p in TRAJ_ROOT.rglob("seg_*.pdb")})


def register_of(repdir: Path) -> str | None:
    name = repdir.name
    if name.startswith("par_"):
        return "parallel"
    if name.startswith("anti_"):
        return "antiparallel"
    return None


def featurize(register: str):
    """Featurize all replicas of one register. Returns (X, traj_slices, frame_index)."""
    ref = load_pdb(REF[register])
    X, slices, index = [], [], []
    for rd in replica_dirs():
        if register_of(rd) != register:
            continue
        segs = sorted(rd.glob("seg_*.pdb"), key=lambda p: int(re.search(r"seg_(\d+)", p.name).group(1)))
        start = len(X)
        for seg in segs:
            for mi, block in enumerate(split_models(seg)):
                tmp = OUT / "_frame.pdb"
                tmp.write_text(block)
                try:
                    feats = frame_features(load_pdb(tmp), ref)
                except Exception as e:  # a corrupt/short frame must not sink the run
                    print(f"  skip {seg.name}#{mi} in {rd.name}: {e}")
                    continue
                X.append(feats)
                index.append({"replica": rd.name, "seg": seg.name, "model": mi})
        if len(X) > start:
            slices.append((start, len(X)))
        print(f"  {rd.name}: {len(X) - start} frames")
    return np.array(X), slices, index


def build_msm(X, slices, n_micro, lag):
    """k-means microstates, within-trajectory transition counts at `lag`,
    Laplace-smoothed row-stochastic T, stationary π."""
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    k = min(n_micro, max(2, len(X) // 50))
    centroids, labels = kmeans2(Z, k, minit="++", seed=0, missing="warn")
    used = np.unique(labels)
    remap = {c: i for i, c in enumerate(used)}
    labels = np.array([remap[c] for c in labels])
    k = len(used)
    C = np.zeros((k, k))
    for a, b in slices:
        seq = labels[a:b]
        for t in range(len(seq) - lag):
            C[seq[t], seq[t + lag]] += 1
    C += 1.0 / k  # Laplace prior -> ergodic, π well-defined
    T = C / C.sum(1, keepdims=True)
    w, v = np.linalg.eig(T.T)
    pi = np.real(v[:, np.argmin(np.abs(w - 1.0))])
    pi = np.abs(pi) / np.abs(pi).sum()
    micro_centroid = np.array([X[labels == i].mean(0) for i in range(k)])
    return labels, pi, micro_centroid, (mu, sd)


def macrostates(micro_centroid, pi, std, n_macro):
    mu, sd = std
    Zc = (micro_centroid - mu) / sd
    n_macro = min(n_macro, len(Zc))
    if len(Zc) <= n_macro:
        macro = np.arange(len(Zc))
    else:
        macro = fcluster(linkage(Zc, method="ward"), t=n_macro, criterion="maxclust") - 1
    out = []
    for m in np.unique(macro):
        members = np.where(macro == m)[0]
        occ = float(pi[members].sum())
        # representative microstate = highest-π member
        rep_micro = members[np.argmax(pi[members])]
        out.append({"macro": int(m), "occupancy": occ, "members": members.tolist(),
                    "rep_micro": int(rep_micro)})
    out.sort(key=lambda d: -d["occupancy"])
    return macro, out


def score_pdb(pdb: Path) -> float | None:
    try:
        r = subprocess.run([sys.executable, str(SCORE), str(pdb)],
                           capture_output=True, text=True, timeout=600)
        m = re.search(r"activity score:\s*([+-]?\d+\.\d+)", r.stdout)
        return float(m.group(1)) if m else None
    except Exception as e:
        print(f"  score failed for {pdb.name}: {e}")
        return None


def extract_medoid_pdb(X, labels, index, rep_micro, dest: Path):
    """Write the frame closest to the rep microstate's centroid as a PDB."""
    members = np.where(labels == rep_micro)[0]
    centroid = X[members].mean(0)
    best = members[np.argmin(((X[members] - centroid) ** 2).sum(1))]
    info = index[best]
    seg = next(TRAJ_ROOT.rglob(f"{info['replica']}/**/{info['seg']}"), None) \
        or next(TRAJ_ROOT.rglob(f"**/{info['replica']}/{info['seg']}"), None)
    blocks = split_models(seg)
    dest.write_text(blocks[info["model"]])
    return info, X[best]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lag", type=int, default=4, help="MSM lag in frames (×50 ps)")
    ap.add_argument("--micro", type=int, default=80, help="target microstates")
    ap.add_argument("--macro", type=int, default=5, help="macrostates per register")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    summary = {}
    for register in ("parallel", "antiparallel"):
        print(f"\n=== {register} ===")
        X, slices, index = featurize(register)
        if len(X) < 100:
            print(f"  too few frames ({len(X)}); skipping {register}")
            continue
        labels, pi, micro_centroid, std = build_msm(X, slices, args.micro, args.lag)
        macro, macros = macrostates(micro_centroid, pi, std, args.macro)
        rows = []
        for d in macros:
            frames = np.isin(labels, d["members"])
            feat_mean = X[frames].mean(0)
            medoid = OUT / f"{register}_macro{d['macro']}_occ{d['occupancy']:.2f}_medoid.pdb"
            info, fvec = extract_medoid_pdb(X, labels, index, d["rep_micro"], medoid)
            activity = score_pdb(medoid)
            row = {"register": register, "macro": d["macro"], "occupancy": round(d["occupancy"], 4),
                   "n_frames": int(frames.sum()),
                   "beta_core_rmsd": round(float(feat_mean[0]), 2),
                   "contact_jaccard": round(float(feat_mean[1]), 3),
                   "interchain_contacts": round(float(feat_mean[2]), 1),
                   "beta_core_rg": round(float(feat_mean[3]), 2),
                   "activity": activity, "medoid_pdb": medoid.name,
                   "medoid_from": f"{info['replica']}/{info['seg']}#{info['model']}"}
            rows.append(row)
            print(f"  macro{d['macro']}: occ={d['occupancy']:.3f} rmsd={feat_mean[0]:.2f} "
                  f"jac={feat_mean[1]:.3f} act={activity}")
        # CSV
        csv = OUT / f"{register}_macrostates.csv"
        keys = list(rows[0].keys())
        csv.write_text(",".join(keys) + "\n" +
                       "\n".join(",".join(str(r[k]) for k in keys) for r in rows) + "\n")
        summary[register] = {"n_frames": len(X), "n_replicas": len(slices), "macrostates": rows}

    (OUT / "occupancy_summary.json").write_text(json.dumps(summary, indent=2))
    if (OUT / "_frame.pdb").exists():
        (OUT / "_frame.pdb").unlink()
    print(f"\nwrote {OUT/'occupancy_summary.json'}")


if __name__ == "__main__":
    main()
