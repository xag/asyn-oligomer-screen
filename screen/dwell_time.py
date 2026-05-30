"""Shape-stability (dwell-time) channel — issue #14.

The static Δactivity channel (stage3.py) is sign-bound: any docked pose
occludes SASA and adds contacts, so it can only ever see *destabilisers*
and never *stabilisers* (anti-targets, #30). This channel removes that
bound by asking a different question:

    Across many short velocity-seeded MD replicas, does the ligand keep
    the oligomer in its toxic shape (β-core intact, chains packed) or push
    it out of that basin?

Per (shape, ligand) pair we run N replicas of the *apo* oligomer and N of
the *docked complex*, each a few ns of NPT MD with a different velocity
seed, dumping a trajectory frame every Δt. Each frame is scored with
shape_metrics (β-core Cα RMSD + inter-chain contact Jaccard vs the toxic
reference) and labelled in-basin / out-of-basin. The per-replica
*dwell fraction* = fraction of frames still in the toxic basin. We then
bootstrap the shift:

    dwell_shift = mean(complex dwell) − mean(apo dwell)

    shift < 0  → the ligand spends *less* time toxic → destabiliser
    shift > 0  → the ligand spends *more* time toxic → stabiliser / anti-target

The sign is free both ways — that is the whole point, and it is what #30
post-processes into the "what to avoid" list.

A per-replica binding-occupancy check (shape_metrics.ligand_bound) guards
the complex side: a ligand that diffuses off the β-core is flagged (low
occupancy) rather than silently scored as "no effect".

Two execution surfaces:

  * Analysis (this pip venv, no MD): `score` and `selftest` subcommands
    operate on trajectory PDBs that already exist. Fully runnable here.

  * MD (`pilot` subcommand): runs the replicas via md_relax.py in the conda
    MD env (`environment-md.yml`, found automatically by md_env.py). This is
    the cluster-scale step the issue flags; the pilot is the cheap decision
    gate before the full sweep.

Usage:
    # score trajectories that already exist (pip venv, no GPU):
    python screen/dwell_time.py score \\
        --reference results/oligomers/<shape>_relaxed.pdb \\
        --apo results/dwell/<shape>/apo_rep*.pdb \\
        --complex results/dwell/<shape>/silibinin_rep*.pdb

    # self-test the bootstrap on synthetic dwell fractions (no inputs):
    python screen/dwell_time.py selftest

    # full pilot (runs replicas in the conda MD env, found automatically):
    python screen/dwell_time.py pilot
"""
from __future__ import annotations

import argparse
import glob
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from shape_metrics import (  # noqa: E402
    BETA_CORE_RANGE,
    CONTACT_CUTOFF,
    NAC_RANGE,
    OCCUPANCY_CUTOFF,
    TOXIC_JACCARD_MIN,
    TOXIC_RMSD_MAX,
    frame_metrics,
    in_toxic_basin,
    ligand_bound,
    load_pdb,
)

ROOT = Path(__file__).resolve().parents[1]
DWELL_DIR = ROOT / "results" / "dwell"

# ---------------------------------------------------------------------------
# Pilot configuration — the cheap "decide whether to commit" run from #14:
# 4 ligands × 2 shapes × N replicas × a few ns. silibinin (top hit) and DHEA
# (top novel) should read as destabilisers; trehalose is a chemical chaperone
# with a different mechanism; caffeine is the negative-control decoy (the
# static screen ranks it ~neutral, Δgat≈-0.006, so its dwell shift should
# straddle zero). Two structurally distinct shapes (parallel vs antiparallel
# β-core) give the multi-shape consistency check #30 requires.
# ---------------------------------------------------------------------------

PILOT_SHAPES = [
    "fusco_parallel_3mer_core70-88_relaxed",
    "fusco_antiparallel_3mer_core70-88_relaxed",
]
PILOT_LIGANDS = ["silibinin", "dhea", "trehalose", "caffeine"]

N_REPLICAS = 10
EQUIL_PS = 100.0
PROD_NS = 2.0
FRAME_INTERVAL_PS = 20.0       # 2 ns / 20 ps = 100 frames per replica
MIN_OCCUPANCY = 0.5            # complex replica is "reliable" if ligand stays bound ≥50% of frames
N_BOOT = 10000


# ---------------------------------------------------------------------------
# Trajectory scoring (pip venv — no MD needed).
# ---------------------------------------------------------------------------

def iter_frames(traj_pdb: Path):
    """Yield each MODEL of a multi-MODEL trajectory PDB as a Biopython Model."""
    structure = load_pdb(traj_pdb)
    for model in structure:
        yield model


def score_trajectory(
    traj_pdb: Path,
    reference,
    core_range: tuple[int, int] = BETA_CORE_RANGE,
    contact_range: tuple[int, int] | None = NAC_RANGE,
    cutoff: float = CONTACT_CUTOFF,
    rmsd_max: float = TOXIC_RMSD_MAX,
    jaccard_min: float = TOXIC_JACCARD_MIN,
    occupancy_cutoff: float = OCCUPANCY_CUTOFF,
) -> dict:
    """Score every frame of one replica trajectory against the toxic
    reference shape. Returns per-frame metrics plus the replica's dwell
    fraction and (for complexes) occupancy fraction.

    For complexes the dwell fraction is computed over *bound* frames only,
    so a ligand that diffuses off does not get credited with the apo
    dynamics of the frames it spent away from the site. ``occupancy``
    surfaces how often it stayed bound."""
    rmsds: list[float] = []
    jaccards: list[float] = []
    basin: list[bool] = []
    bound: list[bool | None] = []
    for model in iter_frames(traj_pdb):
        m = frame_metrics(model, reference, core_range, contact_range, cutoff)
        rmsds.append(m["beta_core_rmsd"])
        jaccards.append(m["contact_jaccard"])
        basin.append(in_toxic_basin(m, rmsd_max, jaccard_min))
        bound.append(ligand_bound(model, core_range, occupancy_cutoff))

    n = len(basin)
    has_ligand = any(b is not None for b in bound)
    if has_ligand:
        bound_mask = [bool(b) for b in bound]
        n_bound = sum(bound_mask)
        occupancy = n_bound / n if n else 0.0
        # Dwell over bound frames only; if the ligand never stays, dwell is
        # undefined (nan) and the replica is flagged via low occupancy.
        in_basin_bound = [bs for bs, bd in zip(basin, bound_mask) if bd]
        dwell = float(np.mean(in_basin_bound)) if in_basin_bound else float("nan")
    else:
        occupancy = float("nan")
        n_bound = n
        dwell = float(np.mean(basin)) if n else float("nan")

    return {
        "traj": traj_pdb.name,
        "n_frames": n,
        "dwell_fraction": dwell,
        "occupancy": occupancy,
        "n_bound": n_bound,
        "mean_rmsd": float(np.mean(rmsds)) if rmsds else float("nan"),
        "mean_jaccard": float(np.mean(jaccards)) if jaccards else float("nan"),
        "per_frame": {
            "beta_core_rmsd": rmsds,
            "contact_jaccard": jaccards,
            "in_basin": [bool(x) for x in basin],
            "bound": [None if b is None else bool(b) for b in bound],
        },
    }


# ---------------------------------------------------------------------------
# Bootstrap of the dwell shift (pip venv — pure numpy, unit-testable).
# ---------------------------------------------------------------------------

def bootstrap_dwell_shift(
    apo_dwell: list[float],
    complex_dwell: list[float],
    n_boot: int = N_BOOT,
    seed: int = 0,
    alpha: float = 0.05,
) -> dict:
    """Bootstrap the dwell shift = mean(complex) − mean(apo) by resampling
    replicas with replacement. Returns the point estimate, a (1−alpha) CI,
    the probabilities that the ligand is a destabiliser (shift<0) or
    stabiliser (shift>0), and a three-way classification:

        destabiliser  — CI entirely below 0  (spends less time toxic)
        stabiliser    — CI entirely above 0  (spends more time toxic; anti-target)
        inconclusive  — CI straddles 0

    NaN dwell fractions (replicas where the ligand never stayed bound) are
    dropped before resampling and counted in ``n_apo`` / ``n_complex``."""
    apo = np.array([x for x in apo_dwell if np.isfinite(x)], dtype=float)
    cpx = np.array([x for x in complex_dwell if np.isfinite(x)], dtype=float)
    if len(apo) == 0 or len(cpx) == 0:
        return {
            "shift": float("nan"), "ci_low": float("nan"), "ci_high": float("nan"),
            "prob_destabiliser": float("nan"), "prob_stabiliser": float("nan"),
            "classification": "no-data", "n_apo": int(len(apo)), "n_complex": int(len(cpx)),
            "n_dropped_apo": int(np.sum(~np.isfinite(apo_dwell))) if len(apo_dwell) else 0,
            "n_dropped_complex": int(np.sum(~np.isfinite(complex_dwell))) if len(complex_dwell) else 0,
        }
    rng = np.random.default_rng(seed)
    a_idx = rng.integers(0, len(apo), size=(n_boot, len(apo)))
    c_idx = rng.integers(0, len(cpx), size=(n_boot, len(cpx)))
    shifts = cpx[c_idx].mean(axis=1) - apo[a_idx].mean(axis=1)
    shift = float(cpx.mean() - apo.mean())
    ci_low = float(np.quantile(shifts, alpha / 2))
    ci_high = float(np.quantile(shifts, 1 - alpha / 2))
    prob_destab = float(np.mean(shifts < 0))
    prob_stab = float(np.mean(shifts > 0))
    if ci_high < 0:
        cls = "destabiliser"
    elif ci_low > 0:
        cls = "stabiliser"
    else:
        cls = "inconclusive"
    return {
        "shift": shift, "ci_low": ci_low, "ci_high": ci_high,
        "prob_destabiliser": prob_destab, "prob_stabiliser": prob_stab,
        "classification": cls, "n_apo": int(len(apo)), "n_complex": int(len(cpx)),
        "n_dropped_apo": int(np.sum(~np.isfinite(np.asarray(apo_dwell, dtype=float)))),
        "n_dropped_complex": int(np.sum(~np.isfinite(np.asarray(complex_dwell, dtype=float)))),
    }


def summarise_pair(
    reference_pdb: Path,
    apo_trajs: list[Path],
    complex_trajs: list[Path],
    seed: int = 0,
    **score_kw,
) -> dict:
    """Score apo + complex replica trajectories and bootstrap the shift."""
    reference = load_pdb(reference_pdb)
    apo_scores = [score_trajectory(t, reference, **score_kw) for t in apo_trajs]
    cpx_scores = [score_trajectory(t, reference, **score_kw) for t in complex_trajs]
    apo_dwell = [s["dwell_fraction"] for s in apo_scores]
    cpx_dwell = [s["dwell_fraction"] for s in cpx_scores]
    boot = bootstrap_dwell_shift(apo_dwell, cpx_dwell, seed=seed)
    # Average occupancy over the complex replicas that actually carried a
    # ligand. A nan occupancy means "no LIG residue in this trajectory"
    # (setup issue), which is distinct from a finite-but-low occupancy
    # ("ligand diffused off the site"). Keep them apart so the report does
    # not blame a missing ligand on diffusion — and so np.mean is never
    # handed an all-nan slice.
    finite_occ = [s["occupancy"] for s in cpx_scores if np.isfinite(s["occupancy"])]
    boot["ligand_present"] = bool(finite_occ)
    mean_occ = float(np.mean(finite_occ)) if finite_occ else float("nan")
    boot["mean_complex_occupancy"] = mean_occ
    boot["occupancy_ok"] = bool(np.isfinite(mean_occ) and mean_occ >= MIN_OCCUPANCY)
    return {
        "reference": reference_pdb.name,
        "apo_dwell": apo_dwell,
        "complex_dwell": cpx_dwell,
        "apo_scores": apo_scores,
        "complex_scores": cpx_scores,
        "bootstrap": boot,
    }


# ---------------------------------------------------------------------------
# MD replica orchestration (runs md_relax.py in the conda MD env).
# ---------------------------------------------------------------------------

def _md_python() -> Path:
    from md_env import md_python  # noqa: E402  (resolves the conda MD env)
    return md_python()


def _run_md(cmd: list[str], out_pdb: Path, tag: str) -> None:
    print("  $ " + " ".join(cmd), flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.stdout:
        for line in res.stdout.splitlines():
            print(f"  md > {line}", flush=True)
    if res.returncode != 0 or not out_pdb.exists():
        if res.stderr:
            print(res.stderr, file=sys.stderr)
        raise RuntimeError(f"md_relax failed (rc={res.returncode}) for {tag}")


def run_apo_replicas(
    apo_pdb: Path, out_dir: Path, n_replicas: int = N_REPLICAS,
    equil_ps: float = EQUIL_PS, prod_ns: float = PROD_NS,
    frame_interval_ps: float = FRAME_INTERVAL_PS, base_seed: int = 1000,
    skip_existing: bool = True,
) -> list[Path]:
    md_python = _md_python()
    md_relax = Path(__file__).parent / "md_relax.py"
    out_dir.mkdir(parents=True, exist_ok=True)
    trajs: list[Path] = []
    for i in range(n_replicas):
        traj = out_dir / f"apo_rep{i:02d}.pdb"
        if skip_existing and traj.exists():
            print(f"  apo replica {i} cache hit: {traj.name}", flush=True)
            trajs.append(traj)
            continue
        cmd = [
            str(md_python), str(md_relax),
            "--apo-pdb", str(apo_pdb),
            "--out-pdb", str(out_dir / f"apo_rep{i:02d}_final.pdb"),
            "--equil-ps", str(equil_ps), "--prod-ps", str(prod_ns * 1000.0),
            "--seed", str(base_seed + i),
            "--traj-out", str(traj), "--traj-interval-ps", str(frame_interval_ps),
        ]
        _run_md(cmd, traj, f"apo_rep{i}")
        trajs.append(traj)
    return trajs


def run_complex_replicas(
    complex_pdb: Path, ligand_smiles: str, out_dir: Path, tag: str,
    n_replicas: int = N_REPLICAS, equil_ps: float = EQUIL_PS, prod_ns: float = PROD_NS,
    frame_interval_ps: float = FRAME_INTERVAL_PS, base_seed: int = 2000,
    skip_existing: bool = True,
) -> list[Path]:
    md_python = _md_python()
    md_relax = Path(__file__).parent / "md_relax.py"
    out_dir.mkdir(parents=True, exist_ok=True)
    trajs: list[Path] = []
    for i in range(n_replicas):
        traj = out_dir / f"{tag}_rep{i:02d}.pdb"
        if skip_existing and traj.exists():
            print(f"  {tag} replica {i} cache hit: {traj.name}", flush=True)
            trajs.append(traj)
            continue
        cmd = [
            str(md_python), str(md_relax),
            "--complex-pdb", str(complex_pdb), "--ligand-smiles", ligand_smiles,
            "--out-pdb", str(out_dir / f"{tag}_rep{i:02d}_final.pdb"),
            "--equil-ps", str(equil_ps), "--prod-ps", str(prod_ns * 1000.0),
            "--seed", str(base_seed + i),
            "--traj-out", str(traj), "--traj-interval-ps", str(frame_interval_ps),
        ]
        _run_md(cmd, traj, f"{tag}_rep{i}")
        trajs.append(traj)
    return trajs


def _ensure_complex(mol_id: str, shape_stem: str) -> tuple[Path, str]:
    """Dock mol_id onto the shape (via stage3) if the complex PDB is missing.
    Returns (complex_pdb, smiles)."""
    import stage3  # local import: pulls in rdkit/meeko, only needed for docking

    shape_pdb = ROOT / "results" / "oligomers" / f"{shape_stem}.pdb"
    pair_tag = f"{mol_id}_{shape_stem}"
    complex_pdb = ROOT / "results" / "stage3" / f"{pair_tag}_complex.pdb"
    smiles = stage3.load_vicinity_molecule(mol_id).get("smiles")
    if not smiles:
        raise ValueError(f"vicinity molecule {mol_id!r} has no SMILES; cannot dock")
    if not complex_pdb.exists():
        print(f"  docking {mol_id} onto {shape_stem} (stage3) ...", flush=True)
        stage3.perturb_oligomer(mol_id, shape_pdb)
    return complex_pdb, smiles


def run_pilot(
    shapes: list[str] = PILOT_SHAPES, ligands: list[str] = PILOT_LIGANDS,
    n_replicas: int = N_REPLICAS, prod_ns: float = PROD_NS,
    skip_existing: bool = True,
) -> dict:
    """The #14 pilot: shapes × ligands × replicas, apo + complex, then
    score + bootstrap. Writes results/dwell/pilot_summary.{json,csv}."""
    summary: dict = {"config": {
        "shapes": shapes, "ligands": ligands, "n_replicas": n_replicas,
        "prod_ns": prod_ns, "frame_interval_ps": FRAME_INTERVAL_PS,
        "rmsd_max": TOXIC_RMSD_MAX, "jaccard_min": TOXIC_JACCARD_MIN,
        "min_occupancy": MIN_OCCUPANCY,
    }, "pairs": []}

    for shape_stem in shapes:
        shape_pdb = ROOT / "results" / "oligomers" / f"{shape_stem}.pdb"
        if not shape_pdb.exists():
            raise FileNotFoundError(f"shape PDB not found: {shape_pdb}")
        shape_dir = DWELL_DIR / shape_stem
        print(f"\n=== shape {shape_stem}: {n_replicas} apo replicas ===", flush=True)
        apo_trajs = run_apo_replicas(shape_pdb, shape_dir, n_replicas, prod_ns=prod_ns,
                                     skip_existing=skip_existing)
        for mol_id in ligands:
            print(f"\n=== {mol_id} × {shape_stem} ===", flush=True)
            complex_pdb, smiles = _ensure_complex(mol_id, shape_stem)
            cpx_trajs = run_complex_replicas(complex_pdb, smiles, shape_dir, mol_id,
                                             n_replicas, prod_ns=prod_ns,
                                             skip_existing=skip_existing)
            pair = summarise_pair(shape_pdb, apo_trajs, cpx_trajs,
                                  seed=hash((shape_stem, mol_id)) & 0xFFFF)
            pair["shape"] = shape_stem
            pair["ligand"] = mol_id
            # Drop the bulky per-frame arrays from the persisted summary.
            for s in pair["apo_scores"] + pair["complex_scores"]:
                s.pop("per_frame", None)
            summary["pairs"].append(pair)
            _print_pair_line(pair)

    DWELL_DIR.mkdir(parents=True, exist_ok=True)
    (DWELL_DIR / "pilot_summary.json").write_text(json.dumps(summary, indent=2))
    _write_summary_csv(summary, DWELL_DIR / "pilot_summary.csv")
    print(f"\nwrote {DWELL_DIR / 'pilot_summary.json'} and pilot_summary.csv")
    return summary


# ---------------------------------------------------------------------------
# Reporting.
# ---------------------------------------------------------------------------

def _print_pair_line(pair: dict) -> None:
    b = pair["bootstrap"]
    occ = b.get("mean_complex_occupancy", float("nan"))
    if not b.get("ligand_present", True):
        flag = "  [NO LIGAND IN TRAJ]"
    elif not b.get("occupancy_ok", False):
        flag = "  [LOW OCCUPANCY]"
    else:
        flag = ""
    print(
        f"  {pair.get('ligand',''):14} shift={b['shift']:+.3f} "
        f"CI[{b['ci_low']:+.3f},{b['ci_high']:+.3f}] "
        f"{b['classification']:13} occ={occ:.2f}{flag}",
        flush=True,
    )


def _write_summary_csv(summary: dict, path: Path) -> None:
    import csv
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "shape", "ligand", "dwell_shift", "ci_low", "ci_high",
            "classification", "prob_destabiliser", "prob_stabiliser",
            "mean_complex_occupancy", "occupancy_ok", "n_apo", "n_complex",
        ])
        for p in summary["pairs"]:
            b = p["bootstrap"]
            w.writerow([
                p.get("shape", ""), p.get("ligand", ""),
                f"{b['shift']:.4f}", f"{b['ci_low']:.4f}", f"{b['ci_high']:.4f}",
                b["classification"], f"{b['prob_destabiliser']:.4f}",
                f"{b['prob_stabiliser']:.4f}", f"{b.get('mean_complex_occupancy', float('nan')):.4f}",
                b.get("occupancy_ok", False), b["n_apo"], b["n_complex"],
            ])


# ---------------------------------------------------------------------------
# CLI.
# ---------------------------------------------------------------------------

def _expand(globs: list[str]) -> list[Path]:
    out: list[Path] = []
    for g in globs:
        matched = [Path(p) for p in sorted(glob.glob(g))]
        out.extend(matched if matched else [Path(g)])
    return out


def _cmd_score(args) -> None:
    apo = _expand(args.apo)
    cpx = _expand(args.complex)
    pair = summarise_pair(Path(args.reference), apo, cpx, seed=args.seed)
    b = pair["bootstrap"]
    print(f"reference: {pair['reference']}")
    print(f"apo dwell fractions:     {[round(x,3) for x in pair['apo_dwell']]}")
    print(f"complex dwell fractions: {[round(x,3) for x in pair['complex_dwell']]}")
    if not b.get("ligand_present", True):
        occ_note = "no LIG residue in complex trajectories"
    elif b.get("occupancy_ok"):
        occ_note = "OK"
    else:
        occ_note = "LOW — ligand left the site"
    print(f"mean complex occupancy:  {b.get('mean_complex_occupancy', float('nan')):.3f} ({occ_note})")
    print(f"dwell shift: {b['shift']:+.4f}  CI95[{b['ci_low']:+.4f}, {b['ci_high']:+.4f}]")
    print(f"P(destabiliser)={b['prob_destabiliser']:.3f}  P(stabiliser)={b['prob_stabiliser']:.3f}")
    print(f"==> {b['classification'].upper()}")


def _cmd_selftest(args) -> None:
    """Bootstrap sanity check on synthetic dwell fractions (no MD)."""
    rng = np.random.default_rng(0)
    # Destabiliser: complex dwells much less than apo.
    apo = list(rng.uniform(0.8, 1.0, 10))
    destab = list(rng.uniform(0.1, 0.4, 10))
    stab_apo = list(rng.uniform(0.4, 0.6, 10))
    stab = list(rng.uniform(0.85, 1.0, 10))          # stabiliser: dwells more
    neutral = list(rng.uniform(0.4, 0.6, 10))         # decoy: same as apo

    cases = {
        "destabiliser": (apo, destab, "destabiliser"),
        "stabiliser": (stab_apo, stab, "stabiliser"),
        "neutral-decoy": (stab_apo, neutral, "inconclusive"),
    }
    ok = True
    for name, (a, c, expect) in cases.items():
        b = bootstrap_dwell_shift(a, c, seed=1)
        got = b["classification"]
        status = "ok" if got == expect else "FAIL"
        if got != expect:
            ok = False
        print(f"  {name:14} shift={b['shift']:+.3f} "
              f"CI[{b['ci_low']:+.3f},{b['ci_high']:+.3f}] -> {got:13} "
              f"(expect {expect}) [{status}]")
    print("selftest:", "PASS" if ok else "FAIL")
    if not ok:
        sys.exit(1)


def _cmd_pilot(args) -> None:
    run_pilot(n_replicas=args.replicas, prod_ns=args.prod_ns,
              skip_existing=not args.no_skip)


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    ps = sub.add_parser("score", help="score existing replica trajectories (pip venv)")
    ps.add_argument("--reference", required=True, help="toxic-shape reference PDB")
    ps.add_argument("--apo", nargs="+", required=True, help="apo replica trajectory PDB(s)/glob(s)")
    ps.add_argument("--complex", nargs="+", required=True, help="complex replica trajectory PDB(s)/glob(s)")
    ps.add_argument("--seed", type=int, default=0)
    ps.set_defaults(func=_cmd_score)

    pt = sub.add_parser("selftest", help="bootstrap sanity check on synthetic data (no MD)")
    pt.set_defaults(func=_cmd_selftest)

    pp = sub.add_parser("pilot", help="run the full #14 pilot (uses the conda MD env)")
    pp.add_argument("--replicas", type=int, default=N_REPLICAS)
    pp.add_argument("--prod-ns", type=float, default=PROD_NS)
    pp.add_argument("--no-skip", action="store_true", help="re-run replicas even if cached")
    pp.set_defaults(func=_cmd_pilot)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
