"""Author + run distributable dwell-time experiments as resumable chunks (#34).

This is the dwell-time-specific layer on top of the generic ``chunk_store``:

* ``create`` turns a spec (shape, ligands, seeds, lengths) into the chunk DAG and
  publishes it. Docking + NAC-core truncation happen here (central, CPU) so every
  GPU chunk is a uniform, self-contained work unit.
* ``work`` is the contributor loop: pull a runnable chunk, run the matching
  ``md_relax`` step (build / equilibrate / segment), push the outputs. ``--n``
  bounds how many chunks one invocation runs so each call stays short; many workers
  can ``work`` the same published experiment concurrently (leasing prevents clashes).
* ``status`` shows progress; ``score`` merges each replica's segment frames into one
  trajectory and runs the existing GPU-free dwell scorer (``dwell_time``).

Chunk DAG per (shape, ligand) pair:  build → {equilibrate(seed) → segment(seed, i)…}.
Only the per-replica segment chain has an internal dependency; seeds and pairs are
independent — the workload fans out wide. See the approved plan and #34.
"""
from __future__ import annotations

import argparse
import collections
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import chunk_store as store

ROOT = Path(__file__).resolve().parent.parent
MD_RELAX = Path(__file__).resolve().parent / "md_relax.py"

# Sensible dwell defaults (overridable on `create`).
DEFAULT_EQUIL_PS = 100.0
DEFAULT_PROD_PS = 2000.0
DEFAULT_SEGMENT_PS = 100.0
DEFAULT_TRAJ_INTERVAL_PS = 20.0
DEFAULT_TEMPERATURE_K = 300.0
APO = "apo"   # the ligand label for the apo baseline arm


# ---------------------------------------------------------------------------
# Seed derivation — distinct but deterministic per chunk (reproducible runs).
# ---------------------------------------------------------------------------

def _segment_seed(base_seed: int, index: int) -> int:
    # Langevin noise is memoryless, so each segment may use a fresh stream; make
    # it deterministic so a re-run of the same chunk reproduces bit-for-bit.
    return (base_seed * 100003 + index + 1) & 0x7FFFFFFF


# ---------------------------------------------------------------------------
# DAG construction.
# ---------------------------------------------------------------------------

def _pair_tag(shape: str, ligand: str) -> str:
    return f"{shape}__{ligand}"


def enumerate_chunks(spec: dict) -> list[dict]:
    """Build the chunk DAG (ids, consumes/produces, params) from a spec. Pure —
    no IO. ``initial_artifacts`` (the core PDBs) are registered separately by
    :func:`cmd_create`."""
    shape = spec["shape"]
    seeds = spec["seeds"]
    prod_ps = spec["prod_ps"]
    segment_ps = spec["segment_ps"]
    n_seg = max(1, math.ceil(prod_ps / segment_ps))
    chunks: list[dict] = []

    for ligand in spec["ligands"]:
        pair = _pair_tag(shape, ligand)
        core = f"{pair}/core.pdb"
        sys_xml = f"{pair}/system.xml"
        solv = f"{pair}/solvated.pdb"
        smiles = spec["ligand_smiles"].get(ligand)  # None for apo

        chunks.append({
            "id": f"build__{pair}", "kind": "build",
            "consumes": [core], "produces": [sys_xml, solv],
            "params": {"ligand": ligand, "smiles": smiles, "rect_box": spec["rect_box"],
                       "is_complex": ligand != APO},
            "meta": {"shape": shape, "ligand": ligand},
        })

        for seed in seeds:
            sd = f"{pair}/s{seed}"
            state0 = f"{sd}/state_0.xml"
            chunks.append({
                "id": f"equil__{pair}__s{seed}", "kind": "equilibrate",
                "consumes": [sys_xml, solv], "produces": [state0],
                "params": {"equil_ps": spec["equil_ps"], "seed": seed,
                           "temperature_k": spec["temperature_k"]},
                "meta": {"shape": shape, "ligand": ligand, "seed": seed},
            })
            for i in range(n_seg):
                this_ps = segment_ps if (i + 1) * segment_ps <= prod_ps else (prod_ps - i * segment_ps)
                state_in = f"{sd}/state_{i}.xml"
                state_out = f"{sd}/state_{i+1}.xml"
                seg_out = f"{sd}/seg_{i}.pdb"
                chunks.append({
                    "id": f"seg__{pair}__s{seed}__{i}", "kind": "segment",
                    "consumes": [sys_xml, solv, state_in], "produces": [state_out, seg_out],
                    "params": {"segment_ps": round(this_ps, 4), "index": i,
                               "traj_interval_ps": spec["traj_interval_ps"],
                               "seed": _segment_seed(seed, i),
                               "temperature_k": spec["temperature_k"]},
                    "meta": {"shape": shape, "ligand": ligand, "seed": seed, "index": i},
                })
    return chunks


# ---------------------------------------------------------------------------
# create
# ---------------------------------------------------------------------------

def cmd_create(args) -> None:
    import dwell_time  # local: pulls Bio.PDB; truncation + (complex) docking
    shape = args.shape
    shape_pdb = ROOT / "results" / "oligomers" / f"{shape}.pdb"
    if not shape_pdb.exists():
        raise FileNotFoundError(f"shape PDB not found: {shape_pdb}")

    ligands = list(args.ligands)
    if args.apo and APO not in ligands:
        ligands = [APO] + ligands
    if not ligands:
        ligands = [APO]

    seeds = [args.base_seed + i for i in range(args.replicas)]
    lo, hi = dwell_time.CHUNK_RANGE

    spec = {
        "shape": shape, "arm": "mixed", "ligands": ligands,
        "ligand_smiles": {}, "seeds": seeds,
        "equil_ps": args.equil_ps, "prod_ps": args.prod_ps, "segment_ps": args.segment_ps,
        "traj_interval_ps": args.traj_interval_ps, "temperature_k": args.temperature_k,
        "rect_box": True, "core_range": [lo, hi],
    }

    # Central, CPU prep: produce each pair's NAC-core PDB (the build chunk's input
    # and the toxic reference for scoring). Apo truncates the shape; a real ligand
    # docks first, then truncates the complex.
    scratch = store.experiment_dir(args.exp_id) / "_prep"
    scratch.mkdir(parents=True, exist_ok=True)
    initial: dict[str, Path] = {}
    for ligand in ligands:
        pair = _pair_tag(shape, ligand)
        if ligand == APO:
            core = dwell_time._truncate_chunk(shape_pdb, scratch / f"{pair}_core.pdb", lo, hi)
        else:
            complex_pdb, smiles = dwell_time._ensure_complex(ligand, shape)
            spec["ligand_smiles"][ligand] = smiles
            core = dwell_time._truncate_chunk(complex_pdb, scratch / f"{pair}_core.pdb", lo, hi)
        initial[f"{pair}/core.pdb"] = core

    chunks = enumerate_chunks(spec)
    store.create_experiment(args.exp_id, spec, chunks, initial)
    s = store.status_summary(args.exp_id)
    kinds = ", ".join(f"{k}×{v['total']}" for k, v in s["by_kind"].items())
    print(f"created experiment {args.exp_id}: {s['n_chunks']} chunks ({kinds}) "
          f"over {len(ligands)} pair(s) × {len(seeds)} seed(s)", flush=True)
    print(f"  {args.prod_ps:.0f} ps prod / {args.segment_ps:.0f} ps segments "
          f"→ {math.ceil(args.prod_ps/args.segment_ps)} segments per replica", flush=True)


# ---------------------------------------------------------------------------
# work — pull / execute / push one chunk at a time
# ---------------------------------------------------------------------------

def _venv_python() -> str:
    return sys.executable  # run_chunks runs in the pip venv


def _conda_python() -> str:
    from md_env import md_python
    return str(md_python())


def _run(cmd: list[str], tag: str) -> None:
    # Stream the child's output live (unbuffered) so a contributor sees the MD
    # heartbeat — step, speed, elapsed — as it happens, instead of silence until
    # the chunk ends. Keep a tail of recent lines for the error message.
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    tail = collections.deque(maxlen=20)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, encoding="utf-8", errors="replace", bufsize=1, env=env)
    for line in proc.stdout:
        line = line.rstrip()
        tail.append(line)
        print(f"  · {line}", flush=True)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"{tag} failed (exit {proc.returncode}):\n" + "\n".join(tail))


def execute_chunk(ch: dict, infile, scratch: Path) -> dict[str, Path]:
    """Run the md_relax step for one chunk; return {artifact_id: produced file}.

    ``infile(aid)`` resolves a consumed artifact id to a readable local path and
    ``scratch`` is where outputs are written — both supplied by the caller, so the
    same execution path serves the local worker (inputs from the on-disk store) and
    a remote contributor (inputs downloaded from the coordinator). No store access
    here; ``push`` / the coordinator content-addresses the returned files."""
    scratch = Path(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    p = ch["params"]
    kind = ch["kind"]

    if kind == "build":
        core = infile(ch["consumes"][0])
        prefix = scratch / ch["id"]
        cmd = [_venv_python(), str(MD_RELAX), "--rect-box", "--prepare-only", str(prefix)]
        if p["is_complex"]:
            cmd = [_conda_python(), str(MD_RELAX), "--complex-pdb", str(core),
                   "--ligand-smiles", p["smiles"], "--rect-box", "--prepare-only", str(prefix)]
        else:
            cmd = [_venv_python(), str(MD_RELAX), "--apo-pdb", str(core),
                   "--rect-box", "--prepare-only", str(prefix)]
        _run(cmd, ch["id"])
        sys_xml, solv = ch["produces"]
        return {sys_xml: Path(f"{prefix}_system.xml"), solv: Path(f"{prefix}_solvated.pdb")}

    if kind == "equilibrate":
        sys_xml, solv = ch["consumes"]
        state0 = ch["produces"][0]
        out = scratch / f"{ch['id']}_state.xml"
        cmd = [_venv_python(), str(MD_RELAX),
               "--system-xml", str(infile(sys_xml)), "--solvated-pdb", str(infile(solv)),
               "--equilibrate", str(out), "--equil-ps", str(p["equil_ps"]),
               "--temperature-k", str(p["temperature_k"]), "--seed", str(p["seed"]),
               "--report-interval-ps", "2"]
        _run(cmd, ch["id"])
        return {state0: out}

    if kind == "segment":
        sys_xml, solv, state_in = ch["consumes"]
        state_out, seg_out = ch["produces"]
        out_state = scratch / f"{ch['id']}_state.xml"
        out_seg = scratch / f"{ch['id']}_seg.pdb"
        cmd = [_venv_python(), str(MD_RELAX), "--segment",
               "--system-xml", str(infile(sys_xml)), "--solvated-pdb", str(infile(solv)),
               "--state-in", str(infile(state_in)), "--state-out", str(out_state),
               "--seg-out", str(out_seg), "--segment-ps", str(p["segment_ps"]),
               "--traj-interval-ps", str(p["traj_interval_ps"]),
               "--temperature-k", str(p["temperature_k"]), "--seed", str(p["seed"]),
               "--report-interval-ps", "2"]
        _run(cmd, ch["id"])
        return {state_out: out_state, seg_out: out_seg}

    raise ValueError(f"unknown chunk kind {kind!r}")


def cmd_work(args) -> None:
    worker = args.worker or f"local-{int(time.time())}"
    done = 0
    while args.n <= 0 or done < args.n:
        ch = store.pull(args.exp_id, worker, lease_seconds=args.lease_seconds)
        if ch is None:
            print("no runnable chunk (experiment complete or all leased)", flush=True)
            break
        manifest = store.load_manifest(args.exp_id)
        print(f"\n=== chunk {ch['id']} [{ch['kind']}] ===", flush=True)
        t0 = time.time()
        scratch = store.experiment_dir(args.exp_id) / "_scratch"
        infile = lambda aid: store.artifact_file(args.exp_id, manifest, aid)  # noqa: E731
        try:
            outputs = execute_chunk(ch, infile, scratch)
        except Exception as e:  # noqa: BLE001 — record + move on
            store.mark_failed(args.exp_id, ch["id"], repr(e))
            print(f"  FAILED: {e}", flush=True)
            if args.stop_on_fail:
                raise
            done += 1
            continue
        wall = time.time() - t0
        store.push(args.exp_id, ch["id"], outputs, wall)
        print(f"  pushed {ch['id']} ({wall/60:.1f} min)", flush=True)
        done += 1

    s = store.status_summary(args.exp_id)
    print(f"\nworker {worker}: ran {done} chunk(s). "
          f"done {s['counts']['done']}/{s['n_chunks']}"
          f"{' — COMPLETE' if s['complete'] else ''}", flush=True)


# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------

def cmd_status(args) -> None:
    s = store.status_summary(args.exp_id)
    c = s["counts"]
    print(f"experiment {args.exp_id}: {c['done']}/{s['n_chunks']} done"
          f"{' — COMPLETE' if s['complete'] else ''}", flush=True)
    print(f"  pending={c['pending']} runnable={c['runnable']} leased={c['leased']} "
          f"failed={c['failed']}", flush=True)
    for kind, k in s["by_kind"].items():
        print(f"  {kind:12} {k['done']}/{k['total']}", flush=True)
    print(f"  GPU wall-clock so far: {s['wall_seconds_total']/60:.1f} min", flush=True)


# ---------------------------------------------------------------------------
# score — merge each replica's segments, run the GPU-free dwell scorer
# ---------------------------------------------------------------------------

def merge_segments(seg_paths: list[Path], out_path: Path) -> Path:
    """Concatenate per-segment frame PDBs into one multi-MODEL replica trajectory
    (renumbered MODEL serials), so the existing per-file dwell scorer sees one
    continuous replica. Keeps only ATOM/HETATM/TER records."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mi = 0
    with out_path.open("w") as out:
        for sp in seg_paths:
            for line in Path(sp).read_text().splitlines():
                rec = line[:6]
                if rec == "MODEL ":
                    mi += 1
                    out.write(f"MODEL     {mi:>4}\n")
                elif rec == "ENDMDL":
                    out.write("ENDMDL\n")
                elif rec in ("ATOM  ", "HETATM", "TER   "):
                    out.write(line + "\n")
        out.write("END\n")
    return out_path


def cmd_score(args) -> None:
    import dwell_time
    manifest = store.load_manifest(args.exp_id)
    spec = manifest["spec"]
    shape = spec["shape"]
    n_seg = max(1, math.ceil(spec["prod_ps"] / spec["segment_ps"]))
    merged_dir = store.experiment_dir(args.exp_id) / "_scored"
    merged_dir.mkdir(parents=True, exist_ok=True)

    # Reference = the apo NAC-core (the toxic shape); fall back to first ligand.
    ref_ligand = APO if APO in spec["ligands"] else spec["ligands"][0]
    ref = store.artifact_file(args.exp_id, manifest, f"{_pair_tag(shape, ref_ligand)}/core.pdb")

    # Merge each (ligand, seed) replica's segments into one trajectory.
    per_ligand_trajs: dict[str, list[Path]] = {}
    for ligand in spec["ligands"]:
        pair = _pair_tag(shape, ligand)
        trajs: list[Path] = []
        for seed in spec["seeds"]:
            seg_ids = [f"{pair}/s{seed}/seg_{i}.pdb" for i in range(n_seg)]
            if not all(manifest["artifacts"].get(a, {}).get("present") for a in seg_ids):
                continue  # replica not finished yet
            seg_files = [store.artifact_file(args.exp_id, manifest, a) for a in seg_ids]
            traj = merge_segments(seg_files, merged_dir / f"{pair}_s{seed}.pdb")
            trajs.append(traj)
        per_ligand_trajs[ligand] = trajs

    # Per-replica dwell fractions (always — validates merge + scorer).
    reference = dwell_time.load_pdb(ref)
    print(f"reference: {ref.name}\n", flush=True)
    for ligand, trajs in per_ligand_trajs.items():
        if not trajs:
            print(f"  {ligand:14} — no completed replicas yet", flush=True)
            continue
        fracs = [dwell_time.score_trajectory(t, reference)["dwell_fraction"] for t in trajs]
        print(f"  {ligand:14} dwell fractions {[round(x, 3) for x in fracs]} "
              f"(n={len(fracs)})", flush=True)

    # Pairwise shift vs apo if both arms present.
    apo_trajs = per_ligand_trajs.get(APO, [])
    others = [l for l in spec["ligands"] if l != APO and per_ligand_trajs.get(l)]
    if apo_trajs and others:
        print("\ndwell shift vs apo:", flush=True)
        for ligand in others:
            pair = dwell_time.summarise_pair(ref, apo_trajs, per_ligand_trajs[ligand],
                                             seed=hash((shape, ligand)) & 0xFFFF)
            b = pair["bootstrap"]
            print(f"  {ligand:14} shift={b['shift']:+.3f} "
                  f"CI[{b['ci_low']:+.3f},{b['ci_high']:+.3f}] {b['classification']}", flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    pc = sub.add_parser("create", help="author + publish a chunk experiment")
    pc.add_argument("exp_id")
    pc.add_argument("--shape", required=True, help="oligomer stem under results/oligomers/")
    pc.add_argument("--ligands", nargs="*", default=[], help="vicinity-molecule ids (docked)")
    pc.add_argument("--apo", action="store_true", help="include the apo baseline arm")
    pc.add_argument("--replicas", type=int, default=2)
    pc.add_argument("--base-seed", type=int, default=1000)
    pc.add_argument("--equil-ps", type=float, default=DEFAULT_EQUIL_PS)
    pc.add_argument("--prod-ps", type=float, default=DEFAULT_PROD_PS)
    pc.add_argument("--segment-ps", type=float, default=DEFAULT_SEGMENT_PS)
    pc.add_argument("--traj-interval-ps", type=float, default=DEFAULT_TRAJ_INTERVAL_PS)
    pc.add_argument("--temperature-k", type=float, default=DEFAULT_TEMPERATURE_K)
    pc.set_defaults(func=cmd_create)

    pw = sub.add_parser("work", help="pull/run/push runnable chunks one at a time")
    pw.add_argument("exp_id")
    pw.add_argument("--n", type=int, default=1, help="max chunks to run (<=0 = drain)")
    pw.add_argument("--worker", default=None)
    pw.add_argument("--lease-seconds", type=float, default=store.DEFAULT_LEASE_SECONDS)
    pw.add_argument("--stop-on-fail", action="store_true")
    pw.set_defaults(func=cmd_work)

    ps = sub.add_parser("status", help="progress summary")
    ps.add_argument("exp_id")
    ps.set_defaults(func=cmd_status)

    pscore = sub.add_parser("score", help="merge replica segments + run the dwell scorer")
    pscore.add_argument("exp_id")
    pscore.set_defaults(func=cmd_score)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
