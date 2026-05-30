"""Oligomer ensemble sweep: build → relax → score → summary.

Sweeps 11 topologies across seed / arrangement / oligomer-size / β-core-range
axes, plus two negative controls (all-coil trimer and all-coil monomer) and
one truncated-scoring control (reference trimer scored on residues 60–100 only
to isolate the tail vs β-core contribution to the activity score).

Run from the pipeline/ directory with the pipeline venv:
    .venv\\Scripts\\python.exe oligomers/run_ensemble.py

Flags:
    --skip-existing   skip build+relax for topologies whose relaxed PDB already
                      exists (default: on — idempotent re-runs just re-score)
    --no-skip         force rebuild and re-relax even if relaxed PDB exists
    --dry-run         print the command plan and exit
    --summary-only    re-score all existing relaxed PDBs and rewrite the summary
                      without running any build/relax steps

Two Python binaries:
    VENV_PY — this script's own interpreter (pipeline venv; has PeptideBuilder,
               pydssp, biopython, pandas)
    MD_PY   — the conda MD env (environment-md.yml), located automatically by
               screen/md_env.py; used for md_relax.py.

MD relaxation parameters (same as the reference trimer, 2026-05-27):
    --vacuum-min-iter 5000   resolve inter-chain clashes before solvation
    --collapse-ps 500        OBC2 implicit-solvent MD; tails fold up
    --no-explicit            write heavy-atom PDB, skip TIP3P box
    --restrain-residues <range>   lock β-core Cα so topology survives collapse
    --restrain-chains <chains>    only restrain chains in the oligomer
    --restrain-k 1000        kJ/mol/nm² spring constant
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import time
from pathlib import Path

# Docstrings and help strings contain non-ASCII (α, β, →, Å, …) which would
# otherwise crash --help on Windows consoles defaulting to cp1252.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

PIPELINE   = Path(__file__).resolve().parents[1]           # pipeline/
OLIGO_DIR  = PIPELINE / "oligomers"
RESULTS    = PIPELINE / "results" / "oligomers"
RESULTS.mkdir(parents=True, exist_ok=True)

VENV_PY = sys.executable                                    # this interpreter
# The conda MD env is only needed for the build/relax subprocesses.
# `--summary-only` re-scores existing relaxed PDBs and never invokes MD, so
# resolving it is deferred into _require_md_py(), called from main() only when
# MD work will actually be performed.


def _require_md_py() -> Path:
    sys.path.insert(0, str(PIPELINE / "screen"))
    from md_env import md_python  # locates the conda MD env (environment-md.yml)
    py = md_python()
    if not py.exists():
        sys.exit(f"ERROR: MD conda python {py} does not exist.")
    return py

# ---------------------------------------------------------------------------
# Topology grid
# ---------------------------------------------------------------------------
# Each entry:
#   tag            — output filename stem (no extension)
#   n_mers         — 1 / 2 / 3 / 4
#   core           — (start, end) residue numbers; set start > end for all-coil
#   arrangement    — "parallel" | "antiparallel"
#   seed           — int; seeds the coil φ/ψ draws
#   label          — short human label for the summary table
# ---------------------------------------------------------------------------
TOPOLOGIES: list[dict] = [
    # ── reference (already relaxed; will be skipped when --skip-existing) ──
    dict(tag="fusco_parallel_3mer_core70-88",
         n_mers=3, core=(70, 88), arrangement="parallel", seed=42,
         label="parallel 3mer 70-88 s42 [ref]"),

    # ── seed variation: same topology, different coil conformations ──
    dict(tag="fusco_parallel_3mer_core70-88_s123",
         n_mers=3, core=(70, 88), arrangement="parallel", seed=123,
         label="parallel 3mer 70-88 s123"),
    dict(tag="fusco_parallel_3mer_core70-88_s777",
         n_mers=3, core=(70, 88), arrangement="parallel", seed=777,
         label="parallel 3mer 70-88 s777"),

    # ── arrangement variation ──
    dict(tag="fusco_antiparallel_3mer_core70-88",
         n_mers=3, core=(70, 88), arrangement="antiparallel", seed=42,
         label="antiparallel 3mer 70-88 s42"),
    dict(tag="fusco_antiparallel_3mer_core70-88_s123",
         n_mers=3, core=(70, 88), arrangement="antiparallel", seed=123,
         label="antiparallel 3mer 70-88 s123"),

    # ── oligomer size variation ──
    dict(tag="fusco_parallel_2mer_core70-88",
         n_mers=2, core=(70, 88), arrangement="parallel", seed=42,
         label="parallel 2mer 70-88 s42"),
    dict(tag="fusco_parallel_4mer_core70-88",
         n_mers=4, core=(70, 88), arrangement="parallel", seed=42,
         label="parallel 4mer 70-88 s42"),

    # ── β-core range variation ──
    dict(tag="fusco_parallel_3mer_core65-83",
         n_mers=3, core=(65, 83), arrangement="parallel", seed=42,
         label="parallel 3mer 65-83 s42"),
    dict(tag="fusco_parallel_3mer_core73-91",
         n_mers=3, core=(73, 91), arrangement="parallel", seed=42,
         label="parallel 3mer 73-91 s42"),

    # ── negative controls: all-coil (core_start > core_end → no β-strand) ──
    dict(tag="ctrl_coil_3mer",
         n_mers=3, core=(200, 0), arrangement="parallel", seed=42,
         label="[ctrl] all-coil 3mer"),
    dict(tag="ctrl_coil_monomer",
         n_mers=1, core=(200, 0), arrangement="parallel", seed=42,
         label="[ctrl] all-coil monomer"),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def chain_ids_for(n_mers: int) -> str:
    return ",".join(list("ABCDEFGH")[:n_mers])


def has_core(core: tuple[int, int]) -> bool:
    """False for all-coil controls (core_start > core_end)."""
    return core[0] <= core[1]


def build_cmd(topo: dict) -> list[str]:
    cs, ce = topo["core"]
    return [
        str(VENV_PY),
        str(OLIGO_DIR / "build_fusco_trimer.py"),
        "--n-mers",     str(topo["n_mers"]),
        "--core-start", str(cs),
        "--core-end",   str(ce),
        "--arrangement", topo["arrangement"],
        "--seed",       str(topo["seed"]),
        "--tag",        topo["tag"],
        "--out-dir",    str(RESULTS),
    ]


def relax_cmd(topo: dict, starting_pdb: Path, relaxed_pdb: Path) -> list[str]:
    cs, ce = topo["core"]
    cmd = [
        str(_require_md_py()),
        str(PIPELINE / "md_relax.py"),
        "--apo-pdb",         str(starting_pdb),
        "--out-pdb",         str(relaxed_pdb),
        "--vacuum-min-iter", "5000",
        "--collapse-ps",     "500",
        "--no-explicit",
        "--restrain-k",      "1000",
    ]
    if has_core(topo["core"]):
        cmd += ["--restrain-residues", f"{cs}-{ce}"]
        cmd += ["--restrain-chains",   chain_ids_for(topo["n_mers"])]
    return cmd


def score_structure(pdb: Path) -> float | None:
    """Run score_oligomer.py on a relaxed PDB; parse and return activity."""
    result = subprocess.run(
        [str(VENV_PY), str(OLIGO_DIR / "score_oligomer.py"), str(pdb)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"  [score] FAILED:\n{result.stderr[-800:]}", flush=True)
        return None
    m = re.search(r"activity score:\s*([+-]?\d+\.\d+)", result.stdout)
    if m:
        return float(m.group(1))
    print(f"  [score] could not parse activity from output:\n{result.stdout[-400:]}", flush=True)
    return None


def make_truncated_pdb(source_pdb: Path, res_start: int, res_end: int, out_pdb: Path) -> None:
    """Extract residues [res_start, res_end] from all chains and write a new PDB."""
    from Bio.PDB import PDBParser, PDBIO, Select

    class ResRange(Select):
        def accept_residue(self, residue):
            return res_start <= residue.get_id()[1] <= res_end

    parser = PDBParser(QUIET=True)
    s = parser.get_structure(source_pdb.stem, str(source_pdb))
    io = PDBIO()
    io.set_structure(s)
    io.save(str(out_pdb), ResRange())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-skip",      action="store_true",
                    help="force rebuild and re-relax even if output exists")
    ap.add_argument("--dry-run",      action="store_true",
                    help="print planned commands and exit")
    ap.add_argument("--summary-only", action="store_true",
                    help="re-score existing relaxed PDBs and rewrite the summary")
    args = ap.parse_args()
    skip_existing = not args.no_skip

    if not args.summary_only:
        _require_md_py()

    results: list[dict] = []
    total_t0 = time.time()

    # ── 1. Main topology sweep ─────────────────────────────────────────────
    for topo in TOPOLOGIES:
        tag          = topo["tag"]
        starting_pdb = RESULTS / f"{tag}.pdb"
        relaxed_pdb  = RESULTS / f"{tag}_relaxed.pdb"
        log_pdb      = RESULTS / f"_relax_{tag}.log"

        print(f"\n{'-'*60}", flush=True)
        print(f"  {topo['label']}", flush=True)
        print(f"{'-'*60}", flush=True)

        status = "ok"

        # ── build ──────────────────────────────────────────────────────────
        if args.summary_only:
            if not relaxed_pdb.exists():
                print(f"  [skip] no relaxed PDB at {relaxed_pdb}", flush=True)
                results.append({**topo, "activity": None, "status": "missing"})
                continue
        else:
            if skip_existing and starting_pdb.exists():
                print(f"  [build] already exists: {starting_pdb.name}", flush=True)
            else:
                cmd = build_cmd(topo)
                if args.dry_run:
                    print("  [build] " + " ".join(cmd), flush=True)
                else:
                    print(f"  [build] {starting_pdb.name} ...", flush=True)
                    r = subprocess.run(cmd, capture_output=True, text=True)
                    if r.returncode != 0:
                        print(f"  [build] FAILED:\n{r.stderr[-600:]}", flush=True)
                        status = "build-failed"
                    else:
                        print(r.stdout.strip(), flush=True)

            if status != "ok":
                results.append({**topo, "activity": None, "status": status})
                continue

            # ── relax ──────────────────────────────────────────────────────
            if skip_existing and relaxed_pdb.exists():
                print(f"  [relax] already exists: {relaxed_pdb.name}", flush=True)
            else:
                cmd = relax_cmd(topo, starting_pdb, relaxed_pdb)
                if args.dry_run:
                    print("  [relax] " + " ".join(cmd), flush=True)
                else:
                    print(f"  [relax] {relaxed_pdb.name} ... (log → {log_pdb.name})", flush=True)
                    t0 = time.time()
                    with open(log_pdb, "w") as lf:
                        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
                    elapsed = (time.time() - t0) / 60
                    if r.returncode != 0:
                        print(f"  [relax] FAILED after {elapsed:.1f} min; see {log_pdb.name}", flush=True)
                        status = "relax-failed"
                    else:
                        print(f"  [relax] done in {elapsed:.1f} min", flush=True)

            if status != "ok":
                results.append({**topo, "activity": None, "status": status})
                continue

        if args.dry_run:
            print(f"  [score] (dry-run)", flush=True)
            continue

        # ── score ──────────────────────────────────────────────────────────
        print(f"  [score] {relaxed_pdb.name} ...", flush=True)
        activity = score_structure(relaxed_pdb)
        if activity is None:
            status = "score-failed"
        else:
            print(f"  [score] activity = {activity:+.3f}", flush=True)
        results.append({**topo, "activity": activity, "status": status})

    # ── 2. Truncated scoring control ───────────────────────────────────────
    ref_tag     = "fusco_parallel_3mer_core70-88"
    ref_relaxed = RESULTS / f"{ref_tag}_relaxed.pdb"
    trunc_tag   = "ctrl_trunc_3mer_res60-100"
    trunc_pdb   = RESULTS / f"{trunc_tag}.pdb"

    print(f"\n{'-'*60}", flush=True)
    print(f"  [ctrl] truncated scoring: residues 60-100 of reference trimer", flush=True)
    print(f"{'-'*60}", flush=True)

    if not ref_relaxed.exists():
        print(f"  [skip] reference relaxed PDB not found: {ref_relaxed.name}", flush=True)
    elif args.dry_run:
        print(f"  [trunc] would extract res 60-100 from {ref_relaxed.name}", flush=True)
        print(f"  [score] (dry-run)", flush=True)
    else:
        print(f"  [trunc] extracting residues 60-100 → {trunc_pdb.name} ...", flush=True)
        make_truncated_pdb(ref_relaxed, 60, 100, trunc_pdb)
        print(f"  [score] {trunc_pdb.name} ...", flush=True)
        activity = score_structure(trunc_pdb)
        if activity is not None:
            print(f"  [score] activity = {activity:+.3f}", flush=True)
        results.append(dict(
            tag=trunc_tag, n_mers=3, core=(70, 88), arrangement="parallel", seed=42,
            label="[ctrl] truncated trimer res60-100",
            activity=activity,
            status="ok" if activity is not None else "score-failed",
        ))

    # ── 3. Summary ──────────────────────────────────────────────────────────
    if args.dry_run:
        return

    summary_csv = RESULTS / "ensemble_summary.csv"
    fields = ["tag", "label", "n_mers", "core", "arrangement", "seed", "activity", "status"]
    with open(summary_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in results:
            row = dict(r)
            row["core"] = f"{r['core'][0]}-{r['core'][1]}"
            w.writerow(row)
    print(f"\nWrote {summary_csv}", flush=True)

    # ── ranking table ───────────────────────────────────────────────────────
    scored = [(r["activity"], r["label"]) for r in results if r["activity"] is not None]
    scored.sort(reverse=True)

    # Load deposited anchor scores for context
    anchor_scores_csv = PIPELINE / "results" / "anchor_scores.csv"
    anchors: list[tuple[float, str]] = []
    if anchor_scores_csv.exists():
        import csv as _csv
        with open(anchor_scores_csv) as f:
            reader = _csv.DictReader(f)
            for row in reader:
                try:
                    anchors.append((float(row["activity"]), row["pdb_id"]))
                except (KeyError, ValueError):
                    pass
        anchors.sort(reverse=True)

    total_elapsed = (time.time() - total_t0) / 60
    print(f"\n{'='*66}", flush=True)
    print(f"  Ensemble results  ({total_elapsed:.1f} min total)", flush=True)
    print("="*66, flush=True)

    # Merge generated + anchors into one ranked list
    all_rows: list[tuple[float, str, str]] = []  # (activity, label, source)
    for act, lbl in scored:
        all_rows.append((act, lbl, "generated"))
    for act, pdb_id in anchors:
        all_rows.append((act, pdb_id, "deposited"))
    all_rows.sort(key=lambda x: -x[0])

    max_lbl = max((len(lbl) for _, lbl, _ in all_rows), default=30)
    for act, lbl, src in all_rows:
        marker = ">>>" if src == "generated" else "   "
        print(f"  {marker}  {lbl:{max_lbl}}  {act:+.3f}  [{src}]", flush=True)

    # ── warn about failures ─────────────────────────────────────────────────
    failed = [r for r in results if r["status"] != "ok"]
    if failed:
        print(f"\n  {len(failed)} topology/ies did not complete:", flush=True)
        for r in failed:
            print(f"    {r['tag']}  status={r['status']}", flush=True)


if __name__ == "__main__":
    main()
