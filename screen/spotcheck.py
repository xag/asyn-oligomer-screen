"""Coordinator spot-check: independently re-run a sample of candidate chunks and
emit ``{chunk_id: dwell}`` for ``hf_store ingest --spotcheck-file`` (#43).

The ground-truth anchor of the acceptance gate. Without it, acceptance rests on
contributors agreeing *with each other*; a colluding cluster of fabrications
that happen to agree would pass. Spot-checking makes a random sample also have
to agree with a run the coordinator did itself (the coordinator holds the full,
trusted environment) — and a chunk whose contributor cluster excludes the
coordinator's own value is quarantined and its contributors flagged.

Only ``segment`` chunks carry the scored observable (dwell fraction), so only
they are spot-checked. Re-running reuses ``run_chunks.execute_chunk`` exactly as
a contributor would, pulling the chunk's inputs from the dataset's accepted
artifacts; the seed is baked into the chunk, so an honest re-run lands in the
same dwell cluster (within float non-determinism).

    python screen/spotcheck.py --repo user/asyn-dwell-results --sample 3 --out spotcheck.json
    python screen/hf_store.py ingest <exp> --repo user/asyn-dwell-results --spotcheck-file spotcheck.json
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import tempfile
from pathlib import Path

import hf_store as H
import run_chunks


def select_targets(candidates: list[str], sample_n: int, rng: random.Random) -> list[str]:
    """Pure: choose up to ``sample_n`` chunk ids to re-run (all of them when
    ``sample_n`` <= 0 or exceeds the pool). Sorted for deterministic output."""
    cands = sorted(set(candidates))
    if sample_n <= 0 or sample_n >= len(cands):
        return cands
    return sorted(rng.sample(cands, sample_n))


def _candidates(manifest: dict, submitted: set[str]) -> list[str]:
    """Segment chunks that aren't done yet, have at least one submission, and
    whose inputs are all present — i.e. exactly what's about to face quorum."""
    arts = manifest["artifacts"]
    out = []
    for c in manifest["chunks"]:
        if (c["kind"] == "segment" and c["status"] != "done" and c["id"] in submitted
                and all(arts.get(a, {}).get("present") for a in c["consumes"])):
            out.append(c["id"])
    return out


def cmd_run(args) -> None:
    from huggingface_hub import HfApi, hf_hub_download
    import dwell_time
    api = HfApi(token=args.token)
    manifest = H._download_manifest(api, args.repo, args.token)
    chunks_by_id = {c["id"]: c for c in manifest["chunks"]}

    files = api.list_repo_files(args.repo, repo_type=H.DATASET, token=args.token)
    submitted = {f.split("/")[1] for f in files
                 if f.startswith("submissions/") and f.endswith("/meta.json")}

    targets = select_targets(_candidates(manifest, submitted), args.sample, random.Random(args.seed))
    if not targets:
        Path(args.out).write_text("{}", encoding="utf-8")
        print("no spot-check candidates", flush=True)
        return

    shape = manifest["spec"]["shape"]
    ref_aid = f"{run_chunks._pair_tag(shape, run_chunks.APO)}/core.pdb"
    reference = dwell_time.load_pdb(hf_hub_download(
        args.repo, manifest["artifacts"][ref_aid]["path"], repo_type=H.DATASET, token=args.token))

    out: dict[str, float] = {}
    for cid in targets:
        ch = chunks_by_id[cid]
        scratch = Path(tempfile.mkdtemp(prefix=f"spot_{cid}_"))
        local = {a: Path(hf_hub_download(args.repo, manifest["artifacts"][a]["path"],
                                         repo_type=H.DATASET, token=args.token))
                 for a in ch["consumes"]}
        produced = run_chunks.execute_chunk(ch, lambda a: local[a], scratch)
        pdb = next(p for a, p in produced.items() if a.endswith(".pdb"))
        out[cid] = dwell_time.score_trajectory(pdb, reference)["dwell_fraction"]
        print(f"  spot-checked {cid}: dwell={out[cid]:.3f}", flush=True)

    Path(args.out).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {len(out)} spot-check value(s) -> {args.out}", flush=True)


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repo", required=True, help="HF dataset repo")
    ap.add_argument("--token", default=None, help="HF token (default: cached login / $HF_TOKEN)")
    ap.add_argument("--sample", type=int, default=3, help="how many candidate chunks to re-run (<=0 = all)")
    ap.add_argument("--seed", type=int, default=None, help="RNG seed for reproducible sampling")
    ap.add_argument("--out", default="spotcheck.json", help="output JSON {chunk_id: dwell}")
    ap.set_defaults(func=cmd_run)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
