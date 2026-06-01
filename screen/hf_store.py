"""Hugging Face Hub as the free, zero-ops result store for dwell-time chunks (#34).

The public HF *Dataset* repo **is** the authoritative store — no server, DNS, or
always-on host. Every contribution arrives the same way, with no privileged path:

  publish      seed the repo from a locally-authored experiment (manifest reset to
               pending + the initial inputs only). One-time, by the coordinator.
  work         contributor (anywhere): pull a runnable chunk + its inputs, run the MD,
               upload the outputs under submissions/ (optionally as a PR). The local
               machine runs this too — it is a contributor like any other.
  submit-local  backfill: submit already-computed local chunk outputs as contributions
               (no recompute), e.g. after a run that predates the HF store.
  ingest       coordinator: pull submissions, verify, and write *accepted* results into
               the repo's authoritative artifacts/ + manifest. Two gates — SHA-256
               integrity and >=K observable agreement (dwell-fraction consensus, since
               cross-platform MD is not bit-identical); outliers quarantined.
  status       progress from the repo's manifest.

Layout: manifest.json | artifacts/<sha256>.<ext> (accepted) | submissions/<chunk>/<worker>/{<file>, meta.json}.
Auth via a cached `hf auth login` token (or --token / $HF_TOKEN).
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
import tempfile
import time
from pathlib import Path

import chunk_store as store
import run_chunks

DATASET = "dataset"


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _runnable(manifest: dict, ch: dict) -> bool:
    if ch["status"] != "pending":
        return False
    arts = manifest["artifacts"]
    return all(arts.get(aid, {}).get("present") for aid in ch["consumes"])


def _largest_agreeing(values: list[float], tol: float) -> list[int]:
    """Indices of the largest subset whose values all lie within ``tol`` (max-min<=tol)
    — the consensus cluster; anything outside is an outlier. Pure + unit-testable."""
    if not values:
        return []
    order = sorted(range(len(values)), key=lambda i: values[i])
    best: list[int] = []
    lo = 0
    for hi in range(len(order)):
        while values[order[hi]] - values[order[lo]] > tol:
            lo += 1
        if hi - lo + 1 > len(best):
            best = order[lo:hi + 1]
    return sorted(best)


def _download_manifest(api, repo: str, token) -> dict:
    from huggingface_hub import hf_hub_download
    p = hf_hub_download(repo, "manifest.json", repo_type=DATASET, token=token, force_download=True)
    return json.loads(Path(p).read_text(encoding="utf-8"))


def _submit(api, repo: str, chunk_id: str, worker: str, outputs: dict, wall: float,
            token, pr: bool) -> None:
    """Upload one chunk's outputs as a contribution under submissions/ (uniform path —
    used by both the live `work` loop and the `submit-local` backfill)."""
    from huggingface_hub import CommitOperationAdd
    base = f"submissions/{chunk_id}/{worker}"
    ops, meta_outputs = [], {}
    for aid, path in outputs.items():
        fname = aid.replace("/", "__")
        ops.append(CommitOperationAdd(f"{base}/{fname}", str(path)))
        meta_outputs[aid] = {"sha256": _sha256(path), "file": fname}
    meta = {"chunk_id": chunk_id, "worker": worker, "wall_s": wall, "outputs": meta_outputs}
    ops.append(CommitOperationAdd(f"{base}/meta.json", json.dumps(meta, indent=2).encode("utf-8")))
    api.create_commit(repo, ops, repo_type=DATASET, create_pr=pr, token=token,
                      commit_message=f"submit {chunk_id} by {worker}")


# ---------------------------------------------------------------------------
# publish — seed the authoritative manifest + initial inputs (coordinator)
# ---------------------------------------------------------------------------

def cmd_publish(args) -> None:
    from huggingface_hub import HfApi, CommitOperationAdd
    api = HfApi(token=args.token)
    api.create_repo(args.repo, repo_type=DATASET, exist_ok=True)
    local = store.load_manifest(args.exp_id)
    produced = {aid for c in local["chunks"] for aid in c["produces"]}
    seed = copy.deepcopy(local)
    for c in seed["chunks"]:
        c["status"], c["lease"], c["error"] = "pending", None, None
        c.pop("wall_s", None)
    initial = {}
    for aid, rec in seed["artifacts"].items():
        if aid in produced:
            seed["artifacts"][aid] = {"present": False, "path": None, "sha256": None, "bytes": None}
        else:
            initial[aid] = rec
    edir = store.experiment_dir(args.exp_id)
    ops = [CommitOperationAdd("manifest.json", json.dumps(seed, indent=2).encode("utf-8"))]
    ops += [CommitOperationAdd(rec["path"], str(edir / rec["path"])) for rec in initial.values()]
    api.create_commit(args.repo, ops, repo_type=DATASET, token=args.token,
                      commit_message=f"publish (seed) {args.exp_id}")
    print(f"published seed {args.exp_id} -> {args.repo}: {len(seed['chunks'])} chunks pending, "
          f"{len(initial)} initial input(s)", flush=True)


# ---------------------------------------------------------------------------
# work — contributor: pull a runnable chunk, run it, submit outputs
# ---------------------------------------------------------------------------

def cmd_work(args) -> None:
    from huggingface_hub import HfApi, hf_hub_download
    api = HfApi(token=args.token)
    worker = args.worker or f"hf-{int(time.time())}"
    attempted: set[str] = set()
    done = 0
    while args.n <= 0 or done < args.n:
        manifest = _download_manifest(api, args.repo, args.token)
        ch = next((c for c in manifest["chunks"]
                   if c["id"] not in attempted and _runnable(manifest, c)), None)
        if ch is None:
            print("no runnable chunk to claim (all done/awaiting-ingest, or attempted)", flush=True)
            break
        attempted.add(ch["id"])
        print(f"\n=== chunk {ch['id']} [{ch['kind']}] ===", flush=True)
        scratch = Path(tempfile.mkdtemp(prefix=f"hf_work_{ch['id']}_"))
        local: dict[str, Path] = {}
        for aid in ch["consumes"]:
            rec = manifest["artifacts"][aid]
            local[aid] = Path(hf_hub_download(args.repo, rec["path"], repo_type=DATASET, token=args.token))
        infile = lambda aid: local[aid]  # noqa: E731
        t0 = time.time()
        try:
            outputs = run_chunks.execute_chunk(ch, infile, scratch)
        except Exception as e:  # noqa: BLE001
            print(f"  FAILED: {e}", flush=True)
            done += 1
            if args.stop_on_fail:
                raise
            continue
        _submit(api, args.repo, ch["id"], worker, outputs, time.time() - t0, args.token, args.pr)
        print(f"  submitted {ch['id']}{' as PR' if args.pr else ''}", flush=True)
        done += 1
    print(f"\nworker {worker}: ran {done} chunk(s). run `ingest` to accept them.", flush=True)


# ---------------------------------------------------------------------------
# submit-local — backfill already-computed local chunks as contributions
# ---------------------------------------------------------------------------

def cmd_submit_local(args) -> None:
    from huggingface_hub import HfApi
    api = HfApi(token=args.token)
    manifest = store.load_manifest(args.exp_id)
    worker = args.worker or "local-laptop"
    n = 0
    for c in manifest["chunks"]:
        if c["status"] != "done":
            continue
        outputs = {aid: store.artifact_file(args.exp_id, manifest, aid) for aid in c["produces"]}
        _submit(api, args.repo, c["id"], worker, outputs, float(c.get("wall_s", 0.0)),
                args.token, args.pr)
        n += 1
        print(f"  submitted {c['id']}", flush=True)
    print(f"\nsubmitted {n} locally-done chunk(s) as contributions from {worker}. "
          f"run `ingest` to accept them.", flush=True)


# ---------------------------------------------------------------------------
# ingest — coordinator: verify + write accepted results into the repo
# ---------------------------------------------------------------------------

def cmd_ingest(args) -> None:
    from huggingface_hub import HfApi, hf_hub_download, CommitOperationAdd, CommitOperationDelete
    import dwell_time
    api = HfApi(token=args.token)
    manifest = _download_manifest(api, args.repo, args.token)   # repo = authoritative
    chunks_by_id = {c["id"]: c for c in manifest["chunks"]}
    done_ids = {cid for cid, c in chunks_by_id.items() if c["status"] == "done"}
    stale: set[str] = set()   # submissions for already-done chunks -> delete, don't reprocess
    shape = manifest["spec"]["shape"]
    ref_aid = f"{run_chunks._pair_tag(shape, run_chunks.APO)}/core.pdb"

    files = api.list_repo_files(args.repo, repo_type=DATASET, token=args.token)
    metas = [f for f in files if f.startswith("submissions/") and f.endswith("/meta.json")]

    subs: dict[str, list] = {}
    rejected = 0
    for mp in metas:
        meta = json.loads(Path(hf_hub_download(args.repo, mp, repo_type=DATASET,
                                               token=args.token)).read_text(encoding="utf-8"))
        cid = meta["chunk_id"]
        if cid in done_ids:
            stale.add(cid)
            continue
        base = f"submissions/{cid}/{meta['worker']}"
        outs, ok = {}, True
        for aid, info in meta["outputs"].items():
            fp = Path(hf_hub_download(args.repo, f"{base}/{info['file']}", repo_type=DATASET,
                                      token=args.token))
            if _sha256(fp) != info["sha256"]:
                print(f"  REJECT {cid} by {meta['worker']}: sha256 mismatch on {aid}", flush=True)
                ok = False
                break
            outs[aid] = fp
        if ok:
            subs.setdefault(cid, []).append((meta["worker"], outs, float(meta.get("wall_s", 0.0))))

    all_ops, accepted, waiting, quarantined = [], 0, 0, 0
    accepted_cids: set[str] = set()
    reference = None
    for cid, sublist in subs.items():
        ch = chunks_by_id[cid]
        if len(sublist) < args.min_agree:
            print(f"  awaiting quorum {cid}: {len(sublist)}/{args.min_agree}", flush=True)
            waiting += 1
            continue
        if ch["kind"] == "segment":
            if reference is None:
                ref_rec = manifest["artifacts"][ref_aid]
                reference = dwell_time.load_pdb(hf_hub_download(args.repo, ref_rec["path"],
                                                               repo_type=DATASET, token=args.token))
            dwell = [dwell_time.score_trajectory(
                next(p for aid, p in outs.items() if aid.endswith(".pdb")), reference)["dwell_fraction"]
                for _w, outs, _wl in sublist]
            cluster = _largest_agreeing(dwell, args.agree_tol)
            if len(cluster) < args.min_agree:
                print(f"  QUARANTINE {cid}: no {args.min_agree}-way agreement within "
                      f"{args.agree_tol} (dwell {[round(x, 3) for x in dwell]})", flush=True)
                quarantined += 1
                continue
            mean = sum(dwell[i] for i in cluster) / len(cluster)
            rep = min(cluster, key=lambda i: abs(dwell[i] - mean))
            note = f"{len(cluster)}/{len(sublist)} agree (dwell {[round(x, 3) for x in dwell]})"
        else:
            rep, note = 0, f"{ch['kind']}, {len(sublist)} submission(s), integrity ok"
        worker, outs, wall = sublist[rep]
        for aid, path in outs.items():
            sha, ext = _sha256(path), (Path(aid).suffix or ".bin")
            dest = f"artifacts/{sha}{ext}"
            all_ops.append(CommitOperationAdd(dest, str(path)))
            manifest["artifacts"][aid] = {"present": True, "path": dest, "sha256": sha,
                                          "bytes": Path(path).stat().st_size}
        ch["status"], ch["wall_s"], ch["lease"], ch["error"] = "done", wall, None, None
        done_ids.add(cid)
        accepted_cids.add(cid)
        accepted += 1
        print(f"  accepted {cid}: {note}; from {worker}", flush=True)

    # Clean up: drop submissions now folded into artifacts/ (accepted here) or already
    # done (stale) — keeps the inbox lean so the scheduled ingester doesn't reprocess.
    for cid in accepted_cids | stale:
        all_ops.append(CommitOperationDelete(f"submissions/{cid}", is_folder=True))
    if accepted:
        all_ops.append(CommitOperationAdd("manifest.json", json.dumps(manifest, indent=2).encode("utf-8")))
    if all_ops:
        api.create_commit(args.repo, all_ops, repo_type=DATASET, token=args.token,
                          commit_message=f"ingest +{accepted} accepted, {len(accepted_cids | stale)} cleaned")
    n_done = sum(1 for c in manifest["chunks"] if c["status"] == "done")
    print(f"\ningest: +{accepted} accepted, {waiting} awaiting-quorum, {quarantined} quarantined, "
          f"{rejected} integrity-rejected, {len(stale)} stale-cleaned. "
          f"now {n_done}/{len(manifest['chunks'])} done in the repo.", flush=True)


def cmd_status(args) -> None:
    from huggingface_hub import HfApi
    manifest = _download_manifest(HfApi(token=args.token), args.repo, args.token)
    by_kind: dict[str, list] = {}
    for c in manifest["chunks"]:
        k = by_kind.setdefault(c["kind"], [0, 0])
        k[0] += int(c["status"] == "done")
        k[1] += 1
    n_done = sum(k[0] for k in by_kind.values())
    n_tot = sum(k[1] for k in by_kind.values())
    print(f"{args.repo} [{args.exp_id}]: {n_done}/{n_tot} done", flush=True)
    for kind, (d, t) in by_kind.items():
        print(f"  {kind:12} {d}/{t}", flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("exp_id")
        sp.add_argument("--repo", required=True, help="HF dataset repo, e.g. user/asyn-dwell-results")
        sp.add_argument("--token", default=None, help="HF token (default: cached `hf auth login` / $HF_TOKEN)")

    pp = sub.add_parser("publish", help="seed the repo from a local experiment (manifest pending + inputs)")
    common(pp)
    pp.set_defaults(func=cmd_publish)

    pw = sub.add_parser("work", help="contributor: pull/run/submit runnable chunks")
    common(pw)
    pw.add_argument("--n", type=int, default=1, help="max chunks to run (<=0 = drain currently-runnable)")
    pw.add_argument("--worker", default=None)
    pw.add_argument("--pr", action="store_true", help="submit as a pull request (untrusted contributors)")
    pw.add_argument("--stop-on-fail", action="store_true")
    pw.set_defaults(func=cmd_work)

    psl = sub.add_parser("submit-local", help="backfill already-computed local chunks as contributions")
    common(psl)
    psl.add_argument("--worker", default=None)
    psl.add_argument("--pr", action="store_true")
    psl.set_defaults(func=cmd_submit_local)

    pi = sub.add_parser("ingest", help="coordinator: verify + write accepted results into the repo")
    common(pi)
    pi.add_argument("--min-agree", type=int, default=1,
                    help="min independent submissions that must agree before accepting (>=2 = redundancy)")
    pi.add_argument("--agree-tol", type=float, default=0.2,
                    help="max dwell-fraction spread within the consensus cluster (research param)")
    pi.set_defaults(func=cmd_ingest)

    pst = sub.add_parser("status", help="progress from the repo manifest")
    common(pst)
    pst.set_defaults(func=cmd_status)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
