"""Hugging Face Hub as the free, zero-ops result store for dwell-time chunks (#34).

No server, DNS, TLS, or always-on host: a public HF *Dataset* repo **is** the store.
The maintainer authors an experiment locally (``run_chunks create`` -> local
``chunk_store``) and PUBLISHES it to the repo. Contributors anywhere PULL a runnable
chunk + its inputs from the repo, run the MD, and UPLOAD the outputs back under
``submissions/`` (optionally as a PR, so untrusted contributors never touch the
canonical data). The maintainer INGESTS submissions into the local store
(SHA-256-verified), which marks chunks done + unlocks downstream, then re-publishes.

Strict leasing is dropped in favour of the project's required >=2x redundancy
(duplicate work is wanted, not a bug); the agreement check belongs in ``ingest`` and
is the next build on top of this hash-verified v1.

Layout in the repo:
  manifest.json                                  the chunk DAG + artifact registry
  artifacts/<sha256>.<ext>                       content-addressed inputs (+ accepted outputs)
  submissions/<chunk_id>/<worker>/<file>         a contributor's raw uploaded outputs
  submissions/<chunk_id>/<worker>/meta.json      {chunk_id, worker, wall_s, outputs:{aid:{sha256,file}}}

Commands (auth via a cached ``hf auth login`` token or --token / $HF_TOKEN):
  publish EXP --repo R    upload manifest.json + present artifacts to the HF repo
  work    EXP --repo R    contributor: pull a runnable chunk, run it, upload its outputs
  ingest  EXP --repo R    maintainer: pull submitted outputs, verify, push into local store
  status  EXP             local chunk_store progress (same as run_chunks status)
"""
from __future__ import annotations

import argparse
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
    h = hashlib.sha256()
    h.update(Path(path).read_bytes())
    return h.hexdigest()


def _runnable(manifest: dict, ch: dict) -> bool:
    """Pure copy of chunk_store's runnable predicate, minus leasing — a chunk is
    fair game if it is pending and every artifact it consumes is present."""
    if ch["status"] != "pending":
        return False
    arts = manifest["artifacts"]
    return all(arts.get(aid, {}).get("present") for aid in ch["consumes"])


# ---------------------------------------------------------------------------
# publish — maintainer pushes the experiment (manifest + present artifacts)
# ---------------------------------------------------------------------------

def cmd_publish(args) -> None:
    from huggingface_hub import HfApi
    api = HfApi(token=args.token)
    api.create_repo(args.repo, repo_type=DATASET, exist_ok=True)
    edir = store.experiment_dir(args.exp_id)
    if not (edir / "manifest.json").exists():
        raise FileNotFoundError(f"no local experiment {args.exp_id} at {edir}")
    api.upload_folder(
        repo_id=args.repo, repo_type=DATASET, folder_path=str(edir),
        allow_patterns=["manifest.json", "artifacts/*"],
        commit_message=f"publish {args.exp_id}",
    )
    s = store.status_summary(args.exp_id)
    print(f"published {args.exp_id} -> {args.repo}: {s['counts']['done']}/{s['n_chunks']} done, "
          f"{s['counts']['runnable']} runnable for contributors", flush=True)


# ---------------------------------------------------------------------------
# work — contributor pulls a runnable chunk, runs it, uploads outputs
# ---------------------------------------------------------------------------

def cmd_work(args) -> None:
    from huggingface_hub import HfApi, hf_hub_download, CommitOperationAdd
    api = HfApi(token=args.token)
    worker = args.worker or f"hf-{int(time.time())}"
    attempted: set[str] = set()
    done = 0
    while args.n <= 0 or done < args.n:
        mpath = hf_hub_download(args.repo, "manifest.json", repo_type=DATASET,
                                token=args.token, force_download=True)
        manifest = json.loads(Path(mpath).read_text(encoding="utf-8"))
        ch = next((c for c in manifest["chunks"]
                   if c["id"] not in attempted and _runnable(manifest, c)), None)
        if ch is None:
            print("no runnable chunk to claim (all done, leased elsewhere, or attempted)", flush=True)
            break
        attempted.add(ch["id"])
        print(f"\n=== chunk {ch['id']} [{ch['kind']}] ===", flush=True)

        scratch = Path(tempfile.mkdtemp(prefix=f"hf_work_{ch['id']}_"))
        local: dict[str, Path] = {}
        for aid in ch["consumes"]:
            rec = manifest["artifacts"][aid]
            local[aid] = Path(hf_hub_download(args.repo, rec["path"], repo_type=DATASET,
                                              token=args.token))
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
        wall = time.time() - t0

        ops, meta_outputs = [], {}
        base = f"submissions/{ch['id']}/{worker}"
        for aid, path in outputs.items():
            fname = aid.replace("/", "__")
            ops.append(CommitOperationAdd(f"{base}/{fname}", str(path)))
            meta_outputs[aid] = {"sha256": _sha256(path), "file": fname}
        meta = {"chunk_id": ch["id"], "worker": worker, "wall_s": wall, "outputs": meta_outputs}
        ops.append(CommitOperationAdd(f"{base}/meta.json",
                                      json.dumps(meta, indent=2).encode("utf-8")))
        api.create_commit(args.repo, ops, repo_type=DATASET,
                          commit_message=f"result {ch['id']} by {worker}",
                          create_pr=args.pr, token=args.token)
        print(f"  uploaded {ch['id']} ({wall/60:.1f} min){' as PR' if args.pr else ''}", flush=True)
        done += 1

    print(f"\nworker {worker}: ran {done} chunk(s).", flush=True)


# ---------------------------------------------------------------------------
# ingest — maintainer pulls submissions, verifies, pushes into the local store
# ---------------------------------------------------------------------------

def _largest_agreeing(values: list[float], tol: float) -> list[int]:
    """Indices of the largest subset whose values all lie within ``tol`` of each
    other (max-min <= tol). A sliding window over the sorted values — that's the
    consensus cluster; anything outside it is an outlier to quarantine. Pure +
    deterministic so the accept/quarantine decision is unit-testable."""
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


def cmd_ingest(args) -> None:
    """Pull contributor submissions, verify, and push accepted results into the
    local store. Two gates, since a single stranger's trajectory cannot be trusted:
    (1) **integrity** — recomputed SHA-256 must match the worker's claim;
    (2) **>=K agreement** — for scorable (segment) chunks, score each submission's
    trajectory on the dwell observable (cross-platform MD is *not* bit-identical, so
    we compare the science, not bytes) and accept the consensus only if >=``min_agree``
    submissions agree within ``agree_tol``; outliers are quarantined. Prep chunks
    (build/equilibrate) have no trajectory, so they pass on integrity + quorum count."""
    from huggingface_hub import HfApi, hf_hub_download
    import dwell_time
    api = HfApi(token=args.token)
    manifest = store.load_manifest(args.exp_id)          # local = authoritative
    chunks_by_id = {c["id"]: c for c in manifest["chunks"]}
    done_ids = {cid for cid, c in chunks_by_id.items() if c["status"] == "done"}
    shape = manifest["spec"]["shape"]
    ref_aid = f"{run_chunks._pair_tag(shape, run_chunks.APO)}/core.pdb"

    files = api.list_repo_files(args.repo, repo_type=DATASET, token=args.token)
    metas = [f for f in files if f.startswith("submissions/") and f.endswith("/meta.json")]

    # Group integrity-valid submissions by chunk: {cid: [(worker, {aid: path}, wall_s)]}
    subs: dict[str, list] = {}
    rejected = 0
    for mp in metas:
        meta = json.loads(Path(hf_hub_download(args.repo, mp, repo_type=DATASET,
                                               token=args.token)).read_text(encoding="utf-8"))
        cid = meta["chunk_id"]
        if cid in done_ids:
            continue
        base = f"submissions/{cid}/{meta['worker']}"
        outputs, ok = {}, True
        for aid, info in meta["outputs"].items():
            fp = Path(hf_hub_download(args.repo, f"{base}/{info['file']}", repo_type=DATASET,
                                      token=args.token))
            if _sha256(fp) != info["sha256"]:
                print(f"  REJECT {cid} by {meta['worker']}: sha256 mismatch on {aid}", flush=True)
                ok = False
                break
            outputs[aid] = fp
        if ok:
            subs.setdefault(cid, []).append((meta["worker"], outputs, float(meta.get("wall_s", 0.0))))

    accepted = waiting = quarantined = 0
    reference = None
    for cid, sublist in subs.items():
        ch = chunks_by_id[cid]
        if len(sublist) < args.min_agree:
            print(f"  awaiting quorum {cid}: {len(sublist)}/{args.min_agree} submissions", flush=True)
            waiting += 1
            continue
        if ch["kind"] == "segment":
            if reference is None:
                reference = dwell_time.load_pdb(store.artifact_file(args.exp_id, manifest, ref_aid))
            dwell = []
            for _worker, outs, _wall in sublist:
                seg = next(p for aid, p in outs.items() if aid.endswith(".pdb"))
                dwell.append(dwell_time.score_trajectory(seg, reference)["dwell_fraction"])
            cluster = _largest_agreeing(dwell, args.agree_tol)
            if len(cluster) < args.min_agree:
                print(f"  QUARANTINE {cid}: no {args.min_agree}-way agreement within "
                      f"{args.agree_tol} (dwell {[round(x, 3) for x in dwell]})", flush=True)
                quarantined += 1
                continue
            mean = sum(dwell[i] for i in cluster) / len(cluster)
            rep = min(cluster, key=lambda i: abs(dwell[i] - mean))
            worker, outs, wall = sublist[rep]
            store.push(args.exp_id, cid, outs, wall)
            print(f"  accepted {cid}: {len(cluster)}/{len(sublist)} agree "
                  f"(dwell {[round(x, 3) for x in dwell]}); consensus from {worker}", flush=True)
        else:
            worker, outs, wall = sublist[0]
            store.push(args.exp_id, cid, outs, wall)
            print(f"  accepted {cid} ({ch['kind']}, {len(sublist)} submission(s), integrity ok; "
                  f"from {worker})", flush=True)
        done_ids.add(cid)
        accepted += 1

    s = store.status_summary(args.exp_id)
    print(f"\ningest: +{accepted} accepted, {waiting} awaiting-quorum, {quarantined} quarantined, "
          f"{rejected} integrity-rejected. now {s['counts']['done']}/{s['n_chunks']} done.", flush=True)
    print("  re-run `publish` to push the newly-unlocked chunks' inputs to contributors.", flush=True)


def cmd_status(args) -> None:
    run_chunks.cmd_status(args)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_common(sp):
        sp.add_argument("exp_id")
        sp.add_argument("--repo", required=True, help="HF dataset repo, e.g. user/asyn-dwell-results")
        sp.add_argument("--token", default=None, help="HF token (default: cached `hf auth login` / $HF_TOKEN)")

    pp = sub.add_parser("publish", help="push manifest + present artifacts to the HF repo")
    add_common(pp)
    pp.set_defaults(func=cmd_publish)

    pw = sub.add_parser("work", help="contributor: pull/run/upload one runnable chunk")
    add_common(pw)
    pw.add_argument("--n", type=int, default=1, help="max chunks to run (<=0 = drain runnable)")
    pw.add_argument("--worker", default=None)
    pw.add_argument("--pr", action="store_true", help="upload as a pull request (untrusted contributors)")
    pw.add_argument("--stop-on-fail", action="store_true")
    pw.set_defaults(func=cmd_work)

    pi = sub.add_parser("ingest", help="maintainer: verify + push submitted outputs into the local store")
    add_common(pi)
    pi.add_argument("--min-agree", type=int, default=1,
                    help="min independent submissions that must agree before accepting (>=2 = redundancy)")
    pi.add_argument("--agree-tol", type=float, default=0.2,
                    help="max dwell-fraction spread within the consensus cluster (research param)")
    pi.set_defaults(func=cmd_ingest)

    pst = sub.add_parser("status", help="local chunk_store progress")
    pst.add_argument("exp_id")
    pst.set_defaults(func=cmd_status)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
