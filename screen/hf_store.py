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
               the repo's authoritative artifacts/ + manifest. Gates (see contrib_gate):
               SHA-256 integrity, then a distinct-pseudonym, reputation-weighted
               dwell-fraction consensus (cross-platform MD is not bit-identical), with
               an optional coordinator spot-check as ground truth. Per-contributor
               outcomes are written to outcomes/ for health to fold into reputation.
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
import contrib_gate
import run_chunks
from contrib_gate import Reputation, Submission

DATASET = "dataset"


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _runnable(manifest: dict, ch: dict) -> bool:
    if ch["status"] != "pending":
        return False
    arts = manifest["artifacts"]
    return all(arts.get(aid, {}).get("present") for aid in ch["consumes"])


def _submission_key_sha(outputs_meta: dict) -> str:
    """The SHA-256 that identifies a submission for de-dup + archiving: the
    scored trajectory (``.pdb``) if present, else the first output. Pure."""
    for aid, info in outputs_meta.items():
        if aid.endswith(".pdb"):
            return info["sha256"]
    return next(iter(outputs_meta.values()))["sha256"]


def _parse_reputations(obj: dict) -> dict:
    """Parse health's published ``reputations.json`` into ``{pseudonym:
    Reputation}``. Pure + unit-testable; missing fields default to fresh."""
    out: dict[str, Reputation] = {}
    for pseudo, rec in (obj or {}).items():
        rec = rec or {}
        out[pseudo] = Reputation(
            agreed=int(rec.get("agreed", 0)),
            outlier=int(rec.get("outlier", 0)),
            spot_pass=int(rec.get("spot_pass", 0)),
            spot_fail=int(rec.get("spot_fail", 0)),
            allowlist_bonus=float(rec.get("allowlist_bonus", 0.0)),
        )
    return out


def _load_reputations(api, repo: str, token) -> dict:
    """Download health's ``reputations.json`` from the repo if present, else {}.
    health publishes it; the dataset stays the single rendezvous point — no
    direct service-to-service call, so ingest runs offline-of-health."""
    from huggingface_hub import hf_hub_download
    try:
        p = hf_hub_download(repo, "reputations.json", repo_type=DATASET,
                            token=token, force_download=True)
    except Exception:  # noqa: BLE001 — absent file = everyone fresh
        return {}
    try:
        return _parse_reputations(json.loads(Path(p).read_text(encoding="utf-8")))
    except Exception:  # noqa: BLE001
        return {}


def _load_spotchecks(path) -> dict:
    """Optional ``{chunk_id: coordinator_dwell}`` the coordinator produced by
    re-running sampled chunks itself; absent → no spot-checks this run."""
    if not path:
        return {}
    return {k: float(v) for k, v in
            json.loads(Path(path).read_text(encoding="utf-8")).items()}


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
    # Guard the single manifest slot: re-publishing seeds the manifest to
    # all-pending and orphans every artifact a different live experiment had
    # accepted. A `demo`/test publish silently wiping a real run is exactly how
    # 48 h of accepted results became unreferenced. So refuse to overwrite a
    # manifest whose exp_id differs from this one unless --force is given.
    try:
        existing = _download_manifest(api, args.repo, args.token)
    except Exception:  # noqa: BLE001 — no manifest yet = first publish, nothing to clobber
        existing = None
    if existing and existing.get("exp_id") not in (None, args.exp_id) and not args.force:
        raise SystemExit(
            f"refusing to publish '{args.exp_id}': the repo already holds a different "
            f"live experiment '{existing.get('exp_id')}'. Publishing would reset its "
            f"manifest and orphan its accepted artifacts. Re-run with --force to override."
        )
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
# publish-molecules — the registry health's GET /molecules reads
# ---------------------------------------------------------------------------
# molecules.json is the candidate-molecule registry: the primary seed list
# (data/vicinity_molecules.js, emitted as JSON by scripts/build_molecules_json.mjs)
# plus contributed proposals the broker appends at runtime. Republishing the
# primary set MERGES — it preserves every source!="primary" entry already in the
# file — so refreshing the seed never clobbers contributors' proposals.

def cmd_publish_molecules(args) -> None:
    from huggingface_hub import HfApi, CommitOperationAdd, hf_hub_download
    api = HfApi(token=args.token)
    primary = json.loads(Path(args.file).read_text(encoding="utf-8"))
    if not isinstance(primary, list):
        raise SystemExit(f"{args.file} must be a JSON array of molecule records")
    contributed = []
    try:
        p = hf_hub_download(args.repo, "molecules.json", repo_type=DATASET, token=args.token)
        existing = json.loads(Path(p).read_text(encoding="utf-8"))
        contributed = [m for m in existing if isinstance(m, dict) and m.get("source") != "primary"]
    except Exception:  # noqa: BLE001 — no registry yet = first publish
        pass
    merged = list(primary) + contributed
    api.create_commit(
        args.repo,
        [CommitOperationAdd("molecules.json", json.dumps(merged, indent=2).encode("utf-8"))],
        repo_type=DATASET, token=args.token,
        commit_message="publish molecules registry",
    )
    print(f"published molecules.json -> {args.repo}: {len(primary)} primary "
          f"+ {len(contributed)} contributed", flush=True)


# ---------------------------------------------------------------------------
# enqueue-awaiting — turn contributed submissions into runnable work (coordinator)
# ---------------------------------------------------------------------------

def _norm_chunk(ch: dict) -> dict:
    """The pending-chunk shape chunk_store.create_experiment writes — mirrored here
    so appended chunks match the rest of the manifest exactly."""
    return {
        "id": ch["id"], "kind": ch["kind"],
        "consumes": list(ch.get("consumes", [])),
        "produces": list(ch.get("produces", [])),
        "params": dict(ch.get("params", {})),
        "meta": dict(ch.get("meta", {})),
        "status": "pending", "lease": None, "wall_s": None,
        "error": None, "updated_at": None,
    }


# A contributed submission earns a cheap SCREENING pass, not the full curated
# treatment: one pose, a couple of seeds, short trajectories — a few dozen chunks
# instead of thousands, so an unvetted molecule can't blow up the work list. It's
# also the right science: screen cheaply, and promote a promising hit to the full
# budget as a moderation step (re-create it as a primary molecule).
#
# INVARIANT: these reduce only the *count* of chunks (poses, seeds, total prod_ps).
# They never touch equil_ps / segment_ps — those are inherited from the curated
# spec — so a contributed chunk's wall time always equals (or is less than) one of
# ours, never more. Do not add a SCREEN_SEGMENT_PS/SCREEN_EQUIL_PS here.
SCREEN_N_POSES = 1
SCREEN_N_SEEDS = 2
SCREEN_PROD_PS = 1000.0


def cmd_enqueue_awaiting(args) -> None:
    """Enqueue contributed molecules still `prep:awaiting` into the live experiment:
    generate each one's dock->build->equilibrate->segment chunks (docking stays a
    contributor `dock` chunk) against the manifest's spec + shared apo core, on a
    reduced screening budget, and append them. Marks the molecule `prep:ready`.
    Idempotent — molecules already in the spec, or chunk ids already present, are
    skipped. Adds nothing to our prioritisation: the new chunks are contributed-tier,
    run only via the opt-out."""
    from huggingface_hub import HfApi, CommitOperationAdd, hf_hub_download
    from run_chunks import enumerate_chunks, APO

    api = HfApi(token=args.token)
    manifest = _download_manifest(api, args.repo, args.token)
    if manifest.get("exp_id") not in (None, args.exp_id):
        raise SystemExit(f"repo holds experiment {manifest.get('exp_id')!r}, not {args.exp_id!r}")

    try:
        p = hf_hub_download(args.repo, "molecules.json", repo_type=DATASET,
                            token=args.token, force_download=True)
        registry = json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001 — no registry yet = nothing to enqueue
        registry = []

    spec = manifest["spec"]
    have = set(spec.get("ligands", []))
    artifacts = manifest["artifacts"]
    seen_ids = {c["id"] for c in manifest["chunks"]}

    added_mols, added_chunks = [], 0
    for m in registry:
        if not isinstance(m, dict) or m.get("source") != "contributed" or m.get("prep") != "awaiting":
            continue
        mid, smiles = m.get("id"), m.get("smiles")
        if not mid or not smiles or mid == APO or mid in have:
            continue
        lig_spec = {
            **spec, "ligands": [mid], "ligand_smiles": {mid: smiles},
            "n_poses": SCREEN_N_POSES,
            "seeds": list(spec.get("seeds", []))[:SCREEN_N_SEEDS] or [1],
            "prod_ps": min(float(spec.get("prod_ps", SCREEN_PROD_PS)), SCREEN_PROD_PS),
        }
        # Enforce the invariant: a contributed chunk must never run longer than ours.
        assert lig_spec.get("segment_ps") == spec.get("segment_ps") \
            and lig_spec.get("equil_ps") == spec.get("equil_ps"), \
            "screening budget must not change per-chunk MD length (equil_ps/segment_ps)"
        try:
            chunks = enumerate_chunks(lig_spec)
        except Exception as e:  # noqa: BLE001 — one bad submission shouldn't block the rest
            print(f"  skip {mid}: {e}", flush=True)
            continue
        for ch in chunks:
            if ch["id"] in seen_ids:
                continue
            manifest["chunks"].append(_norm_chunk(ch))
            seen_ids.add(ch["id"])
            for aid in ch.get("produces", []):
                artifacts.setdefault(aid, {"present": False, "path": None, "sha256": None, "bytes": None})
            added_chunks += 1
        spec.setdefault("ligands", []).append(mid)
        spec.setdefault("ligand_smiles", {})[mid] = smiles
        m["prep"] = "ready"
        added_mols.append(mid)

    if not added_mols:
        print("nothing to enqueue (no contributed molecules awaiting prep)", flush=True)
        return

    api.create_commit(
        args.repo,
        [CommitOperationAdd("manifest.json", json.dumps(manifest, indent=2).encode("utf-8")),
         CommitOperationAdd("molecules.json", json.dumps(registry, indent=2).encode("utf-8"))],
        repo_type=DATASET, token=args.token,
        commit_message=f"enqueue {len(added_mols)} contributed molecule(s): {', '.join(added_mols)}",
    )
    print(f"enqueued {len(added_mols)} molecule(s), +{added_chunks} chunks: "
          f"{', '.join(added_mols)}", flush=True)


# ---------------------------------------------------------------------------
# seed-dag — seed the whole experiment DAG (dock/build/replica) into the store
# ---------------------------------------------------------------------------

def cmd_seed_dag(args) -> None:
    """Seed the whole experiment DAG into the work store: the MD spec, the initial
    artifacts (the apo core that dock/build consume), and every job — dock, build and
    replica cursor — each wired by needs/produces. The store gates dispatch on
    dependencies and unblocks downstream work as artifacts register, so the entire
    pipeline flows through one dispatch path (no manifest dispatch, no build→replica
    bridge). Idempotent — re-seeding leaves committed progress untouched. POSTs to
    health /seed (cron-authenticated)."""
    import urllib.request
    from huggingface_hub import HfApi, hf_hub_download

    api = HfApi(token=args.token)
    manifest = _download_manifest(api, args.repo, args.token)
    spec = manifest.get("spec", {})
    arts = manifest.get("artifacts", {})
    chunks = manifest.get("chunks", [])
    prod_ps = spec.get("prod_ps")

    try:
        p = hf_hub_download(args.repo, "molecules.json", repo_type=DATASET,
                            token=args.token, force_download=True)
        registry = json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        registry = []
    prio = {}
    for m in registry:
        if isinstance(m, dict) and m.get("id"):
            prio[m["id"]] = 10 ** 9 if m.get("source") == "contributed" else (m.get("priority") or 100)
    priority_of = lambda lig: 0 if lig == "apo" else prio.get(lig, 100)  # noqa: E731

    # Initial artifacts = those no chunk produces (the apo core every dock/build
    # consumes); register them with their dataset paths so dependents can unblock.
    produced = {aid for c in chunks for aid in c.get("produces", [])}
    artifacts = {aid: rec["path"] for aid, rec in arts.items()
                 if aid not in produced and rec and rec.get("path")}

    dock, build, replicas = [], [], []
    for ch in chunks:
        kind, meta = ch.get("kind"), ch.get("meta", {})
        lig = meta.get("ligand", "apo")
        if kind in ("dock", "build"):
            job = {
                "id": ch["id"], "ligand": lig,
                "needs": {"core.pdb": ch["consumes"][0]},
                "produces": ch["produces"], "params": ch.get("params", {}),
                "priority": priority_of(lig),
            }
            (dock if kind == "dock" else build).append(job)
        elif kind == "equilibrate":
            # One replica cursor per equilibrate; its segments are subsumed by the
            # cursor, so segment chunks are skipped. It needs the build's system.
            sys_xml, solv = ch["consumes"][0], ch["consumes"][1]
            replicas.append({
                "rid": ch["id"][len("equil__"):],   # equil__<prefix>__s<seed> → <prefix>__s<seed>
                "ligand": lig, "pose": meta.get("pose"), "seed": meta.get("seed"),
                "target_ps": prod_ps,
                "needs": {"system.xml": sys_xml, "solvated.pdb": solv},
                "priority": priority_of(lig),
            })

    body = json.dumps({
        "spec": {
            "equil_ps": spec.get("equil_ps"), "segment_ps": spec.get("segment_ps"),
            "temperature_k": spec.get("temperature_k"), "traj_interval_ps": spec.get("traj_interval_ps"),
            "checkpoint_s": args.checkpoint_s,
        },
        "artifacts": artifacts, "dock": dock, "build": build, "replicas": replicas,
    }).encode("utf-8")
    req = urllib.request.Request(
        args.site.rstrip("/") + "/api/screen/v1/seed", data=body,
        headers={"Content-Type": "application/json", "Authorization": f"Bearer {args.secret}"})
    with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
        print(f"seeded DAG: {len(dock)} dock, {len(build)} build, {len(replicas)} replica(s) — "
              f"{r.read().decode('utf-8')}", flush=True)


# ---------------------------------------------------------------------------
# project — mirror the work store's artifacts into manifest.json (scoring view)
# ---------------------------------------------------------------------------

def cmd_project(args) -> None:
    """Project the work store back into manifest.json: pull the registered artifacts
    from health /state and mark them present (with their dataset paths), and mark a
    chunk done once all its produced artifacts exist. The store stays the source of
    truth for dispatch + progress; the manifest is the generated view that scoring,
    spot-check and archival read. Cron-authenticated."""
    import urllib.request
    from huggingface_hub import HfApi, CommitOperationAdd

    api = HfApi(token=args.token)
    manifest = _download_manifest(api, args.repo, args.token)

    req = urllib.request.Request(
        args.site.rstrip("/") + "/api/screen/v1/state",
        headers={"Authorization": f"Bearer {args.secret}"})
    with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
        state = json.loads(r.read().decode("utf-8"))
    artifacts = state.get("artifacts", {}) or {}

    arts = manifest.setdefault("artifacts", {})
    for aid, path in artifacts.items():
        arts[aid] = {"present": True, "path": path, "sha256": None, "bytes": None}
    for ch in manifest.get("chunks", []):
        prod = ch.get("produces", [])
        if prod and all(arts.get(a, {}).get("present") for a in prod):
            ch["status"] = "done"

    api.create_commit(
        args.repo, [CommitOperationAdd("manifest.json", json.dumps(manifest, indent=2).encode("utf-8"))],
        repo_type=DATASET, token=args.token, commit_message="project work store -> manifest")
    print(f"projected {len(artifacts)} artifact(s) into manifest.json", flush=True)


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

    reputations = _load_reputations(api, args.repo, args.token)
    spotchecks = _load_spotchecks(args.spotcheck_file)

    files = api.list_repo_files(args.repo, repo_type=DATASET, token=args.token)
    metas = [f for f in files if f.startswith("submissions/") and f.endswith("/meta.json")]

    # Gather integrity-valid submissions per chunk. Each record keeps the full
    # output set (to archive on accept) plus the identifying SHA the gate de-dups
    # on and the contributor pseudonym (the `worker` id minted by the dispatcher).
    subs: dict[str, list[dict]] = {}
    rejected = 0
    for mp in metas:
        meta = json.loads(Path(hf_hub_download(args.repo, mp, repo_type=DATASET,
                                               token=args.token)).read_text(encoding="utf-8"))
        cid = meta["chunk_id"]
        if cid in done_ids:
            stale.add(cid)
            continue
        worker = meta["worker"]
        base = f"submissions/{cid}/{worker}"
        outs, ok = {}, True
        for aid, info in meta["outputs"].items():
            fp = Path(hf_hub_download(args.repo, f"{base}/{info['file']}", repo_type=DATASET,
                                      token=args.token))
            if _sha256(fp) != info["sha256"]:
                print(f"  REJECT {cid} by {worker}: sha256 mismatch on {aid}", flush=True)
                ok = False
                rejected += 1
                break
            outs[aid] = fp
        if ok:
            subs.setdefault(cid, []).append({
                "worker": worker, "outs": outs,
                "wall": float(meta.get("wall_s", 0.0)),
                "key_sha": _submission_key_sha(meta["outputs"]),
                "ts": float(meta.get("iat", 0.0)),
            })

    all_ops, accepted, waiting, quarantined = [], 0, 0, 0
    accepted_cids: set[str] = set()
    outcomes: list[dict] = []   # per-pseudonym, for health to fold into reputation
    reference = None
    now = time.time()
    for cid, recs in subs.items():
        ch = chunks_by_id[cid]
        observable = ch["kind"] == "segment"   # only segments carry a dwell value

        # Score the observable for segment chunks (the existing GPU-free scorer).
        if observable:
            if reference is None:
                ref_rec = manifest["artifacts"][ref_aid]
                reference = dwell_time.load_pdb(hf_hub_download(args.repo, ref_rec["path"],
                                                               repo_type=DATASET, token=args.token))
            for r in recs:
                pdb = next(p for aid, p in r["outs"].items() if aid.endswith(".pdb"))
                r["dwell"] = dwell_time.score_trajectory(pdb, reference)["dwell_fraction"]

        by_key = {r["key_sha"]: r for r in recs}
        sub_objs = [Submission(pseudonym=r["worker"], sha256=r["key_sha"],
                               dwell=r.get("dwell"), ts=r["ts"]) for r in recs]
        d = contrib_gate.decide(sub_objs, reputations, tol=args.agree_tol,
                                quorum_weight=args.quorum_weight,
                                coordinator_dwell=spotchecks.get(cid),
                                observable=observable)
        contributors = {r["worker"] for r in recs}
        cluster = set(d.cluster)

        if d.status == "accept":
            rep = by_key[d.representative]
            for aid, path in rep["outs"].items():
                sha, ext = _sha256(path), (Path(aid).suffix or ".bin")
                dest = f"artifacts/{sha}{ext}"
                all_ops.append(CommitOperationAdd(dest, str(path)))
                manifest["artifacts"][aid] = {"present": True, "path": dest, "sha256": sha,
                                              "bytes": Path(path).stat().st_size}
            ch["status"], ch["wall_s"], ch["lease"], ch["error"] = "done", rep["wall"], None, None
            done_ids.add(cid)
            accepted_cids.add(cid)
            accepted += 1
            good = "spot_pass" if cid in spotchecks else "agreed"
            for w in contributors:
                outcomes.append({"ts": now, "chunk_id": cid, "pseudonym": w,
                                 "outcome": good if w in cluster else "outlier"})
            print(f"  accepted {cid}: {d.note}", flush=True)
        elif d.status == "spotcheck_fail":
            quarantined += 1
            for w in contributors:
                outcomes.append({"ts": now, "chunk_id": cid, "pseudonym": w,
                                 "outcome": "spot_fail" if w in cluster else "outlier"})
            print(f"  QUARANTINE {cid}: {d.note}", flush=True)
        else:   # awaiting — undecided, emit nothing (never penalise an incomplete)
            waiting += 1
            print(f"  awaiting {cid}: {d.note}", flush=True)

    # Clean up: drop submissions now folded into artifacts/ (accepted here) or already
    # done (stale) — keeps the inbox lean so the scheduled ingester doesn't reprocess.
    for cid in accepted_cids | stale:
        all_ops.append(CommitOperationDelete(f"submissions/{cid}", is_folder=True))
    if accepted:
        all_ops.append(CommitOperationAdd("manifest.json", json.dumps(manifest, indent=2).encode("utf-8")))
    # Outcomes go to a fresh timestamped file (append-only by convention — HF
    # commits can't append), which health reads to update reputation.
    if outcomes:
        all_ops.append(CommitOperationAdd(
            f"outcomes/{int(now)}.jsonl",
            ("\n".join(json.dumps(o) for o in outcomes) + "\n").encode("utf-8")))
    if all_ops:
        api.create_commit(args.repo, all_ops, repo_type=DATASET, token=args.token,
                          commit_message=f"ingest +{accepted} accepted, {len(accepted_cids | stale)} cleaned")
    n_done = sum(1 for c in manifest["chunks"] if c["status"] == "done")
    print(f"\ningest: +{accepted} accepted, {waiting} awaiting, {quarantined} quarantined, "
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
    pp.add_argument("--force", action="store_true",
                    help="overwrite even if the repo holds a different live experiment "
                         "(orphans its accepted artifacts)")
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
    pi.add_argument("--quorum-weight", type=float, default=1.0,
                    help="summed reputation weight a chunk's agreeing cluster must reach "
                         "to be accepted (a fresh contributor = 0.1; trust-scaled successor "
                         "to the old --min-agree count)")
    pi.add_argument("--agree-tol", type=float, default=0.2,
                    help="max dwell-fraction spread within the consensus cluster (research param)")
    pi.add_argument("--spotcheck-file", default=None,
                    help="optional JSON {chunk_id: coordinator_dwell} from coordinator "
                         "re-runs; a chunk whose cluster excludes it is quarantined")
    pi.set_defaults(func=cmd_ingest)

    pst = sub.add_parser("status", help="progress from the repo manifest")
    common(pst)
    pst.set_defaults(func=cmd_status)

    pm = sub.add_parser("publish-molecules",
                        help="publish/refresh molecules.json (primary registry; preserves contributed)")
    pm.add_argument("--repo", required=True, help="HF dataset repo, e.g. user/asyn-dwell-results")
    pm.add_argument("--token", default=None, help="HF token (default: cached `hf auth login` / $HF_TOKEN)")
    pm.add_argument("--file", default="molecules.json",
                    help="primary molecules JSON, e.g. from `node scripts/build_molecules_json.mjs`")
    pm.set_defaults(func=cmd_publish_molecules)

    pe = sub.add_parser("enqueue-awaiting",
                        help="generate chunks for contributed molecules awaiting prep and append to the live manifest")
    common(pe)
    pe.set_defaults(func=cmd_enqueue_awaiting)

    psr = sub.add_parser("seed-dag",
                         help="seed the whole experiment DAG (dock/build/replica jobs) into the work store (POST health /seed)")
    common(psr)
    psr.add_argument("--site", required=True, help="health site root, e.g. https://<site>")
    psr.add_argument("--secret", required=True, help="CRON_SECRET shared with health (Bearer auth)")
    psr.add_argument("--checkpoint-s", type=float, default=10.0,
                     help="wall-seconds between contributor checkpoints stamped onto units")
    psr.set_defaults(func=cmd_seed_dag)

    ppj = sub.add_parser("project",
                         help="mirror work-store artifacts into manifest.json for scoring (GET health /state)")
    common(ppj)
    ppj.add_argument("--site", required=True, help="health site root, e.g. https://<site>")
    ppj.add_argument("--secret", required=True, help="CRON_SECRET shared with health (Bearer auth)")
    ppj.set_defaults(func=cmd_project)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
