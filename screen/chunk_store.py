"""File-based experiment store for distributable MD chunks (issue #34, local scope).

This is the local stand-in for the future crowdsource coordinator: it holds an
*experiment* — a published spec plus a DAG of self-contained work *chunks* — and
exposes the pull/lease/push API a network service would later wrap. Nothing here
is GPU- or MD-specific; the dwell-time semantics (how a spec becomes chunks, how a
chunk is executed) live in ``run_chunks.py``. Keeping that split means the same
store serves any future workload the project wants to crowdsource.

Model
-----
* An **artifact** is a content-addressed blob (a built system, a serialised State,
  a frame PDB, an initial structure). Each has a stable logical *id*; once produced
  its bytes are copied under ``artifacts/<sha256>.<ext>`` and marked ``present``.
* A **chunk** ``consumes`` some artifact ids and ``produces`` others. It is
  **runnable** iff it is ``pending``, unleased (or its lease expired), and every
  artifact it consumes is ``present``. That single predicate gives the DAG:
  ``build`` (consumes the initial structure) unlocks each independent ``equilibrate``
  (per seed), which unlocks that replica's ``segment`` chain. Replicas never depend
  on each other — the workload is embarrassingly parallel across seeds and pairs.
* ``pull`` atomically leases one runnable chunk; ``push`` ingests its outputs, marks
  it done, and thereby flips the next chunk runnable. Many workers can pull the same
  published experiment concurrently; leasing (with an expiry, so a dead worker's
  chunk is reclaimable) prevents double work.

Concurrency is guarded by an atomic per-experiment lock (``os.mkdir`` is atomic on
both Windows and POSIX); the JSON manifest is written via a temp file + ``os.replace``
so a reader never sees a half-written file.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import time
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXPERIMENTS_DIR = ROOT / "experiments"

DEFAULT_LEASE_SECONDS = 3600.0
_LOCK_STALE_SECONDS = 120.0
_LOCK_RETRY_SECONDS = 0.05


# ---------------------------------------------------------------------------
# Paths.
# ---------------------------------------------------------------------------

def experiment_dir(exp_id: str) -> Path:
    return EXPERIMENTS_DIR / exp_id


def _manifest_path(exp_id: str) -> Path:
    return experiment_dir(exp_id) / "manifest.json"


def _artifacts_dir(exp_id: str) -> Path:
    return experiment_dir(exp_id) / "artifacts"


# ---------------------------------------------------------------------------
# Atomic lock + manifest IO.
# ---------------------------------------------------------------------------

@contextmanager
def _locked(exp_id: str):
    """Per-experiment mutex via atomic directory creation. Breaks a lock older
    than _LOCK_STALE_SECONDS (a crashed worker)."""
    lock = experiment_dir(exp_id) / ".lock"
    while True:
        try:
            os.mkdir(lock)
            break
        except FileExistsError:
            try:
                age = time.time() - lock.stat().st_mtime
            except FileNotFoundError:
                continue
            if age > _LOCK_STALE_SECONDS:
                # Stale: steal it. Best-effort; a racing stealer just retries.
                try:
                    os.utime(lock, None)
                except FileNotFoundError:
                    pass
                break
            time.sleep(_LOCK_RETRY_SECONDS)
    try:
        yield
    finally:
        try:
            os.rmdir(lock)
        except FileNotFoundError:
            pass


def load_manifest(exp_id: str) -> dict:
    return json.loads(_manifest_path(exp_id).read_text(encoding="utf-8"))


def _save_manifest(exp_id: str, manifest: dict) -> None:
    path = _manifest_path(exp_id)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# Construction.
# ---------------------------------------------------------------------------

def create_experiment(
    exp_id: str,
    spec: dict,
    chunks: list[dict],
    initial_artifacts: dict[str, Path],
) -> dict:
    """Create + publish an experiment. ``chunks`` is the prebuilt DAG (each with
    ``id``/``kind``/``consumes``/``produces``/``params`` and any descriptive
    fields). ``initial_artifacts`` maps artifact-id → source file present from the
    start (e.g. the NAC-core PDB a build chunk consumes). Idempotent-safe: refuses
    to clobber an existing experiment."""
    edir = experiment_dir(exp_id)
    if _manifest_path(exp_id).exists():
        raise FileExistsError(f"experiment {exp_id!r} already exists at {edir}")
    _artifacts_dir(exp_id).mkdir(parents=True, exist_ok=True)

    artifacts: dict[str, dict] = {}
    for aid, src in initial_artifacts.items():
        artifacts[aid] = _ingest_artifact(exp_id, src)

    # Register every produced artifact id as not-yet-present.
    for ch in chunks:
        for aid in ch.get("produces", []):
            artifacts.setdefault(aid, {"present": False, "path": None, "sha256": None, "bytes": None})

    norm_chunks = []
    for ch in chunks:
        norm_chunks.append({
            "id": ch["id"],
            "kind": ch["kind"],
            "consumes": list(ch.get("consumes", [])),
            "produces": list(ch.get("produces", [])),
            "params": dict(ch.get("params", {})),
            "meta": dict(ch.get("meta", {})),
            "status": "pending",
            "lease": None,
            "wall_s": None,
            "error": None,
            "updated_at": None,
        })

    manifest = {
        "exp_id": exp_id,
        "created_at": _now_iso(),
        "published": True,
        "spec": spec,
        "artifacts": artifacts,
        "chunks": norm_chunks,
    }
    _save_manifest(exp_id, manifest)
    return manifest


# ---------------------------------------------------------------------------
# Artifact helpers (content-addressed).
# ---------------------------------------------------------------------------

def _ingest_artifact(exp_id: str, src: Path) -> dict:
    """Copy ``src`` into the experiment's artifacts/ under a content hash; return
    its registry record."""
    src = Path(src)
    data = src.read_bytes()
    sha = hashlib.sha256(data).hexdigest()
    ext = src.suffix or ".bin"
    dest = _artifacts_dir(exp_id) / f"{sha}{ext}"
    if not dest.exists():
        tmp = dest.with_suffix(dest.suffix + ".tmp")
        tmp.write_bytes(data)
        os.replace(tmp, dest)
    return {
        "present": True,
        "path": str(dest.relative_to(experiment_dir(exp_id))).replace("\\", "/"),
        "sha256": sha,
        "bytes": len(data),
    }


def artifact_file(exp_id: str, manifest: dict, artifact_id: str) -> Path:
    """Absolute path to a present artifact's blob (for a worker to read as input)."""
    rec = manifest["artifacts"].get(artifact_id)
    if not rec or not rec.get("present"):
        raise KeyError(f"artifact {artifact_id!r} not present in {exp_id}")
    return experiment_dir(exp_id) / rec["path"]


# ---------------------------------------------------------------------------
# Runnability + lease lifecycle.
# ---------------------------------------------------------------------------

def _is_runnable(manifest: dict, chunk: dict, now: float) -> bool:
    if chunk["status"] != "pending":
        return False
    lease = chunk.get("lease")
    if lease is not None and lease.get("expires_at", 0) > now:
        return False
    arts = manifest["artifacts"]
    return all(arts.get(aid, {}).get("present") for aid in chunk["consumes"])


def _reclaim_expired(manifest: dict, now: float) -> None:
    for ch in manifest["chunks"]:
        lease = ch.get("lease")
        if ch["status"] == "pending" and lease is not None and lease.get("expires_at", 0) <= now:
            ch["lease"] = None


def pull(exp_id: str, worker_id: str, lease_seconds: float = DEFAULT_LEASE_SECONDS) -> dict | None:
    """Atomically lease the next runnable chunk for ``worker_id``. Returns the
    chunk dict (with resolved consume/produce ids) or None if nothing is runnable
    right now. Resolve input file paths via :func:`artifact_file`."""
    with _locked(exp_id):
        manifest = load_manifest(exp_id)
        now = time.time()
        _reclaim_expired(manifest, now)
        for ch in manifest["chunks"]:
            if _is_runnable(manifest, ch, now):
                ch["lease"] = {"worker": worker_id, "expires_at": now + lease_seconds,
                               "leased_at": now}
                ch["updated_at"] = _now_iso()
                _save_manifest(exp_id, manifest)
                return ch
        return None


def push(exp_id: str, chunk_id: str, outputs: dict[str, Path], wall_s: float) -> dict:
    """Ingest a finished chunk's output files (artifact-id → produced file), mark
    the chunk done and its produced artifacts present. Flips downstream chunks
    runnable. ``outputs`` must cover exactly the chunk's ``produces``."""
    with _locked(exp_id):
        manifest = load_manifest(exp_id)
        ch = _find_chunk(manifest, chunk_id)
        missing = set(ch["produces"]) - set(outputs)
        if missing:
            raise ValueError(f"push for {chunk_id} missing outputs: {sorted(missing)}")
        for aid, src in outputs.items():
            if aid not in ch["produces"]:
                raise ValueError(f"{aid!r} is not produced by {chunk_id}")
            manifest["artifacts"][aid] = _ingest_artifact(exp_id, src)
        ch["status"] = "done"
        ch["lease"] = None
        ch["wall_s"] = float(wall_s)
        ch["error"] = None
        ch["updated_at"] = _now_iso()
        _save_manifest(exp_id, manifest)
        return ch


def mark_failed(exp_id: str, chunk_id: str, error: str) -> dict:
    with _locked(exp_id):
        manifest = load_manifest(exp_id)
        ch = _find_chunk(manifest, chunk_id)
        ch["status"] = "failed"
        ch["lease"] = None
        ch["error"] = str(error)[:2000]
        ch["updated_at"] = _now_iso()
        _save_manifest(exp_id, manifest)
        return ch


def reset_chunk(exp_id: str, chunk_id: str) -> dict:
    """Return a failed/leased chunk to pending so it can be retried."""
    with _locked(exp_id):
        manifest = load_manifest(exp_id)
        ch = _find_chunk(manifest, chunk_id)
        ch["status"] = "pending"
        ch["lease"] = None
        ch["error"] = None
        ch["updated_at"] = _now_iso()
        _save_manifest(exp_id, manifest)
        return ch


# ---------------------------------------------------------------------------
# Reporting.
# ---------------------------------------------------------------------------

def status_summary(exp_id: str) -> dict:
    manifest = load_manifest(exp_id)
    now = time.time()
    by_kind: dict[str, dict[str, int]] = {}
    counts = {"pending": 0, "runnable": 0, "leased": 0, "done": 0, "failed": 0}
    wall_total = 0.0
    for ch in manifest["chunks"]:
        st = ch["status"]
        counts[st] = counts.get(st, 0) + 1
        if st == "pending":
            lease = ch.get("lease")
            if lease is not None and lease.get("expires_at", 0) > now:
                counts["leased"] += 1
            elif _is_runnable(manifest, ch, now):
                counts["runnable"] += 1
        if ch.get("wall_s"):
            wall_total += ch["wall_s"]
        k = by_kind.setdefault(ch["kind"], {"total": 0, "done": 0})
        k["total"] += 1
        k["done"] += int(st == "done")
    return {
        "exp_id": exp_id,
        "counts": counts,
        "by_kind": by_kind,
        "wall_seconds_total": wall_total,
        "n_chunks": len(manifest["chunks"]),
        "complete": counts["done"] == len(manifest["chunks"]),
    }


# ---------------------------------------------------------------------------
# Internals.
# ---------------------------------------------------------------------------

def _find_chunk(manifest: dict, chunk_id: str) -> dict:
    for ch in manifest["chunks"]:
        if ch["id"] == chunk_id:
            return ch
    raise KeyError(f"chunk {chunk_id!r} not in experiment {manifest['exp_id']}")


def _now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime())
