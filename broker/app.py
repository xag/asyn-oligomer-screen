"""Upload broker for crowdsourced dwell-time chunks (#43).

Runs as a free Hugging Face *Space*. It is the only writer to the dataset's
``submissions/`` inbox and the only holder of the dataset write token, so
contributors need no Hugging Face account: their notebook just POSTs a finished
chunk here with the signed lease health issued. The broker

  1. verifies the lease (shared SCREEN_LEASE_KEY) — identity = the ``pseudonym``
     inside it; the email never reaches here;
  2. size-caps + SHA-256s each uploaded artifact;
  3. writes ``submissions/<chunk_id>/<pseudonym>/{<file>, meta.json}`` in exactly
     the shape ``hf_store.cmd_ingest`` already reads.

It is deliberately dumb: no scoring, no acceptance. The reputation-weighted,
distinct-pseudonym quorum + spot-check all happen later in ``ingest``. The worst
an attacker with a leaked token can do is overwrite their *own*
chunk×pseudonym slot — a candidate that still has to survive ingest.

Secrets (set as Space secrets): HF_TOKEN (dataset write), SCREEN_LEASE_KEY
(shared with health), SCREEN_DATASET_REPO (e.g. user/asyn-dwell-results).
"""
from __future__ import annotations

import hashlib
import io
import json
import os

from fastapi import FastAPI, Form, HTTPException, UploadFile
from huggingface_hub import CommitOperationAdd, HfApi

import lease as L  # bundled copy of screen/lease.py

MAX_FILE_BYTES = int(os.environ.get("BROKER_MAX_FILE_BYTES", 64 * 1024 * 1024))
MAX_TOTAL_BYTES = int(os.environ.get("BROKER_MAX_TOTAL_BYTES", 256 * 1024 * 1024))
DATASET = "dataset"

app = FastAPI(title="asyn dwell-time upload broker")


def _repo() -> str:
    repo = os.environ.get("SCREEN_DATASET_REPO")
    if not repo:
        raise HTTPException(503, "SCREEN_DATASET_REPO not configured")
    return repo


@app.get("/")
def info():
    return {"service": "asyn dwell-time upload broker", "ok": True}


@app.post("/submit")
async def submit(
    lease: str = Form(...),
    wall_s: float = Form(0.0),
    manifest: str = Form(...),   # JSON {artifact_id: uploaded_filename}
    files: list[UploadFile] = None,
):
    payload = L.verify_lease(lease)
    if payload is None:
        raise HTTPException(401, "invalid or expired lease")
    chunk_id = payload["chunk_id"]
    pseudonym = payload["pseudonym"]

    try:
        out_map: dict[str, str] = json.loads(manifest)
        assert isinstance(out_map, dict)
    except Exception:
        raise HTTPException(400, "manifest must be a JSON object {artifact_id: filename}")

    by_name = {f.filename: f for f in (files or [])}
    ops, meta_outputs, total = [], {}, 0
    base = f"submissions/{chunk_id}/{pseudonym}"
    for aid, fname in out_map.items():
        up = by_name.get(fname)
        if up is None:
            raise HTTPException(400, f"missing uploaded file for {aid!r} ({fname!r})")
        data = await up.read()
        if len(data) > MAX_FILE_BYTES:
            raise HTTPException(413, f"{aid} exceeds per-file cap ({MAX_FILE_BYTES} bytes)")
        total += len(data)
        if total > MAX_TOTAL_BYTES:
            raise HTTPException(413, "submission exceeds total size cap")
        dest = aid.replace("/", "__")
        ops.append(CommitOperationAdd(f"{base}/{dest}", io.BytesIO(data)))
        meta_outputs[aid] = {"sha256": hashlib.sha256(data).hexdigest(), "file": dest}

    if not meta_outputs:
        raise HTTPException(400, "no outputs submitted")

    meta = {
        "chunk_id": chunk_id,
        "worker": pseudonym,           # the ingest's contributor id
        "iat": payload.get("iat", 0),
        "wall_s": float(wall_s),
        "outputs": meta_outputs,
    }
    ops.append(CommitOperationAdd(f"{base}/meta.json",
                                  io.BytesIO(json.dumps(meta, indent=2).encode("utf-8"))))

    HfApi(token=os.environ["HF_TOKEN"]).create_commit(
        _repo(), ops, repo_type=DATASET,
        commit_message=f"submit {chunk_id} by {pseudonym}",
    )
    return {"status": "received", "chunk_id": chunk_id, "outputs": list(meta_outputs)}
