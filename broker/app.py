"""Write broker for the crowdsourced dwell-time screen.

Runs as a free Hugging Face *Space*. It is the only writer to the dataset and
the only holder of the dataset write token, so contributors need no Hugging Face
account. health is the read + identity + dispatch surface (no token, no GPU);
every endpoint here *mutates the dataset*, which is why it lives on the
token-holder:

  POST /api/screen/v1/results    a finished chunk's outputs → submissions/ inbox
                                 (auth = the signed lease health issued).
  POST /api/screen/v1/molecules  a candidate molecule proposal → appended to
                                 molecules.json as source="contributed",
                                 prep="awaiting". Validated (SMILES parses) but
                                 not trusted: the metadata is advisory and the
                                 coordinator's prep pass docks it into runnable
                                 chunks; moderation governs promotion to primary.

`/submit` is kept as a deprecated alias of the results endpoint for one release
(old clients). The broker is deliberately dumb about results: no scoring, no
acceptance — the reputation-weighted, distinct-pseudonym quorum + spot-check all
happen later in `ingest`. The worst an attacker with a leaked lease can do is
overwrite their *own* chunk×pseudonym slot, a candidate that still has to
survive ingest.

Secrets (Space secrets): HF_TOKEN (dataset write), SCREEN_LEASE_KEY (shared with
health), SCREEN_DATASET_REPO (e.g. user/asyn-dwell-results).
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import re
import urllib.request

from fastapi import FastAPI, Form, HTTPException, UploadFile
from huggingface_hub import CommitOperationAdd, HfApi
from pydantic import BaseModel

import lease as L  # bundled copy of screen/lease.py

MAX_FILE_BYTES = int(os.environ.get("BROKER_MAX_FILE_BYTES", 64 * 1024 * 1024))
MAX_TOTAL_BYTES = int(os.environ.get("BROKER_MAX_TOTAL_BYTES", 256 * 1024 * 1024))
DATASET = "dataset"

# Mirror of data/vicinity_molecules.js enums — keep in sync. A contributed
# proposal must use these; everything else is free text the moderator curates.
GROUPS = {"endogenous-metabolites", "dietary", "neurotransmitter", "metal",
          "lipid", "gut-derived", "environmental"}
ROLES = {"anchor", "candidate", "both"}
ROUTES = {"endogenous", "diet", "supplement", "precursor", "microbiome",
          "environmental", "injection", "none-known"}
FEASIBILITIES = {"native", "achievable", "low-bioavailability", "invasive-only", "unknown"}

app = FastAPI(title="asyn dwell-time write broker")


def _repo() -> str:
    repo = os.environ.get("SCREEN_DATASET_REPO")
    if not repo:
        raise HTTPException(503, "SCREEN_DATASET_REPO not configured")
    return repo


def _api() -> HfApi:
    return HfApi(token=os.environ["HF_TOKEN"])


@app.get("/")
def info():
    return {"service": "asyn dwell-time write broker", "ok": True}


# --- results ----------------------------------------------------------------

async def _ingest_results(lease: str, wall_s: float, manifest: str,
                          files: list[UploadFile] | None, token: str | None) -> dict:
    payload = L.verify_lease(lease)
    if payload is None:
        raise HTTPException(401, "invalid or expired lease")
    chunk_id = payload["chunk_id"]
    pseudonym = payload["pseudonym"]

    # Optional identity token the client carries everywhere: a consistency guard,
    # not the authority. Identity is bound at dispatch (the lease's pseudonym);
    # if a token is also sent it must agree, so a token can never silently
    # re-credit someone else's lease.
    if token:
        tp = L.verify_lease(token)
        if not tp or tp.get("pseudonym") != pseudonym:
            raise HTTPException(400, "token does not match the lease")

    try:
        out_map: dict[str, str] = json.loads(manifest)
        assert isinstance(out_map, dict)
    except Exception:
        raise HTTPException(400, "manifest must be a JSON object {artifact_id: filename}")

    by_name = {f.filename: f for f in (files or [])}
    ops, meta_outputs, total = [], {}, 0
    # Unit ids carry '#'/'@' (e.g. rid#segment@5). They land in dataset paths that
    # are later fetched as resolve_base + path, where '#' would be read as a URL
    # fragment — so make the storage path URL-safe. The real id stays in meta.json.
    safe_chunk = re.sub(r"[^A-Za-z0-9_.-]", "_", chunk_id)
    base = f"submissions/{safe_chunk}/{pseudonym}"
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

    _api().create_commit(_repo(), ops, repo_type=DATASET,
                         commit_message=f"results {chunk_id} by {pseudonym}")
    # Return the stored dataset paths so a work-store client can point the replica
    # cursor at the state it just uploaded (resolvable via the dataset resolve base).
    paths = {aid: f"{base}/{rec['file']}" for aid, rec in meta_outputs.items()}
    return {"status": "accepted", "chunk_id": chunk_id, "outputs": list(meta_outputs), "paths": paths}


@app.post("/api/screen/v1/results")
async def results(
    lease: str = Form(...),
    wall_s: float = Form(0.0),
    manifest: str = Form(...),   # JSON {artifact_id: uploaded_filename}
    token: str | None = Form(None),
    files: list[UploadFile] = None,
):
    return await _ingest_results(lease, wall_s, manifest, files, token)


@app.post("/submit")
async def submit(
    lease: str = Form(...),
    wall_s: float = Form(0.0),
    manifest: str = Form(...),
    token: str | None = Form(None),
    files: list[UploadFile] = None,
):
    # Deprecated alias of /api/screen/v1/results. Kept for one release so older
    # clients keep working; new clients use the versioned path. Returns the same
    # body, with the legacy "received" status preserved for old parsers.
    out = await _ingest_results(lease, wall_s, manifest, files, token)
    out["status"] = "received"
    return out


# --- molecules --------------------------------------------------------------

class MoleculeProposal(BaseModel):
    name: str
    smiles: str
    group: str | None = None
    role: str | None = "candidate"
    pdb_ligand: str | None = None
    mw_da: float | None = None
    cns_conc: dict | None = None
    delivery: dict | None = None
    evidence: str | None = None
    refs: list[str] | None = None


MOLECULES_FILE = "molecules.json"


def _slug(name: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "-", name.strip().lower()).strip("-")
    return s or "molecule"


def _read_molecules(repo: str) -> list[dict]:
    url = f"https://huggingface.co/datasets/{repo}/raw/main/{MOLECULES_FILE}"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:  # noqa: S310
            data = json.loads(r.read().decode("utf-8"))
        return data if isinstance(data, list) else []
    except Exception:
        return []  # absent before the coordinator first publishes it


@app.post("/api/screen/v1/molecules")
def post_molecule(mol: MoleculeProposal):
    repo = _repo()

    name = (mol.name or "").strip()
    smiles = (mol.smiles or "").strip()
    if not name or not smiles:
        raise HTTPException(400, "name and smiles are required")

    # Structure must parse — the only field acted on automatically (docked, then
    # simulated). Lazy import so the results path never pays for rdkit.
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        raise HTTPException(503, "molecule validation unavailable (rdkit not installed)")
    if Chem.MolFromSmiles(smiles) is None:
        raise HTTPException(400, f"could not parse SMILES {smiles!r}")

    if mol.group is not None and mol.group not in GROUPS:
        raise HTTPException(400, f"group must be one of {sorted(GROUPS)}")
    role = mol.role or "candidate"
    if role not in ROLES:
        raise HTTPException(400, f"role must be one of {sorted(ROLES)}")
    if mol.delivery and mol.delivery.get("route") and mol.delivery["route"] not in ROUTES:
        raise HTTPException(400, f"delivery.route must be one of {sorted(ROUTES)}")
    if mol.delivery and mol.delivery.get("feasibility") and mol.delivery["feasibility"] not in FEASIBILITIES:
        raise HTTPException(400, f"delivery.feasibility must be one of {sorted(FEASIBILITIES)}")

    molecules = _read_molecules(repo)
    existing = {m.get("id") for m in molecules if isinstance(m, dict)}
    mid = _slug(name)
    if mid in existing:
        raise HTTPException(409, f"a molecule with id {mid!r} already exists")

    record = {
        "id": mid,
        "name": name,
        "source": "contributed",   # primary | contributed
        "prep": "awaiting",        # awaiting | ready (set by the coordinator prep pass)
        "group": mol.group,
        "role": role,
        "smiles": smiles,
        "pdb_ligand": mol.pdb_ligand,
        "mw_da": mol.mw_da,
        "cns_conc": mol.cns_conc,
        "delivery": mol.delivery,
        "evidence": mol.evidence,
        "refs": mol.refs or [],
    }
    molecules.append(record)
    _api().create_commit(
        repo,
        [CommitOperationAdd(MOLECULES_FILE,
                            io.BytesIO(json.dumps(molecules, indent=2).encode("utf-8")))],
        repo_type=DATASET,
        commit_message=f"molecule proposal {mid}",
    )
    return {"id": mid, "name": name, "source": "contributed", "prep": "awaiting"}
