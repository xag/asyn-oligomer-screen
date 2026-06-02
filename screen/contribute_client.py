"""Account-free contributor loop for the dwell-time screen (#43).

This is what runs on a volunteer's GPU (a Colab/Kaggle free GPU, or their own
machine). It needs no Hugging Face account and no write token — only an email,
the same identity the health site already uses:

  1. ask health for a work assignment (a signed lease + the chunk + where to
     fetch its inputs from the *public* dataset);
  2. download the inputs (plain HTTPS GET — public dataset, no auth);
  3. run the MD step with the existing ``run_chunks.execute_chunk``;
  4. POST the outputs + the lease to the broker, which verifies the lease and
     writes them into the dataset's ``submissions/`` inbox.

Health dispatch decides *what* to run (never the contributor), so a given email
is never handed the same chunk twice and distinct contributors spread across
distinct chunks. Acceptance happens later in ``hf_store ingest``.

Usage:
    python contribute_client.py --health-url https://<site>/screen \
        --broker-url https://<space>.hf.space --email you@example.com --n 3
"""
from __future__ import annotations

import argparse
import sys
import tempfile
import time
from pathlib import Path

import requests

import run_chunks


def dispatch(health_url: str, email: str) -> dict:
    """Request one assignment. Returns the JSON ('assigned' / 'idle' / 'busy')."""
    r = requests.post(health_url, params={"action": "dispatch", "email": email}, timeout=60)
    if r.status_code in (429, 503) or not r.ok:
        try:
            return r.json()
        except Exception:  # noqa: BLE001
            r.raise_for_status()
    return r.json()


def download_inputs(resolve_base: str, inputs: dict, scratch: Path) -> dict:
    """Fetch each consumed artifact to ``scratch``; return {artifact_id: Path}."""
    local: dict[str, Path] = {}
    for aid, path in inputs.items():
        dest = scratch / aid.replace("/", "__")
        dest.parent.mkdir(parents=True, exist_ok=True)
        with requests.get(resolve_base + path, stream=True, timeout=300) as r:
            r.raise_for_status()
            with dest.open("wb") as fh:
                for block in r.iter_content(chunk_size=1 << 20):
                    fh.write(block)
        local[aid] = dest
    return local


def submit(broker_url: str, lease: str, outputs: dict, wall_s: float) -> dict:
    """Upload produced artifacts + the lease to the broker."""
    manifest = {aid: Path(p).name for aid, p in outputs.items()}
    files = [("files", (Path(p).name, open(p, "rb"))) for p in outputs.values()]
    try:
        import json as _json
        r = requests.post(
            broker_url.rstrip("/") + "/submit",
            data={"lease": lease, "wall_s": str(wall_s), "manifest": _json.dumps(manifest)},
            files=files, timeout=600,
        )
    finally:
        for _, (_, fh) in files:
            fh.close()
    r.raise_for_status()
    return r.json()


def run_once(health_url: str, broker_url: str, email: str) -> str:
    """One pull → run → submit cycle. Returns a short status string."""
    asg = dispatch(health_url, email)
    status = asg.get("status")
    if status != "assigned":
        return asg.get("message") or status or "no assignment"

    chunk = asg["chunk"]
    print(f"  assigned {chunk['id']} [{chunk['kind']}]", flush=True)
    scratch = Path(tempfile.mkdtemp(prefix=f"contrib_{chunk['id']}_"))
    local = download_inputs(asg["resolve_base"], asg.get("inputs", {}), scratch)

    t0 = time.time()
    outputs = run_chunks.execute_chunk(chunk, lambda aid: local[aid], scratch)
    wall = time.time() - t0

    res = submit(broker_url, asg["lease"], outputs, wall)
    return f"submitted {chunk['id']} ({wall / 60:.1f} min) -> {res.get('status')}"


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--health-url", required=True, help="https://<site>/screen")
    ap.add_argument("--broker-url", required=True, help="https://<space>.hf.space")
    ap.add_argument("--email", required=True)
    ap.add_argument("--n", type=int, default=1, help="how many chunks to run")
    args = ap.parse_args()

    done = 0
    for _ in range(max(1, args.n)):
        try:
            msg = run_once(args.health_url, args.broker_url, args.email)
        except Exception as e:  # noqa: BLE001 — report + stop the loop
            print(f"  error: {e}", flush=True)
            break
        print(f"  {msg}", flush=True)
        if not msg.startswith("submitted"):
            break   # idle / busy — nothing more to do right now
        done += 1
    print(f"\nthank you — ran {done} chunk(s) this session.", flush=True)


if __name__ == "__main__":
    main()
