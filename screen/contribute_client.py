"""Account-free contributor loop for the dwell-time screen (#43).

Runs on a volunteer's GPU (a free Colab/Kaggle GPU, or their own machine). It
needs no Hugging Face account, no write token, and — crucially — no email typed
into the notebook. Identity is handed off *from the app*: the runner pairs the
session (a short code + a link the signed-in contributor opens on the site),
receives a session token carrying only their pseudonym, then:

  1. asks health for an assignment with that token (health decides what to run);
  2. downloads the chunk's inputs over plain HTTPS (public dataset, no auth);
  3. runs the MD step via the existing ``run_chunks.execute_chunk``;
  4. POSTs the outputs + lease to the broker (URL comes back in the dispatch
     response), which writes them into the dataset's ``submissions/`` inbox.

Acceptance happens later in ``hf_store ingest``.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
from pathlib import Path

import requests

import run_chunks


# --- pairing: get a session token without typing an email --------------------

def pair(health_url: str, *, poll_s: float = 3.0, timeout_s: float = 600.0) -> tuple[str, str]:
    """Hand identity off to the app. Prints a link for the (already signed-in)
    contributor to open, polls until linked, returns (session_token, broker_url)."""
    r = requests.post(health_url, params={"action": "pair_start"}, timeout=60)
    r.raise_for_status()
    d = r.json()
    if "verification_url" not in d:
        raise RuntimeError(d.get("error", "pairing not available"))
    print("To link this session to your account, open:\n"
          f"    {d['verification_url']}\n"
          "(you're already signed in there — just confirm). Waiting…", flush=True)
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        pr = requests.get(health_url, params={"action": "pair_poll", "code": d["code"]}, timeout=30).json()
        if pr.get("status") == "linked":
            print("  linked — thank you!", flush=True)
            return pr["token"], pr.get("broker_url", "")
        if pr.get("status") == "expired":
            raise RuntimeError("pairing code expired — re-run to get a fresh one")
        time.sleep(poll_s)
    raise RuntimeError("pairing timed out — re-run to try again")


# --- one pull → run → submit cycle ------------------------------------------

def dispatch(health_url: str, token: str) -> dict:
    r = requests.post(health_url, params={"action": "dispatch", "token": token}, timeout=60)
    if not r.ok and r.status_code not in (429, 503):
        r.raise_for_status()
    return r.json()


def download_inputs(resolve_base: str, inputs: dict, scratch: Path) -> dict:
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
    manifest = {aid: Path(p).name for aid, p in outputs.items()}
    files = [("files", (Path(p).name, open(p, "rb"))) for p in outputs.values()]
    try:
        r = requests.post(
            broker_url.rstrip("/") + "/submit",
            data={"lease": lease, "wall_s": str(wall_s), "manifest": json.dumps(manifest)},
            files=files, timeout=600,
        )
    finally:
        for _, (_, fh) in files:
            fh.close()
    r.raise_for_status()
    return r.json()


def run_once(health_url: str, token: str) -> str:
    """One pull → run → submit cycle. Returns a short status string."""
    asg = dispatch(health_url, token)
    if asg.get("status") != "assigned":
        return asg.get("message") or asg.get("error") or asg.get("status") or "no assignment"

    broker_url = asg.get("broker_url")
    if not broker_url:
        return "the site has no broker configured yet — results can't be submitted"

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
    ap.add_argument("--token", default=None, help="session token (else pair interactively)")
    ap.add_argument("--n", type=int, default=5, help="how many chunks to run")
    args = ap.parse_args()

    token = args.token or pair(args.health_url)[0]
    done = 0
    for _ in range(max(1, args.n)):
        try:
            msg = run_once(args.health_url, token)
        except Exception as e:  # noqa: BLE001
            print(f"  error: {e}", flush=True)
            break
        print(f"  {msg}", flush=True)
        if not msg.startswith("submitted"):
            break
        done += 1
    print(f"\nthank you — ran {done} chunk(s) this session.", flush=True)


if __name__ == "__main__":
    main()
