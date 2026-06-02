"""Contributor runner for the dwell-time screen (#43).

Runs on a volunteer's GPU (a free Colab/Kaggle GPU, or their own machine). No
account, no email, no sign-in: it starts immediately with an anonymous session.
The contributor sets a time budget; the runner pulls the simulation the screen
needs next, runs it, sends the result back, and reports each step and the time
spent so they can stop whenever they like.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import requests

import run_chunks


def _fmt(seconds: float) -> str:
    s = int(seconds)
    return f"{s // 60}:{s % 60:02d} min" if s >= 60 else f"{s}s"


def _gpu_line() -> str:
    try:
        out = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                             capture_output=True, text=True, timeout=10)
        if out.returncode == 0 and out.stdout.strip():
            return f"GPU: {out.stdout.strip().splitlines()[0]}"
    except Exception:  # noqa: BLE001
        pass
    return "No GPU detected — set Runtime > Change runtime type > GPU (runs will be slow otherwise)."


def describe(chunk: dict) -> str:
    """A plain-language line for what a chunk computes."""
    p = chunk.get("params", {})
    kind = chunk["kind"]
    if kind == "build":
        return "preparing the simulation system"
    if kind == "equilibrate":
        return f"warming up a replica (seed {p.get('seed', '?')})"
    if kind == "segment":
        return f"simulating dynamics — segment {p.get('index', '?')}, seed {p.get('seed', '?')}"
    return kind


def start_session(health_url: str) -> str:
    """Get an anonymous session token — instant, no account or email."""
    r = requests.post(health_url, params={"action": "anon_session"}, timeout=60)
    r.raise_for_status()
    d = r.json()
    if "token" not in d:
        raise RuntimeError(d.get("error", "could not start a session"))
    return d["token"]


def show_claim_link(health_url: str, token: str) -> None:
    """Show a one-click link to credit this (anonymous) session to the
    contributor. Touches nothing in the kernel, so the run keeps going whether
    or not they open it. Optional — runs count either way."""
    url = f"{health_url}?action=claim&token={token}"
    try:
        from IPython.display import display, HTML
        display(HTML(
            '<div style="padding:10px;margin:6px 0;border:1px solid #6c6;border-radius:6px;background:#f3fff3">'
            '🏅 <b>Optional — credit these runs to you:</b> '
            f'<a href="{url}" target="_blank" rel="noopener">claim them &amp; build reputation</a></div>'))
    except Exception:  # noqa: BLE001 — not in a notebook
        print(f"Optional — credit these runs to you:\n    {url}\n", flush=True)


# --- one pull → run → submit cycle ------------------------------------------

def dispatch(health_url: str, token: str, done: str | None = None) -> dict:
    params = {"action": "dispatch", "token": token}
    if done:
        params["done"] = done   # release the just-finished lease so the next one can be assigned
    r = requests.post(health_url, params=params, timeout=60)
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


def run_once(health_url: str, token: str, done: str | None = None) -> dict:
    """Pull → run → submit one chunk. Returns {stop|chunk_id, seconds, lease, ...}
    and prints what's happening as it goes."""
    asg = dispatch(health_url, token, done)
    if asg.get("status") != "assigned":
        return {"stop": True, "msg": asg.get("message") or asg.get("error") or "no work available"}
    broker = asg.get("broker_url")
    if not broker:
        return {"stop": True, "msg": "this site has no result broker configured yet"}

    chunk = asg["chunk"]
    print(f"  ▶ {describe(chunk)}", flush=True)
    scratch = Path(tempfile.mkdtemp(prefix=f"contrib_{chunk['id']}_"))
    local = download_inputs(asg["resolve_base"], asg.get("inputs", {}), scratch)
    t0 = time.time()
    outputs = run_chunks.execute_chunk(chunk, lambda aid: local[aid], scratch)
    secs = time.time() - t0
    submit(broker, asg["lease"], outputs, secs)
    print(f"    done in {_fmt(secs)} · sent back ✓", flush=True)
    return {"stop": False, "chunk_id": chunk["id"], "seconds": secs, "lease": asg["lease"]}


# --- the contributor's whole session ----------------------------------------

def contribute(health_url: str, minutes: float = 30) -> None:
    """Pair, then run simulations for up to `minutes`, reporting progress.
    Stops on the time budget, when there's no work, or when interrupted —
    always with a clear message, never an opaque loop."""
    print(_gpu_line(), flush=True)
    token = start_session(health_url)
    show_claim_link(health_url, token)

    budget = max(1.0, float(minutes)) * 60.0
    start = time.time()
    done_count, compute = 0, 0.0
    last_lease = None
    print(f"Running for up to {int(minutes)} min. Stop whenever you like — an "
          "unfinished task is simply reassigned, so nothing is wasted.\n", flush=True)

    # Steady liveness line every 10 s. Two jobs: shows the run is alive between
    # the MD heartbeats, and — since Colab redraws a cleared output only when new
    # output arrives — makes the output reappear within seconds of coming back to
    # the notebook, so it never looks cancelled.
    stop_hb = threading.Event()

    def _heartbeat():
        while not stop_hb.wait(10):
            print(f"  · still running — {_fmt(time.time() - start)} elapsed, {done_count} done", flush=True)
    threading.Thread(target=_heartbeat, daemon=True).start()

    while True:
        remaining = budget - (time.time() - start)
        if remaining <= 0:
            print("\nTime budget reached — wrapping up.", flush=True)
            break
        print(f"[{_fmt(time.time() - start)} in · {_fmt(remaining)} left]", flush=True)
        try:
            r = run_once(health_url, token, done=last_lease)
        except KeyboardInterrupt:
            stop_hb.set()
            print("\nStopped. The current task will be reassigned.", flush=True)
            return
        except Exception as e:  # noqa: BLE001
            print(f"    task failed ({e}) — it will be reassigned to someone else\n", flush=True)
            last_lease = None
            continue
        if r.get("stop"):
            print(f"\n{r['msg']}.", flush=True)
            break
        done_count += 1
        compute += r["seconds"]
        last_lease = r["lease"]
        print(f"    so far this session: {done_count} simulation(s), {_fmt(compute)} of compute\n", flush=True)

    stop_hb.set()
    print(f"\nThank you — you ran {done_count} simulation(s) ({_fmt(compute)} of compute).", flush=True)


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--health-url", required=True, help="https://<site>/screen")
    ap.add_argument("--minutes", type=float, default=30, help="how long to run")
    args = ap.parse_args()
    contribute(args.health_url, minutes=args.minutes)


if __name__ == "__main__":
    main()
