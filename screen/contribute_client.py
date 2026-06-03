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


# Single output sink so the 10-second liveness dots never collide with real
# log lines (the contributor-facing messages from describe() and the session).
# A dot is written without a newline; real text starts on a fresh line if dots
# were pending, so the log reads as "....<line>" — not a line buried after a
# long dot trail.
_out_lock = threading.Lock()
_dots_pending = False


def _say(msg: str = "") -> None:
    global _dots_pending
    with _out_lock:
        sys.stdout.write(("\n" if _dots_pending else "") + str(msg) + "\n")
        sys.stdout.flush()
        _dots_pending = False


def _tick() -> None:
    global _dots_pending
    with _out_lock:
        sys.stdout.write(".")
        sys.stdout.flush()
        _dots_pending = True


def _gpu_line() -> str:
    try:
        out = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                             capture_output=True, text=True, timeout=10)
        if out.returncode == 0 and out.stdout.strip():
            return f"GPU: {out.stdout.strip().splitlines()[0]}"
    except Exception:  # noqa: BLE001
        pass
    return "No GPU detected — set Runtime > Change runtime type > GPU (runs will be slow otherwise)."


def _molecule(meta: dict) -> str:
    """The candidate-molecule id (e.g. 'dopamine', 'l-dopa'), or '' for the apo
    baseline arm — the id is already human-readable, so it's shown as-is."""
    lig = meta.get("ligand", "") or ""
    return "" if lig in ("", "apo") else lig


def describe(chunk: dict) -> str:
    """A plain-language line saying, concretely, what this run simulates: the
    toxic α-synuclein core on its own (the baseline) or with a candidate molecule
    docked onto it, and which stage of the simulation it is."""
    meta = chunk.get("meta", {})
    mol = _molecule(meta)
    # The system under simulation, named the same way in every stage.
    system = (f"α-synuclein with {mol} docked onto it" if mol
              else "the α-synuclein baseline (no molecule)")
    kind = chunk["kind"]
    if kind == "build":
        return (f"Setting up: docking {mol} onto the toxic α-synuclein core and "
                "surrounding it with water." if mol
                else "Setting up the baseline: the toxic α-synuclein core alone, "
                     "surrounded with water.")
    if kind == "equilibrate":
        return f"Warming up to body temperature: {system}."
    if kind == "segment":
        part = int(chunk.get("params", {}).get("index", 0)) + 1
        return (f"Simulating {system} — watching whether the toxic shape holds "
                f"or loosens (part {part}).")
    return kind


def start_session(health_url: str) -> str:
    """Get an anonymous session token — instant, no account or email. Used only
    when the notebook was opened without a personal identity token (the no-sign-in
    path); the per-user notebook from the site already carries its own token."""
    r = requests.post(health_url, params={"action": "anon_session"}, timeout=60)
    r.raise_for_status()
    d = r.json()
    if "token" not in d:
        raise RuntimeError(d.get("error", "could not start a session"))
    return d["token"]


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
    _say(f"  ▶ {describe(chunk)}")
    scratch = Path(tempfile.mkdtemp(prefix=f"contrib_{chunk['id']}_"))
    local = download_inputs(asg["resolve_base"], asg.get("inputs", {}), scratch)
    t0 = time.time()
    outputs = run_chunks.execute_chunk(chunk, lambda aid: local[aid], scratch)
    secs = time.time() - t0
    submit(broker, asg["lease"], outputs, secs)
    _say(f"    done in {_fmt(secs)} · sent back ✓")
    return {"stop": False, "chunk_id": chunk["id"], "seconds": secs, "lease": asg["lease"]}


# --- the contributor's whole session ----------------------------------------

def contribute(health_url: str, minutes: float = 30, token: str | None = None) -> None:
    """Run simulations for up to `minutes`, reporting progress. Stops on the time
    budget, when there's no work, or when interrupted — always with a clear
    message, never an opaque loop.

    `token` is the personal identity token baked into the notebook the site
    generated, so every run is credited to the signed-in contributor with no
    link to click. Opened without one (the plain repo notebook), it falls back to
    an anonymous session — runs still count, they're just not credited to anyone.
    """
    # The MD engine's own per-line output is deliberately not surfaced: it's
    # low-level (atom counts, energies, scratch file names) and obscures what the
    # contributor is actually doing. describe() says that in plain words and the
    # heartbeat dots show it's alive; a failing chunk still reports a log tail
    # (see run_chunks._run), so nothing diagnostic is lost.
    run_chunks._emit = lambda *_: None
    _say(_gpu_line())
    if not token:
        token = start_session(health_url)
        _say("Running anonymously (no sign-in) — your runs count but aren't "
             "credited to you. Open the notebook from the site to get credit.\n")

    budget = max(1.0, float(minutes)) * 60.0
    start = time.time()
    done_count, compute = 0, 0.0
    last_lease = None
    _say(f"Running for up to {int(minutes)} min. Stop whenever you like — an "
         "unfinished task is simply reassigned, so nothing is wasted.\n")

    # A dot every 10 s during quiet stretches keeps output flowing — Colab
    # redraws cleared output only when new output arrives, so this makes the log
    # reappear within ~10 s of coming back to the notebook, never looking
    # cancelled. Real lines reset the trail (see _say), so it's not endless dots.
    stop_hb = threading.Event()

    def _heartbeat():
        while not stop_hb.wait(10):
            _tick()
    threading.Thread(target=_heartbeat, daemon=True).start()

    while True:
        remaining = budget - (time.time() - start)
        if remaining <= 0:
            _say("\nTime budget reached — wrapping up.")
            break
        _say(f"[{_fmt(time.time() - start)} in · {_fmt(remaining)} left]")
        try:
            r = run_once(health_url, token, done=last_lease)
        except KeyboardInterrupt:
            stop_hb.set()
            _say("\nStopped. The current task will be reassigned.")
            return
        except Exception as e:  # noqa: BLE001
            _say(f"    task failed ({e}) — it will be reassigned to someone else\n")
            last_lease = None
            continue
        if r.get("stop"):
            _say(f"\n{r['msg']}.")
            break
        done_count += 1
        compute += r["seconds"]
        last_lease = r["lease"]
        _say(f"    so far this session: {done_count} simulation(s), {_fmt(compute)} of compute\n")

    stop_hb.set()
    _say(f"\nThank you — you ran {done_count} simulation(s) ({_fmt(compute)} of compute).")


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--health-url", required=True, help="https://<site>/screen")
    ap.add_argument("--minutes", type=float, default=30, help="how long to run")
    ap.add_argument("--token", default=None,
                    help="personal identity token (the site's notebook bakes this in; "
                         "omit to run anonymously)")
    args = ap.parse_args()
    contribute(args.health_url, minutes=args.minutes, token=args.token)


if __name__ == "__main__":
    main()
