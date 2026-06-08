"""Reference contributor runner for the dwell-time screen (#43).

Runs on a volunteer's GPU (a free Colab/Kaggle GPU, or their own machine) and
speaks the public contributor API:

  GET  {base}/api/screen/v1/work     lease the next simulation
  POST {results_url}                 send its outputs back (URL comes from /work)

No sign-in is required: omit the token to run anonymously (runs still count, just
uncredited). To get credit, pass the email-verified token from {base}/screen —
either ``--token`` or the ``ASYN_CONTRIB_TOKEN`` env var. There is no session to
create; the token is a long-lived credential carried on every call.

This is just *one* client. Any program that speaks the API works — the point of
the API is that the runner is replaceable (a pip CLI, a curl one-liner, a Colab
cell). It reports each step and the time spent so a contributor can stop whenever
they like.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import requests

import run_chunks

DEFAULT_BASE_URL = "https://health-two-iota.vercel.app"


def _fmt(seconds: float) -> str:
    s = int(seconds)
    return f"{s // 60}:{s % 60:02d} min" if s >= 60 else f"{s}s"


# Single output sink so the 10-second liveness dots never collide with real
# log lines (the contributor-facing messages from the session). A dot is written
# without a newline; real text starts on a fresh line if dots were pending, so
# the log reads as "....<line>" — not a line buried after a long dot trail.
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


def runnable_kinds() -> list[str]:
    """Which unit kinds this machine can actually run, so /work only hands us those.
    Every box runs the MD itself (pip OpenMM); docking needs RDKit/Meeko/Vina and
    force-field prep needs the OpenFF conda env — a pip-only GPU client (e.g. Colab)
    has neither, so it sticks to the simulations."""
    import importlib.util as iu
    have = lambda m: iu.find_spec(m) is not None  # noqa: E731
    kinds = ["equilibrate", "segment"]
    if have("rdkit") and have("meeko") and have("vina"):
        kinds.append("dock")
    try:
        import md_env
        md_env.md_python()              # raises if the OpenFF conda env is absent
        kinds.append("build")
    except Exception:  # noqa: BLE001
        pass
    return kinds


# --- API endpoints ----------------------------------------------------------

def _work_url(base: str) -> str:
    return base.rstrip("/") + "/api/screen/v1/work"


def get_work(base: str, token: str | None, molecules: str | None,
             done: str | None, kinds: str | None = None) -> dict:
    """GET the next simulation. `molecules` is a comma-separated preference list
    in decreasing interest; `done` releases the lease just finished; `kinds`
    restricts dispatch to the unit kinds this machine can run."""
    params: dict[str, str] = {}
    if molecules:
        params["molecules"] = molecules
    if done:
        params["done"] = done
    if kinds:
        params["kinds"] = kinds
    headers = {"Authorization": f"Bearer {token}"} if token else {}
    r = requests.get(_work_url(base), params=params, headers=headers, timeout=60)
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


def submit_results(results_url: str, lease: str, token: str | None,
                   outputs: dict, wall_s: float) -> dict:
    manifest = {aid: Path(p).name for aid, p in outputs.items()}
    data = {"lease": lease, "wall_s": str(wall_s), "manifest": json.dumps(manifest)}
    if token:
        data["token"] = token
    files = [("files", (Path(p).name, open(p, "rb"))) for p in outputs.values()]
    try:
        r = requests.post(results_url, data=data, files=files, timeout=600)
    finally:
        for _, (_, fh) in files:
            fh.close()
    r.raise_for_status()
    return r.json()


# --- run one unit: download inputs → run by kind → upload → report ----------
# Every unit (dock, build, equilibrate, segment) arrives in one shape: inputs to
# download, params, and an `outputs` map (local filename → artifact id). We run the
# kind-appropriate step, upload the produced files to the broker, and report the
# resulting dataset paths to /checkpoint, which registers the produced artifacts and
# advances the job. The MD steps checkpoint every `checkpoint_s` seconds, so on an
# interrupt the latest state is already saved and only the last few seconds are lost.

_CKPT_RE = re.compile(r"^CHECKPOINT\s+([0-9.]+)\s*$")


def _run_md_streaming(cmd: list[str]) -> float | None:
    """Run md_relax, ticking for liveness and tracking the last `CHECKPOINT <ps>`
    it prints. On KeyboardInterrupt, stop the child (its latest checkpoint is
    already on disk) and re-raise. Returns the last checkpoint ps, or None."""
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, encoding="utf-8", errors="replace", bufsize=1, env=env)
    last_ps = None
    try:
        for line in proc.stdout:
            m = _CKPT_RE.match(line.strip())
            if m:
                last_ps = float(m.group(1))
            else:
                _tick()
        if proc.wait() != 0:
            raise RuntimeError(f"md step exited {proc.returncode}")
        return last_ps
    except KeyboardInterrupt:
        proc.terminate()
        try:
            proc.wait(timeout=20)
        except Exception:  # noqa: BLE001
            proc.kill()
        raise


def run_unit(base: str, w: dict, token: str | None, scratch: Path) -> dict:
    chunk = w["chunk"]
    p = chunk.get("params", {})
    kind = chunk.get("kind")
    local = download_inputs(w["resolve_base"], w.get("inputs", {}), scratch)
    t0 = time.time()
    interrupted, last_ps = False, None
    produced: dict[str, Path] = {}   # local filename → produced file

    if kind in ("equilibrate", "segment"):
        state_out = scratch / "state.xml"
        seg_out = scratch / "seg.pdb"
        py, mdr = sys.executable, str(run_chunks.MD_RELAX)
        if kind == "equilibrate":
            cmd = [py, mdr, "--equilibrate", str(state_out),
                   "--system-xml", str(local["system.xml"]), "--solvated-pdb", str(local["solvated.pdb"]),
                   "--equil-ps", str(p["equil_ps"]), "--temperature-k", str(p["temperature_k"]),
                   "--seed", str(p["seed"]), "--report-interval-ps", "2"]
        else:
            cmd = [py, mdr, "--segment",
                   "--system-xml", str(local["system.xml"]), "--solvated-pdb", str(local["solvated.pdb"]),
                   "--state-in", str(local["state_in.xml"]), "--state-out", str(state_out),
                   "--seg-out", str(seg_out), "--segment-ps", str(p["step_ps"]),
                   "--traj-interval-ps", str(p["traj_interval_ps"]), "--temperature-k", str(p["temperature_k"]),
                   "--seed", str(p["seed"]), "--report-interval-ps", "2",
                   "--checkpoint-s", str(p.get("checkpoint_s", 10))]
        try:
            last_ps = _run_md_streaming(cmd)
        except KeyboardInterrupt:
            interrupted = True
        if not state_out.exists():
            # nothing reached even one checkpoint — let it be reassigned from the cursor.
            if interrupted:
                raise KeyboardInterrupt
            return {"stop": True, "msg": "the simulation produced no checkpoint"}
        produced["state.xml"] = state_out
        if kind == "segment" and seg_out.exists():
            produced["seg.pdb"] = seg_out
    elif kind in ("dock", "build"):
        produced = run_chunks.execute_unit(kind, p, local, scratch)
    else:
        return {"stop": True, "msg": f"unknown work kind {kind!r}"}

    resp = submit_results(w["results_url"], w["lease"], token, produced, time.time() - t0)
    paths = resp.get("paths") or {}

    # Report to /checkpoint: produced artifacts (by id) for the DAG, plus the cursor
    # state + how far an MD segment actually got. The store registers the artifacts
    # and advances the job; a partial segment advances partially (nothing wasted).
    body: dict = {"lease": w["lease"]}
    out_aids = {aid: paths[name] for name, aid in (w.get("outputs") or {}).items() if name in paths}
    if out_aids:
        body["outputs"] = out_aids
    if kind in ("equilibrate", "segment"):
        body["state"] = paths.get("state.xml")
    if kind == "segment":
        from_ps = float(p["from_ps"])
        body["ps_reached"] = from_ps + (last_ps if last_ps is not None else float(p["step_ps"]))
    requests.post(base.rstrip("/") + "/api/screen/v1/checkpoint", json=body, timeout=60)

    secs = time.time() - t0
    _say(f"    {'stopped early — ' if interrupted else ''}advanced {chunk.get('molecule') or kind} "
         f"· sent back ✓ ({_fmt(secs)})")
    if interrupted:
        raise KeyboardInterrupt
    return {"stop": False, "chunk_id": chunk["id"], "seconds": secs, "lease": None}


# --- one work → run → submit cycle ------------------------------------------

def run_once(base: str, token: str | None, molecules: str | None,
             done: str | None = None, kinds: str | None = None) -> dict:
    """Pull → run → submit one chunk. Returns {stop|chunk_id, seconds, lease, ...}
    and prints what's happening as it goes."""
    w = get_work(base, token, molecules, done, kinds)
    if w.get("status") in ("idle", "busy"):
        return {"stop": True, "msg": w.get("message") or w.get("status"),
                "scope": w.get("scope")}
    if "chunk" not in w:
        return {"stop": True, "msg": w.get("error") or "no work available"}
    if not w.get("results_url"):
        return {"stop": True, "msg": "this site has no result broker configured yet"}

    chunk = w["chunk"]
    _say(f"  ▶ {chunk.get('label') or chunk.get('kind')}")
    safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", str(chunk["id"]))
    scratch = Path(tempfile.mkdtemp(prefix=f"contrib_{safe_id}_"))
    return run_unit(base, w, token, scratch)


# --- propose a molecule to test ---------------------------------------------

def submit_molecule(name: str, smiles: str, *, base: str = DEFAULT_BASE_URL,
                    group: str | None = None, evidence: str | None = None,
                    route: str | None = None, feasibility: str | None = None,
                    refs: str | list[str] | None = None) -> dict:
    """Propose a molecule for the screen to dock and test. `name` and `smiles` are
    required (the structure is what gets docked); the rest help it score and show
    up well. The server validates the SMILES and stores it for the next prep pass.
    Returns the stored record, or {} on a rejected/failed submission."""
    name, smiles = (name or "").strip(), (smiles or "").strip()
    if not name or not smiles:
        _say("A name and a SMILES string are both required to submit a molecule.")
        return {}

    payload: dict = {"name": name, "smiles": smiles}
    if group:
        payload["group"] = group
    if evidence:
        payload["evidence"] = evidence
    delivery = {k: v for k, v in (("route", route), ("feasibility", feasibility)) if v}
    if delivery:
        payload["delivery"] = delivery
    if refs:
        items = refs if isinstance(refs, list) else str(refs).split(";")
        cleaned = [s.strip() for s in items if s and s.strip()]
        if cleaned:
            payload["refs"] = cleaned

    try:
        r = requests.post(base.rstrip("/") + "/api/screen/v1/molecules",
                          json=payload, timeout=60)
    except requests.RequestException as e:  # noqa: BLE001
        _say(f"Couldn't reach the screen to submit ({e}).")
        return {}

    if r.status_code == 200:
        rec = r.json()
        _say(f"Submitted “{rec.get('name')}” (id {rec.get('id')}). It's queued for "
             "docking and testing — its results will appear on the results page.")
        return rec
    try:
        msg = r.json().get("error") or r.text
    except Exception:  # noqa: BLE001
        msg = r.text
    _say(f"Not accepted ({r.status_code}): {msg}")
    return {}


# --- the contributor's whole session ----------------------------------------

def contribute(base: str = DEFAULT_BASE_URL, minutes: float | None = None,
               token: str | None = None, molecules: str | None = None) -> None:
    """Run simulations, reporting progress. Runs until the contributor presses
    stop (or work runs out); pass `minutes` for a fixed time budget instead.
    Always ends with a clear message, never an opaque loop.

    `token` credits every run to the signed-in contributor — get it from
    {base}/screen. The personal notebook / launcher injects it via the
    ``ASYN_CONTRIB_TOKEN`` env var (picked up below). With no token at all, runs
    are anonymous — they still count, just uncredited. `molecules` restricts work
    to a comma-separated preference list (decreasing interest); omit it to run
    whatever the screen needs next.
    """
    # The MD engine's own per-line output is deliberately not surfaced: it's
    # low-level (atom counts, energies, scratch file names) and obscures what the
    # contributor is actually doing. The chunk label says that in plain words and
    # the heartbeat dots show it's alive; a failing chunk still reports a log tail
    # (see run_chunks._run), so nothing diagnostic is lost.
    run_chunks._emit = lambda *_: None
    _say(_gpu_line())
    # Tell the server only the unit kinds this machine can run, so it never hands
    # us docking/prep work a pip-only client can't do (which would just bounce).
    kinds = ",".join(runnable_kinds())
    if "dock" not in kinds and "build" not in kinds:
        _say("This machine runs the simulations; docking and prep need extra tools, so they're left to others.")
    token = token or os.environ.get("ASYN_CONTRIB_TOKEN")
    if not token:
        _say("Running anonymously (no token) — your runs count but aren't "
             "credited to you. Get a token at the site's /screen page for credit.\n")

    budget = None if minutes is None else max(1.0, float(minutes)) * 60.0
    start = time.time()
    done_count, compute = 0, 0.0
    last_lease = None
    _say("Running until you stop it (interrupt the cell / Ctrl-C) — an "
         "unfinished task is simply reassigned, so nothing is wasted.\n"
         if budget is None else
         f"Running for up to {int(minutes)} min. Stop whenever you like — an "
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
        elapsed = time.time() - start
        if budget is not None and budget - elapsed <= 0:
            _say("\nTime budget reached — wrapping up.")
            break
        _say(f"[{_fmt(elapsed)} in · {_fmt(budget - elapsed)} left]"
             if budget is not None else f"[{_fmt(elapsed)} in]")
        try:
            r = run_once(base, token, molecules, done=last_lease, kinds=kinds)
        except KeyboardInterrupt:
            stop_hb.set()
            _say("\nStopped. The current task will be reassigned.")
            return
        except Exception as e:  # noqa: BLE001
            _say(f"    task failed ({e}) — it will be reassigned to someone else\n")
            last_lease = None
            continue
        if r.get("stop") and molecules and r.get("scope") == "molecule":
            _say("\nNothing queued for your molecule pick right now — running "
                 "whatever the screen needs next so your GPU isn't idle. Your "
                 "picks run as the screen reaches them.\n")
            molecules = None
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
    ap.add_argument("--base-url", default=DEFAULT_BASE_URL,
                    help="site root, e.g. https://<site> (the API is at <base>/api/screen/v1)")
    ap.add_argument("--minutes", type=float, default=None,
                    help="optional fixed time budget; default runs until stopped or out of work")
    ap.add_argument("--token", default=None,
                    help="email-verified identity token from <base>/screen; omit to run anonymously")
    ap.add_argument("--molecules", default=None,
                    help="comma-separated molecule ids in decreasing priority; omit to trust the queue")
    args = ap.parse_args()
    contribute(args.base_url, minutes=args.minutes, token=args.token, molecules=args.molecules)


if __name__ == "__main__":
    main()
