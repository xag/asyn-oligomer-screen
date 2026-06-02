"""Fold ingest outcomes into reputations.json (#43).

``hf_store ingest`` writes per-contributor outcomes to ``outcomes/<ts>.jsonl``
(one of agreed | spot_pass | outlier | spot_fail). This job accumulates them
into a single ``reputations.json`` that BOTH readers consume: the gate
(``ingest`` → ``contrib_gate.weight_for``) and the dispatcher (health
``/screen``). Reputation is therefore a published artifact in the dataset — one
writer (the coordinator, who holds the token), two readers — so health needs no
Hugging Face token and no Redis reputation store.

Folding is incremental: it adds the new outcome files onto the running totals
already in ``reputations.json`` and then deletes the files it consumed (like
ingest's own inbox cleanup), so a file is never counted twice. The allowlist
bonus is refreshed each run from an optional health endpoint that maps
already-trusted contributors' pseudonyms to a starting weight.

    python screen/reputation.py fold --repo user/asyn-dwell-results \
        [--bonus-url "https://site/screen?action=bonus&secret=..."]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

DATASET = "dataset"
FIELDS = ("agreed", "outlier", "spot_pass", "spot_fail")


def aggregate(records: list[dict], base: dict | None = None, bonus: dict | None = None) -> dict:
    """Pure: fold ``records`` (each {pseudonym, outcome}) onto ``base`` running
    totals, then stamp each contributor's ``allowlist_bonus`` from ``bonus``.

    spot_pass counts as an agreement too (it is an agreement the coordinator
    re-verified); spot_fail is tracked on its own and zeroes weight downstream.
    """
    out: dict[str, dict] = {}
    for pseudo, rec in (base or {}).items():
        out[pseudo] = {f: int((rec or {}).get(f, 0)) for f in FIELDS}

    def slot(p: str) -> dict:
        return out.setdefault(p, {f: 0 for f in FIELDS})

    for r in records:
        p, o = r.get("pseudonym"), r.get("outcome")
        if not p:
            continue
        s = slot(p)
        if o == "agreed":
            s["agreed"] += 1
        elif o == "spot_pass":
            s["agreed"] += 1
            s["spot_pass"] += 1
        elif o == "outlier":
            s["outlier"] += 1
        elif o == "spot_fail":
            s["spot_fail"] += 1

    for p in (bonus or {}):
        slot(p)
    for p, s in out.items():
        s["allowlist_bonus"] = float((bonus or {}).get(p, 0.0))
    return out


def _read_jsonl(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if line:
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return rows


def cmd_fold(args) -> None:
    from huggingface_hub import HfApi, hf_hub_download, CommitOperationAdd, CommitOperationDelete
    api = HfApi(token=args.token)

    files = api.list_repo_files(args.repo, repo_type=DATASET, token=args.token)
    outcome_files = sorted(f for f in files if f.startswith("outcomes/") and f.endswith(".jsonl"))
    if not outcome_files:
        print("no new outcomes to fold", flush=True)
        return

    # Existing running totals (base) — absent on the first fold.
    base = {}
    if "reputations.json" in files:
        p = hf_hub_download(args.repo, "reputations.json", repo_type=DATASET,
                            token=args.token, force_download=True)
        base = json.loads(Path(p).read_text(encoding="utf-8"))

    records = []
    for of in outcome_files:
        p = hf_hub_download(args.repo, of, repo_type=DATASET, token=args.token, force_download=True)
        records += _read_jsonl(Path(p).read_text(encoding="utf-8"))

    bonus = _fetch_bonus(args.bonus_url) if args.bonus_url else _existing_bonus(base)
    reps = aggregate(records, base=base, bonus=bonus)

    ops = [CommitOperationAdd("reputations.json",
                              json.dumps(reps, indent=2, sort_keys=True).encode("utf-8"))]
    ops += [CommitOperationDelete(of) for of in outcome_files]   # consumed — don't double-count
    api.create_commit(args.repo, ops, repo_type=DATASET, token=args.token,
                      commit_message=f"fold {len(outcome_files)} outcome file(s) -> reputations.json")
    print(f"folded {len(records)} outcome(s) from {len(outcome_files)} file(s); "
          f"{len(reps)} contributor(s) tracked.", flush=True)


def _existing_bonus(base: dict) -> dict:
    """Preserve the allowlist bonuses already recorded when no fresh map is given."""
    return {p: float(rec.get("allowlist_bonus", 0.0)) for p, rec in (base or {}).items()
            if rec.get("allowlist_bonus")}


def _fetch_bonus(url: str) -> dict:
    import urllib.request
    try:
        with urllib.request.urlopen(url, timeout=30) as r:   # noqa: S310 — coordinator-supplied URL
            return {k: float(v) for k, v in json.loads(r.read().decode("utf-8")).items()}
    except Exception as e:  # noqa: BLE001
        print(f"  warning: could not fetch bonus map ({e}); keeping existing", flush=True)
        return {}


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    pf = sub.add_parser("fold", help="accumulate outcomes/*.jsonl into reputations.json")
    pf.add_argument("--repo", required=True, help="HF dataset repo")
    pf.add_argument("--token", default=None, help="HF token (default: cached login / $HF_TOKEN)")
    pf.add_argument("--bonus-url", default=None,
                    help="health endpoint mapping allowlisted pseudonyms -> starting bonus")
    pf.set_defaults(func=cmd_fold)
    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
