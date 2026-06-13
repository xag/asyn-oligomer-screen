"""Tests for the pure plumbing helpers in hf_store ingest.

No network, no HF — run directly with the pipeline venv:

    .venv/bin/python tests/test_ingest_plumbing.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))

import hf_store as h  # noqa: E402


def test_submission_key_sha_prefers_pdb():
    outputs = {
        "pair/s1/state_1.xml": {"sha256": "STATE", "file": "a"},
        "pair/s1/seg_0.pdb": {"sha256": "PDB", "file": "b"},
    }
    assert h._submission_key_sha(outputs) == "PDB"


def test_submission_key_sha_falls_back_to_first():
    outputs = {"pair/system.xml": {"sha256": "SYS", "file": "a"}}
    assert h._submission_key_sha(outputs) == "SYS"


def test_parse_reputations_defaults_missing_fields():
    reps = h._parse_reputations({
        "alice": {"agreed": 5, "allowlist_bonus": 1.0},
        "bob": {},
    })
    assert reps["alice"].agreed == 5
    assert reps["alice"].allowlist_bonus == 1.0
    assert reps["alice"].spot_fail == 0
    assert reps["bob"].agreed == 0
    # Fresh contributor sits at the weight floor; an allowlisted one rises above it.
    import contrib_gate as g
    assert g.weight_for(reps["bob"]) == g.WEIGHT_FLOOR
    assert g.weight_for(reps["alice"]) > g.WEIGHT_FLOOR


def test_parse_reputations_empty():
    assert h._parse_reputations({}) == {}
    assert h._parse_reputations(None) == {}


def test_load_spotchecks_absent_is_empty():
    assert h._load_spotchecks(None) == {}


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
