"""Tests for spot-check target selection. Pure; no MD, no network.

Run:  .venv/bin/python tests/test_spotcheck.py
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))

import spotcheck as S  # noqa: E402


def test_returns_all_when_sample_covers_pool():
    cands = ["c3", "c1", "c2"]
    assert S.select_targets(cands, 5, random.Random(0)) == ["c1", "c2", "c3"]
    assert S.select_targets(cands, 0, random.Random(0)) == ["c1", "c2", "c3"]


def test_samples_a_subset_deterministically_with_seed():
    cands = [f"c{i}" for i in range(10)]
    a = S.select_targets(cands, 3, random.Random(42))
    b = S.select_targets(cands, 3, random.Random(42))
    assert a == b and len(a) == 3
    assert set(a).issubset(set(cands))


def test_dedups_and_sorts():
    assert S.select_targets(["c2", "c2", "c1"], 0, random.Random(0)) == ["c1", "c2"]


def test_candidates_picks_only_runnable_segments_with_submissions():
    manifest = {
        "artifacts": {"x/state_0.xml": {"present": True}, "x/sys.xml": {"present": True},
                      "y/state_0.xml": {"present": False}},
        "chunks": [
            {"id": "seg_ready", "kind": "segment", "status": "pending",
             "consumes": ["x/state_0.xml", "x/sys.xml"]},
            {"id": "seg_no_sub", "kind": "segment", "status": "pending",
             "consumes": ["x/state_0.xml"]},
            {"id": "seg_missing_input", "kind": "segment", "status": "pending",
             "consumes": ["y/state_0.xml"]},
            {"id": "seg_done", "kind": "segment", "status": "done",
             "consumes": ["x/state_0.xml"]},
            {"id": "build_ready", "kind": "build", "status": "pending",
             "consumes": ["x/sys.xml"]},
        ],
    }
    submitted = {"seg_ready", "seg_missing_input", "seg_done", "build_ready"}
    assert S._candidates(manifest, submitted) == ["seg_ready"]


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
