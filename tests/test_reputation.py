"""Tests for the reputation-folding aggregation.

Pure logic, no network. Run:  .venv/bin/python tests/test_reputation.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))

import reputation as R  # noqa: E402
import contrib_gate as g  # noqa: E402


def test_counts_outcomes_per_pseudonym():
    recs = [
        {"pseudonym": "a", "outcome": "agreed"},
        {"pseudonym": "a", "outcome": "agreed"},
        {"pseudonym": "a", "outcome": "outlier"},
        {"pseudonym": "b", "outcome": "spot_fail"},
    ]
    reps = R.aggregate(recs)
    assert reps["a"]["agreed"] == 2 and reps["a"]["outlier"] == 1
    assert reps["b"]["spot_fail"] == 1


def test_spot_pass_counts_as_agreement_and_is_tracked():
    reps = R.aggregate([{"pseudonym": "a", "outcome": "spot_pass"}])
    assert reps["a"]["agreed"] == 1
    assert reps["a"]["spot_pass"] == 1


def test_incremental_folding_adds_to_base():
    base = {"a": {"agreed": 5, "outlier": 1, "spot_pass": 0, "spot_fail": 0}}
    reps = R.aggregate([{"pseudonym": "a", "outcome": "agreed"}], base=base)
    assert reps["a"]["agreed"] == 6
    assert reps["a"]["outlier"] == 1


def test_bonus_is_stamped_and_defaults_zero():
    reps = R.aggregate(
        [{"pseudonym": "a", "outcome": "agreed"}],
        bonus={"a": 1.0, "vip": 2.0},
    )
    assert reps["a"]["allowlist_bonus"] == 1.0
    assert reps["vip"]["allowlist_bonus"] == 2.0   # bonus-only contributor created
    # default zero for anyone without a bonus
    reps2 = R.aggregate([{"pseudonym": "b", "outcome": "agreed"}])
    assert reps2["b"]["allowlist_bonus"] == 0.0


def test_output_is_consumable_by_the_gate_weight():
    # The folded record must drive contrib_gate.weight_for end to end.
    reps = R.aggregate([{"pseudonym": "a", "outcome": "agreed"}] * 20, bonus={"a": 1.0})
    rep = g.Reputation(**reps["a"])
    assert g.weight_for(rep) > g.WEIGHT_FLOOR
    # a spot_fail anywhere zeroes it
    bad = R.aggregate([{"pseudonym": "x", "outcome": "spot_fail"}])
    assert g.weight_for(g.Reputation(**bad["x"])) == 0.0


def test_jsonl_parsing_skips_blank_and_bad_lines():
    text = '{"pseudonym":"a","outcome":"agreed"}\n\nnot json\n{"pseudonym":"b","outcome":"outlier"}\n'
    rows = R._read_jsonl(text)
    assert len(rows) == 2


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
