"""Tests for the crowdsourced-submission acceptance gate (#43).

Pure logic, no MD and no network — run directly with the pipeline venv:

    .venv/bin/python tests/test_contrib_gate.py

The functions are also named ``test_*`` so pytest discovers them.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))

import contrib_gate as g  # noqa: E402
from contrib_gate import Submission, Reputation  # noqa: E402


# --- the consensus primitive (unchanged behaviour, moved from hf_store) ------

def test_largest_agreeing_picks_consensus_cluster():
    # 0.50/0.52/0.51 agree within 0.05; 0.90 is the outlier.
    assert g.largest_agreeing([0.50, 0.90, 0.52, 0.51], 0.05) == [0, 2, 3]


def test_largest_agreeing_empty():
    assert g.largest_agreeing([], 0.1) == []


# --- de-duplication: the core Sybil defence -----------------------------------

def test_identical_bytes_under_many_pseudonyms_is_one_vote():
    # The lazy attack: one real output, re-posted under five fake identities.
    subs = [Submission(pseudonym=f"sock{i}", sha256="DEAD", dwell=0.5) for i in range(5)]
    assert len(g.dedup(subs)) == 1


def test_one_pseudonym_many_uploads_is_one_vote():
    subs = [
        Submission("alice", sha256="a1", dwell=0.50, ts=1.0),
        Submission("alice", sha256="a2", dwell=0.51, ts=2.0),
        Submission("alice", sha256="a3", dwell=0.52, ts=3.0),
    ]
    deduped = g.dedup(subs)
    assert len(deduped) == 1
    assert deduped[0].sha256 == "a3"  # latest ts wins


# --- distinct-identity quorum, reputation-weighted ----------------------------

def _fresh(n, dwell=0.50):
    # n distinct, never-seen pseudonyms, each a byte-distinct but agreeing run.
    return [Submission(pseudonym=f"new{i}", sha256=f"h{i}", dwell=dwell) for i in range(n)]


def test_fresh_crowd_needs_numbers_to_reach_quorum():
    reps: dict[str, Reputation] = {}  # everyone at the weight floor (0.1)
    # Two fresh contributors = weight 0.2, below a quorum target of 1.0.
    assert g.decide(_fresh(2), reps, tol=0.05, quorum_weight=1.0).status == "awaiting"
    # Ten fresh contributors = weight 1.0, clears it.
    assert g.decide(_fresh(10), reps, tol=0.05, quorum_weight=1.0).status == "accept"


def test_trusted_contributors_clear_quorum_with_fewer_bodies():
    reps = {
        "vet1": Reputation(agreed=20, allowlist_bonus=1.0),
        "vet2": Reputation(agreed=20, allowlist_bonus=1.0),
    }
    subs = [
        Submission("vet1", sha256="x1", dwell=0.50),
        Submission("vet2", sha256="x2", dwell=0.51),
    ]
    d = g.decide(subs, reps, tol=0.05, quorum_weight=1.0)
    assert d.status == "accept"          # two high-weight votes (~0.85 each) > 1.0
    assert set(d.cluster) == {"vet1", "vet2"}


def test_outlier_excluded_from_cluster_and_representative_is_central():
    reps = {p: Reputation(agreed=20, allowlist_bonus=1.0) for p in ("a", "b", "c")}
    subs = [
        Submission("a", sha256="ha", dwell=0.50),
        Submission("b", sha256="hb", dwell=0.52),
        Submission("c", sha256="hc", dwell=0.95),  # outlier
    ]
    d = g.decide(subs, reps, tol=0.05, quorum_weight=1.0)
    assert d.status == "accept"
    assert "c" not in d.cluster
    assert d.representative in {"ha", "hb"}


def test_spotcheck_failure_rejects_an_agreeing_fabrication():
    # A cluster of agreeing values that the coordinator's own re-run contradicts.
    reps = {p: Reputation(agreed=20, allowlist_bonus=1.0) for p in ("a", "b")}
    subs = [
        Submission("a", sha256="ha", dwell=0.20),
        Submission("b", sha256="hb", dwell=0.21),
    ]
    d = g.decide(subs, reps, tol=0.05, quorum_weight=1.0, coordinator_dwell=0.60)
    assert d.status == "spotcheck_fail"


def test_spotcheck_pass_accepts():
    reps = {p: Reputation(agreed=20, allowlist_bonus=1.0) for p in ("a", "b")}
    subs = [
        Submission("a", sha256="ha", dwell=0.50),
        Submission("b", sha256="hb", dwell=0.51),
    ]
    d = g.decide(subs, reps, tol=0.05, quorum_weight=1.0, coordinator_dwell=0.49)
    assert d.status == "accept"


def test_spot_fail_history_zeroes_weight():
    assert g.weight_for(Reputation(agreed=50, spot_fail=1)) == 0.0
    # A caught fabricator contributes nothing, so a chunk resting only on it waits.
    reps = {"faker": Reputation(agreed=50, spot_fail=1)}
    subs = [Submission("faker", sha256="h", dwell=0.5)]
    assert g.decide(subs, reps, tol=0.05, quorum_weight=0.5).status == "awaiting"


# --- non-observable kinds (build / equilibrate) -------------------------------

def test_non_observable_kind_gated_on_distinct_weight():
    reps = {"vet1": Reputation(agreed=20, allowlist_bonus=1.0),
            "vet2": Reputation(agreed=20, allowlist_bonus=1.0)}
    subs = [Submission("vet1", sha256="s1"), Submission("vet2", sha256="s2")]
    d = g.decide(subs, reps, tol=0.05, quorum_weight=1.0, observable=False)
    assert d.status == "accept"
    assert d.representative in {"s1", "s2"}


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
