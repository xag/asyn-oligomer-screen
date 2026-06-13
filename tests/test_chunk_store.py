"""Chunk store DAG + lease semantics (local scope).

These pin the contract that makes the workload distributable and resumable, with
no MD/OpenMM involved (pure store logic, runs in the bare pip venv):

  * the runnable predicate *is* the DAG — ``build`` first, then both seeds'
    ``equilibrate`` chunks in parallel (no ordering between independent replicas),
    then each replica's ``segment`` chain in order;
  * ``push`` ingesting a chunk's outputs flips its dependents runnable;
  * a chunk is never leased to two workers at once, and an expired lease is
    reclaimed so a dead worker doesn't strand a chunk;
  * ``push`` rejects an output set that doesn't match the chunk's ``produces``.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "screen"))
import chunk_store as store  # noqa: E402


def _touch(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture()
def exp(tmp_path, monkeypatch):
    """A 1-pair / 2-seed / 2-segment experiment mirroring the dwell DAG."""
    monkeypatch.setattr(store, "EXPERIMENTS_DIR", tmp_path / "experiments")
    pair = "shapeX__apo"
    core, sysx, solv = f"{pair}/core.pdb", f"{pair}/system.xml", f"{pair}/solvated.pdb"
    seeds = [1000, 1001]
    chunks = [{"id": f"build__{pair}", "kind": "build",
               "consumes": [core], "produces": [sysx, solv], "params": {}}]
    for sd in seeds:
        s0 = f"{pair}/s{sd}/state_0.xml"
        chunks.append({"id": f"equil__{pair}__s{sd}", "kind": "equilibrate",
                       "consumes": [sysx, solv], "produces": [s0], "params": {"seed": sd}})
        for i in range(2):
            chunks.append({"id": f"seg__{pair}__s{sd}__{i}", "kind": "segment",
                           "consumes": [sysx, solv, f"{pair}/s{sd}/state_{i}.xml"],
                           "produces": [f"{pair}/s{sd}/state_{i+1}.xml", f"{pair}/s{sd}/seg_{i}.pdb"],
                           "params": {"index": i}})
    core_file = _touch(tmp_path / "core.pdb", "ATOM core\n")
    store.create_experiment("exp1", {"shape": "shapeX", "seeds": seeds}, chunks,
                            {core: core_file})
    return tmp_path, "exp1", pair, seeds


def _outfile(tmp_path, name, text="x") -> Path:
    p = tmp_path / "out" / name
    p.parent.mkdir(parents=True, exist_ok=True)
    return _touch(p, text)


def test_build_is_the_only_first_runnable(exp):
    tmp_path, exp_id, pair, seeds = exp
    ch = store.pull(exp_id, "w1")
    assert ch["id"] == f"build__{pair}"
    # build is now leased; nothing else consumes a present artifact yet.
    assert store.pull(exp_id, "w2") is None


def test_build_push_unlocks_both_seeds_in_parallel(exp):
    tmp_path, exp_id, pair, seeds = exp
    build = store.pull(exp_id, "w1")
    sysx, solv = build["produces"]
    store.push(exp_id, build["id"],
               {sysx: _outfile(tmp_path, "system.xml"), solv: _outfile(tmp_path, "solvated.pdb")}, 1.0)

    # Both equilibrate chunks (independent seeds) are now runnable — order-free.
    a = store.pull(exp_id, "w1")
    b = store.pull(exp_id, "w2")
    got = {a["id"], b["id"]}
    assert got == {f"equil__{pair}__s{seeds[0]}", f"equil__{pair}__s{seeds[1]}"}
    # No third equilibrate/segment is runnable yet (segments need state_0).
    assert store.pull(exp_id, "w3") is None


def test_segment_chain_is_ordered_within_a_replica(exp):
    tmp_path, exp_id, pair, seeds = exp
    sd = seeds[0]
    build = store.pull(exp_id, "w1")
    sysx, solv = build["produces"]
    store.push(exp_id, build["id"],
               {sysx: _outfile(tmp_path, "sys"), solv: _outfile(tmp_path, "solv")}, 1.0)

    # Drain the other seed's equilibrate so it can't be confused for this chain.
    while True:
        ch = store.pull(exp_id, "w9", lease_seconds=9999)
        if ch is None:
            break
        if ch["id"] == f"equil__{pair}__s{sd}":
            target = ch
        # leave others leased (parked) so only our seed's chain is exercised
    state0 = target["produces"][0]
    store.push(exp_id, target["id"], {state0: _outfile(tmp_path, "state0")}, 1.0)

    # seg_0 runnable; seg_1 NOT (needs state_1 from seg_0).
    seg0 = store.pull(exp_id, "w1")
    assert seg0["id"] == f"seg__{pair}__s{sd}__0"
    assert store.pull(exp_id, "w1") is None
    st1, _seg = seg0["produces"]
    store.push(exp_id, seg0["id"],
               {st1: _outfile(tmp_path, "s1"), seg0["produces"][1]: _outfile(tmp_path, "f0")}, 1.0)
    seg1 = store.pull(exp_id, "w1")
    assert seg1["id"] == f"seg__{pair}__s{sd}__1"


def test_expired_lease_is_reclaimed(exp):
    tmp_path, exp_id, pair, seeds = exp
    ch = store.pull(exp_id, "dead-worker", lease_seconds=0.0)  # expires immediately
    assert ch["id"] == f"build__{pair}"
    time.sleep(0.01)
    again = store.pull(exp_id, "live-worker")
    assert again is not None and again["id"] == ch["id"]  # reclaimed, not stranded


def test_push_rejects_wrong_outputs(exp):
    tmp_path, exp_id, pair, seeds = exp
    build = store.pull(exp_id, "w1")
    with pytest.raises(ValueError, match="missing outputs"):
        store.push(exp_id, build["id"], {build["produces"][0]: _outfile(tmp_path, "only-one")}, 1.0)


def test_full_drain_reports_complete(exp):
    tmp_path, exp_id, pair, seeds = exp
    n = 0
    while True:
        ch = store.pull(exp_id, "w1")
        if ch is None:
            break
        outs = {aid: _outfile(tmp_path, aid.replace("/", "_")) for aid in ch["produces"]}
        store.push(exp_id, ch["id"], outs, 0.1)
        n += 1
    s = store.status_summary(exp_id)
    assert s["complete"] is True
    assert s["counts"]["done"] == s["n_chunks"] == n
    # 1 build + 2 equilibrate + 2 seeds × 2 segments = 7 chunks.
    assert n == 7
