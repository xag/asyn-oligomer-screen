"""Tests for the dwell-time chunk DAG (pose-sampling redesign).

Pure structure, no MD/docking — run with the pipeline venv:

    .venv/bin/python tests/test_enumerate_chunks.py

Named test_* so pytest discovers them too.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))

import run_chunks as rc  # noqa: E402


def _spec(**over):
    spec = {
        "shape": "shp", "seeds": [1, 2], "prod_ps": 100.0, "segment_ps": 50.0,
        "equil_ps": 10.0, "traj_interval_ps": 2.0, "temperature_k": 300.0,
        "rect_box": True, "ligands": ["apo", "egcg"],
        "ligand_smiles": {"egcg": "O=C"}, "n_poses": 3,
    }
    spec.update(over)
    return spec


def _by_kind(chunks):
    out = {}
    for c in chunks:
        out.setdefault(c["kind"], []).append(c)
    return out


def test_apo_arm_has_no_dock_and_builds_from_apo_core():
    chunks = rc.enumerate_chunks(_spec(ligands=["apo"]))
    kinds = _by_kind(chunks)
    assert "dock" not in kinds
    build = kinds["build"]
    assert len(build) == 1
    assert build[0]["consumes"] == ["shp__apo/core.pdb"]
    # 2 seeds × (1 equil + ceil(100/50)=2 segments)
    assert len(kinds["equilibrate"]) == 2
    assert len(kinds["segment"]) == 4


def test_complex_docks_an_ensemble_onto_the_shared_receptor():
    chunks = rc.enumerate_chunks(_spec(ligands=["egcg"], n_poses=3))
    kinds = _by_kind(chunks)
    # one dock per ligand, producing all poses; docking is the only setup step
    assert len(kinds["dock"]) == 1
    assert "param" not in kinds
    dock = kinds["dock"][0]
    assert dock["consumes"] == ["shp__apo/core.pdb"]   # docks onto the shared receptor core
    assert dock["params"]["n_poses"] == 3
    assert dock["produces"] == [f"shp__egcg/p{j}/core.pdb" for j in range(3)]


def test_complex_fans_dynamics_over_pose_times_seed():
    chunks = rc.enumerate_chunks(_spec(ligands=["egcg"], n_poses=3))
    kinds = _by_kind(chunks)
    # build per pose; equil per (pose, seed); segment per (pose, seed, index)
    assert len(kinds["build"]) == 3
    assert len(kinds["equilibrate"]) == 3 * 2
    assert len(kinds["segment"]) == 3 * 2 * 2
    # every complex dynamics chunk carries its pose in meta (so scoring can
    # aggregate across poses); seeds/index too.
    for c in kinds["segment"]:
        assert c["meta"]["pose"] in (0, 1, 2)
        assert "seed" in c["meta"] and "index" in c["meta"]


def test_build_consumes_its_pose_core():
    chunks = rc.enumerate_chunks(_spec(ligands=["egcg"], n_poses=2))
    builds = [c for c in chunks if c["kind"] == "build"]
    for j, b in enumerate(sorted(builds, key=lambda c: c["params"]["pose"])):
        assert b["consumes"] == [f"shp__egcg/p{j}/core.pdb"]
        assert b["params"]["is_complex"] is True


def test_chunk_ids_are_unique():
    chunks = rc.enumerate_chunks(_spec())
    ids = [c["id"] for c in chunks]
    assert len(ids) == len(set(ids)), "duplicate chunk ids"


def test_consumed_artifacts_are_produced_or_initial():
    # Every consumed artifact must be produced by some chunk or be the apo core
    # (the single registered initial input). Guards the DAG against dangling deps.
    chunks = rc.enumerate_chunks(_spec())
    produced = {a for c in chunks for a in c["produces"]}
    initial = {"shp__apo/core.pdb"}
    for c in chunks:
        for a in c["consumes"]:
            assert a in produced or a in initial, f"{c['id']} consumes unproduced {a}"


def test_pose_count_scales_chunks():
    five = rc.enumerate_chunks(_spec(ligands=["egcg"], n_poses=5))
    one = rc.enumerate_chunks(_spec(ligands=["egcg"], n_poses=1))
    assert len([c for c in five if c["kind"] == "build"]) == 5
    assert len([c for c in one if c["kind"] == "build"]) == 1
    # one dock chunk regardless of pose count (it produces all poses)
    assert len([c for c in five if c["kind"] == "dock"]) == 1
    assert len([c for c in one if c["kind"] == "dock"]) == 1


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)}/{len(fns)} passed")
