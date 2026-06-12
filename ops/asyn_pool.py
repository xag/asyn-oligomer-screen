#!/usr/bin/env python3
"""Pool all accumulated blocks and apply the pre-registered sequential decision.

For each block (= one conformer dir under results/blocks/<conformer>/):
  * score every test/decoy replica trajectory -> primary per-trajectory statistic
    jaccard_decay_per_ns = -(slope of contact Jaccard vs time);
  * matched within-block contrast C = mean_seed[ decay(test) - decay(decoy) ],
    paired by seed where possible (cancels structural + thermal variance).
Across blocks, bootstrap the mean of the per-block contrasts and apply prereg.py:
stop-H1 / stop-null / continue. CPU-only; reads whatever is on disk.

Env: ASYN_REPO (required).
"""
from __future__ import annotations

import glob
import os
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(os.environ["ASYN_REPO"]).resolve()
sys.path.insert(0, str(REPO / "screen"))
sys.path.insert(0, str(Path(__file__).resolve().parent))  # for prereg.py
os.chdir(REPO)

import prereg  # noqa: E402
import dwell_time as dt  # noqa: E402
from dwell_time import score_trajectory  # noqa: E402
from shape_metrics import load_pdb  # noqa: E402

LO, HI = dt.CHUNK_RANGE
BLOCKS_ROOT = REPO / "results" / "blocks"
DT_NS = prereg.TRAJ_INTERVAL_PS / 1000.0
_SEED_RE = re.compile(r"_s(\d+)\.pdb$")


def _decay(traj: Path, ref) -> float:
    """Primary statistic: contact-Jaccard decay rate (1/ns); higher = faster loss."""
    sc = score_trajectory(traj, ref)
    j = np.asarray(sc["per_frame"]["contact_jaccard"], dtype=float)
    if j.size < 3:
        return float("nan")
    t = np.arange(j.size) * DT_NS
    return -float(np.polyfit(t, j, 1)[0])


def _by_seed(block: Path, label: str, ref) -> dict[int, float]:
    out: dict[int, float] = {}
    for p in sorted(glob.glob(str(block / f"{label}_s[0-9]*.pdb"))):
        if p.endswith("_final.pdb"):
            continue
        m = _SEED_RE.search(Path(p).name)
        if not m:
            continue
        try:
            out[int(m.group(1))] = _decay(Path(p), ref)
        except Exception as e:  # noqa: BLE001
            print(f"  ! {Path(p).name}: {type(e).__name__}: {e}")
    return out


def block_contrast(block: Path) -> float | None:
    """Matched within-block contrast: mean over shared seeds of decay(test)-decay(decoy)."""
    ref_pdb = block / f"apo_core{LO}-{HI}.pdb"
    if not ref_pdb.exists():
        return None
    ref = load_pdb(ref_pdb)
    test = _by_seed(block, prereg.TEST, ref)
    decoy = _by_seed(block, prereg.DECOY, ref)
    shared = sorted(set(test) & set(decoy))
    pairs = [test[s] - decoy[s] for s in shared
             if np.isfinite(test[s]) and np.isfinite(decoy[s])]
    if pairs:
        return float(np.mean(pairs))
    # fall back to group means if seeds don't line up
    tv = [v for v in test.values() if np.isfinite(v)]
    dv = [v for v in decoy.values() if np.isfinite(v)]
    if tv and dv:
        return float(np.mean(tv) - np.mean(dv))
    return None


def main() -> None:
    print(f"pre-registration locked {prereg.LOCKED}: primary={prereg.PRIMARY}, "
          f"H1: {prereg.TEST} decay > {prereg.DECOY} (one-sided)")
    if not BLOCKS_ROOT.exists():
        print("no blocks yet — run ops/block.ps1 -Action run first."); return

    blocks = [d for d in sorted(BLOCKS_ROOT.iterdir()) if d.is_dir()]
    contrasts: list[float] = []
    print("\n=== per-block matched contrast (test - decoy, jaccard_decay 1/ns; >0 favours H1) ===")
    for b in blocks:
        c = block_contrast(b)
        if c is None:
            print(f"   {b.name:48} (incomplete)")
            continue
        contrasts.append(c)
        print(f"   {b.name:48} {c:+.4f}")

    n = len(contrasts)
    if n == 0:
        print("\nno complete blocks to pool yet."); return

    arr = np.array(contrasts, dtype=float)
    rng = np.random.default_rng(0)
    boot = arr[rng.integers(0, n, size=(prereg.N_BOOT, n))].mean(axis=1)
    mean = float(arr.mean())
    lo = float(np.quantile(boot, 0.025))
    hi = float(np.quantile(boot, 0.975))
    p_h1 = float(np.mean(boot > 0))

    print(f"\n=== pooled over {n} block(s) ===")
    print(f"   mean contrast = {mean:+.4f} /ns   95% CI [{lo:+.4f}, {hi:+.4f}]")
    print(f"   P({prereg.TEST} destabilises faster than {prereg.DECOY}) = {p_h1:.3f}")

    nullish = (lo > -prereg.NULL_HALFWIDTH) and (hi < prereg.NULL_HALFWIDTH)
    if p_h1 >= prereg.POSTERIOR_THRESHOLD:
        verdict = (f"STOP — H1 SUPPORTED. The channel separates the known destabiliser "
                   f"from the decoy; a powered sweep is justified.")
    elif n >= prereg.MIN_BLOCKS and nullish:
        verdict = ("STOP — PRACTICAL NULL. CI is tight around zero; the channel does not "
                   "see the effect on this model. Do not sweep; the anti-target axis stands.")
    elif n >= prereg.MAX_BLOCKS:
        verdict = (f"STOP — BUDGET REACHED ({n} blocks) still inconclusive. Treat as unresolved.")
    else:
        need = max(prereg.MIN_BLOCKS - n, 1)
        verdict = (f"CONTINUE — collect more blocks (have {n}; >= {prereg.MIN_BLOCKS} before a "
                   f"null call, up to {prereg.MAX_BLOCKS}). ~{need} more night-block(s) minimum.")
    print(f"\n   DECISION: {verdict}")


if __name__ == "__main__":
    main()
