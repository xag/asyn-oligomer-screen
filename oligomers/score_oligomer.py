"""Score a relaxed oligomer PDB against the Stage 2 deposited-anchor calibration.

Reads a relaxed PDB (output of `md_relax.py --apo-pdb ... --collapse-ps ...`),
runs `features.py` on each chain, averages across chains, and reports:

  - per-feature mean values (averaged over all chains by default)
  - z-scores against `results/anchor_features.csv` mean/SD
  - the weighted activity score under the Stage 2 classifier
  - a β-structure gate: structures with no detectable β-SASA and no
    NAC-β activity are capped at 0 regardless of disordered-tail or
    contact-density signal (they are not β-rich oligomers)
  - a small ranking table showing where the oligomer falls relative to
    the deposited inert + graded-active anchors

Default chain behaviour — all chains averaged:
    The auto-inner single-chain default was sensitive to which chain the
    OBC2 collapse happened to expose; a lucky tail conformation on the
    most-buried chain could inflate the score by 2× (see STATUS.md
    "s777 outlier"). The new default averages features over every protein
    chain in the structure before scoring, giving a more representative
    and reproducible number.

Usage:
    python score_oligomer.py results/oligomers/fusco_parallel_3mer_core70-88_relaxed.pdb
    python score_oligomer.py <pdb> --chain B          (single chain)
    python score_oligomer.py <pdb> --chains A,B,C     (explicit list, same as default for A/B/C)
    python score_oligomer.py <pdb> --auto-inner        (old default: most-buried chain)
    python score_oligomer.py <pdb> --no-beta-gate      (report raw score, bypass gate)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE))

from features import FEATURES, ordered_core_full_ids
from classifier import WEIGHTS

ANCHOR_FEATURES_CSV = PIPELINE / "results" / "anchor_features.csv"

# β-structure gate thresholds.
# A structure scoring below BOTH floors on the two β-features is not a
# β-rich NAC oligomer and receives activity = min(raw_activity, 0).
# Set just above floating-point zero to catch genuine "no-β" cases while
# tolerating single-residue DSSP fluctuations.
BETA_SASA_FLOOR = 0.5   # Å² per residue
NAC_FLOOR       = 0.5   # Å² per residue


def features_on_chain(structure, chain_ids: list[str]) -> dict[str, float]:
    core_mask = ordered_core_full_ids(structure)
    out: dict[str, float] = {}
    for name, fn in FEATURES.items():
        try:
            out[name] = float(fn(structure, chain_ids=chain_ids, core_mask=core_mask))
        except Exception as exc:
            print(f"  feature {name} failed: {exc}", file=sys.stderr)
            out[name] = float("nan")
    return out


def score(feats: dict, anchor_df: pd.DataFrame) -> tuple[float, dict[str, float], dict[str, float]]:
    zs: dict[str, float] = {}
    contribs: dict[str, float] = {}
    activity = 0.0
    for col, w in WEIGHTS.items():
        mu = anchor_df[col].mean()
        sd = anchor_df[col].std(ddof=0)
        z = 0.0 if sd == 0 else (feats[col] - mu) / sd
        zs[col] = float(z)
        contribs[col] = float(z * w)
        activity += z * w
    return float(activity), zs, contribs


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("pdb", type=Path)

    chain_grp = p.add_mutually_exclusive_group()
    chain_grp.add_argument("--chain",      type=str, default=None,
                           help="score a single named chain")
    chain_grp.add_argument("--chains",     type=str, default=None,
                           help="comma-separated chain list; features are accumulated "
                                "jointly over the listed chains then scored once")
    chain_grp.add_argument("--auto-inner", action="store_true",
                           help="pick the most buried chain (old default; sensitive to "
                                "single-chain tail artifacts — use all-chain mean instead)")

    p.add_argument("--no-beta-gate", action="store_true",
                   help="report raw activity even when β-structure is absent "
                        "(bypasses the 0-cap applied to coil controls)")
    args = p.parse_args()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(args.pdb.stem, str(args.pdb))
    model     = next(iter(structure))
    anchor_df = pd.read_csv(ANCHOR_FEATURES_CSV)

    # ── Resolve which chains to score ─────────────────────────────────────
    if args.chains:
        per_chain_list = [[c.strip() for c in args.chains.split(",")]]
        chain_label    = "+".join(per_chain_list[0])
        score_mode     = "joint"
    elif args.chain:
        per_chain_list = [[args.chain]]
        chain_label    = args.chain
        score_mode     = "single"
    elif args.auto_inner:
        from assembly import inner_chain_ids
        cids           = inner_chain_ids(structure, top_k=1)
        per_chain_list = [cids]
        chain_label    = cids[0] + " (auto-inner)"
        score_mode     = "single"
    else:
        # Default: each chain scored individually, then the feature vectors
        # are averaged before the final z-score / activity computation.
        all_cids       = [ch.id for ch in model]
        per_chain_list = [[cid] for cid in all_cids]
        chain_label    = f"all {len(all_cids)} chains (per-chain mean)"
        score_mode     = "mean"

    print(f"PDB:   {args.pdb}")
    print(f"chain: {chain_label}")

    # ── Compute features ──────────────────────────────────────────────────
    all_feats: list[dict] = []
    for cids in per_chain_list:
        all_feats.append(features_on_chain(structure, cids))

    # Average feature vectors across all scored chains.
    feat_keys = list(all_feats[0].keys())
    mean_feats = {k: float(np.mean([f[k] for f in all_feats])) for k in feat_keys}

    print("\nfeatures (mean over scored chains):")
    for k, v in mean_feats.items():
        print(f"  {k:35} {v:+.4f}")

    # ── Score ──────────────────────────────────────────────────────────────
    raw_activity, zs, contribs = score(mean_feats, anchor_df)

    print("\nz-scores against deposited-anchor calibration:")
    for k in WEIGHTS:
        print(f"  {k:35} z={zs[k]:+.3f}  w={WEIGHTS[k]:+.2f}  contrib={contribs[k]:+.3f}")

    # ── β-structure gate ──────────────────────────────────────────────────
    beta_absent = (
        mean_feats["exposed_hydrophobic_beta_sasa"] < BETA_SASA_FLOOR
        and mean_feats["nac_active_score"] < NAC_FLOOR
    )
    if beta_absent and not args.no_beta_gate:
        activity = min(raw_activity, 0.0)
        print(f"\n  [beta-gate] β-structure absent "
              f"(ehbs={mean_feats['exposed_hydrophobic_beta_sasa']:.3f} < {BETA_SASA_FLOOR}, "
              f"nac={mean_feats['nac_active_score']:.3f} < {NAC_FLOOR}); "
              f"raw={raw_activity:+.3f} capped at {activity:+.3f}")
    else:
        activity = raw_activity

    print(f"\nactivity score: {activity:+.3f}")

    # Per-chain breakdown when scoring more than one chain independently.
    if score_mode == "mean" and len(per_chain_list) > 1:
        acts = []
        for cids, f in zip(per_chain_list, all_feats):
            a, _, _ = score(f, anchor_df)
            acts.append((cids[0], a))
        print("per-chain raw activities: " +
              "  ".join(f"{cid}={a:+.3f}" for cid, a in acts))

    # ── Ranking table ──────────────────────────────────────────────────────
    anchor_scores_csv = PIPELINE / "results" / "anchor_scores.csv"
    if anchor_scores_csv.exists():
        df  = pd.read_csv(anchor_scores_csv).sort_values("activity", ascending=False)
        own = {"pdb_id": args.pdb.stem, "label": "(oligomer hypothesis)", "activity": activity}
        rows = df.to_dict(orient="records") + [own]
        rows.sort(key=lambda r: -r["activity"])
        print("\nranking vs deposited anchors:")
        for r in rows:
            marker = "<<<" if r is own else "   "
            cls    = r.get("label") or ""
            print(f"  {marker} {r['pdb_id']:>34}  {cls:25s}  {r['activity']:+.3f}")


if __name__ == "__main__":
    main()
