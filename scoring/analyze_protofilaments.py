"""Test the within-class hypothesis from STATUS.md (2026-05-26):

  "single-protofilament fibrils score higher than paired-protofilament
   fibrils within each class"

Two protofilament-count notions are tested independently:
- n_protofilaments          — literature-curated, biological count
                              (what the polymorph actually is)
- n_protofilaments_deposited — geometric count from the chains in the
                              biological assembly file (what the
                              feature pipeline actually sees)

These can differ: several cryo-EM α-syn entries deposit only one
protofilament in their REMARK 350 biological assembly even though the
fibril is paired (the second protofilament is implicit by symmetry).
That distinction matters: the *deposited* count drives features; the
*biological* count is the ground truth we ultimately want to predict.

Run `python validate.py` first to refresh the scored CSV."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def _report(df: pd.DataFrame, pf_col: str) -> None:
    print(f"\n--- correlation against {pf_col} ---")
    for label in ["inert", "graded-active"]:
        sub = df[df.label == label].copy()
        # Monomer 1XQ8 has no meaningful protofilament count — drop it
        # from the inert correlation only.
        if label == "inert":
            sub = sub[sub.pdb_id != "1XQ8"]
        if len(sub) < 2 or sub[pf_col].nunique() < 2:
            print(f"  {label:15} n={len(sub):2}  rho=undefined (single class or constant pf)")
            continue
        rho = spearman(sub[pf_col].to_numpy(), sub["activity"].to_numpy())
        print(
            f"  {label:15} n={len(sub):2}  rho={rho:+.3f}  "
            f"(negative -> more protofilaments -> lower activity, supporting the burial hypothesis)"
        )

    print(f"  group means by {pf_col}:")
    for label in df["label"].unique():
        sub = df[df.label == label]
        if label == "inert":
            sub = sub[sub.pdb_id != "1XQ8"]
        for pf, g in sub.groupby(pf_col):
            mean_a = g["activity"].mean()
            ids = ",".join(g["pdb_id"].tolist())
            print(f"    {label:15} pf={pf}  n={len(g):2}  mean(activity)={mean_a:+.3f}  [{ids}]")


def main() -> None:
    scored_path = RESULTS_DIR / "anchor_scores.csv"
    df = pd.read_csv(scored_path)
    for needed in ("n_protofilaments", "n_protofilaments_deposited"):
        if needed not in df.columns:
            raise SystemExit(f"{needed} column missing — rerun `python validate.py`")

    print("=== anchor protofilament counts (lit vs deposited) ===")
    cols = [
        "pdb_id", "label", "n_chains",
        "n_protofilaments", "n_protofilaments_deposited", "activity",
    ]
    print(
        df[cols]
        .sort_values(["label", "activity"], ascending=[True, False])
        .to_string(index=False)
    )

    _report(df, "n_protofilaments")
    _report(df, "n_protofilaments_deposited")


if __name__ == "__main__":
    main()
