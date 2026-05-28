"""Run the anchor structures through Stage 2 features + the weighted score
and plot separation. This is the project's first go/no-go gate.

By default, structures are expanded to their REMARK 350 BIOMOLECULE 1
biological assembly and feature accumulation runs on the most buried
chain of that assembly. This addresses the asymmetric-unit-size confound
where structures with smaller deposited units (e.g. 8A9L) get inflated
SASA-based features (see STATUS.md, problem 1).

Usage:
    python validate.py            # assembly + inner chain (default)
    python validate.py --au       # legacy: raw asymmetric unit, all chains
    python validate.py --compare  # print both side-by-side
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from anchors import ANCHORS, load
from assembly import inner_chain_ids, load_assembly
from classifier import WEIGHTS, score_table
from features import FEATURES, ordered_core_full_ids
from protofilaments import count_protofilaments

RESULTS_DIR = Path(__file__).parent / "results"


def _chain_ids_for(mode: str, structure) -> list[str] | None:
    if mode == "au":
        return None
    if mode == "assembly_all":
        return None
    return inner_chain_ids(structure, top_k=1)


def compute_features(mode: str = "assembly_inner", use_core_mask: bool = True) -> pd.DataFrame:
    """Compute features under one of three modes.

    - assembly_inner: REMARK 350 assembly, features on the most buried chain
    - assembly_all:   REMARK 350 assembly, features on all chains
    - au:             raw asymmetric unit, all chains (the original behaviour)

    `use_core_mask`: when True, every per-residue feature accumulates only
    over the structurally ordered core (Cα with ≥6 non-sequential
    neighbours within 8 Å). Equalises NMR full-length and cryo-EM
    core-only structures (STATUS.md problem 1).
    """
    rows = []
    for anchor in ANCHORS:
        try:
            if mode == "au":
                structure = load(anchor.pdb_id)
            else:
                structure = load_assembly(anchor.pdb_id)
        except Exception as exc:
            print(f"  {anchor.pdb_id} fetch/parse failed: {exc}")
            continue
        n_chains = sum(1 for _ in next(iter(structure)))
        chain_ids = _chain_ids_for(mode, structure)
        chain_tag = ",".join(chain_ids) if chain_ids else "all"
        core_mask = ordered_core_full_ids(structure) if use_core_mask else None
        if use_core_mask and chain_ids is not None:
            core_in_chain = sum(
                1 for fid in (core_mask or set()) if fid[2] in chain_ids
            )
        else:
            core_in_chain = len(core_mask) if core_mask else 0
        mask_tag = (
            f"core={core_in_chain}" if core_mask is not None
            else "core=fallback(none)" if use_core_mask
            else "core=off"
        )
        try:
            n_pf_geom = count_protofilaments(structure)
        except Exception as exc:
            print(f"    protofilament count failed: {exc}")
            n_pf_geom = 0
        print(
            f"  {anchor.pdb_id} ({anchor.label}): {n_chains} chains, "
            f"using {chain_tag}, {mask_tag}, pf_lit={anchor.n_protofilaments}, pf_dep={n_pf_geom}"
        )
        row = {
            "pdb_id": anchor.pdb_id,
            "label": anchor.label,
            "description": anchor.description,
            "n_chains": n_chains,
            "n_protofilaments": anchor.n_protofilaments,
            "n_protofilaments_deposited": n_pf_geom,
            "chains_used": chain_tag,
            "core_in_chain": core_in_chain,
        }
        for name, fn in FEATURES.items():
            try:
                row[name] = fn(structure, chain_ids=chain_ids, core_mask=core_mask)
            except Exception as exc:
                print(f"    feature {name} failed: {exc}")
                row[name] = float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def plot_features(df: pd.DataFrame, suffix: str = "") -> Path:
    feats = list(FEATURES.keys())
    fig, axes = plt.subplots(1, len(feats), figsize=(3.2 * len(feats), 4.2))
    colors = {"inert": "tab:blue", "active": "tab:red", "graded-active": "tab:orange"}
    for ax, feat in zip(axes, feats):
        for label, color in colors.items():
            sub = df[df.label == label]
            if sub.empty:
                continue
            ax.scatter([0] * len(sub), sub[feat], c=color, label=label, s=80, alpha=0.7)
            for _, row in sub.iterrows():
                ax.annotate(
                    row.pdb_id,
                    (0, row[feat]),
                    xytext=(6, 0),
                    textcoords="offset points",
                    fontsize=8,
                )
        ax.set_title(feat, fontsize=9)
        ax.set_xticks([])
    axes[-1].legend(loc="best", fontsize=8)
    fig.suptitle(f"Stage 2 features across anchors{suffix}", fontsize=11)
    fig.tight_layout()
    path = RESULTS_DIR / f"anchor_features{suffix.replace(' ', '_')}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_activity(df: pd.DataFrame, suffix: str = "") -> Path:
    fig, ax = plt.subplots(figsize=(8, 4.2))
    colors = {"inert": "tab:blue", "active": "tab:red", "graded-active": "tab:orange"}
    for label, color in colors.items():
        sub = df[df.label == label]
        if sub.empty:
            continue
        ax.scatter(sub["activity"], [label] * len(sub), c=color, s=120, alpha=0.7)
        for _, row in sub.iterrows():
            ax.annotate(row.pdb_id, (row["activity"], label), xytext=(4, 6),
                        textcoords="offset points", fontsize=8)
    ax.set_xlabel("activity score (weighted z-sum)")
    ax.axvline(0, color="grey", linewidth=0.5)
    ax.set_title(f"Anchor separation along the activity axis{suffix}")
    fig.tight_layout()
    path = RESULTS_DIR / f"anchor_activity{suffix.replace(' ', '_')}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _print_scored(scored: pd.DataFrame, title: str) -> None:
    print(f"\n=== {title} ===")
    cols = [
        "pdb_id", "label", "n_chains", "n_protofilaments", "n_protofilaments_deposited",
        "chains_used", "core_in_chain", "activity",
    ] + list(FEATURES.keys())
    available = [c for c in cols if c in scored.columns]
    print(scored[available].to_string(index=False))


def run_mode(mode: str, suffix: str, use_core_mask: bool = True) -> pd.DataFrame:
    mask_tag = "core_on" if use_core_mask else "core_off"
    print(f"\ncomputing features for {len(ANCHORS)} anchors [mode={mode}, {mask_tag}]:")
    raw = compute_features(mode=mode, use_core_mask=use_core_mask)
    if raw.empty:
        print("no anchors loaded — aborting")
        return raw
    raw.to_csv(RESULTS_DIR / f"anchor_features{suffix}.csv", index=False)
    scored = score_table(raw)
    scored.to_csv(RESULTS_DIR / f"anchor_scores{suffix}.csv", index=False)
    plot_features(raw, suffix=suffix)
    plot_activity(scored, suffix=suffix)
    _print_scored(scored, f"scored anchors (mode={mode}, {mask_tag}, descending activity)")
    return scored


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    args = sys.argv[1:]
    use_core_mask = "--no-core-mask" not in args

    print("weights:")
    for name, w in WEIGHTS.items():
        print(f"  {name:42} {w:+.2f}")
    print(f"core mask: {'ON' if use_core_mask else 'OFF'}")

    if "--compare" in args:
        run_mode("au", "_au", use_core_mask=use_core_mask)
        run_mode("assembly_all", "_asm_all", use_core_mask=use_core_mask)
        run_mode("assembly_inner", "", use_core_mask=use_core_mask)
    elif "--au" in args:
        run_mode("au", "_au", use_core_mask=use_core_mask)
    elif "--assembly-all" in args:
        run_mode("assembly_all", "_asm_all", use_core_mask=use_core_mask)
    else:
        run_mode("assembly_inner", "", use_core_mask=use_core_mask)


if __name__ == "__main__":
    main()
