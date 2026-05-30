"""Run the anchor structures through Stage 2 features + the weighted score
and plot separation. This is the project's first go/no-go gate.

By default, structures are expanded to their REMARK 350 BIOMOLECULE 1
biological assembly and feature accumulation runs on the most buried
chain of that assembly. This addresses the asymmetric-unit-size confound
where structures with smaller deposited units (e.g. 8A9L) get inflated
SASA-based features.

Usage:
    python validate.py            # assembly + inner chain (default)
    python validate.py --au       # legacy: raw asymmetric unit, all chains
    python validate.py --compare  # print both side-by-side
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from anchors import ANCHORS, load
from assembly import inner_chain_ids, load_assembly
from classifier import WEIGHTS, score_table
from features import FEATURES, ordered_core_full_ids
from protofilaments import count_protofilaments

RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"

PALETTE = {
    "inert":         "#2f6f96",
    "graded-active": "#e08134",
    "active":        "#c1432d",
}

PRETTY_FEATURE = {
    "exposed_hydrophobic_beta_sasa": "exposed hydrophobic β-SASA",
    "membrane_insertion_propensity": "membrane-insertion propensity",
    "nac_active_score":              "NAC β-accessibility",
    "contact_density":               "contact density (signed −)",
    "disordered_hydrophobic_exposure": "disordered hydrophobic exposure",
}


def _pairwise_auc(active_scores, inert_scores) -> float:
    a = np.asarray(active_scores, dtype=float)
    i = np.asarray(inert_scores, dtype=float)
    if a.size == 0 or i.size == 0:
        return float("nan")
    wins = (a[:, None] > i[None, :]).sum() + 0.5 * (a[:, None] == i[None, :]).sum()
    return float(wins) / (a.size * i.size)


def _stagger_label_offsets(values, gap: float = 0.35):
    """Vertical pt offsets for crowded scatter labels along x.
    Adjacent points within `gap` cycle through above/below levels."""
    order = sorted(range(len(values)), key=lambda k: values[k])
    cycle = (24, -24, 42, -42)
    offsets = [0] * len(values)
    cursor = 0
    last_v = float("-inf")
    for rank, k in enumerate(order):
        v = values[k]
        if v - last_v < gap:
            cursor = (cursor + 1) % len(cycle)
        else:
            cursor = 0
        offsets[k] = cycle[cursor]
        last_v = v
    return offsets


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
    core-only structures.
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
    n_cols = 2
    n_rows = (len(feats) + n_cols) // n_cols  # reserve last cell for legend/summary
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(13.0, 4.6 * n_rows), dpi=160)
    axes_flat = axes.flatten() if n_rows * n_cols > 1 else [axes]

    inert_mask = df.label == "inert"
    active_mask = df.label.isin(["active", "graded-active"])

    for ax, feat in zip(axes_flat, feats):
        sub = df[df[feat].notna()].sort_values(feat).reset_index(drop=True)
        y_pos = np.arange(len(sub))
        point_colors = [PALETTE.get(c, "#888") for c in sub.label]

        # Faint per-class mean lines
        for cls in ("inert", "graded-active", "active"):
            m_vals = df.loc[df.label == cls, feat].dropna()
            if not m_vals.empty:
                ax.axvline(m_vals.mean(), color=PALETTE[cls],
                           linewidth=1.2, alpha=0.35, zorder=1)

        ax.scatter(sub[feat], y_pos, c=point_colors, s=78, alpha=0.95,
                   edgecolors="white", linewidths=1.1, zorder=3)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(sub.pdb_id, fontsize=8.5)
        ax.set_title(PRETTY_FEATURE.get(feat, feat),
                     fontsize=10.5, fontweight="semibold", pad=6, loc="left")
        ax.grid(True, axis="x", alpha=0.22, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", length=3, color="#999")

        auc_f = _pairwise_auc(
            df.loc[active_mask, feat].dropna().tolist(),
            df.loc[inert_mask, feat].dropna().tolist(),
        )
        if not np.isnan(auc_f):
            ax.text(0.985, 0.04, f"AUC  {auc_f:.2f}", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=9, color="#333",
                    bbox=dict(boxstyle="round,pad=0.32", facecolor="white",
                              edgecolor="#d0d0d0", linewidth=0.8))

    # Legend / summary in the trailing empty cell(s)
    for empty_ax in axes_flat[len(feats):]:
        empty_ax.axis("off")
    if len(axes_flat) > len(feats):
        legend_ax = axes_flat[len(feats)]
        used = [c for c in ("inert", "graded-active", "active") if c in df.label.values]
        handles = [plt.Line2D([0], [0], marker="o", color="w",
                              markerfacecolor=PALETTE[c], markeredgecolor="white",
                              markeredgewidth=1.1, markersize=12,
                              label=c.replace("-", " "))
                   for c in used]
        leg = legend_ax.legend(handles=handles, loc="upper center", fontsize=11,
                               frameon=False, title="anchor class",
                               title_fontsize=11, borderpad=1)
        leg.get_title().set_fontweight("semibold")
        legend_ax.text(
            0.5, 0.42,
            "points sorted by feature value\nvertical lines: per-class mean\n"
            "AUC: pairwise Mann–Whitney\n(graded-active vs inert)",
            ha="center", va="top", fontsize=9.5, color="#444",
            transform=legend_ax.transAxes,
        )

    fig.suptitle(f"Per-feature anchor profiles{suffix}",
                 fontsize=13.5, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    path = RESULTS_DIR / f"anchor_features{suffix.replace(' ', '_')}.png"
    fig.savefig(path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def plot_activity(df: pd.DataFrame, suffix: str = "") -> Path:
    class_order = [c for c in ("inert", "graded-active", "active")
                   if not df[df.label == c].empty]
    y_base = {label: i for i, label in enumerate(class_order)}

    fig, ax = plt.subplots(figsize=(12.5, 5.4), dpi=160)

    # Distribution underlay
    for cls in class_order:
        sub = df[df.label == cls]
        if sub.empty:
            continue
        y = y_base[cls]
        bp = ax.boxplot(
            sub.activity, vert=False, positions=[y], widths=0.50,
            whis=(0, 100), showfliers=False, patch_artist=True,
            manage_ticks=False, zorder=1,
        )
        for box in bp["boxes"]:
            box.set(facecolor=PALETTE[cls], alpha=0.13,
                    edgecolor=PALETTE[cls], linewidth=1.0)
        for w in bp["whiskers"] + bp["caps"]:
            w.set(color=PALETTE[cls], linewidth=0.9, alpha=0.55)
        for m in bp["medians"]:
            m.set(color=PALETTE[cls], linewidth=2.0)

    # Class-mean marker + label below the row
    for cls in class_order:
        sub = df[df.label == cls]
        if sub.empty:
            continue
        m = sub.activity.mean()
        y = y_base[cls]
        ax.scatter([m], [y - 0.34], marker="v", s=90, color=PALETTE[cls],
                   edgecolor="white", linewidth=1.0, zorder=4)
        ax.text(m, y - 0.46, f"μ = {m:+.2f}", ha="center", va="top",
                fontsize=9, color=PALETTE[cls], fontweight="semibold")

    # Points and staggered labels
    for cls in class_order:
        sub = df[df.label == cls].sort_values("activity").reset_index(drop=True)
        if sub.empty:
            continue
        y_c = y_base[cls]
        jitters = [((i % 2) * 2 - 1) * 0.045 for i in range(len(sub))]
        ys = [y_c + j for j in jitters]
        offsets = _stagger_label_offsets(sub.activity.tolist(), gap=0.40)
        for i, (_, row) in enumerate(sub.iterrows()):
            ax.scatter(row.activity, ys[i], c=PALETTE[cls], s=145,
                       alpha=0.96, edgecolor="white", linewidth=1.4, zorder=5)
            ax.annotate(
                row.pdb_id, (row.activity, ys[i]),
                xytext=(0, offsets[i]), textcoords="offset points",
                ha="center", va="bottom" if offsets[i] > 0 else "top",
                fontsize=9, color="#1f1f1f", zorder=6,
                arrowprops=dict(arrowstyle="-", color="#bbbbbb", linewidth=0.7,
                                shrinkA=2, shrinkB=4)
                if abs(offsets[i]) > 25 else None,
            )

    # Classifier zero
    ax.axvline(0, color="#9a9a9a", linewidth=0.9, linestyle=(0, (4, 4)), zorder=2)

    # AUC callout
    inert_vals = df.loc[df.label == "inert", "activity"].tolist()
    active_vals = df.loc[df.label.isin(["graded-active", "active"]), "activity"].tolist()
    auc = _pairwise_auc(active_vals, inert_vals)
    if not np.isnan(auc):
        ax.text(0.987, 0.965,
                f"pairwise AUC  {auc:.2f}\n"
                f"n = {len(active_vals)} active · {len(inert_vals)} inert",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=10.5, fontweight="semibold", color="#222",
                bbox=dict(boxstyle="round,pad=0.55", facecolor="white",
                          edgecolor="#cccccc", linewidth=1.0))

    ax.set_yticks(list(y_base.values()))
    ax.set_yticklabels([c.replace("-", " ") for c in class_order],
                       fontsize=11.5, fontweight="medium")
    ax.set_xlabel("activity score   (weighted z-sum of five surface-biophysics features)",
                  fontsize=10.5)
    ax.set_title("Anchor separation along the activity axis",
                 fontsize=13.5, fontweight="bold", loc="left", pad=10)
    ax.grid(True, axis="x", alpha=0.22, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", length=3, color="#999")
    ax.set_ylim(-0.95, len(class_order) - 0.15)

    fig.tight_layout()
    path = RESULTS_DIR / f"anchor_activity{suffix.replace(' ', '_')}.png"
    fig.savefig(path, dpi=160, bbox_inches="tight", facecolor="white")
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
