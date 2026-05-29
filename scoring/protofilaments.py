"""Protofilament count from chain geometry.

Stage 2 (2026-05-26) suggested that within-class activity-score ordering
might track protofilament count: single-protofilament fibrils (more
solvent-exposed NAC, fewer inter-protofilament contacts) score higher
than paired-protofilament fibrils. To test that, we need a per-anchor
protofilament count. This module computes it geometrically from the
biological assembly, so we don't have to hand-curate from polymorph
papers.

Algorithm.
- Compute one centroid per chain.
- Fibril axis = principal direction of the *short* pair-vectors between
  chain centroids. In a fibril each chain's nearest neighbours along
  the stack sit ~4.8 Å away along the fibril axis; chains in other
  protofilaments are 25-40 Å away laterally. Weighting toward the
  shortest pairs picks the within-PF rung step rather than the inter-PF
  separation or the chain's own intrinsic long axis (which is what
  all-Cα PCA would pick on ssNMR deposits like 2N0A, where each chain
  is ~150 Å long but the chain centroids stack over ~60 Å).
- Project chain centroids onto the plane perpendicular to the fibril
  axis.
- Connected components on the graph where lateral distance between two
  centroids is below `lateral_threshold`. Number of components = number
  of protofilaments.

The lateral threshold (12 Å default) separates chains stacked within
one protofilament (lateral spread < ~5 Å, even with helical wobble)
from chains in distinct protofilaments (inter-protofilament gap ~25-40
Å at the interface for α-syn polymorphs).

Returns 1 for single-chain assemblies (monomer 1XQ8) — semantically
"one stack", with `n_chains` distinguishing monomer from fibril.
"""
from __future__ import annotations

import numpy as np
from Bio.PDB.Structure import Structure

from features import THREE_TO_ONE, _first_model


def _chain_centroids(structure: Structure) -> tuple[list[str], np.ndarray, list[np.ndarray]]:
    """Return (chain_ids, chain_centroids, per_chain_ca_coords). Skips
    chains with fewer than 3 Cα atoms — short ligand-like fragments
    would otherwise pollute the principal-axis fit."""
    per_chain_ca: list = []
    chain_ids: list[str] = []
    chain_centroids: list = []
    for chain in _first_model(structure):
        ca_coords = []
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if residue.get_resname() not in THREE_TO_ONE:
                continue
            if "CA" not in {atom.get_name() for atom in residue}:
                continue
            ca_coords.append(residue["CA"].coord)
        if len(ca_coords) < 3:
            continue
        arr = np.asarray(ca_coords)
        chain_ids.append(chain.id)
        chain_centroids.append(arr.mean(axis=0))
        per_chain_ca.append(arr)
    if not per_chain_ca:
        return [], np.empty((0, 3)), []
    return chain_ids, np.asarray(chain_centroids), per_chain_ca


def _fibril_axis_from_short_pairs(centroids: np.ndarray) -> np.ndarray:
    """Estimate the fibril axis from the shortest centroid-to-centroid
    vectors. The shortest n_chains pairs are taken as a proxy for
    nearest-neighbour rungs within a single protofilament; their
    common direction is the fibril growth axis."""
    n = len(centroids)
    diffs = centroids[:, None, :] - centroids[None, :, :]
    dists = np.linalg.norm(diffs, axis=-1)
    iu, ju = np.triu_indices(n, k=1)
    lengths = dists[iu, ju]
    vectors = diffs[iu, ju]
    # Take the shortest ~n pairs (each chain contributes ~1 nearest
    # neighbour). Clamp to ≥3 for the smallest assemblies.
    k = min(len(lengths), max(n, 3))
    order = np.argsort(lengths)[:k]
    short = vectors[order]
    # Second-moment tensor (sign-invariant); top eigenvector = mean
    # direction of the chosen pair vectors.
    units = short / np.linalg.norm(short, axis=-1, keepdims=True).clip(min=1e-9)
    M = units.T @ units
    eigvals, eigvecs = np.linalg.eigh(M)
    return eigvecs[:, -1]


def _mean_chain_lateral_spread(
    per_chain_ca: list, fibril_axis: np.ndarray,
) -> float:
    """Mean across chains of each chain's lateral spread (max Cα distance
    from chain centroid, measured in the plane perpendicular to the
    fibril axis). Sets the natural lateral scale of one protofilament:
    chains in the same PF should sit within ~half this distance of each
    other laterally; chains in distinct PFs sit well beyond it."""
    spreads = []
    for arr in per_chain_ca:
        ax = arr @ fibril_axis
        perp_arr = arr - np.outer(ax, fibril_axis)
        centroid_perp = perp_arr.mean(axis=0)
        spreads.append(float(np.linalg.norm(perp_arr - centroid_perp, axis=1).max()))
    return float(np.mean(spreads)) if spreads else 0.0


def count_protofilaments(
    structure: Structure,
    lateral_threshold: float = 12.0,
    coincident_centroid_spread: float = 1.0,
) -> int:
    """Number of distinct protofilament stacks in the assembly.

    Single-chain assemblies return 1. Multi-chain assemblies cluster
    chain centroids by lateral distance perpendicular to the fibril
    axis (recovered from the shortest centroid-to-centroid pair
    vectors).

    `lateral_threshold` (Å) — minimum lateral distance for two chain
    centroids to count as belonging to distinct protofilaments. The
    effective threshold is `max(lateral_threshold, 0.5 *
    mean_per_chain_lateral_spread)`. The adaptive component handles
    ssNMR ensembles (notably 2N0A) where each chain conformer is
    larger and noisier than a typical cryo-EM rung, so cryo-EM-tuned
    thresholds would over-segment a single PF.

    `coincident_centroid_spread` (Å) — if chain centroids are nearly
    coincident (mean distance from their mean below this), return 1
    rather than trying to extract a fibril axis from noise. Only
    triggers on degenerate inputs."""
    chain_ids, centroids, per_chain_ca = _chain_centroids(structure)
    n = len(chain_ids)
    if n <= 1:
        return n

    centered_c = centroids - centroids.mean(axis=0)
    spread = float(np.linalg.norm(centered_c, axis=1).mean())
    if spread < coincident_centroid_spread:
        return 1

    fibril_axis = _fibril_axis_from_short_pairs(centroids)
    lateral_scale = _mean_chain_lateral_spread(per_chain_ca, fibril_axis)
    effective_threshold = max(lateral_threshold, 0.5 * lateral_scale)

    # Project chain centroids onto the plane perpendicular to fibril_axis.
    axial = centroids @ fibril_axis
    perp = centroids - np.outer(axial, fibril_axis)

    # Connected components: edge iff lateral distance < threshold.
    diffs = perp[:, None, :] - perp[None, :, :]
    lat_dist = np.linalg.norm(diffs, axis=-1)
    adj = lat_dist < effective_threshold

    visited = [False] * n
    components = 0
    for i in range(n):
        if visited[i]:
            continue
        components += 1
        stack = [i]
        while stack:
            k = stack.pop()
            if visited[k]:
                continue
            visited[k] = True
            for j in range(n):
                if not visited[j] and adj[k, j]:
                    stack.append(j)
    return components
