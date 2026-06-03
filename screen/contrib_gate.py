"""Acceptance gate for crowdsourced chunk submissions (#43).

Pure, IO-free decision logic, unit-tested independently of the Hugging Face
round-trip in ``hf_store.py``. It replaces the old count-based ``>=min-agree``
rule with the Sybil-resistant model designed for anonymous, email-keyed
contribution through the ``health`` front door.

Why count is not enough. Cross-platform MD is not bit-identical, so the
original gate compared a continuous observable (the dwell fraction) and
accepted a chunk on ``>=N`` agreeing *uploads*. The moment contribution is
anonymous, one actor can mint N uploads, so "N agreeing uploads" stops meaning
"N independent results." This module moves the gate from *counting uploads* to
*weighing distinct, reputation-bearing contributors*:

  * de-duplication      identical bytes (a replayed file) collapse to one, and
                        at most one submission per contributor pseudonym counts
                        — so the same output re-posted under many pseudonyms is
                        a single vote, not a quorum.
  * distinct-identity   the consensus cluster is measured over distinct
    quorum              pseudonyms, never over upload count.
  * reputation weight   each pseudonym contributes ``weight(rep)`` in
                        ``[WEIGHT_FLOOR, WEIGHT_CEIL]``; a chunk accepts when the
                        cluster's summed weight reaches ``W``. Trusted
                        contributors clear quorum with few bodies; a crowd of
                        fresh identities clears it only in numbers.
  * spot-check anchor   when the coordinator re-ran the chunk itself, its own
                        observable must fall inside the accepted cluster — the
                        ground truth that bootstraps reputation and catches a
                        cluster of agreeing fabrications.

The seed is part of the chunk (deterministic per replica/segment — see
``run_chunks._segment_seed``), so honest re-runs of one chunk share a seed and
differ only by float non-determinism; they are *not* distinguished here. The
unit of independence is the contributor pseudonym, supplied by the dispatcher.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# Per-contributor weight bounds. A never-before-seen pseudonym sits at the
# floor (a real but quiet vote); a long track record approaches the ceiling.
WEIGHT_FLOOR = 0.1
WEIGHT_CEIL = 1.0
# Smoothing constant in the reputation ratio — how much corroborated history it
# takes to move the weight appreciably off the floor.
REP_SMOOTHING = 2.0


def largest_agreeing(values: list[float], tol: float) -> list[int]:
    """Indices of the largest subset whose values all lie within ``tol``
    (max-min <= tol) — the consensus cluster; anything outside is an outlier.

    Pure + unit-testable. Moved here from ``hf_store`` so the gate owns the
    consensus primitive and ``hf_store`` imports it back.
    """
    if not values:
        return []
    order = sorted(range(len(values)), key=lambda i: values[i])
    best: list[int] = []
    lo = 0
    for hi in range(len(order)):
        while values[order[hi]] - values[order[lo]] > tol:
            lo += 1
        if hi - lo + 1 > len(best):
            best = order[lo:hi + 1]
    return sorted(best)


@dataclass
class Submission:
    """One contribution as the gate sees it, after integrity (SHA-256) checks.

    ``dwell`` is the scored observable for ``segment`` chunks (the dwell
    fraction); it is ``None`` for non-observable kinds (build / equilibrate),
    which are gated on distinct-pseudonym weight alone. ``ts`` breaks ties when
    one pseudonym has submitted more than once (the latest wins).
    """
    pseudonym: str
    sha256: str
    dwell: float | None = None
    ts: float = 0.0


@dataclass
class Reputation:
    """Per-pseudonym history, maintained by ``health`` and read by the gate.

    ``allowlist_bonus`` lets an identity already trusted in ``health`` (an
    allowlisted email, or a good summary / digest record) start above the floor
    — the cross-credit between contribution kinds.
    """
    agreed: int = 0
    outlier: int = 0
    spot_pass: int = 0
    spot_fail: int = 0
    allowlist_bonus: float = 0.0


def weight_for(rep: Reputation | None) -> float:
    """Map a contributor's history to a quorum weight in ``[FLOOR, CEIL]``.

    A single confirmed spot-check failure zeroes the weight: a caught fabricator
    contributes nothing to any quorum (and is separately dispatch-banned by
    ``health``). Otherwise weight rises with corroborated agreements and falls
    with outliers; a fresh identity (no history, no bonus) sits at the floor.
    """
    if rep is None:
        return WEIGHT_FLOOR
    if rep.spot_fail > 0:
        return 0.0
    num = rep.agreed - 2 * rep.outlier + rep.allowlist_bonus
    denom = rep.agreed + rep.outlier + REP_SMOOTHING
    score = 0.0 if denom <= 0 else num / denom
    score = max(0.0, min(1.0, score))
    return WEIGHT_FLOOR + (WEIGHT_CEIL - WEIGHT_FLOOR) * score


def dedup(subs: list[Submission]) -> list[Submission]:
    """Collapse to one vote per *distinct result* and per *distinct pseudonym*.

    1. Identical bytes (same SHA-256) collapse to a single submission — a
       replayed file, even under many pseudonyms, is one vote. This is what
       makes "submit the same output N times" worthless.
    2. Of what remains, keep at most one submission per pseudonym (latest by
       ``ts``) — so a single identity cannot occupy several cluster slots.
    """
    by_sha: dict[str, Submission] = {}
    for s in subs:
        by_sha.setdefault(s.sha256, s)  # first occurrence of these exact bytes
    by_pseudo: dict[str, Submission] = {}
    for s in by_sha.values():
        cur = by_pseudo.get(s.pseudonym)
        if cur is None or s.ts >= cur.ts:
            by_pseudo[s.pseudonym] = s
    return list(by_pseudo.values())


@dataclass
class Decision:
    """Outcome for one chunk. ``status`` is one of:

      accept          quorum weight reached (and spot-check, if any, passed);
                      ``representative`` is the SHA-256 to archive.
      awaiting        not enough distinct-pseudonym weight yet — wait for more.
      spotcheck_fail  a cluster formed but the coordinator's own re-run sits
                      outside it: reject and flag the members.
    """
    status: str
    weight: float
    cluster: list[str] = field(default_factory=list)   # pseudonyms accepted
    representative: str | None = None                   # SHA-256 to archive
    note: str = ""


def _cluster_weight(members: list[Submission], reps: dict[str, Reputation]) -> float:
    return sum(weight_for(reps.get(m.pseudonym)) for m in members)


def decide(
    subs: list[Submission],
    reputations: dict[str, Reputation],
    *,
    tol: float,
    quorum_weight: float,
    coordinator_dwell: float | None = None,
    observable: bool = True,
) -> Decision:
    """Decide one chunk's fate from its submissions and contributor history.

    ``observable=False`` is the build / equilibrate path: deterministic *setup*,
    not a measurement, so there is nothing to corroborate — accept on the first
    integrity-valid result. ``coordinator_dwell`` is supplied only when this
    chunk was sampled for a spot-check re-run.
    """
    deduped = dedup(subs)
    reps = reputations

    if not observable:
        # Setup chunks (build / equilibrate) construct and warm up the system;
        # the output is preparation, not a value to agree on. A bad setup is
        # caught downstream — the chunks that consume it crash or land outside
        # the segment consensus — so waiting for a distinct-pseudonym quorum here
        # buys no safety and, for a small or early crowd, is never reached: the
        # chunk stays pending and is re-dispatched forever. So accept the first
        # integrity-valid upload. Quorum is kept only for the observable segment
        # path below, where the dwell value genuinely needs corroboration.
        if not deduped:
            return Decision("awaiting", 0.0, note="no integrity-valid submissions yet")
        rep = max(deduped, key=lambda s: weight_for(reps.get(s.pseudonym)))
        wsum = _cluster_weight(deduped, reps)
        return Decision("accept", wsum, [s.pseudonym for s in deduped], rep.sha256,
                        note=f"setup chunk accepted on first valid result "
                             f"({len(deduped)} contributor(s))")

    scored = [s for s in deduped if s.dwell is not None]
    if not scored:
        return Decision("awaiting", 0.0, note="no scored submissions yet")
    vals = [s.dwell for s in scored]
    idx = largest_agreeing(vals, tol)
    cluster = [scored[i] for i in idx]
    wsum = _cluster_weight(cluster, reps)
    if wsum < quorum_weight:
        return Decision("awaiting", wsum,
                        cluster=[s.pseudonym for s in cluster],
                        note=f"cluster {len(cluster)}/{len(scored)} "
                             f"weight {wsum:.2f}/{quorum_weight:.2f} "
                             f"(dwell {[round(v, 3) for v in vals]})")

    mean = sum(s.dwell for s in cluster) / len(cluster)
    if coordinator_dwell is not None and abs(coordinator_dwell - mean) > tol:
        return Decision("spotcheck_fail", wsum,
                        cluster=[s.pseudonym for s in cluster],
                        note=f"coordinator dwell {coordinator_dwell:.3f} outside "
                             f"cluster mean {mean:.3f} +/- {tol}")

    representative = min(cluster, key=lambda s: abs(s.dwell - mean))
    return Decision("accept", wsum, [s.pseudonym for s in cluster], representative.sha256,
                    note=f"{len(cluster)}/{len(scored)} agree, weight {wsum:.2f} "
                         f"(dwell {[round(v, 3) for v in vals]})")
