#!/usr/bin/env python3
"""How much MD must be simulated to read a ligand effect in the dwell channel
beyond chance — grounded in the conformational space the channel actually samples.

The dwell observable (screen/dwell_time.py) is, per replica,
  dwell = fraction of trajectory frames with beta-core RMSD < TOXIC_RMSD_MAX
          AND inter-chain contact Jaccard > TOXIC_JACCARD_MIN.
The contrast we want is dwell(complex) - dwell(apo): negative = destabiliser,
positive = stabiliser/anti-target.

For "how much data" to be a number, three properties of the sampled space have to
hold; the diagnostics (EXPERIMENTS.md E1/E2) say two of them do not:

1. SHAPE OF THE NOISE. The reference is not metastable: the beta-core RMSD drifts
   at v ~ 0.8 A/ns with no plateau in 2 ns (E2). So a single replica's dwell is a
   first-passage measurement (when does RMSD cross the cutoff), not an equilibrium
   occupancy. For drift-diffusion dX = v dt + sqrt(2D) dW, first passage to L has
   mean T = L/v and CV^2 = 2D/(L v) — a coefficient of variation of order 1. The
   E1 apo spread (mean 0.20, replicas 0.01-0.56) gives CV ~ 0.9 directly, and
   pins D. The per-replica noise is therefore ~ the size of the mean and cannot be
   reduced by anything except more replicas (averaging) — it is intrinsic.

2. DYNAMIC RANGE. apo dwell ~ 0.20 means the system is already ~80% out of the
   basin within 2 ns. The destabiliser direction (downward) has a ceiling of only
   0.20, while the per-replica SD is ~0.18 — so the *best possible* destabiliser
   signal is about one replica-SD. The observable is floored.

3. SIGN. Any bound ligand mechanically restrains the core (occupancy raises dwell),
   so a true destabiliser is pushed toward looking like a stabiliser — a systematic
   bias on top of the noise (E1: every arm shifted positive).

Given (1), the data requirement is the standard two-arm test
  n_per_arm = (z_alpha + z_power)^2 * 2 * sigma^2 / delta^2
with sigma the per-replica dwell SD and delta the dwell shift to detect; alpha
matches the channel's 95% CI (dwell_time.bootstrap_dwell_shift). But (2)+(3) cap
delta well below sigma in the direction of interest, so on THIS observable no
finite n gives a clean beyond-chance destabiliser call. The budget only becomes
finite after the observable is fixed: establish a metastable basin (apo dwell -> ~1,
full 0..1 range, low CV) or switch to the first-passage RATE. This script reports
both the intrinsic noise and what each target delta would then cost.

CPU-only, stdlib. Defaults are the recorded E1/E2 numbers.
"""
from __future__ import annotations

import argparse
import math
from statistics import NormalDist

N = NormalDist()
# Expected range / sigma for a normal sample (d2 control-chart constant), n=10.
_D2 = {5: 2.326, 8: 2.847, 10: 3.078, 12: 3.258, 15: 3.472, 20: 3.735}


def sigma_from_range(lo: float, hi: float, n: int) -> float:
    """Estimate per-replica SD from the observed min/max of n replicas."""
    return (hi - lo) / _D2.get(n, 3.078)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # Recorded apo dwell (EXPERIMENTS.md E1): mean 0.20, per-replica 0.01-0.56, n=10.
    ap.add_argument("--apo-mean", type=float, default=0.20)
    ap.add_argument("--apo-lo", type=float, default=0.01)
    ap.add_argument("--apo-hi", type=float, default=0.56)
    ap.add_argument("--n-rep", type=int, default=10)
    # Recorded landscape (E2): drift v, start and cutoff RMSD.
    ap.add_argument("--drift", type=float, default=0.8, help="beta-core RMSD drift, A/ns")
    ap.add_argument("--rmsd-start", type=float, default=1.38, help="A")
    ap.add_argument("--rmsd-cutoff", type=float, default=3.0, help="A (TOXIC_RMSD_MAX)")
    # Per-replica MD cost and ensemble width.
    ap.add_argument("--ns-per-rep", type=float, default=2.1, help="prod+equil ns per replica")
    ap.add_argument("--n-shapes", type=int, default=9)
    # Test stringency: channel uses a 95% CI (two-sided alpha 0.05) + 80% power.
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--power", type=float, default=0.80)
    args = ap.parse_args()

    sigma = sigma_from_range(args.apo_lo, args.apo_hi, args.n_rep)
    cv = sigma / args.apo_mean
    L = args.rmsd_cutoff - args.rmsd_start
    T_fp = L / args.drift                       # mean first-passage time, ns
    D = cv * cv * L * args.drift / 2.0          # implied RMSD-coordinate diffusion, A^2/ns
    ceiling = args.apo_mean                     # destabiliser can only push dwell 0.20 -> 0
    snr1 = ceiling / sigma                      # best-case signal / per-replica noise

    z = N.inv_cdf(1 - args.alpha / 2) + N.inv_cdf(args.power)
    coef = z * z * 2.0 * sigma * sigma          # n_per_arm = coef / delta^2

    print("=== the sampled space (from E1/E2) ===")
    print(f"  per-replica dwell SD sigma  = {sigma:.3f}   (apo mean {args.apo_mean:.2f}, "
          f"range {args.apo_lo:.2f}-{args.apo_hi:.2f}, n={args.n_rep})")
    print(f"  coefficient of variation CV = {cv:.2f}  -> per-replica noise ~ the mean itself")
    print(f"  first-passage picture: cross L={L:.2f} A at v={args.drift} A/ns -> mean T={T_fp:.2f} ns")
    print(f"    consistent CV^2=2D/(L v) gives RMSD-coordinate D ~ {D:.2f} A^2/ns (no barrier)")
    print(f"  destabiliser signal ceiling = {ceiling:.2f} (floored baseline); "
          f"best-case SNR/replica = {snr1:.2f}")

    print(f"\n=== IF the observable is fixed (metastable basin or rate readout), "
          f"replicas/arm to detect a shift delta ===")
    print(f"  (two-arm test, {int((1-args.alpha)*100)}% CI + {int(args.power*100)}% power; "
          f"z-sum={z:.2f})")
    print(f"  {'delta':>7} {'n/arm':>7} {'ns/pair':>9} {'ns x{0} shapes'.format(args.n_shapes):>16}")
    for delta in (0.20, 0.10, 0.05, 0.02):
        n_arm = math.ceil(coef / (delta * delta))
        ns_pair = n_arm * 2 * args.ns_per_rep
        print(f"  {delta:>7.2f} {n_arm:>7d} {ns_pair:>9.0f} {ns_pair * args.n_shapes:>16.0f}")

    print(f"\n  reference: E1 ran {args.n_rep} replicas x 2 ns/arm "
          f"(~{args.n_rep * args.ns_per_rep * 4:.0f} ns total over 4 arms).")
    print("  On the current (floored, non-metastable, occupancy-confounded) observable,")
    print("  the destabiliser delta is capped near sigma, so no finite n clears chance.")


if __name__ == "__main__":
    main()
