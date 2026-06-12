"""Pre-registration for the dwell-time channel validation. LOCKED 2026-06-12,
before any blocked run was collected, so the analysis cannot be tuned to pass.
Imported by asyn_pool.py and enforced there. Do NOT edit after the first block
is collected — that would void the pre-registration.

Hypothesis (one-sided): a known alpha-synuclein destabiliser (silibinin) drives
the toxic-oligomer inter-chain contacts apart FASTER than an inert decoy
(caffeine). The pilot showed contact-Jaccard is far lower-variance than RMSD,
so it is the primary coordinate.
"""

LOCKED = "2026-06-12"

# Primary per-trajectory statistic: contact-Jaccard decay rate (1/ns) =
# -(least-squares slope of inter-chain contact Jaccard vs time). Larger =>
# contacts lost faster => more destabilising.
PRIMARY = "jaccard_decay_per_ns"
SECONDARY = "rmsd_rise_per_ns"

TEST = "silibinin"   # known destabiliser (positive control)
DECOY = "caffeine"   # inert decoy

# Design: block on starting conformer (the dominant initial-condition variance).
# Within a block, apo / test / decoy share matched velocity seeds; the contrast
# is computed WITHIN block (cancels structural variance), then averaged across
# blocks. Trajectory params fixed to the pilot's so blocks are comparable.
REPLICAS_PER_ARM_PER_BLOCK = 5
PROD_NS = 2.0
TRAJ_INTERVAL_PS = 20.0

# Sequential decision over accumulated blocks, bootstrapping the per-block
# contrast  C_block = mean_seed[ jaccard_decay(test) - jaccard_decay(decoy) ]:
#   stop-H1   : P(mean C > 0) >= POSTERIOR_THRESHOLD            (channel validated)
#   stop-null : >= MIN_BLOCKS collected AND 95% CI within +/- NULL_HALFWIDTH of 0
#   continue  : otherwise, up to MAX_BLOCKS
POSTERIOR_THRESHOLD = 0.97
NULL_HALFWIDTH = 0.005   # 1/ns; a CI this tight around 0 is a practical null
MIN_BLOCKS = 6
MAX_BLOCKS = 30
N_BOOT = 20000
