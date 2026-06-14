# EXPERIMENTS — lab notebook

Append-only record of experiments performed: **hypothesis**, **method** (with a
pointer to the script that ran it), where the **data** lives, the **result**, and
the **conclusion**. Data files are the source of truth; entries link rather than
restate. Pre-registered designs are locked in [`prereg.py`](prereg.py). Science
only — no hardware, host, account, or path details.

---

## E1 — Dwell-time pilot: does the channel recover a known destabiliser? (2026-06-12)

**Hypothesis.** A short-MD "toxic-basin dwell" readout separates a known α-syn
destabiliser (silibinin) from an inert decoy (caffeine): silibinin spends less
time in the toxic basin than apo, and clearly less than caffeine.

**Method.** NAC-core chunk (~55k-atom) of `fusco_parallel_3mer_core70-88_relaxed`;
arms apo + silibinin, dhea, caffeine (trehalose dropped — an RDKit ring-perception
failure in `md_relax` ligand prep). 10 velocity-seeded replicas/arm × 2 ns, CUDA
OpenMM. Dwell = fraction of frames with β-core RMSD < 3.0 Å **and** contact-Jaccard
> 0.5 vs the relaxed reference; shift = mean(complex) − mean(apo), bootstrap 95% CI.
`ops/run_pilot.ps1` + `ops/asyn_launch.py`. Data
`results/dwell/fusco_parallel_3mer_core70-88_relaxed/`.

**Result.** apo mean dwell 0.20 (per-replica 0.01–0.56). Gate all *inconclusive*:
silibinin +0.012 [−0.110,+0.120]; dhea +0.047 [−0.100,+0.196]; caffeine +0.098
[−0.071,+0.290]. The known destabiliser did not separate; the decoy shifted most.

**Conclusion.** The channel is not validated at this design/power.

---

## E2 — Diagnostic + ligand-blind rescore: why null? (2026-06-12)

**Method.** CPU re-score of the E1 trajectories (no new MD): per-replica occupancy;
apo per-frame RMSD/Jaccard distribution vs the thresholds; apo-vs-control frame
time-courses; and a tuning-free rate observable (β-core RMSD-rise and
contact-Jaccard-decay slopes), contrast bootstrapped. `ops/analyze_dwell.ps1` +
`ops/asyn_diag.py`.

**Result.** Occupancy ~0.92–1.00 → ligands stayed bound. apo β-core RMSD
min/median/max = 1.38 / 3.73 / 5.76 Å against a 3.0 Å basin cutoff → the model
relaxes above the basin and keeps drifting (~0.8 Å/ns, no plateau in 2 ns) while
inter-chain contacts stay intact (Jaccard ~0.77). Rate rescore: silibinin
inconclusive on both axes (faint stabilising lean); caffeine the lone crossing
(faster RMSD-rise) — one false positive across 6 tests at n=10.

**Conclusion.** The null is unresolved — underpowered, with a mis-set basin
threshold and a reference that is not metastable on this timescale — not a proven
negative.

---

## E3 — Blocked validation, pre-registered (2026-06-12)

**Hypothesis (one-sided).** silibinin's contact-Jaccard decays faster than
caffeine's.

**Design.** Locked in [`ops/prereg.py`](ops/prereg.py) (blocked on conformer;
matched seeds; within-block contrast; sequential stop). Primary statistic
`jaccard_decay_per_ns` = −(slope of inter-chain contact Jaccard vs time).
Runner `ops/asyn_block.py`; pooled decision `ops/asyn_pool.py`.

**Method.** 9 blocks = 9 distinct Fusco-topology conformers (parallel/antiparallel,
2/3/4-mer, β-core 65–83/70–88/73–91, plus seed variants). Each block: apo + silibinin
+ caffeine, 5 matched velocity seeds/arm × 2 ns; within-block contrast
`mean_seed[decay(silibinin) − decay(caffeine)]`. Data
`results/blocks/<conformer>/`.

**Result.** Per-block matched contrast ranged −0.0256 to +0.0027 /ns; negative in
7/9 blocks. Pooled mean **−0.0072 /ns, 95% CI [−0.0131, −0.0019], P(H1) = 0.002** —
the destabiliser's contacts decay *no faster* than the decoy's (slightly slower).
Two design notes surfaced: (i) the sequential H1 branch trips trivially at one
block, where the single-value bootstrap gives a zero-width CI and P = 1.0, so only
multi-block pooling is interpretable; (ii) the practical-null band (±0.005 /ns) is
tighter than the between-conformer spread, so the null branch cannot fire even as
the mean stabilises — the formal verdict stays *continue* despite a decisive
against-H1 signal.

**Conclusion.** H1 not supported. The shape-stability (dwell/rate) channel does not
separate a known destabiliser from an inert decoy on this model, consistent with E1
and E2 — the channel is not validated and the contributor sweep stays paused. The
docking Δact_gated ranking and the covalent anti-target channel are independent of
this readout and unaffected.

---

## E4 — How much MD clears the noise? The data requirement is set by the space, not the budget (2026-06-13)

**Question.** How much simulation does the dwell channel need to read a ligand
effect beyond chance — given what E1/E2 showed about the conformational space it
samples?

**Method.** Reduce the channel to its sampled coordinate (β-core RMSD; contacts
stay intact at Jaccard ~0.77, so RMSD is the binding constraint) and ask what the
landscape, not the sample size, allows. Closed form on the recorded E1 dwell spread
and E2 drift; no new MD. Calculator [`ops/power.py`](ops/power.py).

**Data.** E1 apo dwell mean 0.20, per-replica 0.01–0.56 (n=10). E2 drift ~0.8 Å/ns,
no plateau in 2 ns; β-core RMSD start 1.38 Å vs 3.0 Å basin cutoff.

**Result.** Three properties of the space decide this before any budget does:
(i) **noise shape** — the reference is not metastable, so a replica's dwell is a
first-passage measurement (drift-diffusion crossing of the cutoff), not an
equilibrium occupancy. Its coefficient of variation is order 1 (E1 gives CV ≈ 0.9,
σ ≈ 0.18; mean first-passage L/v ≈ 2 ns, implied RMSD-coordinate D ≈ 0.5 Å²/ns, no
barrier). The per-replica noise is ~the mean itself and only averages down with
replicas. (ii) **floored range** — apo at 0.20 means the system is already ~80% out
of the basin in 2 ns, so the destabiliser direction has a ceiling of 0.20 against
σ ≈ 0.18: best-case SNR ≈ 1 per replica. (iii) **sign confound** — a bound ligand
mechanically restrains the core (occupancy raises dwell), pushing true destabilisers
toward looking like stabilisers (E1: every arm shifted positive). On this observable
no finite replica count clears chance. *If* the observable is fixed, the cost is the
ordinary two-arm test n/arm = (z_α+z_power)²·2σ²/δ²: ~50 replicas/arm/shape for a
δ=0.10 shift (~1.9 µs across the 9-shape ensemble), ~200 for δ=0.05 (~7.6 µs) —
versus the ~84 ns E1 actually spent.

**Conclusion.** "How much data clears the noise" is the wrong axis: data is not the
binding constraint. The dwell observable is floored, non-metastable, and
occupancy-confounded, so more replicas sharpen a number that has no clean signal in
the direction of interest. The fix is the observable, not the budget — establish a
metastable basin (longer equilibration / restrained or better reference, so apo
dwell → ~1 and the full 0–1 range opens at low CV) or score the first-passage rate
directly. Only then does the ~50–200 replica/arm/shape estimate apply.

---

## E5 — The hand-built shape is not metastable, but toxic *character* survives relaxation (2026-06-13)

**Question.** #56: every channel scores against one hand-built oligomer shape that
E2 showed is not metastable. Released under unbiased dynamics, which conformations
does the oligomer actually occupy, and does any *stable* basin remain toxic-looking?

**Method.** Unbiased MD — β-core restraints **off** — from the relaxed Fusco trimer
cores, parallel `fusco_parallel_3mer_core70-88` and antiparallel
`fusco_antiparallel_3mer_core70-88`, on the NAC-core chunk (residues 58–102,
~55k-atom rect box). Reuses `screen/md_relax.py`'s resumable segment path (which
applies no position restraints) for the dynamics and `screen/shape_metrics.py` for
per-frame β-core Cα RMSD + inter-chain contact Jaccard vs each shape's own relaxed
reference. 3 velocity-seeded replicas/shape × (200 ps equil + 50 ns production),
frames every 50 ps. Post-equilibration frames (second half/replica) pooled and
clustered in standardised (RMSD, Jaccard) space (Ward; k by silhouette); each basin
medoid scored with `oligomers/score_oligomer.py`. Baselines: each relaxed-reference
chunk and the coil-trimer chunk (non-toxic floor). Driver + analysis
[`ops/stable_states.py`](ops/stable_states.py). Data `results/stable_states/`.

**Result.** β-core RMSD drift collapses from E2's ~0.8 Å/ns to ~0.01–0.06 Å/ns over
the last 20 ns — the cores leave the hand-built shape and settle 6–9 Å away.
Parallel → **1 basin** (mean RMSD 6.2 Å, Jaccard 0.61); antiparallel → **2 basins**
(68% at 6.6 Å, 32% at 8.8 Å; Jaccard 0.51–0.56). **No** frame meets the strict
toxic-basin gate (RMSD ≤ 3 Å *and* Jaccard ≥ 0.5): occupancy of the hand-built
geometry is 0. Yet every stable basin's Stage-2 activity stays high — parallel
**11.4**, antiparallel **10.4 / 11.3** — against relaxed-reference baselines (parallel
16.3, antiparallel 13.7) and the coil floor 3.7. Both registers converge to the same
~10–11 band regardless of starting arrangement. One antiparallel replica still drifts
−0.06 Å/ns, so the 8.8 Å basin is the least settled.

**Conclusion.** Two-part answer to #56. (i) The *specific* conformation the channels
dock against is **not** metastable — vacated within tens of ns to a basin 6–9 Å away,
so any ranking that turns on the precise pose/geometry is scored against a structure
the model does not hold (confirms and extends E2). (ii) But toxic *character* is
robust: every populated, ≈plateaued basin still scores ~65–70% of its reference
activity, far above the coil floor, and parallel and antiparallel relax into the same
activity band. So a stable, still-toxic-looking ensemble does exist — it is just not
the hand-built pose. A defensible target is the relaxed ensemble (or an experimental
structure), not the hand-built conformation; pose-dependent docking should be
re-grounded there. Longer sampling would tighten the least-settled antiparallel basin
but does not change the activity verdict.

---

## E6 — Diverse starts do not funnel to one shape: a register-split, history-dependent landscape (2026-06-14)

**Question.** E5 used a single starting region. Do *independent* starting
conformations converge to a common ensemble (a real attractor), or to distinct
basins (rugged, history-dependent)? Directly tests the "would another shape
converge to it?" gap left open by E5. (Follows #56.)

**Method.** 8 independent starts, each relaxed 50 ns unbiased (no restraints,
`screen/md_relax.py` segment path) on the NAC-core chunk, one replica each:
parallel core70-88 (seeds 42/123/777), antiparallel core70-88 (seeds 42/123),
parallel core65-83, parallel core73-91, and the all-coil trimer (a disordered
start). Per start: β-core RMSD/Jaccard vs its own start (plateau slope, drift)
and end-basin activity (`oligomers/score_oligomer.py`). Cross-start: pairwise
β-core RMSD among the relaxed medoids (+ the hand-built target), hierarchical
cluster counts, and a classical-MDS map. Driver + analysis
[`ops/funnel.py`](ops/funnel.py). Data `results/funnel_states/`.

**Result.** The end-basins do not coincide: mean pairwise β-core RMSD **14.1 Å**
(min 6.0, max 21.5); cluster counts 8 / 4 / **3** at 4 / 8 / 12 Å cuts. **Register
is the dominant split.** All 5 parallel starts — including the shifted cores
65-83/73-91 and different coil draws — relax into one loose ~7–9 Å neighbourhood
that *contains the hand-built target*; the 2 antiparallel starts cluster with each
other but **15–18 Å** from every parallel one; the coil start never organises
(still drifting +0.23 Å/ns, isolated). End-basin activity is itself
start-dependent: parallel 7.5–16.0, antiparallel 8.4 and **4.6**, coil 4.4 (floor
≈ 3.7).

**Conclusion.** No single global attractor. *Within* a register, independent starts
converge — the parallel family reproducibly relaxes to ~the hand-built region, so
the parallel target is not arbitrary *within its register*; *across* registers they
settle into distinct, non-interconverting basins on this timescale, and disorder
does not fold in 50 ns. The landscape is rugged and history-dependent: which basin —
and how toxic — depends on the starting register. This **qualifies E5's "toxic
character is robust" reading**: with broader sampling toxicity is not uniform (the
antiparallel s123 start fell to the coil floor). For #56 the open question is no
longer whether a toxic-scoring basin exists (several do) but *which register/basin is
biologically real* — which these single-replica, single-quench, trimer-scale,
membrane-free simulations cannot decide; only an experimental oligomer structure can.
