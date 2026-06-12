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

**Design.** Locked in [`prereg.py`](prereg.py) (blocked on conformer; matched
seeds; within-block contrast; sequential stop). `ops/block.ps1`.

**Data.** None collected yet.
