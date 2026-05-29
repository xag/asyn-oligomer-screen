# Oligomer generation (Track A)

The original three-step in-silico discovery plan is:

1. **Generate candidate α-syn toxic oligomer conformations** — the
   ones that are not in the PDB.
2. **Score which have toxicity features** (exposed hydrophobic
   β-SASA, membrane-disruption geometry, NAC accessibility) using
   Stage 2 features.
3. **Find molecules that destabilise the toxic shapes** — raise
   their free energy so they dissipate or revert to inert fibrils.

Step 1 is the unsolved part. This directory holds the code for it.

<details>
<summary><b>Plain English</b></summary>

This folder builds an atomic-resolution 3D model of the toxic α-syn
oligomer — the thing that no public database has. The construction
respects what experiments *do* tell us (three protein chains stuck
together by a small β-sheet around residues 70–88, with the rest of
each chain disordered) but invents the parts the experiments don't
pin down, then lets physics-based simulation relax the model into a
plausible shape. We make 11 such models with different choices (sheet
direction, chain count, sheet boundaries) and all of them score
"more toxic" than every deposited mature clump — evidence the
ranking is not an artifact of one lucky build.
</details>

The Stage 2 framework (`features.py`, `classifier.py`) is otherwise
correct, but until step 1 produces structures it has been calibrated
against **deposited fibrils only** — that recovers ordinal mutant
severity but does not perform discovery: docking known α-syn
modulators against fibrils reproduces published science with more
geometric detail. Track A is the pivot to scoring **generated**
oligomers, which the PDB does not provide.

## Approach: topology-prior coarse model + MD refinement

Atomic structures of toxic α-syn oligomers do not exist publicly.
Fusco et al. (Science 2017) characterised the Type B* species by
ssNMR — ~35 structured residues per monomer in a β-rich core around
residues 70-88, the rest disordered — but did not deposit
coordinates. Restrained simulated annealing from the published DARR
correlations is the gold-standard reconstruction, but the constraint
set is sparse and the toolchain (XPLOR-NIH / Rosetta-with-NMR) is
Linux-native and crusty.

We use the Fusco paper as a **topology prior**, not as literal
distance restraints:

- 3 α-syn monomers (full length, residues 1-140)
- residues 70-88 in extended β-strand conformation, one strand per
  monomer (default; the range is a CLI knob)
- strands assembled into a small β-sheet (parallel or antiparallel)
- N-terminal 1-69 and C-terminal 89-140 as extended PPII initially,
  let MD collapse them into disordered conformations
- positional restraints on the β-core Cα during MD so the prescribed
  topology survives equilibration; tails are unrestrained

By sweeping the CLI knobs (n_mers, β-core range, parallel vs
antiparallel, strand spacing), we generate an *ensemble* of plausible
toxic-oligomer starting structures, not a single guess.

## Pipeline

```
build_fusco_trimer.py
  → results/oligomers/<tag>.pdb         (extended starting structure)

md_relax.py (in the conda md env)
  --apo-pdb <tag>.pdb
  --collapse-ps 1000        (1 ns OBC2 implicit MD: tails fold up
                             so the explicit-solvent box stays small)
  --equil-ps 100 --prod-ps 100
  --restrain-residues 70-88
  --restrain-chains A,B,C
  --restrain-k 1000
  → results/oligomers/<tag>_relaxed.pdb (relaxed all-atom structure,
                                         heavy atoms only)

score_oligomer.py (TODO)
  → Stage 2 features on the relaxed trimer, compared against the
    deposited-fibril anchor table — does the trimer rank above the
    fibrils, as the active-side anchor it was built to be?
```

## Files

- `build_fusco_trimer.py` — sequence + φ/ψ-based construction of one
  oligomer hypothesis. Pure-Python, depends on PeptideBuilder. Runs
  in the pipeline pip venv.
- `score_oligomer.py` — *(planned)* runs `features.py` on a relaxed
  oligomer PDB and reports the Stage 2 activity score against the
  deposited-fibril anchor calibration.

## Open questions

- **Is the framework right?** If the relaxed trimer scores high on
  hydrophobic-β-SASA and NAC-active-score, well above the fibril
  anchors, the existing Stage 2 features are right and we have the
  missing anchor. If it doesn't, the features themselves need a
  rethink — possibly the `ordered_core_full_ids` mask (calibrated
  for fibril-grade packing density) excludes too much of the
  oligomer.
- **One topology or an ensemble?** A single trimer hypothesis is a
  guess. The intent is to generate ~10 topologies (parallel vs
  antiparallel × dimer/trimer/tetramer × shifted β-core ranges) and
  see which ones cluster in the high-activity region of Stage 2.
- **Sampling.** 100 ps explicit production after 1 ns implicit
  collapse is enough to relax side chains and tails but not enough
  to find a thermodynamic minimum. The restraint on the β-core
  Cα assumes we *know* the right β-core; longer/restraint-free
  sampling might reveal that we don't.
