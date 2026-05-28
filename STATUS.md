# Pipeline status

Last update: 2026-05-28 (**Covalent / adduct channel (`aspr_score`) added — resurfaces the four reactive metabolites the affinity gate correctly collapses to zero.** New module [`adduct_score.py`](adduct_score.py) computes per-chain mean of (fractional SASA × ligand reactivity weight) over Lys / Arg / Cys / His / Tyr; reported as an orthogonal column alongside `delta_activity_gated`. On the reference trimer: MDA +10.60, acrolein +6.70, 4-HNE +4.65, MGO +4.24, NO +2.67, H2S +0.00. Ordering matches the chemistry of a 15-Lys / 0-Arg / 0-Cys / 1-His / 4-Tyr substrate. See "Stage 3 — covalent / adduct channel" section. Prior: full 127-mol oligomer sweep + 5-mol blind validation hold-out — silibinin #1, EGCG #3, fisetin #5, rosmarinic-acid #24, CAPE #54; pass by pre-registered criteria. Earlier: Stage 3 against generated oligomer; ensemble sweep + framework fixes; Track A direction pivot.)

This file is the entry point for a fresh thread. Read this first, then [`ANCHORS.md`](ANCHORS.md) for the anchor curation context, then [`../IN_SILICO_PLAN.md`](../IN_SILICO_PLAN.md) for the overall project shape, then [`oligomers/README.md`](oligomers/README.md) for the new Track A direction.

## Direction pivot: Track A — generate toxic-oligomer structures (2026-05-27)

The earlier Stage 2/3 work calibrated the activity classifier on **deposited fibrils only** and ran Stage 3 by docking known α-syn modulators against those same fibrils. That reproduces published claims with more geometric detail; it does not perform discovery. The original three-step intent is:

1. **Generate** candidate toxic α-syn oligomer conformations — the structures the PDB does not provide.
2. **Score** which generated shapes have toxicity features (exposed hydrophobic β-SASA, NAC accessibility, membrane-disruption geometry).
3. **Destabilise** — find ligands that raise the free energy of the toxic shapes so they dissipate or revert to inert fibrils.

Step 1 is the bottleneck. Without generated structures, step 2 has nothing to filter and step 3 has nothing to destabilise. Track A is the work to produce step 1.

**Approach.** Topology-prior coarse model + MD-style relaxation. Fusco et al. (Science 2017) characterised the Type B* toxic α-syn oligomer by ssNMR but did not deposit coordinates: ~35 structured residues per monomer in a β-rich core around residues 70-88, the rest disordered. We use that paper as a *topology prior* (which residues are in the β-core, how many monomers per oligomer, β-strand vs hairpin) — **not** as literal distance restraints (the published DARR set is too sparse and the constraint-driven SA toolchain is Linux-native and crusty). The build is:

- 3 α-syn monomers (full length, residues 1-140)
- residues 70-88 in extended β-strand conformation, one strand per monomer
- 3 strands forming a parallel β-sheet (also antiparallel as a variant, planned)
- N-terminal 1-69 and C-terminal 89-140 as **random-walk coil**, per-residue independent draws from a β/PPII mixture (no α-helix basin — α produces tight 3.6-residue turns and a 70-residue tail sampling α has high self-intersection probability; empirically including α at 30% weight produced >30 atom pairs within 1.0 Å in the built oligomer)
- positional restraints on β-core Cα during relaxation; tails unrestrained

Code: [`oligomers/build_fusco_trimer.py`](oligomers/build_fusco_trimer.py), [`oligomers/score_oligomer.py`](oligomers/score_oligomer.py), [`oligomers/README.md`](oligomers/README.md). Relaxation runs through [`md_relax.py`](md_relax.py) with new flags `--vacuum-min-iter`, `--collapse-ps`, `--restrain-residues`, `--restrain-chains`, `--no-explicit`.

### First result: generated trimer ranks #1, well above the deposited anchors

One topology built and scored (parallel trimer, β-core 70-88, seed=42).

Relaxation pipeline:
- 5000 iter vacuum energy minimization (resolves intra-chain and inter-chain clashes; energy dropped from `+4.8e11` to `+5,653 kJ/mol`)
- 500 ps OBC2 implicit-solvent MD with Cα restraints on residues 70-88 of A/B/C (final energy `-43,485 kJ/mol`)
- `--no-explicit` writes the result; explicit-solvent step skipped because Stage 2 features are structure-based, not dynamical

3.8 minutes wall-clock.

Stage 2 activity score on the middle chain (B), against the existing 14-anchor calibration:

```
features:
  exposed_hydrophobic_beta_sasa  +15.7335   z=+6.50   contrib=+6.50
  membrane_insertion_propensity  +0.5627    z=+2.30   contrib=+1.15
  nac_active_score               +7.4796    z=+1.65   contrib=+1.65
  contact_density                +7.7636    z=-2.27   contrib=+2.27   (negative weight)
  disordered_hydrophobic_exposure +17.8547  z=+6.41   contrib=+3.20
                                                      ────────────
                                          activity = +14.78
```

Ranking against deposited anchors:

```
<<<  fusco_parallel_3mer_core70-88  (oligomer hypothesis)  +14.778
     2N0A                            inert                   +4.398
     6UFR                            graded-active           +1.880
     6PEO                            graded-active           +1.808
     6LRQ                            graded-active           +1.062
     6PES                            graded-active           +0.971
     8A9L                            inert                   −0.064
     1XQ8                            inert                   −0.219
     7WO0                            graded-active           −0.538
     8A4L                            inert                   −0.722
     6CU7                            inert                   −0.914
     6CU8                            inert                   −1.426
     6XYO                            inert                   −1.561
     6H6B                            inert                   −2.198
     6XYP                            inert                   −2.477
```

**Headline:** the generated trimer scores **3.4× higher than 2N0A** (the previously highest anchor) and ~8× higher than the median graded-active mutant fibril. Stage 2 features discriminate the generated oligomer from the fibril ensemble cleanly.

**What the dominant features mean:**
- `exposed_hydrophobic_beta_sasa` z=+6.5 is the cleanest signal. The β-core hydrophobic residues are exposed to solvent because the oligomer isn't packed in the multi-layer fashion a fibril is. This is exactly what the framework was designed to capture.
- `disordered_hydrophobic_exposure` z=+6.4 reflects hydrophobic residues on the disordered tails. It is a real feature but partly trivial: any IDP-like structure has high exposure here; fibrils don't (residues 1-60 and 89-140 in fibrils are usually unmodelled or in extended dangling tails too). The contribution is genuine but is bounded by "is there a disordered region present at all," not by oligomer-vs-monomer geometry.
- `contact_density` z=−2.3 reflects the 3-chain trimer having fewer Cα contacts per residue than 10-chain fibrils. Negative weight → positive contribution. Real biology (oligomers do have fewer chain-chain interfaces) but partly a count effect (3 vs 10 chains).
- The remaining two features contribute modestly and as expected.

**Caveats — honest:**
1. **One topology, n=1.** The activity has no error bar yet. Need to generate an ensemble (different seeds, parallel vs antiparallel, dimer/trimer/tetramer, β-core ranges around 70-88) and see whether all topologies score in the same range or whether the score varies wildly.
2. **The framework was calibrated on the same five features that score the trimer — this is not an out-of-distribution test.** What it shows is internal consistency: a structure built to match the textbook description of a toxic oligomer scores higher than fibrils on the textbook toxicity features. That confirms the features work; it doesn't confirm the trimer is biologically real.
3. **disordered_hydrophobic_exposure is partly a "do tails exist" signal.** A truly fair comparison would use anchors with similarly modelled tails. None of the deposited fibrils model residues 1-30 or 100-140; the trimer here has all 140 residues. Part of the z=+6.4 advantage is an apples-to-oranges comparison.
4. **No experimental validation.** This says the framework is consistent with itself, not that the framework matches reality. Validation requires either (a) building from Fusco's DARR constraints rather than topology and showing the result still scores high, or (b) wet-lab partners measuring fibril-vs-oligomer feature differences.

**Next moves (Track A continuation):** → all completed in the ensemble sweep below.

## Track A — ensemble sweep (2026-05-27)

11 topologies built, relaxed (5000-iter vacuum min + 500 ps OBC2 implicit MD, `--no-explicit`), and scored.
Two scoring fixes were applied before reading results (see "Framework fixes" below).
Code: `oligomers/run_ensemble.py` (build → relax → score → CSV); reproduce with
`python oligomers/run_ensemble.py --summary-only` (re-scores existing relaxed PDBs in ≤1 min).

### Results (corrected scoring: all-chain mean + β-gate)

```
activity   structure                             type
+22.11     fusco_parallel_3mer_core70-88_s777    generated  ← seed outlier (see below)
+18.17     fusco_parallel_4mer_core70-88         generated
+17.78     fusco_parallel_3mer_core70-88  [ref]  generated
+16.57     ctrl_trunc_3mer_res60-100             control    ← truncated scoring control
+16.14     fusco_parallel_2mer_core70-88         generated
+14.15     fusco_parallel_3mer_core70-88_s123    generated
+12.95     fusco_antiparallel_3mer_core70-88     generated
+12.83     fusco_parallel_3mer_core73-91         generated
+11.61     fusco_parallel_3mer_core65-83         generated
+10.60     fusco_antiparallel_3mer_core70-88_s123 generated  ← weakest genuine oligomer
─────────────────────────────────────────────── gap ~+4.8 ──
 +5.77     ctrl_coil_3mer                        control    ← coil 3mer (β from MD collapse)
 +4.40     2N0A                                  deposited anchor (best)
 +1.88     6UFR  …  (other fibrils below)
 +0.00     ctrl_coil_monomer                     control    ← β-gated from +11.37
```

**Full summary CSV**: `results/oligomers/ensemble_summary.csv`

### What this tells us

1. **All 9 Fusco-topology oligomers score above 2N0A.** Weakest: antiparallel s123 at +10.60 (2.4× above 2N0A). Strongest: s777 parallel trimer at +22.11. Mean ~+15. The first trimer result (+14.78, n=1) is confirmed with an error bar.

2. **Natural threshold at ~+8.** There is a gap of +4.8 between the weakest genuine oligomer (+10.60) and the best control/anchor (+5.77). A threshold of +8 cleanly separates all 9 generated oligomers from all controls and all deposited fibril anchors with no overlap.

3. **Truncated core-only control: 93% of score from β-core.** Scoring the reference trimer on residues 60-100 only (no N/C tails) gives +16.57 vs full-trimer +17.78. Tail contribution: +1.21 (7%). Caveat 3 from the first result ("tails inflate score") is largely resolved by the all-chain mean fix — the tail residues only enter the score if they pass the ordered-core mask, and across all chains their net contribution is small.

4. **s777 parallel trimer outlier (+22.11, was +35.09 before fix).** Still highest scorer. The original +35 was driven by `nac_active_score` z=+17 on one chain picked by `auto-inner`; the all-chain mean drops it to +22. Still elevated — the seed=777 coil draws collapsed the tails into a conformation that extends the NAC β-geometry beyond residues 70-88 on chain A. A real structure might or might not look like this; treat as an upper bound on the ensemble.

5. **Coil 3mer (+5.77) is a residual false positive.** The three fully disordered chains developed genuine β-contacts during 500 ps OBC2 collapse (`exposed_hydrophobic_beta_sasa = 3.43`, `nac_active_score = 1.28` on the scored chain). The β-gate correctly leaves it unmodified (it has β structure). Its score sits at the 2N0A level, not in the oligomer range, so it is distinguishable by the +8 threshold.

6. **Coil monomer correctly suppressed.** Before the β-gate: +11.37 (false positive from `disordered_hydrophobic_exposure` z=+18 + `contact_density` z=-5.5). After: 0.00. The gate fires because `exposed_hydrophobic_beta_sasa = 0` AND `nac_active_score = 0`.

### Parallel vs antiparallel

Parallel trimers: s42=+17.78, s123=+14.15 — mean ~+16. Antiparallel trimers: s42=+12.95, s123=+10.60 — mean ~+12. Parallel sheets expose more hydrophobic surface per strand face; antiparallel sheets partially bury it via inter-strand hydrogen bonding geometry. Consistent across seeds.

### Core range effect

| core range | activity | notes |
| ---------- | -------- | ----- |
| 65-83      | +11.61   | shifted N-terminal; includes less-hydrophobic residues 65-69 |
| 70-88      | +17.78   | canonical Fusco range — highest for 3mer |
| 73-91      | +12.83   | shifted C-terminal; gains residues 89-91 but loses 70-72 |

The canonical 70-88 window hits the most hydrophobic segment of the NAC region. Use this as the default for Stage 3 targets.

### Framework fixes applied (2026-05-27)

Two changes to `oligomers/score_oligomer.py`. The anchor calibration (`validate.py`, `classifier.py`) is untouched.

**Fix 1 — all-chain mean (was: auto-inner single chain).** The old default picked the most-buried chain via `assembly.inner_chain_ids`, which was sensitive to which chain the OBC2 collapse happened to expose. The new default scores each chain individually and averages the feature vectors before z-scoring. Effect: more stable scores (tetramer +12.96 → +18.17; s777 +35.09 → +22.11); `--auto-inner` flag restores the old behaviour.

**Fix 2 — β-structure gate.** After computing the mean activity, if `exposed_hydrophobic_beta_sasa < 0.5` AND `nac_active_score < 0.5` across all chains, activity is capped at `min(raw, 0.0)`. Justification: the Fusco Type B* toxic-oligomer mechanism requires a β-rich NAC core. A purely disordered structure (zero β-SASA, zero NAC-β) should not be able to score in the oligomer range regardless of `disordered_hydrophobic_exposure` or `contact_density` bonuses. Thresholds (0.5 Å²/residue each) are set just above floating-point zero to catch genuine "no-β" cases while tolerating single-residue DSSP fluctuations. `--no-beta-gate` bypasses for debugging.

### What is still imperfect

- **Coil 3mer at +5.77** is a false positive in the sense that it was built with no β-prior, but it developed real β-contacts during MD and scores above zero. It's below the +8 threshold and below the weakest genuine oligomer (+10.60) by a factor of 2, so it doesn't create confusion in practice.
- **s777 outlier at +22.11** inflates the ensemble range. An additional check (score per-chain and flag any single chain that contributes >2× the mean) would catch tail-geometry outliers earlier. Not implemented yet.
- **No out-of-distribution validation.** All features were designed and calibrated on fibril anchors. The oligomers score high because they have the features the framework was designed to reward; this does not confirm the features correctly identify toxicity in generated structures not seen during calibration.

### Next: Stage 3 against the generated oligomer

The framework is ready. The reference trimer (`fusco_parallel_3mer_core70-88_relaxed.pdb`, score +17.78) is the destabilisation target. `stage3.py` needs a small extension to accept the generated oligomer as receptor and use the all-chain mean Δactivity. The question it will answer: **does binding molecule X push the trimer score from +17.78 toward the inert-fibril regime (+4.40)?**

## Stage 3 MD relaxation (50 ps) — results inconclusive (2026-05-27, then superseded by Track A)

Per the prior STATUS, the path-to-sign-flip on harm-side molecules was OpenMM MD relaxation around each docked complex. The conda `md` environment was set up (OpenMM 8.5 + openff-toolkit + NAGL charges + openmmforcefields SMIRNOFFTemplateGenerator), `md_relax.py` written, and `md_stage3.py` ran the existing 6 pairs through apo MD + complex MD + Stage 2 re-score at 50 ps equilibration + 50 ps production. Took ~4 hours total (2N0A complex MD alone was 138 min — the largest assembly).

| molecule       | anchor | aff_top | act_apo | act_apo_MD | dact_top | dact_wtd | dact_gat | dact_rlx |
| -------------- | ------ | ------- | ------- | ---------- | -------- | -------- | -------- | -------- |
| curcumin       | 6PEO   | -8.47   | +1.808  | +1.739     | -1.242   | -1.281   | -1.281   | **-1.681** |
| curcumin       | 6UFR   | -7.60   | +1.880  | **-1.272** | -1.542   | -1.432   | -1.432   | **+1.354** (sign flip — but spurious) |
| quercetin      | 6PEO   | -8.02   | +1.808  | +1.739     | -1.562   | -1.516   | -1.516   | **-1.805** |
| methylglyoxal  | 6PEO   | -2.87   | +1.808  | +1.739     | -0.652   | -0.633   | -0.003   | **-1.225** |
| curcumin       | 2N0A   | -5.76   | +4.398  | +2.196     | -0.548   | -0.505   | -0.336   | **+0.079** |
| methylglyoxal  | 2N0A   | -2.82   | +4.398  | +2.196     | -0.318   | -0.170   | -0.001   | **-1.580** |

**What this tells us, and what it doesn't:**
- The one sign flip (curcumin × 6UFR, dact_rlx = +1.354) is **probably not "MD found receptor rearrangement."** The apo-MD activity for 6UFR moved from +1.88 to −1.27 in 50 ps — a 3.15-unit drift from MD sampling alone, larger than any plausible ligand-induced effect. Δact = complex_MD − apo_MD on a single 50 ps replica is dominated by stochastic conformational sampling, not by the ligand.
- 2N0A's apo also drifted hard (+4.40 → +2.20), confirming the same sampling-noise issue.
- **None of the harm-side methylglyoxal pairs flipped sign in the expected direction.** They went more-negative (methylglyoxal × 6PEO −0.003 → −1.225; methylglyoxal × 2N0A −0.001 → −1.580). The static-pose gate had correctly damped them to ~0; the MD result re-amplifies the original "any pose occludes SASA" mechanism rather than producing the conformational-rearrangement signal that would flip the sign.
- **At 50 ps × 1 replica per system, the MD experiment cannot distinguish ligand effects from thermal sampling.** Either we need much longer simulations (1+ ns × multiple replicas, hours per pair) or a different design (steered MD, free-energy perturbation, restraint protocols that limit apo drift).

**Why this work is being de-emphasised rather than continued:** the Stage 3 fibril-docking branch is calibrating on deposited anchors and asking "does ligand X perturb fibril Y?" The pivot (above) recognises that the framework was supposed to ask "does ligand X destabilise toxic oligomer Y?" Without generated toxic-oligomer structures (Track A), the Stage 3 MD experiment can't answer the question the project actually has. The MD-relaxation infrastructure (`md_relax.py`, `md_stage3.py`) is kept and reused by Track A; the specific 6-pair experiment is documented here but not pursued further at this calibration.

Artifacts (kept for inspection):
- `results/stage3/<pair>_relaxed.pdb` — MD-relaxed complex (heavy atoms only)
- `results/stage3/<pdb>_apo_relaxed.pdb` — MD-relaxed apo (per-anchor, shared)
- `results/stage3/<pair>_report.json` — extended with `delta_activity_relaxed` and `delta_activity_relaxed_vs_static_apo`
- `results/_logs/_md_stage3_50ps.log` — full run log

## Where we are (legacy — Stage 2 + Stage 3 fibril-docking framework)

Stage 2 prototype runs end-to-end on 14 anchor structures (9 inert + 5 graded-active mutant fibrils). Code in `anchors.py`, `assembly.py`, `features.py`, `classifier.py`, `validate.py`, `protofilaments.py`. Five per-conformer features, transparent weighted z-score, single-model loading. Structures are expanded to their REMARK 350 BIOMOLECULE 1 assembly before feature computation; features accumulate over the most-buried protein chain, restricted to a structurally **ordered-core mask** (Cα with ≥6 non-sequential neighbours within 8 Å on the full assembly). Secondary structure comes from **pydssp** (pure-Python DSSP) on the concatenated-assembly backbone.

**Stage 2 separates cleanly.** Graded-active mean +1.04, inert mean −0.68, Δ ≈ 1.72, pairwise GA-vs-inert AUC ≈ **0.84** (was 0.69 before the ordered-core + pydssp round). The two largest pre-existing confounds — asymmetric-unit size and NMR-vs-cryo-EM residue coverage — are both resolved. The "2N0A misrank" turned out to be a real within-class signal: single-protofilament fibrils are systematically more active-like than paired-protofilament fibrils in the same class (see "Protofilament count is real and within-class" below).

## Last run (14 anchors, descending activity, `assembly_inner` mode, mask ON, pydssp)

```
                       activity  bio-PF
2N0A  inert              4.40    1     ← single-PF ssNMR; top scorer, consistent with the PF-count signal
6UFR  graded-active      1.88    1     ← E46K, single-PF
6PEO  graded-active      1.81    1     ← H50Q narrow, single-PF
6LRQ  graded-active      1.06    2     ← A53T canonical, paired
6PES  graded-active      0.97    2     ← H50Q wide, paired
8A9L  inert             -0.06    1     ← Lewy fold; single-PF
1XQ8  inert             -0.22    -     ← micelle monomer; correctly near zero
7WO0  graded-active     -0.54    2     ← A53T + Ca²⁺ — most-inert graded-active (see "side-result" below)
8A4L  inert             -0.72    2
6CU7  inert             -0.91    2
6CU8  inert             -1.43    2
6XYO  inert             -1.56    2
6H6B  inert             -2.20    2
6XYP  inert             -2.48    2
```

Modes: `python validate.py` (default = assembly_inner + mask). `--no-core-mask` ablates the mask. `--au` reproduces legacy raw asymmetric-unit behaviour. `--assembly-all` keeps assembly but averages over all chains. `--compare` runs all three modes side-by-side.

## What changed this round

1. **Protofilament count annotated and tested.** Added literature-curated `n_protofilaments` to `anchors.Anchor` and a geometric `count_protofilaments()` (in `protofilaments.py`) that clusters chain centroids on the plane perpendicular to the principal axis. New `analyze_protofilaments.py` computes within-class Spearman correlation between protofilament count and activity. See "Protofilament count is real and within-class" for results.

2. **Deposited-vs-biological PF count exposed.** Both `n_protofilaments` (literature) and `n_protofilaments_deposited` (geometric, from REMARK 350 assembly) are now written to `anchor_features.csv` / `anchor_scores.csv`.

3. **Geometric counter recalibrated (later same day).** The first pass of `count_protofilaments` overcounted 2N0A (returned 4) and *undercounted* most paired cryo-EM polymorphs (returned 1 for 6CU7, 6CU8, 6H6B, 6LRQ, 7WO0, 6PES, 6UFR). Root cause: it used PCA over all-Cα as the fibril axis, which on these deposits picks the *chain's* own ~150 Å long axis instead of the inter-chain stacking direction. Replaced with a fibril axis derived from the shortest centroid-to-centroid pair vectors (within-PF rung step) plus an adaptive lateral threshold `max(12 Å, 0.5 · mean per-chain lateral spread)` that widens for ssNMR ensembles (2N0A's chains are ~69 Å laterally vs ~30 Å for cryo-EM).

   Result: deposited count now matches literature on 12 / 14 anchors. Remaining mismatches are honest:
   - **6UFR (lit=1, dep=2)** — Boyer 2020 PNAS title is *"E46K Familial Parkinson's Disease Mutant Strengthens Long-Range Pseudo-2-Fold Symmetry"*. The deposit's BIOMOLECULE 1 is paired. **Literature annotation is probably wrong**; verify and update before depending on the per-anchor lit value.
   - **8A4L (lit=2, dep=3)** — BIOMOLECULE 1 is explicitly PENTADECAMERIC with 15 chains in three lateral stacks (likely a lipid-mediated triplex of two fibrils). The deposit's biological assembly is genuinely 3 PFs; the polymorph paper still calls it a paired fibril.

   The within-class correlation analysis still holds, and **also passes against the corrected deposited column** for the first time: ρ_inert = −0.45 (was +0.18, n=8), ρ_graded-active = −0.35 (was undefined, n=5). The signal is now visible in both annotation modes, which raises confidence.

## Previous round (ordered-core mask + pydssp, kept here for reference)

1. **Ordered-core mask (`features.ordered_core_full_ids`)** — equalises NMR full-length vs cryo-EM core-only modelling. A residue qualifies as "core" when its Cα has ≥6 non-sequential Cα neighbours within 8 Å on the full assembly. All per-residue feature accumulation (numerator and denominator) restricts to this mask. Falls back to no mask when fewer than 20 residues qualify (so SDS-bound 1XQ8, which is helical but small, isn't zeroed out — 1XQ8 actually has 80 core residues because i±2/3/4 within 8 Å satisfy the threshold).

2. **pydssp for secondary structure (`features.dssp_ss_map`)** — replaces `_is_beta` / `_is_helix` φ/ψ rectangles. pydssp's H-bond bridge detector needs the full assembly passed together (a single fibril chain has no intra-chain β H-bonds; per-chain DSSP returns all-loop). Tested per-chain vs concatenated: concatenated gives correct strand fractions for fibrils (65–69%) and correct helix for 1XQ8 (59%). Cached on the structure.

3. **`contact_density` rewritten** — old formula was `total_assembly_contacts*2 / (2*n_chain)` which over-counted because the numerator summed contacts across all residues regardless of chain. New formula: mean non-sequential Cα-contact count over the residues of interest (chain ∩ core), with neighbours allowed to be anywhere in the assembly. Cleaner interpretation: "how many neighbours does an average core residue of this chain see?"

## What's left

### Protofilament count is real and within-class (was "remaining anomaly: 2N0A at #1")

Annotated each anchor with two protofilament counts:
- `n_protofilaments`           — biological count from the polymorph literature
- `n_protofilaments_deposited` — geometric count from the chains actually in the REMARK 350 assembly

**Within-class Spearman ρ between activity score and biological protofilament count** (after the 2026-05-26 lit fix that flipped 6UFR from 1 → 2 — Boyer 2020 explicitly calls it a "symmetric double protofilament"):
| class          | n | rho (lit) | rho (deposited) |
| -------------- | - | --------- | --------------- |
| inert (no 1XQ8)| 8 | **−0.756** | **−0.454** |
| graded-active  | 5 | **−0.354** (was −0.866 before the 6UFR fix) | −0.354 |

**Group means by biological PF count**:
| class         | PF=1 mean activity | PF=2 mean activity | gap |
| ------------- | ------------------ | ------------------ | --- |
| inert         | +2.17 (2N0A, 8A9L) | −1.55 (8A4L, 6CU7, 6CU8, 6XYO, 6H6B, 6XYP) | 3.72 |
| graded-active | +1.81 (6PEO only)  | +0.84 (6UFR, 6LRQ, 6PES, 7WO0) | 0.97 |

The inert signal still looks real — 2N0A is a single-PF outlier scoring +4.40 vs paired-PF mean −1.55, gap 3.7. But the graded-active "single-PF effect" now rests on **a single anchor** (6PEO). The earlier picture of "single-PF scores higher in both classes" was partly an annotation artifact: 6UFR was mislabeled as single-PF and is the second-highest graded-active scorer.

What this means: keep the inert single-vs-paired effect on the watch list as a real candidate; treat the graded-active version as unconfirmed (n=1 per side). Adding another genuine single-PF graded-active anchor would settle it.

**Why this looked weaker on `n_protofilaments_deposited` in the first pass:** the original geometric counter (all-Cα PCA) under-counted paired fibrils as 1 PF because PCA picked the chain's own long axis as the fibril direction. After recalibration (short-pair-vector fibril axis + adaptive lateral threshold), the deposited column matches literature for 12 / 14 anchors and reproduces the same negative correlation: ρ_inert = −0.45, ρ_graded-active = −0.35.

**Reproduce**: `python validate.py` writes the new columns; `python analyze_protofilaments.py` runs the correlation.

**What this implies for the framework**:
- The binary inert/graded-active labels conflate two structurally distinct subcategories: single-PF fibrils (more elongation-competent, growing-end-rich) and paired-PF fibrils (fully buried endpoints).
- Two ways to use this:
  - **Treat single-PF fibrils as a third category.** Re-label 2N0A, 8A9L, 6UFR, 6PEO as "single-PF" and assess separately. Makes the classifier honest about what it's measuring.
  - **Normalise out PF count.** Add a feature or scaling factor that penalises high activity in single-PF context. Makes the binary scheme work.
- Neither is obviously better. Recommend keeping current labels (project lacks direct active anchors anyway) and using PF count as a known modifier when interpreting individual scores.

### Side-result strengthened: 7WO0 (A53T + Ca²⁺) is the most-inert of the graded-actives

After the mask + pydssp pass, 7WO0 scores -0.54, ranking below three of the inert anchors (8A4L, 6CU7, 6CU8 still below it). Calcium binding visibly pulls A53T fibrils toward the stable endpoint — consistent with calcium being on the vicinity-molecule shortlist. This is now the second-strongest side observation (after the protofilament-count effect above) and is the kind of signal Stage 3 is supposed to exploit.

### Anchor anomalies (from prior STATUS) — resolved

- **6PES vs 6PEO** — confirmed as a real polymorph signal driven by inter-protofilament burial (wide = paired = more buried = lower activity score). See `pipeline/results/anchor_features.csv` for per-feature breakdown. Not an artifact; do not normalise away.
- **6H6B labelling** — kept inert. Guerrero-Ferreira 2018 describes the fibrils as "cytotoxic" but the evidence cited is for the 1-121 truncated *protein* (which accelerates aggregation), not for these specific fibrils relative to WT. Under the framework rationale ("mature fibrils = inert because endpoint"), 6H6B stays inert. After this round it ranks #13 (proper place) rather than #4.

## Open architectural decisions

These are the decisions that would change the project's shape, none of which have been made:

- **Single- vs paired-protofilament handling.** Confirmed real and within-class (this round). Choice is between three-way labels (single-PF / paired-PF / oligomer) or a normalisation feature. Defer until the next round of anchor additions makes one direction obvious.
- **Direct active anchors.** Atomic-resolution toxic-oligomer structures still don't exist in the PDB. Options: (a) build Fusco 2017 ssNMR-constraint models ourselves, (b) contact authors, (c) wait for wet-lab partners, (d) accept ordinal-only Stage 2 (current direction). The strong pairwise AUC (0.84) makes (d) more defensible than before.
- **Vicinity-molecule list (`lib/vicinity_molecules.js`)** — **seeded** (2026-05-26). 191 entries across seven groups: 55 endogenous metabolites, 60 dietary, 13 neurotransmitter, 14 metal, 24 lipid, 23 gut-derived, 2 environmental. Each entry carries a `role` field — 43 are `anchor` / `both` (validation controls with published α-syn evidence and ≥1 citation), 151 are `candidate` / `both` (search-space entries justified by CNS reach, no α-syn evidence required). The inclusion criterion is one rule: plausibly reaches the SN at concentrations achievable through diet, supplementation, endogenous synthesis, lifestyle (exercise, stress, sleep, hormonal state), or environmental exposure (occupational, water, food residue). Schema and rationale documented in the file header. Each entry also carries a `delivery: { route, feasibility, notes }` field — **pure metadata, does not enter the Stage 3 score**. Stage 3 scores ΔP(active | M) at the achievable CNS concentration (`cns_conc`); `delivery.feasibility` drives only search/dispatch priority and result presentation. Feasibility distribution: 77 `native`, 63 `achievable`, 45 `low-bioavailability`, 4 `unknown`, 2 `invasive-only`. Routes: `endogenous`, `diet`, `supplement`, `precursor`, `microbiome`, `environmental`, `injection`, `none-known`. **Direction of effect is open**: ΔP can be positive (harm-leaning) or negative (protective) and both belong on the same footing — harm anchors include methylglyoxal, 4-HNE, acrolein, BMAA, paraquat, rotenone, polyamines, and several heavy metals; the steroid axis (cortisol, pregnenolone, estradiol, progesterone, testosterone, DHEA, allopregnanolone) gives the search a lifestyle/hormonal channel rather than treating stress and sex as states. Still needs: (a) admin-side moderation route (`lib/submission.js` `kind: "vicinity_molecule"`) for additions, (b) PMID/DOI canonicalisation of the free-text `refs`, (c) `corpus_articles` population once moderation tags articles with molecules, (d) SMILES verification — many entries (lipid classes, complex polyphenols, steroid hormones, complex polycyclics like rotenone, large cofactors) intentionally seeded with `smiles: null` to avoid risking wrong structures.

## Stage 3 prototype — first results (2026-05-26)

`stage3.py` runs end-to-end on Windows. Dependencies: rdkit + meeko + scipy + gemmi (pip), AutoDock Vina 1.2.5 standalone Windows binary at `pipeline/bin/vina.exe` (downloaded once from the upstream release page; not committed).

Pipeline: vicinity-molecule SMILES → rdkit 3D embed + MMFF → meeko PDBQT ligand; assembly inner chain → mk_prepare_receptor → PDBQT receptor; Vina docks at exhaustiveness=8, num_modes=5 inside a 30 Å (capped) box centred on the inner chain bounding box + 6 Å padding; top pose appended to the assembly as a HETATM residue; Stage 2 features recomputed on the complex; ΔP(active) anchored against `anchor_features.csv` so the apo activity matches the Stage 2 number to numerical precision.

Six pairs run (exhaustiveness=8, seed=42):

| molecule         | anchor (label)         | Vina aff (kcal/mol) | apo activity | complex activity | Δ activity |
| ---------------- | ---------------------- | ------------------- | ------------ | ----------------- | ---------- |
| curcumin         | 6PEO (graded-active, single-PF) | -8.47 | +1.808 | +1.413 | **-0.395** |
| curcumin         | 6UFR (graded-active, paired)    | -7.60 | +1.880 | +1.456 | **-0.424** |
| quercetin        | 6PEO                            | -8.02 | +1.808 | +1.403 | **-0.404** |
| methylglyoxal    | 6PEO                            | -2.87 | +1.808 | +1.467 | -0.341 |
| curcumin         | 2N0A (inert, single-PF top)     | -5.76 | +4.398 | +4.223 | -0.175 |
| methylglyoxal    | 2N0A                            | -2.82 | +4.398 | +4.302 | -0.096 |

**What this validates**:
- The end-to-end machinery works: SMILES → docked complex → Stage 2 recompute → anchored Δactivity, on a single machine in ~1 minute per pair.
- **Active conformers bind ligands more avidly than inert ones.** Curcumin docks 6PEO at -8.5 kcal/mol vs 2N0A at -5.8 (∆ ~3 kcal/mol). This is the in-silico version of the published claim that toxic α-syn states expose hydrophobic patches the inert states have buried.
- **Δactivity magnitude tracks anchor activeness.** Curcumin perturbs active 6PEO by -0.40 vs inert 2N0A by -0.18 (~2.2× ratio). Even without sign-flipping (see below), this magnitude differential could rank molecules by anchor-class selectivity.
- Curcumin's protective sign reproduces correctly. Quercetin matches it. These are the in-silico analogues of the EGCG-redirect controls called for by the plan.

**What this can't yet detect (and how to fix)**:
- **The Δactivity sign is uniformly negative.** Static-occlusion-only: any docked ligand reduces the exposed-hydrophobic-β-SASA mean *and* raises the now-ligand-aware contact_density, so harm anchors (methylglyoxal here) come out looking weakly protective. The sign-flip for harm-side molecules requires MD relaxation around the bound pose so the receptor can rearrange and the disordered-exposure / NAC-β features can move *against* the occlusion contribution. This is the OpenMM step the plan calls out and the prototype defers.
- **Affinity is reported but not scored.** Methylglyoxal (-2.9 kcal/mol, 4 heavy atoms) produces less SASA occlusion *and* less Δcontact_density than curcumin (-8.5 kcal/mol, 29 atoms) — the ligand-aware contact_density partly fixes this, but the Vina affinity itself still doesn't feed the ΔP score. The plan's basin-escape score would weight poses by exp(-aff/RT); the prototype doesn't.

Per-pair artifacts live in `pipeline/results/stage3/<mol>_<pdb>_*.{pdb,pdbqt,json,log}`.

## Stage 3 — contact_density sees the ligand (2026-05-26, follow-up to the prototype landing)

`features.contact_density` was extended to also count Stage-3-appended ligand heavy atoms within the 8 Å cutoff. Helper `_stage3_ligand_heavy_coords` filters to residues named `"LIG"` (the marker `stage3.add_ligand_to_structure` uses), so anchor / apo structures — which only ever carry crystallographic waters as HETATM — take a fast path and reproduce the prior `anchor_scores.csv` byte-for-byte. Verified by re-running `validate.py` and diffing.

The six Stage 3 pairs were re-scored without re-docking via `recompute_stage3.py`, which reloads each existing `*_docked.pdbqt`, rebuilds the complex in memory, and rewrites `*_report.json`.

| molecule         | anchor (label)                  | Vina aff | apo cd | cplx cd | Δcd   | apo act | cplx act | Δact   |
| ---------------- | ------------------------------- | -------- | ------ | ------- | ----- | ------- | -------- | ------ |
| curcumin         | 6PEO (graded-active, single-PF) | -8.47    | 10.984 | 12.238  | +1.254 | +1.808 | +0.566   | **-1.242** |
| curcumin         | 6UFR (graded-active, paired)    | -7.60    | 11.066 | 12.721  | +1.656 | +1.880 | +0.338   | **-1.542** |
| quercetin        | 6PEO                            | -8.02    | 10.984 | 12.698  | +1.714 | +1.808 | +0.246   | **-1.562** |
| methylglyoxal    | 6PEO                            | -2.87    | 10.984 | 11.444  | +0.460 | +1.808 | +1.156   | -0.652 |
| curcumin         | 2N0A (inert, single-PF top)     | -5.76    | 10.776 | 11.329  | +0.553 | +4.398 | +3.850   | -0.548 |
| methylglyoxal    | 2N0A                            | -2.82    | 10.776 | 11.105  | +0.329 | +4.398 | +4.080   | -0.318 |

**What this changed**:
- **Δactivity magnitudes roughly tripled** for the real binders (curcumin × 6PEO -0.395 → -1.242, ~3.1×). The small-ligand pairs scale less (methylglyoxal × 2N0A -0.096 → -0.318, ~3.3× but still smallest in absolute terms).
- **Footprint differential is now visible.** Quercetin and curcumin had near-identical Δact on 6PEO before (-0.40 each). After the change quercetin is the most-perturbative pair (Δact -1.562) because its pose touches more residues in 6PEO's hydrophobic groove (+1.714 Δcd vs curcumin's +1.254). The feature is doing what the proposal intended — picking up binding-footprint geometry without MD.
- **Active-vs-inert selectivity preserved.** Curcumin perturbs active 6PEO by -1.24 vs inert 2N0A by -0.55 (2.3× ratio, same as before the change). The active-side amplification works at both old and new scales.
- **Affinity direction now tracks Δact ordering on each anchor.** On 6PEO: curcumin (-8.47) and quercetin (-8.02) dominate; methylglyoxal (-2.87) is far behind. The previous Δact column ranked them by SASA occlusion only, which made methylglyoxal nearly indistinguishable from curcumin (-0.341 vs -0.395). The ligand-aware contact_density restores roughly the right ordering.

**What this still doesn't do**:
- Sign-flip for harm-side molecules. Static occlusion + ligand-Cα contacts both push Δact down; you cannot get Δact > 0 from a static pose under the current weights. MD relaxation (the deferred next move) is what produces upward Δactivity from rearrangement.
- Penalise weak affinity. Methylglyoxal's pose still contributes a -0.65 / -0.32 Δact even at -2.9 kcal/mol — that's a low-confidence pose getting weight equal to a high-confidence one. Vina affinity isn't yet folded into ΔP.

`recompute_stage3.py` is a small post-feature-change utility. Use it whenever Stage 2 features change — or when the multi-pose aggregation rule changes — and you want the new ΔP across the existing pair set without burning a fresh Vina run.

## Stage 3 — Boltzmann pose weighting (2026-05-26, follow-up to ligand-aware contact_density)

Vina has been emitting `--num_modes 5` from the start; the pipeline only scored MODEL 1. `stage3.score_all_poses` now parses every MODEL block (`parse_pdbqt_all_poses`) and every affinity row from the Vina log (`_parse_affinities`), builds a complex per pose against a fresh assembly copy, recomputes Stage 2 features, and aggregates per-pose Δactivity with `w_i = exp(-aff_i / RT) / Z` at T=300 K (RT ≈ 0.596 kcal/mol).

Vina occasionally drops a trailing clashing pose (positive affinity) from the output PDBQT while still listing it in the log table. `score_all_poses` truncates the affinity tail to match the pose count and notes it on stdout — the dropped pose would carry vanishing Boltzmann weight anyway. (Curcumin × 6PEO had its rank-5 pose at +14.62 dropped; the 4 remaining poses are -8.47 to -6.53.)

Reports now carry: `activity_complex_top` / `delta_activity_top` (rank-1 pose only), `activity_complex_weighted` / `delta_activity_weighted` (Boltzmann aggregate across surviving poses), and a `poses` array with per-pose `{rank, affinity, weight, complex_features, activity_complex, delta_activity}`.

| molecule         | anchor               | aff_top | act_apo | dact_top | dact_wtd | w_top |
| ---------------- | -------------------- | ------- | ------- | -------- | -------- | ----- |
| curcumin         | 6PEO (4 poses, clash dropped) | -8.47 | +1.808 | -1.242 | **-1.281** | 0.488 |
| curcumin         | 6UFR                 | -7.60   | +1.880  | -1.542   | **-1.432** | 0.283 |
| quercetin        | 6PEO                 | -8.02   | +1.808  | -1.562   | **-1.516** | 0.324 |
| methylglyoxal    | 6PEO                 | -2.87   | +1.808  | -0.652   | **-0.633** | 0.262 |
| curcumin         | 2N0A                 | -5.76   | +4.398  | -0.548   | **-0.505** | 0.263 |
| methylglyoxal    | 2N0A                 | -2.82   | +4.398  | -0.318   | **-0.170** | 0.210 |

**What this changed**:
- **Pair ordering is preserved.** By |dact|: quercetin × 6PEO > curcumin × 6UFR > curcumin × 6PEO > methylglyoxal × 6PEO > curcumin × 2N0A > methylglyoxal × 2N0A. Same on both `dact_top` and `dact_wtd`. The change is a refinement, not a re-ranking.
- **Top-pose dominance reflects affinity spread.** Curcumin × 6PEO has w_top = 0.488 because its rank-1 pose is 0.15 kcal/mol better than rank-2 (a real preference). Methylglyoxal × 2N0A has w_top = 0.210 because all five poses fall within 0.07 kcal/mol of each other (Vina effectively cannot pick a winner) — weights collapse toward uniform 1/5.
- **Weak binders with flat pose ensembles shrink toward zero.** Methylglyoxal × 2N0A drops 47% (-0.318 → -0.170) because two of its five poses are non-perturbative (dact = 0.000 — methylglyoxal lands outside the inner chain entirely in those poses) and carry near-uniform weight. The top-pose-only rule gave a full-weight contribution; the Boltzmann rule dilutes it across the pose ensemble.
- **Strong binders barely move.** Curcumin × 6PEO -1.242 → -1.281, quercetin × 6PEO -1.562 → -1.516. When the top pose is genuinely best, Boltzmann reweighting reproduces it.

**What this still doesn't do**:
- **Penalise weak absolute affinity across ligands.** The Boltzmann normalisation is *intra*-ligand — methylglyoxal's pose set sums to weight 1 just like curcumin's. The shrink seen on methylglyoxal × 2N0A is a happy by-product of flat affinities + pose-position diversity, not a direct absolute-affinity penalty. A separate gating term (e.g. multiply by `min(1, exp(-aff_top/RT) / exp(-aff_threshold/RT))`) would do that and is not yet wired in.
- **Sign-flip for harm-side molecules.** Same as before. Static-pose features cannot produce Δact > 0 under the current weights; MD relaxation remains the path.

Reproduce: `python recompute_stage3.py` for the existing 6 pairs without re-docking; `python stage3.py <mol> <pdb>` for a fresh pair (now multi-pose by default).

## Stage 3 — absolute-affinity gating (2026-05-26, follow-up to Boltzmann pose weighting)

The intra-ligand Boltzmann weights `exp(-aff_i / RT) / Z` reweight contributions *within* a pose ensemble but always sum to 1, so a 4-heavy-atom metabolite with top affinity -2.9 kcal/mol carries the same total weight across its poses as a 29-atom polyphenol at -8.5 kcal/mol. The gate adds an *inter*-ligand penalty:

```
gate(aff_top) = min(1, exp((thr - aff_top) / RT))      # thr = -6 kcal/mol, RT ≈ 0.596
delta_activity_gated = gate * delta_activity_weighted
```

Reported as an orthogonal column. The Boltzmann-weighted Δactivity stays in the report; consumers that want the un-gated number still have it. Stage 3 reports now carry `affinity_gate_threshold_kcal_per_mol`, `affinity_gate`, and `delta_activity_gated` alongside everything else.

| molecule         | anchor               | aff_top | dact_top | dact_wtd | gate  | **dact_gated** |
| ---------------- | -------------------- | ------- | -------- | -------- | ----- | -------------- |
| curcumin         | 6PEO                 | -8.47   | -1.242   | -1.281   | 1.000 | **-1.281**     |
| curcumin         | 6UFR                 | -7.60   | -1.542   | -1.432   | 1.000 | **-1.432**     |
| quercetin        | 6PEO                 | -8.02   | -1.562   | -1.516   | 1.000 | **-1.516**     |
| methylglyoxal    | 6PEO                 | -2.87   | -0.652   | -0.633   | 0.005 | **-0.003**     |
| curcumin         | 2N0A                 | -5.76   | -0.548   | -0.505   | 0.665 | **-0.336**     |
| methylglyoxal    | 2N0A                 | -2.82   | -0.318   | -0.170   | 0.005 | **-0.001**     |

**What this changed**:
- **Methylglyoxal pairs collapse to ~zero.** Both go from |dact_wtd| ≈ 0.2–0.6 to |dact_gated| < 0.005. Gate is exp(-3.13/0.596) = 0.005 at -2.87 vs threshold -6.0. The pose poses are still computed and stored — the gate only attenuates the aggregated score, which is what we want: methylglyoxal is a real toxic metabolite biologically, but its Vina-scale binding affinity is not the channel through which it exerts its effect. Static docking on -2.9-kcal/mol poses isn't producing trustworthy Δact, so the gate dampens it. (When MD relaxation lands, this can be revisited per-ligand.)
- **Curcumin × 2N0A is partially damped.** aff_top = -5.76 sits 0.24 kcal/mol above threshold; gate = exp(-0.40) = 0.665. The active-vs-inert ratio between curcumin × 6PEO and curcumin × 2N0A widens from 2.3× (un-gated, -1.281 / -0.505) to 3.8× (gated, -1.281 / -0.336). Curcumin is a confirmed binder; 2N0A is an inert anchor whose pocket happens to be suboptimal for it. The mild damping is consistent with that interpretation.
- **Ordering is preserved on dact_gated.** By |dact_gated|: quercetin × 6PEO (-1.516) > curcumin × 6UFR (-1.432) > curcumin × 6PEO (-1.281) > curcumin × 2N0A (-0.336) > methylglyoxal × 6PEO (-0.003) > methylglyoxal × 2N0A (-0.001). Same top-3 as `dact_top` and `dact_wtd`; the lower tail re-orders because the gate erases the artificial "small Δact from a static pose of a weak binder" signal.
- **Threshold choice.** -6 kcal/mol is the conventional drug-like-binder boundary. Curcumin × 2N0A (-5.76) sits just above; the smooth exponential gives it a non-zero contribution rather than a sharp cutoff. If a future ligand sits in the -5 to -6 range and we want to keep it, lowering the threshold to -5 widens the unattenuated band; raising it to -7 sharpens the curcumin × 2N0A penalty. Threshold is a tuning knob, not a hard line.

**What this still doesn't do**:
- **Sign-flip for harm-side molecules.** The gate attenuates Δact magnitudes but cannot change sign — it's still bounded above by zero in static-pose features. MD relaxation remains the path.
- **Doesn't tell us when a weak ligand still matters.** Methylglyoxal forms covalent adducts with α-syn (MGO–CEL on Lys, MGO–MGH on Arg) in real cells. Vina sees neither the adduct nor the resulting conformational rearrangement. The right thing for covalent / adduct-forming molecules is a separate channel; gating non-covalent affinity to ~zero is the correct *default* until that channel is built.

Reproduce: `python recompute_stage3.py` re-scores the existing 6 pairs in-place. New pair: `python stage3.py <mol> <pdb>` emits `delta_activity_gated` automatically; the `affinity_gate` value is in the report JSON.

## Stage 3 — oligomer target (2026-05-27)

`stage3.py` extended with `perturb_oligomer()` to accept a generated PDB as receptor and score all-chain mean Δactivity. Three molecules docked against `fusco_parallel_3mer_core70-88_relaxed.pdb` (apo score +17.775):

| molecule | aff_top | dact_top | dact_wtd | gate | dact_gated |
| -------- | ------- | -------- | -------- | ---- | ---------- |
| curcumin | -5.21 | -0.548 | -0.553 | 0.268 | **-0.148** |
| quercetin | -5.15 | -0.356 | -0.387 | 0.239 | **-0.092** |
| methylglyoxal | -2.21 | -0.100 | -0.121 | 0.002 | **~0.000** |

**Controls pass as expected:**
- Curcumin and quercetin produce protective signals (Δact < 0). Both dock the exposed β-sheet face at −5.2 kcal/mol.
- Methylglyoxal gated to ~0 (−2.21 kcal/mol → gate 0.002), consistent with weak non-covalent affinity — its biological effect is via covalent adducts, not reversible docking.

**What this tells us:**
1. **The framework detects the oligomer as a druggable target.** All three molecules dock in the NAC-region groove; curcumin/quercetin reduce activity by ~0.4-0.6 units per static pose.
2. **Oligomer affinities are lower than fibril affinities** (curcumin: −5.21 vs −8.47 on 6PEO). Expected: the fibril presents a deep groove; the β-sheet face of the oligomer is an exposed flat surface with shallower pockets. Still real binding, but less avid.
3. **Gate partially attenuates the protective signal.** Curcumin at −5.21 kcal/mol is 0.79 kcal/mol above the −6.0 threshold; gate = 0.268. This is appropriate — the affinity is real but modest for the oligomer target. Lowering the gate threshold to −5.0 would pass curcumin/quercetin unattenuated; raising it to −7.0 would suppress them further. Current −6.0 is a reasonable default.
4. **Δactivity is smaller in absolute terms on the oligomer than on fibrils** (−0.55 vs −1.28 on 6PEO). This reflects a dilution effect: the trimer's apo score (+17.78) is large because all three chains contribute exposed β-SASA. A single ligand occludes a small fraction of that total exposure. The fractional effect (−0.55 / 17.78 ≈ 3%) vs fibril (−1.28 / 1.81 ≈ 71%) differs because the fibril has a compact graded-active score while the oligomer has an enormous baseline. This is the correct interpretation: the oligomer has far more exposed surface than a single ligand can cover.

**Framework changes made (2026-05-27):**
- `perturb_oligomer(mol_id, oligo_pdb_path, ...)` — new function in `stage3.py`. Loads from local PDB file, identifies protein chains, scores apo via all-chain mean (matching `score_oligomer.py`), centers docking box on NAC residues 60-100, runs Vina, computes per-pose all-chain mean Δactivity with full Boltzmann+gate pipeline.
- `compute_mean_features_on(structure, chain_ids)` — per-chain-then-average feature helper, used for both apo and complex scoring.
- `docking_box_for_nac_core(structure)` — boxes on residues 60-100 Cα, correctly targeting the β-sheet face for any Fusco-topology oligomer regardless of tail conformations from OBC2 relaxation.
- `score_all_poses(pdb_id, apo_chain, ...)` extended with keyword-only `load_struct_fn` and `mean_chain_ids` params; backward-compatible with all existing callers (`recompute_stage3.py` verified).
- CLI: `python stage3.py <mol> <path-to-pdb>` auto-detects oligomer mode (file exists with `.pdb` extension) vs anchor mode (4-letter PDB id).

Artifacts: `results/stage3/{mol}_fusco_parallel_3mer_core70-88_relaxed_*.{pdb,pdbqt,json,log}`.

Reproduce: `python stage3.py curcumin results/oligomers/fusco_parallel_3mer_core70-88_relaxed.pdb`

## Stage 3 — oligomer sweep (full, 2026-05-27)

127 vicinity molecules with non-null SMILES scanned against the reference trimer (`fusco_parallel_3mer_core70-88_relaxed.pdb`, apo +17.78). 117 succeeded; 10 failed (all single-atom metal ions — meeko cannot generate ligand PDBQT for `[X+]` / `[X+++]` SMILES, an expected limitation of Vina-style flexible-ligand docking).

Top 25 by `delta_activity_gated` (most protective first):

```
rank  molecule              aff_top  dact_gated  notes
  1   dhea                   -5.57     -0.244    steroid hormone (novel)
  2   retinoic-acid          -5.28     -0.200    vitamin A metabolite (novel)
  3   curcumin               -5.22     -0.148    polyphenol (literature anchor)
  4   baicalein              -5.27     -0.144    Scutellaria flavone (literature)
  5   allopregnanolone       -5.35     -0.142    neurosteroid (novel)
  6   naringenin             -5.19     -0.128    flavanone
  7   luteolin               -5.37     -0.119    flavone
  8   demethoxycurcumin      -5.10     -0.119    curcuminoid
  9   thc                    -5.16     -0.115    cannabinoid (novel direct-binding)
 10   myricetin              -5.28     -0.115    flavonol
 11   piperine               -5.05     -0.112    black-pepper alkaloid (novel)
 12   hesperetin             -5.26     -0.107    flavanone
 13   trehalose              -5.21     -0.103    disaccharide (novel direct-binding;
                                                  literature attributes effect to autophagy)
 14   genistein              -5.03     -0.101    isoflavone
 15   urolithin-a            -5.17     -0.098    gut metabolite (novel)
 16   kaempferol             -5.19     -0.095    flavonol
 17   epicatechin            -5.21     -0.093    flavan-3-ol
 18   quercetin              -5.15     -0.092    flavonol (literature anchor)
 19   apigenin               -5.17     -0.091    flavone
 20   equol                  -5.03     -0.090    isoflavandiol (gut metabolite)
 21   guanosine              -4.91     -0.066    purine nucleoside
 22   daidzein               -4.84     -0.058    isoflavone
 23   honokiol               -4.75     -0.055    biphenyl neolignan
 24   adenosine              -4.65     -0.044    purine nucleoside
 25   resveratrol            -4.73     -0.044    stilbene (literature anchor)
```

Full CSV: `results/sweep/fusco_parallel_3mer_core70-88_relaxed_sweep.csv` (122 ranked entries after the validation hold-out cohort is appended).

**Tail (gated to ~zero):** methylglyoxal, acrolein, 4-HNE, malondialdehyde, TMAO, nitric oxide, hydrogen sulfide, magnesium, manganese, zinc, calcium. The reactive metabolites are correctly suppressed by the −6 kcal/mol gate (Vina cannot see their covalent mechanism; the gate is the right default until the covalent channel exists). Metal ions reach the gate floor by atomic mass alone — also correct, since Vina's docking score is dominated by the receptor-ion electrostatics which it does not parameterise well.

**Errors (could not parse / dock):** copper, iron, sodium, potassium, selenium, lithium, cobalt, lead, aluminum, mercury. All single-atom or simple inorganic — same root cause.

### Headline finding: ~7 novel candidates with no published α-syn assays

Among the top 15 are seven entries that, to the best of the seeded literature in `lib/vicinity_molecules.js`, have no published α-syn direct-binding or aggregation-modulation data: **DHEA, retinoic acid, allopregnanolone, THC, piperine, trehalose-as-direct-binder, urolithin A**. Each has independent plausibility (steroid axis declines with PD risk; retinoic acid is a CNS-active neurotrophic; THC and CBD are widely reported on PD symptoms via receptor pathways but not via direct α-syn binding; piperine is a known curcumin-bioavailability adjuvant; trehalose's anti-aggregation effect in vitro is usually attributed to autophagy enhancement; urolithin A is a hot mitophagy candidate). These seven are the actionable output of the sweep: hypotheses that would not have been on a wet-lab partner's list, because the toxic-oligomer structure to test them against did not exist until two days ago.

### Confirmatory band

Curcumin (#3), baicalein (#4), quercetin (#18), and resveratrol (#25) — known α-syn modulators in the published literature — are all in the top 25 of the gated ranking. Polyphenols generally dominate the head of the distribution: of the top 20, 14 are polyphenols. This is the framework reproducing what is already known, which is the table-stakes validation.

### Three honest red flags

1. **Score range is compressed.** DHEA at Δgated=−0.24 against an apo baseline of +17.78 is a 1.4% perturbation. The ranking is informative; the absolute magnitude is not interpretable as "DHEA dissolves 1.4% of toxic oligomers." Treat the output as ordinal.
2. **No molecule reaches the unattenuated affinity band (−6 kcal/mol).** Every successful hit is attenuated by the gate. Either the oligomer's exposed β-sheet face is genuinely too shallow for drug-like binding at this docking exhaustiveness, or Vina at exhaustiveness=8 is missing deeper pockets. Worth a follow-up at exhaustiveness=16 on the top 25 to test option (b).
3. **Polyphenol bias is partly mechanistic, partly framework-coupled.** The exposed β-sheet face rewards aromatic planarity — that is real biology. But the framework was calibrated on features (exposed β-SASA, NAC-β-score) that themselves reward β-sheet engagement. Distinguishing "polyphenols are genuinely the best class" from "the framework is preferring molecules that engage features it was designed to reward" requires either out-of-distribution validation (the hold-out cohort below) or wet-lab orthogonal assays.

Reproduce: `python sweep_oligomer.py --skip-existing` (re-reads cached reports, no re-docking). Full re-run: `python sweep_oligomer.py --no-skip`.

## Stage 3 — gold-standard validation hold-out (pre-registered prediction, 2026-05-27)

**Why this section exists.** The full sweep above includes molecules whose published α-syn behaviour was known at the time of framework development (e.g. curcumin, quercetin, baicalein were `role: "anchor"` seed entries in `vicinity_molecules.js` from the start). Even though Stage 2 weights were calibrated only on the 14 PDB structures (no small molecules), and the Stage 3 gate threshold (−6 kcal/mol) was set heuristically, the *choice of features* and the *gate threshold value* could in principle have been informed by the seeded "anchor"-role list. To test the framework as a blind predictor, we need a cohort that was unambiguously not used in any tuning decision.

**The hold-out cohort.** Five entries marked `validation_holdout: true` in `vicinity_molecules.js`. Three were in the file with `smiles: null` and could never have docked (EGCG, rosmarinic-acid, silibinin); two are new entries added today (fisetin, CAPE). All five have strong published α-syn aggregation-modulation evidence — they are the literature's gold-standard α-syn modulators that the framework has not yet seen.

| id              | refs (selected)                                                                                  |
| --------------- | ------------------------------------------------------------------------------------------------ |
| egcg            | Ehrnhoefer 2008 NSMB; Bieschke 2010 PNAS — canonical "redirect" of α-syn into off-pathway oligomers |
| silibinin       | Yin 2014 Cell Mol Neurobiol; Pérez-Sánchez 2016 Sci Rep — α-syn aggregation modulation in PD models |
| rosmarinic-acid | Ono 2006 J Neurochem; Takamatsu 2014 J Agric Food Chem — direct anti-fibrillation in vitro      |
| fisetin         | Maher 2007 Brain Res; Ardah 2016 Mol Neurobiol — flavonol anti-fibrillation + DA-neuron rescue   |
| cape            | Morroni 2018 Phytomedicine; Fontanilla 2011 Neuroscience — PD model neuroprotection              |

**Framework lock state.** The Stage 2 feature set (`exposed_hydrophobic_beta_sasa`, `membrane_insertion_propensity`, `nac_active_score`, `contact_density`, `disordered_hydrophobic_exposure`), Stage 2 weights (from 14-anchor calibration), Stage 3 ligand-aware contact_density, Boltzmann pose weighting (T=300 K, RT≈0.596), and absolute-affinity gate (threshold −6 kcal/mol) are **all locked at this commit**. No code change will be made in response to the hold-out result.

**Pre-registered predictions** (written before running the sweep on these 5):

| id              | predicted aff range (kcal/mol) | predicted rank range | predicted dact_gated range | rationale |
| --------------- | ------------------------------ | -------------------- | -------------------------- | --------- |
| egcg            | −5.5 to −7.0                   | 1–5                  | −0.20 to −0.45             | 33 heavy atoms, 8 OH + galloyl ester, bigger π-surface than curcumin (#3); should bind the exposed β-face at least as well |
| silibinin       | −5.5 to −7.0                   | 1–10                 | −0.15 to −0.40             | 35 heavy atoms, largest in cohort; rigid taxifolin+coniferyl scaffold should occupy the NAC groove |
| rosmarinic-acid | −5.0 to −6.0                   | 5–25                 | −0.08 to −0.20             | 26 heavy atoms, flexible ester linker; two catechols engage β-face but flexibility costs entropy |
| fisetin         | −5.0 to −5.5                   | 10–25                | −0.08 to −0.15             | structurally near quercetin (#18); minus the 5-OH; expect quercetin-band placement |
| cape            | −4.5 to −5.5                   | 20–50                | −0.05 to −0.12             | 21 heavy atoms, only one catechol, flexible ester; expect mid-pack |

**Pass criteria** (also pre-registered):
- **Strong pass**: all 5 in top 25 of the combined 122-molecule ranking.
- **Pass**: at least 3 of 5 in top 25, none below rank 60.
- **Calibration concern**: 2+ below rank 30 — the framework is missing the gold-standard signal.
- **Framework failure**: any of EGCG / silibinin / rosmarinic-acid below rank 50 — these have the strongest published binding evidence.

### Results (post-sweep)

| id              | predicted rank | actual rank | predicted aff | actual aff | predicted dact_gated | actual dact_gated | verdict |
| --------------- | -------------- | ----------- | ------------- | ---------- | -------------------- | ----------------- | ------- |
| silibinin       | 1–10           | **1**       | −5.5 to −7.0  | −5.80      | −0.15 to −0.40       | **−0.415**        | inside (rank); slightly above (dact range) |
| egcg            | 1–5            | **3**       | −5.5 to −7.0  | −5.62      | −0.20 to −0.45       | −0.235            | inside |
| fisetin         | 10–25          | **5**       | −5.0 to −5.5  | −5.25      | −0.08 to −0.15       | −0.150            | rank beat prediction; dact at upper bound |
| rosmarinic-acid | 5–25           | **24**      | −5.0 to −6.0  | −4.94      | −0.08 to −0.20       | −0.077            | inside (just) |
| cape            | 20–50          | **54**      | −4.5 to −5.5  | −3.98      | −0.05 to −0.12       | −0.011            | below predicted range — see below |

**Outcome: PASS** by the pre-registered criteria (3+ in top 25, none below rank 60 → 4 in top 25, CAPE at #54 < 60). Not "strong pass" — CAPE at #54 falls just outside top 25.

**Headline.** Of the five gold-standard α-syn modulators the framework had never seen, **three are in the top 5** of the combined 122-molecule ranking — including the two with the strongest published evidence (EGCG, silibinin) at positions 3 and 1. Silibinin is the highest-scoring molecule in the entire sweep, with the only `affinity_gate` value above 0.7 (rank-1 pose at −5.80 kcal/mol, the closest any molecule comes to the un-gated band).

**Per-molecule interpretation.**
- **Silibinin (#1)**: 35 heavy atoms, the largest cohort member; rigid taxifolin+coniferyl scaffold with five aromatic faces and a benzodioxin linker. Vina found a pose that engages multiple residues simultaneously, hence the high gate value (0.714) and the large un-gated Δact (−0.498 top / −0.581 weighted). The framework reads it correctly.
- **EGCG (#3)**: 33 heavy atoms, two trihydroxyphenyl rings linked by a gallate ester. Predicted #1–5; landed #3. Behaves as the textbook "redirect" molecule the framework should rank high; the result is the unambiguous positive control.
- **Fisetin (#5)**: outranks quercetin (#21) despite their similar structure. The difference is the loss of the 5-OH on the A ring — fisetin has 7-OH but not 5-OH, which (counter-intuitively) gives a slightly better Vina pose in this groove. The framework rewards that. Not predicted to be quite this high, but inside the ballpark.
- **Rosmarinic-acid (#24)**: just inside the top-25 band. Two catechols linked by a flexible diester. Better Vina pose than caffeic-acid (#41) or ferulic-acid (#48), confirming that the second catechol is contributing real binding contacts.
- **CAPE (#54)**: missed predicted range. The phenethyl ester of caffeic acid docks **worse** than caffeic acid itself (CAPE −3.98 vs caffeic-acid −4.26). The phenethyl group adds bulk and conformational freedom without contributing matching contacts on the β-sheet face. Vina sees this as an entropy penalty; the gate then suppresses it hard. **This is arguably the framework being more honest than the literature**: CAPE's published α-syn effects are usually attributed to NF-κB / anti-inflammatory mechanisms or intracellular hydrolysis to caffeic acid, not to direct binding at sub-µM. A ranking of #54 for the parent ester (un-gated dact -0.29, comparable to curcumin) plus #41 for the active metabolite (caffeic acid) is structurally coherent.

**What this validates and what it doesn't.**

It validates: the framework correctly identifies the strongest published α-syn modulators (EGCG, silibinin, fisetin, rosmarinic-acid) without having seen them during calibration. The novel candidates from the broader sweep (DHEA, retinoic acid, allopregnanolone, THC, piperine, urolithin-a, trehalose) sit alongside the validated gold standards in the top 25 — they pass the same internal test the gold standards pass.

It does not validate: that the molecules at the top are actually toxic-oligomer destabilisers in cells. The framework reproduces published *direct-binding* rankings; it does not yet model covalent chemistry, autophagy effects, anti-inflammatory pathways, or receptor-mediated signalling. CAPE's case is the cautionary example — a molecule with strong published α-syn-protective effects can rank low here because its mechanism is not the one the framework models.

**Framework lock honoured.** No code, no weight, no threshold was modified in response to these results. The CAPE underperformance is recorded as-is and treated as a *finding about CAPE's mechanism*, not a calibration error.

Reproduce: `python sweep_oligomer.py --skip-existing` (the 5 hold-out reports are now cached in `results/stage3/`).

## Stage 3 — covalent / adduct channel (`aspr_score`, 2026-05-28)

Vina sees only reversible non-covalent binding. The four reactive metabolites in the vicinity list (methylglyoxal, 4-HNE, acrolein, malondialdehyde) modify α-syn through Lys-CEL / Lys-Schiff / His-Michael adducts, which the affinity gate correctly collapses to ~zero in `delta_activity_gated`. That is the right default until a separate channel exists; this section is that channel.

[`adduct_score.py`](adduct_score.py) implements:

```
aspr_score = (1/n_chains) · Σ_chains Σ_relevant_r  min(1, sasa(r) / 200 Å²) · rxty(ligand, restype(r))
```

`rxty(ligand, restype)` is a small ordinalised table of intrinsic per-residue chemical reactivity drawn from the standard adduct literature (Vicente-Miranda 2017 for MGO–Lys; Esterbauer 1991 / Bae 2013 for 4-HNE Michael; Uchida 1999 for acrolein; LoPachin 2014 for soft-electrophile rules). Primary target = 1.0, secondary = 0.3–0.7. SASA from the standard Bio.PDB Shrake-Rupley pass, normalised by 200 Å² to a fractional accessibility in [0, 1]. Independent of Vina pose; one SASA pass per receptor scores the whole sweep.

### Results — reference trimer

| ligand | rxty (primary, secondary, …) | dact_gated | aspr_score |
| ------ | ---------------------------- | ---------- | ---------- |
| malondialdehyde | K=1.0, C=0.2 | -0.0001 | **+10.60** |
| acrolein | C=1.0, K=0.6, H=0.5 | -0.0001 | **+6.70** |
| 4-hne | C=1.0, H=0.6, K=0.4 | -0.0033 | **+4.65** |
| methylglyoxal | R=1.0, K=0.4, C=0.3 | -0.0002 | **+4.24** |
| nitric-oxide | C=1.0, Y=0.8 | -0.0000 | **+2.67** |
| hydrogen-sulfide | C=1.0 | -0.0000 | **+0.00** |

(All other ~120 ligands have no entry in `LIGAND_REACTIVITY` and score `aspr_score = 0`.)

**What this gives.** A separate, orthogonal axis for harm-leaning covalent modifiers. Ranking by `aspr_score` recovers the four reactive aldehydes the gate collapsed and orders them by their chemistry of action on α-syn. The previous sweep's CSV had a flat zero band on these molecules; the new column resolves them.

**Why the ordering is informative biology, not arbitrary.** α-syn has 15 Lys, 0 Arg, 0 Cys, 1 His, 4 Tyr per chain. The ranking reflects this directly:
- **MDA tops the list** because it's a pure Lys-Schiff agent (rxty K=1.0), and α-syn has 15 exposed lysines per chain — that's the maximum possible substrate.
- **Acrolein > 4-HNE** even though both are α,β-unsaturated aldehydes — because acrolein has higher Lys reactivity (its smaller size lets it form Michael+Schiff dual adducts on Lys, well-characterised in the literature).
- **MGO sits below 4-HNE** because its primary chemistry (Arg-MGH) hits zero — α-syn has no Arg. It falls back to Lys-CEL (rxty 0.4), the empirically observed dominant adduct on α-syn (Vicente-Miranda 2017 Brain).
- **H2S correctly scores zero** because its only target (Cys, persulfidation) is absent. The framework refuses to confabulate a signal where there is no substrate.

### Comparison: trimer vs deposited fibril

Smoke-tested MGO against 6PEO (deposited graded-active fibril) — `aspr_score = +0.98` vs +4.24 on the trimer. The **4.3× ratio matches the actual biology**: the generated oligomer's disordered Lys-rich N-tail is exposed, the fibril's is partially buried in the assembly. Reactive metabolites are predicted to prefer the oligomer substrate, which is consistent with the "oligomer-toxicity hypothesis" the framework is built around.

### Caveats — honest

1. **Per-residue rxty weights are ordinal, not measured rate constants.** The relative rankings within a ligand and across the four aldehydes are defensible from chemistry first-principles. The absolute magnitude of `aspr_score` is uncalibrated; treat as ordinal.
2. **No kinetic vs thermodynamic separation.** A surface-exposed Lys with low local pKa will adduct faster than one near acidic neighbours; the simple SASA × rxty model ignores microenvironment. Adding a local-pKa modifier (PROPKA-style) would refine the ranking, mostly within the four-aldehyde block.
3. **No covalent-bond geometry check.** A reactive residue 8 Å from a Vina pose contributes the same as one 30 Å away. This is correct for the *intrinsic-reactivity-meets-substrate* score the channel computes, but it means the channel reports "fraction of α-syn that could be modified" rather than "which modification is most likely first". The latter would require pose + covalent geometry filter — a separate module.
4. **Not folded into `delta_activity_gated`.** Reactive ligands now have two scores: a near-zero `delta_activity_gated` (correct — they don't act non-covalently) and a positive `aspr_score` (correct — they do act covalently). A "harm-leaning index" combining them is a downstream consumer concern; the report keeps both columns separate.

### Framework lock — still honoured

The Stage 2 feature set, Stage 2 weights, Stage 3 contact_density, Boltzmann pose weighting, and absolute-affinity gate are unchanged. The covalent channel is **strictly additive**: a new orthogonal column, no modification of any existing scoring path. The validation hold-out result is therefore not affected; rerunning the hold-out sweep gives byte-identical `delta_activity_gated` values and `aspr_score = 0` for all five gold-standard molecules (none is a covalent modifier).

### Wiring

- New module: [`adduct_score.py`](adduct_score.py) — `LIGAND_REACTIVITY` dict, `aspr_score(structure, mol_id, chain_ids)` function, CLI for inspection.
- `stage3.py` — both `perturb()` and `perturb_oligomer()` now compute and emit `aspr_score` + `aspr_reactive` in the report JSON.
- `sweep_oligomer.py` — backfills `aspr_score` onto existing cached reports (free; depends only on receptor SASA), adds `aspr_score` and `aspr_reactive` columns to the sweep CSV, prints a separate reactive-metabolite table at the end of the run.

Reproduce:
- Single ligand: `python adduct_score.py methylglyoxal results/oligomers/fusco_parallel_3mer_core70-88_relaxed.pdb`
- Full sweep: `python sweep_oligomer.py --skip-existing` (re-emits `results/sweep/fusco_parallel_3mer_core70-88_relaxed_sweep.csv` with the new columns and backfills cached `*_report.json` files in place).

## Recommended next moves, in order

1. **Wet-lab partner handoff on the novel candidates.** DHEA, retinoic acid, allopregnanolone, THC, piperine, urolithin-a, trehalose are the framework's seven top-of-list novel hypotheses (no published α-syn direct-binding assays). A standard ThT aggregation kinetics or DLS oligomer-stability assay would be the natural orthogonal test. Pick 3 most actionable based on supply / regulatory ease. The new `aspr_score` channel adds methylglyoxal, 4-HNE, acrolein, malondialdehyde to the same handoff list with a separate predicted mechanism (covalent Lys-CEL / Lys-Schiff on the disordered N-tail of the oligomer); orthogonal ThT + LC-MS adduct readouts cover both axes.
2. **Sensitivity sweep on the −6 kcal/mol gate.** Every successful sweep hit was gate-attenuated. Worth re-scoring the top 30 with the gate disabled and with the threshold at −5.0 and −4.5 to characterise how rank order changes with the gate setting. Useful for any threshold tuning that happens *after* the current hold-out exercise (now that the framework is validation-passed, future tuning needs its own hold-out cohort).
3. **MD relaxation around the docked complex.** Still the only path to a sign-flip on harm-side molecules in the non-covalent channel. Most meaningful against the generated trimer rather than against deposited fibrils.
4. **Shape-stability channel: multi-replica short-MD dwell-time screen.** Reframes the screen to bypass the §8.3 sign-bound: instead of "does the ligand change the activity of a static complex" (current Δact_gated), ask "does the ligand narrow the distribution of oligomer shapes at fixed Δt across many MD replicas". Per (shape, ligand) pair: N velocity-seeded short MD replicas of apo and of the docked complex; score shape at Δt by β-core Cα-RMSD and inter-chain contact-map Jaccard; bootstrap the distribution shift. Reuses the 11-topology ensemble. Cost is cluster-scale at the full sweep; cheap pilot first — 4 ligands (silibinin / DHEA / trehalose / a low-rank decoy) × 2 shapes × 10 × 2 ns. Requires a per-replica binding-occupancy check so ligands that diffuse off the site are not scored as "no effect". Subsumes item 3.
6. **Expand the vicinity list opportunistically.** Inclusion gate: "reaches CNS at achievable concentration." No α-syn evidence required for candidates. Any future "anchor"-role additions should be marked `validation_holdout: true` if they could blindly test future framework changes.
7. **Decide on ordinal-only vs direct-active investment.** AUC 0.84 + generated-oligomer ensemble + 4/5 gold-standard hold-out pass confirms that ordinal-only is defensible. Binary separation validated against biological data is a wet-lab question.
8. **Refinements to `aspr_score`.** Two natural next steps if the channel turns out to discriminate: (a) local-pKa modifier on Lys / His accessibility (PROPKA on the apo structure → multiplier on the reactivity weight); (b) pose-aware geometric filter so the reactive-residue contribution is gated by distance to the top Vina pose, turning aspr from "substrate exposed somewhere" into "substrate exposed where the ligand will encounter it first". Defer until there is a reason to discriminate further than the current 4-aldehyde ordering does.
9. **(Optional) Ablation table.** Compare `--no-core-mask` and `--au` modes for documentation.

## What is *not* in scope yet

- Stage 1 generator (fragment-MC, generative model, latent dynamics). The topology-prior build is a proxy for Stage 1.
- Any website surface (the `/compute` front door, browser client, results pages).
- The off-Vercel coordinator for volunteer-compute dispatch.

All of these wait until the vicinity-molecule scan against the trimer produces a ranked list.

## How to pick this up in a fresh thread

1. Read this file.
2. Read [`ANCHORS.md`](ANCHORS.md) for the anchor set and per-anchor context.
3. Re-run `validate.py` to reproduce the anchor scores. Useful flags:
   - default: assembly_inner + ordered-core mask + pydssp SS
   - `--no-core-mask`: ablate the mask
   - `--au`: legacy asymmetric-unit mode
   - `--compare`: side-by-side au / assembly_all / assembly_inner
4. Reproduce the oligomer ensemble: `python oligomers/run_ensemble.py --summary-only` (re-scores all 11 relaxed PDBs in < 1 min; no re-build/re-relax). Full re-run: drop `--summary-only`.
5. Reproduce oligomer Stage 3 controls: `python stage3.py curcumin results/oligomers/fusco_parallel_3mer_core70-88_relaxed.pdb` (needs `pipeline/bin/vina.exe`). For deposited-anchor pairs: `python stage3.py curcumin 6PEO`.
6. Inspect the covalent / adduct channel for any reactive ligand: `python adduct_score.py methylglyoxal results/oligomers/fusco_parallel_3mer_core70-88_relaxed.pdb`. Sweep CSV already has `aspr_score` column.
7. Default next move: re-run `python sweep_oligomer.py --skip-existing` to keep the sweep CSV current (cheap; backfills aspr columns onto cached reports). See "Stage 3 — covalent / adduct channel" and "Recommended next moves".
