# Pipeline status

Last update: 2026-05-30.

This file is the entry point for a fresh thread. It is a short index. **The GitHub issue tracker is the source of truth** for results, next moves, open decisions, and caveats. After reading this file, read [`ANCHORS.md`](ANCHORS.md) for anchor curation context, then [`../IN_SILICO_PLAN.md`](../IN_SILICO_PLAN.md) for the overall project shape, then [`oligomers/README.md`](oligomers/README.md) for the Track A direction.

## Project shape

The proposal is a three-step pipeline:

1. **Generate** candidate toxic α-syn oligomer conformations — the structures the PDB does not provide.
2. **Score** which generated shapes have toxicity features (exposed hydrophobic β-SASA, NAC accessibility, membrane-disruption geometry).
3. **Modulate** — find ligands that raise the free energy of the toxic shapes so they dissipate or revert to inert fibrils, and flag ligands that do the opposite (lower the free energy of toxic shapes and stabilise them) as anti-targets to avoid.

Current state: Track A (topology-prior coarse build + MD relaxation) produces step 1; Stage 2 features and ordered-core mask score step 2 with AUC 0.84 on 14 fibril anchors; Stage 3 docks vicinity-list molecules against the generated oligomer for step 3, with three orthogonal channels (`delta_activity_gated` for non-covalent rearrangement, `aspr_score` for covalent adduct propensity, affinity gate for absolute-affinity cutoff).

## Headline results

| Issue | Headline |
| ----- | -------- |
| [#6](https://github.com/xag/asyn-oligomer-screen/issues/6) | Validation hold-out PASS — silibinin #1, EGCG #3, fisetin #5, rosmarinic-acid #24, CAPE #54 |
| [#5](https://github.com/xag/asyn-oligomer-screen/issues/5) | Full 127-mol oligomer sweep — 7 novel non-polyphenol top hits (DHEA, retinoic acid, allopregnanolone, THC, piperine, trehalose, urolithin A) |
| [#7](https://github.com/xag/asyn-oligomer-screen/issues/7) | Covalent/adduct channel `aspr_score` — 4 reactive metabolites recovered in chemistry-consistent order |
| [#2](https://github.com/xag/asyn-oligomer-screen/issues/2) | Track A ensemble — 11 generated topologies all score above the best fibril anchor (threshold +8) |
| [#1](https://github.com/xag/asyn-oligomer-screen/issues/1) | Stage 2 anchor calibration — AUC 0.84 on 14 anchors |
| [#8](https://github.com/xag/asyn-oligomer-screen/issues/8) | Protofilament count is a real within-class signal (ρ_inert = −0.76) |

Full list: [issues with `result` label](https://github.com/xag/asyn-oligomer-screen/issues?q=is%3Aissue+label%3Aresult).

## Top next moves

| Issue | Move |
| ----- | ---- |
| [#11](https://github.com/xag/asyn-oligomer-screen/issues/11) | Wet-lab handoff on the 7 novel candidates + 4 reactive metabolites — **package drafted** ([docs/HANDOFF.md](docs/HANDOFF.md), `node docs/build_handoff.mjs`); awaiting wet-lab partner |
| [#14](https://github.com/xag/asyn-oligomer-screen/issues/14) | Shape-stability channel (multi-replica short-MD dwell-time) — harness in repo; apo chunking validated, complex side pending conda OpenFF |
| [#30](https://github.com/xag/asyn-oligomer-screen/issues/30) | Anti-target flagging — stabilisers of toxic shapes (symmetric output of [#14](https://github.com/xag/asyn-oligomer-screen/issues/14)) |
| [#12](https://github.com/xag/asyn-oligomer-screen/issues/12) | Sensitivity sweep on the −6 kcal/mol affinity gate |
| [#16](https://github.com/xag/asyn-oligomer-screen/issues/16) | Refine `aspr_score`: PROPKA local-pKa + pose-aware geometric filter |

Full list: [issues with `next-step` label](https://github.com/xag/asyn-oligomer-screen/issues?q=is%3Aissue+label%3Anext-step).

## Open decisions

| Issue | Decision |
| ----- | -------- |
| [#19](https://github.com/xag/asyn-oligomer-screen/issues/19) | Single- vs paired-protofilament handling: three-way labels or normalisation feature? |
| [#20](https://github.com/xag/asyn-oligomer-screen/issues/20) | Direct active anchors: build Fusco DARR models, contact authors, or accept ordinal-only? |
| [#22](https://github.com/xag/asyn-oligomer-screen/issues/22) | Ordinal-only vs direct-active investment |
| [#21](https://github.com/xag/asyn-oligomer-screen/issues/21) | Vicinity-list canonicalisation: SMILES verification + PMID/DOI scheme |

Full list: [issues with `decision` label](https://github.com/xag/asyn-oligomer-screen/issues?q=is%3Aissue+label%3Adecision).

## Active caveats

| Issue | Caveat |
| ----- | ------ |
| [#23](https://github.com/xag/asyn-oligomer-screen/issues/23) | Δact is sign-bound by static-pose features — cannot flip positive without MD |
| [#24](https://github.com/xag/asyn-oligomer-screen/issues/24) | Polyphenol bias is partly framework-coupled |
| [#27](https://github.com/xag/asyn-oligomer-screen/issues/27) | Score range compressed — treat outputs as ordinal |
| [#28](https://github.com/xag/asyn-oligomer-screen/issues/28) | No molecule reaches the unattenuated affinity band (-6 kcal/mol) |
| [#29](https://github.com/xag/asyn-oligomer-screen/issues/29) | Hold-out is feature-internal — not orthogonal-mechanism (CAPE underranked) |

Full list: [issues with `caveat` label](https://github.com/xag/asyn-oligomer-screen/issues?q=is%3Aissue+label%3Acaveat).

## What is *not* in scope yet

- Stage 1 generator (fragment-MC, generative model, latent dynamics). The topology-prior build is a proxy for Stage 1.
- Any website surface (the `/compute` front door, browser client, results pages).
- The off-Vercel coordinator for volunteer-compute dispatch. (The per-replica MD *chunk* it would dispatch now exists and is validated on a basic GPU — apo side, pip-only, see [#14](https://github.com/xag/asyn-oligomer-screen/issues/14); the coordinator and front door themselves are still not built.)

## How to pick this up in a fresh thread

1. Read this file.
2. Read [`ANCHORS.md`](ANCHORS.md) for the anchor set and per-anchor context.
3. Browse [GitHub issues](https://github.com/xag/asyn-oligomer-screen/issues) by label (`result`, `next-step`, `decision`, `caveat`) for full context on any specific topic.
4. Re-run `scoring/validate.py` to reproduce the anchor scores. Useful flags:
   - default: assembly_inner + ordered-core mask + pydssp SS
   - `--no-core-mask`: ablate the mask
   - `--au`: legacy asymmetric-unit mode
   - `--compare`: side-by-side au / assembly_all / assembly_inner
5. Reproduce the oligomer ensemble: `python oligomers/run_ensemble.py --summary-only` (re-scores all 11 relaxed PDBs in < 1 min). Full re-run: drop `--summary-only`.
6. Reproduce oligomer Stage 3 controls: `python screen/stage3.py curcumin results/oligomers/fusco_parallel_3mer_core70-88_relaxed.pdb` (needs `bin/vina.exe`). For deposited-anchor pairs: `python screen/stage3.py curcumin 6PEO`.
7. Inspect the covalent / adduct channel for any reactive ligand: `python screen/adduct_score.py methylglyoxal results/oligomers/fusco_parallel_3mer_core70-88_relaxed.pdb`. Sweep CSV has `aspr_score` column.
8. Default next move: re-run `python screen/sweep_oligomer.py --skip-existing` to keep the sweep CSV current (cheap; backfills aspr columns onto cached reports). Pick the highest-priority `next-step` issue from the tracker.
9. Shape-stability / dwell-time channel ([#14](https://github.com/xag/asyn-oligomer-screen/issues/14)): `python screen/dwell_time.py selftest` (bootstrap sanity, no MD) and `score` (aggregate existing replica trajectories) run in the pip venv. The MD replicas run as independent per-GPU chunks: `uv sync --group md`, then `python screen/md_relax.py --apo-pdb <core-truncated>.pdb --rect-box --seed N --traj-out ...` (apo side, pip-only; complex side needs the `ASYN_MD_PYTHON` conda env for OpenFF). See README §9.
