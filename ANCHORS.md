# Anchor structures for Stage 2

These structures define what *active* and *inert* mean operationally for the α-synuclein classifier. Stage 2 must rank known-inert structures low and known-active structures high, with visible separation between the clusters, before any generator is built. If the features don't separate the anchors, the framework needs different features — or a different signal entirely.

This list is a working draft. Entries with confirmed PDB IDs verified against RCSB are marked `verified`. The remaining `verify` entries either need verification or do not exist in the PDB (see notes).

## What "passing" looks like

A weighted-feature score, computed on each anchor structure, that places:

- the high-confidence inert anchors at the low end of the activity axis,
- the high-confidence active anchors at the high end,
- with separation greater than within-group spread,
- and with at least one feature that monotonically tracks the active/inert axis on its own (interpretability check).

The graded-active anchors (mutant fibrils) should fall ordinally where the literature expects them. They are sanity checks, not pass/fail.

## A material gap in the active side

Atomic-resolution structures of toxic α-syn oligomers do not exist in the PDB. Verified during the first pass of curation:

- **Fusco 2017 Type B* oligomers** — characterised by ssNMR DARR correlations but no PDB or BMRB coordinate deposition. The published constraints are usable if we model coordinates ourselves; the atomic structure is not directly available.
- **Annular protofibrils, β-barrel pore models** — exist only as computational models in supplementary materials of various papers, not as deposited entries.
- **Anle138b-bound oligomer** (Antonschmidt 2022, Nat Commun) — the deposited PDB 8A4L is a *lipidic fibril* with anle138b in a cavity, not an oligomer.

This means the active side of the classifier is anchored only by graded-active mutant fibrils for now. Direct active anchors require one of: building Fusco's constraints into a model ourselves, contacting authors for coordinates, or — eventually — wet-lab collaborators contributing new structures.

This is a project-shape result, not a curation failure. We design the classifier knowing the active anchors are weak, and accept that the first pass validates only ordinal correctness across mutant severity, not active/inert separation.

## Inert anchors

| PDB | Method | n_PF (bio) | Description | Citation | Status |
| --- | --- | --- | --- | --- | --- |
| 1XQ8 | solution NMR | n/a | α-syn bound to SDS micelle; canonical helical conformation associated with the membrane-binding physiological state | Ulmer & Bax, J Biol Chem 2005 | verified |
| 2N0A | ssNMR | 1 | first atomic-resolution α-syn fibril; Greek-key topology over the NAC region | Tuttle et al., Nat Struct Mol Biol 2016 | verified |
| 6CU7 | cryo-EM | 2 | "rod" polymorph, full-length recombinant fibril | Li et al., Cell Research 2018 | verified |
| 6CU8 | cryo-EM | 2 | "twister" polymorph, same study | Li et al., Cell Research 2018 | verified |
| 6H6B | cryo-EM | 2 | recombinant 1-121 truncation fibril; paper calls it "cytotoxic" but the cited toxicity is for the truncated *protein* (accelerated aggregation), not these fibrils vs WT | Guerrero-Ferreira et al., eLife 2018 | verified — kept inert per framework rationale (mature fibril = endpoint) |
| 6XYO | cryo-EM | 2 | MSA Type I brain-derived fibril | Schweighauser et al., Nature 2020 | verified |
| 6XYP | cryo-EM | 2 | MSA Type II-1 brain-derived fibril | Schweighauser et al., Nature 2020 | verified |
| 8A9L | cryo-EM | 1 | Lewy fold from PD / PDD / DLB brain, 2.2 Å | Yang et al., Nature 2022 | verified |
| 8A4L | cryo-EM | 2 | lipidic fibril, polymorph L2A (anle138b-bound version exists) | Antonschmidt et al., Nat Commun 2022 | verified |

Rationale for treating mature fibrils as inert: the fibril is the assembly endpoint. Neurons coexist with fibril inclusions for years; the "fast fibril" hypothesis holds that accelerating fibril formation is *neuroprotective* by depleting toxic oligomers. Fibrils form one bound of the activity axis. The micelle-bound helical monomer (1XQ8) forms the other inert bound — it's the physiological conformation.

Note: 1XQ8 (physiological) and the fibril cores (disease endpoint) are both "inert" in the activity-score sense but are otherwise very different. A future, finer classifier may want to separate them; the current binary active/inert framing folds them together.

## Graded-active anchors (familial mutant fibrils)

These are not direct active anchors. The mutation shifts the conformer ensemble toward toxicity and accelerates pathology in carriers, so the fibrils carry a sub-toxic-pathway signature — but they are still fibrils, not oligomers. The classifier should rank them above WT fibrils, ordinally.

| PDB | Mutation | n_PF (bio) | Polymorph / context | Citation | Status |
| --- | --- | --- | --- | --- | --- |
| 6LRQ | A53T | 2 | full-length A53T fibril | Sun et al., Cell Research 2020 | verified |
| 7WO0 | A53T | 2 | A53T fibril induced by **Ca²⁺** | (newer) | verified — calcium is a vicinity molecule, so this doubles as a Stage 3 perturbation reference |
| 6UFR | E46K | 2 | "more stable, pathogenic" fibril, 2.5 Å — symmetric double protofilament (corrected 2026-05-26 from 1; the Boyer 2020 abstract describes it as a "symmetric double protofilament", confirmed by the geometric counter and the RCSB C2 biological assembly summary) | Boyer et al., PNAS 2020 | verified |
| 6PEO | H50Q | 1 | narrow polymorph (single protofilament) | Boyer et al., NSMB 2019 | verified |
| 6PES | H50Q | 2 | wide polymorph (two protofilaments) | Boyer et al., NSMB 2019 | verified — scores below 6PEO; confirmed as real signal driven by inter-protofilament burial, not an artifact |

**A30P** does not have a useful deposited structure. The A30P mutation is the most informative case theoretically — it slows fibril formation while stabilising off-pathway oligomers, dissociating fibril propensity from toxicity. But A30P fibrils adopt a WT-like Greek-key fold, and its toxicity-relevant oligomer has no PDB deposit. Skipped for now.

## Notes on the anchor strategy

- **The active side is currently graded only.** We have no direct active anchors. The classifier can be validated for ordinal correctness on mutant severity (A53T / E46K / H50Q should score ≥ WT fibrils) but not for the binary active/inert separation that the project actually needs.
- **Adding 7WO0 (A53T + Ca²⁺) is a near-future asset.** Calcium is in the candidate vicinity-molecule list. Stage 3 will eventually compare A53T-alone vs A53T-with-Ca²⁺ — and 7WO0 already provides the bound state experimentally. This is the first place where the anchor set and the candidate molecule list intersect.
- **Direct active anchors require new work.** Options: (a) build models from Fusco 2017 ssNMR constraints — possible but project-defining; (b) contact authors for coordinates; (c) wait for wet-lab partners. A combination of (a) and (c) is the realistic path.
- **Lipidic fibril (8A4L) is interesting context.** It captures α-syn in a lipid environment, even if as a fibril rather than an oligomer — useful for Stage 3 membrane-context perturbations.
- **Class imbalance.** Inert is heavily overrepresented. The classifier will need class weighting once we score across classes.
- **Protofilament count: real inert signal, anecdotal on graded-active.** After the 6UFR lit fix (1 → 2, verified against Boyer 2020), the within-class Spearman ρ between `n_protofilaments` and activity is **−0.76 inert** (n=8; 2N0A and 8A9L vs six paired anchors; group-mean gap +2.2 vs −1.6 = 3.7) and **−0.35 graded-active** (n=5; was −0.87 before the fix because 6UFR was wrongly grouped with the single-PFs). The graded-active "single-PF effect" now rests on a single anchor (6PEO), so treat it as unconfirmed until another genuine single-PF graded-active anchor is added. The recalibrated geometric column shows the same direction (ρ = −0.45 inert, −0.35 graded-active). See `scoring/analyze_protofilaments.py` and [issue #8](https://github.com/xag/asyn-oligomer-screen/issues/8).
- **Deposited vs biological protofilament count.** `n_protofilaments_deposited` from `protofilaments.count_protofilaments` reflects the lateral-stack count of the REMARK 350 biological assembly; `n_protofilaments` is the literature-curated biological count. After the 2026-05-26 recalibration of the geometric counter (short-pair-vector fibril axis + adaptive lateral threshold) and the 6UFR lit fix (1 → 2, per Boyer 2020 "symmetric double protofilament"), the two columns agree on 13 / 14 anchors. The remaining mismatch is **8A4L** (lit=2, dep=3): BIOMOLECULE 1 is explicitly PENTADECAMERIC with three lateral stacks, likely a lipid-mediated triplex of two fibrils rather than the canonical paired polymorph the paper centres on.

## Sources beyond the PDB

PDB is the primary source but not the only one. The pipeline should also accept:

- ssNMR-derived models published with coordinates at BMRB or in supplementary materials.
- Computational models published with coordinates in supplementary files.
- Coordinate sets reconstructed from constraints (e.g. Fusco DARR data → model).

A small loader extension will need to handle each source.
