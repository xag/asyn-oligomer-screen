# data

Input files consumed by the pipeline.

<details>
<summary><b>Plain English</b></summary>

This folder holds the inputs: the list of ~190 candidate molecules to
test (in `vicinity_molecules.js`) and a downloaded cache of the 14 known
α-syn structures used for calibration (in `anchors/`). The candidate
list was assembled under one rule: each molecule must plausibly reach
the brain through diet, supplementation, the body's own metabolism,
lifestyle (exercise, sleep, stress, hormones), or environmental
exposure. No prior α-syn evidence is required — that's the *output* of
the screen, not an input. Each entry carries a chemical-structure code
(SMILES), a rough brain concentration estimate, and references.
</details>

| path | content |
| --- | --- |
| `vicinity_molecules.js` | 191-entry candidate-molecule list (id, SMILES, CNS-concentration estimate, references). The single inclusion criterion is that the molecule plausibly reaches the substantia nigra via diet, supplementation, endogenous synthesis, lifestyle, or environmental exposure. |
| `anchors/` | PDB-file cache for the 14 anchor structures. **Gitignored**; auto-populated by `scoring/anchors.py:fetch()` on first run. |

`vicinity_molecules.js` is a hand-authored JS object-literal list (rather
than JSON) so each entry can carry free-text comments alongside the
machine-readable fields. It is parsed by tolerant regex in
[screen/stage3.py](../screen/stage3.py) and
[screen/sweep_oligomer.py](../screen/sweep_oligomer.py); no Node round-trip is required.
