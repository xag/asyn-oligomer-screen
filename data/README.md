# data

Input files consumed by the pipeline.

| path | content |
| --- | --- |
| `vicinity_molecules.js` | 191-entry candidate-molecule list (id, SMILES, CNS-concentration estimate, references). The single inclusion criterion is that the molecule plausibly reaches the substantia nigra via diet, supplementation, endogenous synthesis, lifestyle, or environmental exposure. |
| `anchors/` | PDB-file cache for the 14 anchor structures. **Gitignored**; auto-populated by `scoring/anchors.py:fetch()` on first run. |

`vicinity_molecules.js` is a hand-authored JS object-literal list (rather
than JSON) so each entry can carry free-text comments alongside the
machine-readable fields. It is parsed by tolerant regex in
[screen/stage3.py](../screen/stage3.py) and
[screen/sweep_oligomer.py](../screen/sweep_oligomer.py); no Node round-trip is required.
