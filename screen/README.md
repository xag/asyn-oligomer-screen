# screen

Step 3 of the three-step pipeline: **perturbation screen** of candidate
small molecules against a target conformer (a generated oligomer or a
deposited fibril).

For each (molecule, target) pair, the screen docks the ligand with
AutoDock Vina, recomputes the [scoring/](../scoring/) activity score on
each pose, and reports Δactivity = activity(complex) − activity(apo)
under Boltzmann pose weights and an absolute-affinity gate. A separate
covalent-adduct channel (`adduct_score.py`) handles reactive
electrophiles that docking cannot see.

<details>
<summary><b>Plain English</b></summary>

This folder asks the central question of the project: for each candidate
molecule, does sitting on the toxic α-syn oligomer make it look *less*
toxic? It uses AutoDock Vina (a standard docking tool) to find the
molecule's preferred binding spot, then recomputes the toxicity score on
the resulting complex and compares it to the empty oligomer. Multiple
binding spots are weighted by how energetically favourable they are, and
weak binders are penalised. A second channel handles reactive molecules
(damaging metabolic byproducts) that bind permanently rather than
reversibly — these are anti-targets, not candidates.
</details>

## Reversible-binding channel

| file | role |
| --- | --- |
| `stage3.py` | one (molecule, target) pair: dock + score + report |
| `sweep_oligomer.py` | full vicinity-list sweep against one target → ranked CSV |
| `recompute_stage3.py` | re-score cached docked poses after a feature change (no re-docking) |

## Covalent-adduct channel

| file | role |
| --- | --- |
| `adduct_score.py` | per-residue SASA × intrinsic reactivity; pose-independent |

## Optional MD relaxation

`md_relax.py` and `md_stage3.py` run OpenMM MD around docked complexes.
They require a separate conda env with OpenMM + openff-toolkit +
openmmforcefields (the pip pipeline cannot host them cleanly); point
`$ASYN_MD_PYTHON` at that interpreter before invoking.

## Examples

```bash
.venv/bin/python screen/stage3.py curcumin results/oligomers/fusco_parallel_3mer_core70-88_relaxed.pdb
.venv/bin/python screen/sweep_oligomer.py --skip-existing
.venv/bin/python screen/adduct_score.py methylglyoxal 6PEO
```
