# scoring

Step 2 of the three-step pipeline: a structure-based **activity score** for
α-syn conformers, calibrated on 14 deposited fibrils (graded-active vs
inert, pairwise AUC ≈ 0.84).

The score is a weighted z-combination of five per-conformer features
computed on the inner chain of the REMARK 350 biological assembly. It is
used (a) to calibrate against the deposited anchors, (b) to score
generated oligomer models (see [oligomers/](../oligomers/)), and (c) to
score apo-vs-complex Δactivity during the perturbation screen (see
[screen/](../screen/)).

## Library modules

| file | role |
| --- | --- |
| `anchors.py` | PDB structure loader with disk cache; the 14-anchor curation table |
| `assembly.py` | REMARK 350 biological assembly builder; inner-chain selection |
| `features.py` | the five per-conformer features + ordered-core mask |
| `classifier.py` | feature weights and the weighted-z `score_table` |
| `protofilaments.py` | geometric protofilament counter |

## Runners

```bash
.venv/bin/python scoring/validate.py             # anchor calibration → results/anchor_*.csv + plots
.venv/bin/python scoring/analyze_protofilaments.py   # within-class protofilament-count correlation
```

`validate.py` is the project's first go/no-go gate. Re-run it whenever a
feature, weight, or anchor changes.
