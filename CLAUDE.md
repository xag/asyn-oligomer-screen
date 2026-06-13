# CLAUDE.md — working notes for this repo

`asyn-oligomer-screen` is a three-step in-silico pipeline (toxic α-syn oligomer
shapes → toxicity scoring → ligands that destabilise them, flagging stabilisers
as anti-targets). **Outputs are food/supplement/lifestyle/exposure guidance, not
drug candidates** — keep that scope guardrail. Orient via [`README.md`](README.md)
(the paper) and [`ANCHORS.md`](ANCHORS.md).

## Conventions (things the code can't tell you)

- **Add deps with `uv add`** (`uv add --group <name>` for optional stacks), never
  `uv pip install`.
- **`docs/` is generated — regenerate, never hand-edit** (each generated file names
  its own build command in its header).
- **Prefer subtractive doc edits**: delete resolved/superseded notes rather than
  appending "done/now-implemented" beside them; land trivial doc fixes on `main`.

## Writing issues

Issues are the worklist — one open problem each, written plain. No templates, no
section scaffolding.

- **Title states the problem as a statement** — not a question, and not a neutral
  description of how something works. *"Covalent channel may flag protective
  molecules as harmful anti-targets"*, not *"Does the covalent channel measure
  harm?"* nor *"Covalent channel measures reactivity, not direction"*.
- **Put an action in the title only if it is concretely defined.** A real, scoped
  action can lead (*"…we need to work out how much data clears the noise"* → a
  power calculation). A vague one (*"check it survives"*, *"resolve the direction"*)
  does not — drop it and state only the problem.
- **Body gives the context** — the tension and why it matters, a sentence or two.
  It does not answer the title or re-pose it as a question.
- **Link only real issues.** Never cite an issue number that does not exist.
- **Close from evidence.** When an experiment in [`EXPERIMENTS.md`](EXPERIMENTS.md)
  resolves an issue, close it with a comment citing the experiment (E1, E2, …). If
  an experiment does not address the issue, do not cite it.
