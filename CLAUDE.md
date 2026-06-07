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
