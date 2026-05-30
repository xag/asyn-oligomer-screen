# CLAUDE.md — working notes for this repo

`asyn-oligomer-screen` is a three-step in-silico pipeline (generate toxic α-syn
oligomer shapes → score toxicity features → find ligands that destabilise them,
and flag ones that stabilise them as anti-targets). Outputs are
food/supplement/lifestyle/exposure guidance, **not drug candidates**. Start a
fresh thread by reading [`STATUS.md`](STATUS.md).

## Documentation system — keep the 4 layers in sync

The docs are layered, and a change that lands a feature must update all the
layers it touches **in the same change set** (the `#11` commit `62d9327` is the
model: it touched README + STATUS + docs together):

1. **GitHub issues = source of truth** for every result / next-step / decision /
   caveat. Each is a *labelled* issue (`result` / `next-step` / `decision` /
   `caveat`). New issues **lead with rationale** — *Why this matters* + *How it
   fits the 3-step shape* — not just what-to-do (see `.github/ISSUE_TEMPLATE`).
2. **STATUS.md** — curated index of the tracker: bump **`Last update:`** date,
   keep the mirror tables (Headline results / Top next moves / Open decisions /
   Active caveats) current, and add a "how to pick up" repro step for any new
   runnable surface.
3. **README.md** — the paper. When a channel lands, add a **§9 Reproduction**
   command block (one per channel, `.venv/bin/python ...`), update the
   "Status — work in progress" bullets, and keep the plain-English `<details>`
   mirror in sync.
4. **docs/** — `HANDOFF.md` and `index.html` are **generated** by `build_*.mjs`.
   Never hand-edit them; re-run the build script and commit the regenerated file.

**Docs-sync checklist when a `screen/` channel or next-step lands:**
- [ ] Update / close its **issue** (it's the source of truth).
- [ ] **STATUS.md**: bump the date + update the relevant table + add a how-to-pick-up step.
- [ ] **README**: add a §9 Reproduction block + update the Status bullets + plain-English mirror.
- [ ] Regenerate `docs/` via `build_*.mjs` if the handoff/display changed.
- [ ] Commit message references the **issue #** and names the doc deltas.

**Prefer subtractive edits.** Many doc updates should make the file *shorter*,
not longer. When a next-step, caveat, limitation, or "not in scope" item is
resolved or superseded, **delete it** from the README / STATUS lists rather than
appending a "done / now-implemented / superseded" note beside it — the closed
issue holds the history, so the living docs stay lean. The model is commit
`f7e2285`, whose whole README delta was `README.md | 1 -`: landing the channel
*removed* it from the next-steps list. Resist accreting "update" sentences onto
an existing bullet; rewrite or remove the bullet instead. Next-step / caveat /
not-in-scope lists should shrink as work lands, not grow.

## Conventions

- **Dependencies: `uv add`** (and `uv add --group <name>` for optional stacks),
  not `uv pip install`. The `md` group holds the pip-installable MD deps
  (`openmm`, `pdbfixer`).
- **MD environment is split.** The *apo* dwell-time / relaxation path runs
  pip-only (`uv sync --group md`, OpenCL/CPU/CUDA via `pick_platform`). The
  *docked-complex* path needs OpenFF/SMIRNOFF, which is **not pip-installable**
  (`openff-toolkit` is yanked from PyPI), so it needs a conda env pointed at by
  `$ASYN_MD_PYTHON`. Keep OpenFF imports lazy so the apo path never requires it.
- **Distributed-chunk direction:** the dwell-time MD is meant to split into
  independent per-replica chunks small enough for a volunteer's basic GPU
  (truncate to the NAC core + `md_relax --rect-box` → ~55k atoms), with central
  GPU-free scoring via `dwell_time.py score`. See [`STATUS.md`](STATUS.md)
  "not in scope yet" for the coordinator/website boundary.
- Windows host: scripts that print non-ASCII (β, Å, →) must reconfigure stdout
  to UTF-8 in `main()` (cp1252 console + redirect/capture otherwise crashes).
