# CLAUDE.md — working notes for this repo

`asyn-oligomer-screen` is a three-step in-silico pipeline (generate toxic α-syn
oligomer shapes → score toxicity features → find ligands that destabilise them,
and flag ones that stabilise them as anti-targets). Outputs are
food/supplement/lifestyle/exposure guidance, **not drug candidates**. Start a
fresh thread by browsing the [GitHub issues](https://github.com/xag/asyn-oligomer-screen/issues)
(labelled `result` / `next-step` / `decision` / `caveat` — the source of truth)
and reading [`README.md`](README.md) (the paper) + [`ANCHORS.md`](ANCHORS.md).

## Documentation system — keep the 3 layers in sync

The **GitHub issue tracker is the single source of truth**; the repo docs do
**not** duplicate it (no mirror tables, no parallel status index — that path
leads to noodle soup and drift). A change that lands a feature updates only the
layers it touches, **in the same change set**:

1. **GitHub issues = source of truth** for every result / next-step / decision /
   caveat — *and* for out-of-scope boundaries. Each is a *labelled* issue
   (`result` / `next-step` / `decision` / `caveat`). New issues **lead with
   rationale** — *Why this matters* + *How it fits the 3-step shape* — not just
   what-to-do (see `.github/ISSUE_TEMPLATE`). Don't re-list any of this in the repo.
2. **README.md** — the paper. When a channel lands, add a **§9 Reproduction**
   command block (one per channel, `.venv/bin/python ...`) and keep the
   plain-English `<details>` mirror in sync. Don't maintain a next-steps / caveats
   list here — link to the issue label instead.
3. **docs/** — `HANDOFF.md` and `index.html` are **generated** by `build_*.mjs`.
   Never hand-edit them; re-run the build script and commit the regenerated file.

**Docs-sync checklist when a `screen/` channel or next-step lands:**
- [ ] Update / close its **issue** (it's the source of truth).
- [ ] **README**: add a §9 Reproduction block + keep the plain-English mirror in sync.
- [ ] Regenerate `docs/` via `build_*.mjs` if the handoff/display changed.
- [ ] Commit message references the **issue #** and names the doc deltas.

**Don't duplicate the tracker.** Results, next-steps, decisions, caveats, and
"not in scope yet" boundaries belong in the issues, not in a parallel list in the
repo. The moment such a list starts forming in a doc it drifts and goes stale —
link to the issue label instead.

**Status lives in one mutable slot — never append-only, never copied.** Most
status ("we're leaning 10–20 ns", "X is the gate", "blocked on Y") is judgment,
not something you can probe — so it *must* be written, and the only question is
where. Two habits create every stale-duplicate trace; both are banned:
- **Append-only writes** — issue *comments* and "update: …" notes. Old and new
  coexist with no signal which is current, and a cold reader quotes the
  confident-sounding stale one. **Issue comments are for discussion / decisions-
  in-progress only — never status, progress, or audit logs.** An issue's current
  state lives *only in its body*, overwritten in place. **Open/closed carries one
  meaning across every label: open = an open loop that needs attention, closed =
  settled** — so the open set *is* the worklist. `next-step` open = to-do,
  `decision` open = undecided, `caveat` open = active concern; each is **closed**
  when handled. A **`result` is a record of a completed run — it has no to-do, so
  it is filed *closed* on creation** (an open `result` reads as work you owe and
  pollutes the worklist). The README is the readable synthesis; closed `result`
  issues are the queryable record behind it.
- **Copies** — the same claim in a body *and* a comment, or in two issues, or
  issue + doc. They drift apart. One owner per fact; everywhere else links.

The cross-cutting **"what's next" is derived, not authored**: make blockers
explicit with a `Blocked-by #N` line in the body, then "what's next" = open
`next-step` issues with no open blocker (a `gh` query over the dependency graph).
Don't hand-maintain a critical-path narrative — that's the STATUS.md trap again.
Establish current *mechanical* state (env / artifacts / tests) by **running the
check**, never by reading what a past session wrote about it.

**Prefer subtractive edits.** Many doc updates should make the file *shorter*,
not longer. When a caveat, limitation, or reproduction note in the README is
resolved or superseded, **delete it** rather than appending a "done /
now-implemented / superseded" note beside it — the issue holds the history, so
the living docs stay lean. The model is commit `f7e2285`, whose whole README
delta was `README.md | 1 -`. Resist accreting "update" sentences onto an existing
bullet; rewrite or remove the bullet instead. And land a trivial subtractive doc
fix straight on `main` — don't strand it in an unmerged branch, or `main` (what
the reader sees) stays stale.

## Conventions

- **Dependencies: `uv add`** (and `uv add --group <name>` for optional stacks),
  not `uv pip install`. The `md` group holds the pip-installable MD deps
  (`openmm`, `pdbfixer`).
- **MD environment is split — but only for *building* a complex system.** All
  dwell/relaxation *dynamics* run pip-only (`uv sync --group md`, OpenCL/CPU/CUDA
  via `pick_platform`): the apo path natively, and the docked-complex path via a
  serialised OpenMM `System` (`md_relax --system-xml/--solvated-pdb`, #34/#37).
  Only the one-time ligand parametrisation (`md_relax --prepare-only`) uses
  OpenFF/SMIRNOFF, which runs in the conda MD env (`environment-md.yml`, located
  automatically by `screen/md_env.py`). Keep OpenFF imports lazy so every run
  path except `--prepare-only` stays pip-only.
- **Distributed-chunk direction:** the dwell-time MD is meant to split into
  independent per-replica chunks small enough for a volunteer's basic GPU
  (truncate to the NAC core + `md_relax --rect-box` → ~55k atoms), with central
  GPU-free scoring via `dwell_time.py score`. The coordinator/website boundary
  (volunteer-compute dispatch + contributor runtime) is out of scope — tracked in
  issues [#34](https://github.com/xag/asyn-oligomer-screen/issues/34) / [#35](https://github.com/xag/asyn-oligomer-screen/issues/35).
- Windows host: scripts that print non-ASCII (β, Å, →) must reconfigure stdout
  to UTF-8 in `main()` (cp1252 console + redirect/capture otherwise crashes).
