---
title: Asyn Dwell Broker
emoji: 🧬
colorFrom: indigo
colorTo: green
sdk: docker
app_port: 7860
pinned: false
short_description: Write broker for the crowdsourced dwell-time screen
---

# Asyn dwell-time write broker

The only writer to the [dwell-time result dataset](https://huggingface.co/datasets/xagg/asyn-dwell-results),
and the only holder of its write token — so contributors need no Hugging Face
account. `health` is the read + identity + dispatch surface (no token, no GPU);
every endpoint here *mutates the dataset*, which is why it lives on the
token-holder.

## Endpoints

- `POST /api/screen/v1/results` — a finished chunk's outputs. The client sends
  the signed lease `health` issued (plus `wall_s`, a `manifest` of
  `{artifact_id: filename}`, the files, and optionally its identity `token`).
  The broker verifies the lease (shared `SCREEN_LEASE_KEY`), cross-checks the
  token against the lease if present, size-caps + hashes the files, and writes
  `submissions/<chunk_id>/<pseudonym>/{<file>, meta.json}` in the shape
  `hf_store ingest` reads. No scoring or acceptance — the reputation-weighted,
  distinct-pseudonym quorum + spot-check all happen later in `ingest`.
  `/submit` is a deprecated alias kept for one release.
- `POST /api/screen/v1/molecules` — a candidate-molecule proposal (JSON; `name`
  + `smiles` required). The SMILES must parse (rdkit) — it is the only field
  acted on automatically. The entry is appended to `molecules.json` as
  `source: "contributed"`, `prep: "awaiting"`; the rest of the metadata is
  advisory. The coordinator's prep pass docks awaiting entries into runnable
  chunks; moderation governs promotion to `primary`.

## Secrets (Space settings → Variables and secrets)

| Secret | Value |
| --- | --- |
| `HF_TOKEN` | a token with **write** to the dataset repo |
| `SCREEN_LEASE_KEY` | the exact string also set on the health site (lease verification) |
| `SCREEN_DATASET_REPO` | `xagg/asyn-dwell-results` |

## Files

`app.py` (FastAPI), `lease.py` (copy of `screen/lease.py` — keep in sync;
`tests/test_lease.py` pins the cross-language contract), `Dockerfile`,
`requirements.txt`. Provision/update from the repo with `python broker/provision.py <lease-key-file>`.
