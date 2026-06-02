---
title: Asyn Dwell Broker
emoji: 🧬
colorFrom: indigo
colorTo: green
sdk: docker
app_port: 7860
pinned: false
short_description: Upload broker for the crowdsourced dwell-time screen (#43)
---

# Asyn dwell-time upload broker (#43)

The only writer to the [dwell-time result dataset](https://huggingface.co/datasets/xagg/asyn-dwell-results)'s
`submissions/` inbox, and the only holder of its write token — so contributors
need no Hugging Face account. A volunteer's notebook POSTs a finished chunk here
with the signed lease the `health` site issued; the broker verifies the lease
(shared `SCREEN_LEASE_KEY`), size-caps + hashes the files, and writes
`submissions/<chunk_id>/<pseudonym>/{<file>, meta.json}` in the shape
`hf_store ingest` reads. It does no scoring or acceptance — the
reputation-weighted, distinct-pseudonym quorum + spot-check all happen later in
`ingest`.

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
