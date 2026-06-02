"""Verify the work-assignment leases minted by the health front door (#43).

health (`lib/screen.js`) signs a lease as ``b64url(body).b64url(hmac_sha256(
SCREEN_LEASE_KEY, body))`` where ``body`` is a compact, sorted-key JSON object.
Verification only needs to HMAC the *exact decoded body bytes* — never to
reproduce the JSON canonical form — so this stays a faithful second
implementation with no serialisation drift. The broker (an HF Space) bundles a
copy of this module; ``SCREEN_LEASE_KEY`` is the shared secret between health
and the broker (health never shares its master ENCRYPTION_KEY).

The broker is the only writer to the dataset's ``submissions/``; it identifies
the contributor solely from the ``pseudonym`` inside a verified lease — it never
sees, or needs, the email.
"""
from __future__ import annotations

import base64
import hashlib
import hmac
import json
import os
import time


def lease_key() -> bytes:
    key = os.environ.get("SCREEN_LEASE_KEY")
    if not key:
        raise RuntimeError("SCREEN_LEASE_KEY is not set (shared secret with health).")
    return key.encode("utf-8")


def _b64url_decode(s: str) -> bytes:
    # JS base64url drops padding; restore it for Python's decoder.
    return base64.urlsafe_b64decode(s + "=" * (-len(s) % 4))


def verify_lease(token: str, *, key: bytes | None = None, now_ms: float | None = None) -> dict | None:
    """Return the lease payload if the token is well-formed, untampered and
    unexpired; else ``None``. Constant-time signature compare."""
    if not isinstance(token, str) or token.count(".") != 1:
        return None
    body_b64, sig_b64 = token.split(".")
    try:
        body = _b64url_decode(body_b64)
        sig = _b64url_decode(sig_b64)
    except Exception:  # noqa: BLE001
        return None
    expect = hmac.new(key or lease_key(), body, hashlib.sha256).digest()
    if not hmac.compare_digest(sig, expect):
        return None
    try:
        payload = json.loads(body)
    except Exception:  # noqa: BLE001
        return None
    now = time.time() * 1000 if now_ms is None else now_ms
    exp = payload.get("exp")
    if isinstance(exp, (int, float)) and now > exp:
        return None
    return payload
