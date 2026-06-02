"""Cross-language lease verification (#43).

The golden token below was produced by the *real* health signer
(``lib/screen.js`` ``signLease``) with SCREEN_LEASE_KEY='test-shared-lease-key-123'
and a fixed payload. If this test passes, the Python broker verifies exactly
what the JS front door signs — the contract that lets the two services share a
secret but nothing else. Regenerate with:

    SCREEN_LEASE_KEY=test-shared-lease-key-123 node --input-type=module -e "..."

Run:  .venv/bin/python tests/test_lease.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "screen"))

import lease as L  # noqa: E402

KEY = b"test-shared-lease-key-123"
# JS-signed; payload: v=1, exp_id=exp-demo, chunk_id=seg__a__s1__0, kind=segment,
# pseudonym=abcd1234abcd1234, iat=1700000000000, exp=1700000600000, nonce=fixed-nonce-xyz
GOLDEN = ("eyJjaHVua19pZCI6InNlZ19fYV9fczFfXzAiLCJleHAiOjE3MDAwMDA2MDAwMDAsImV4cF9pZCI6"
          "ImV4cC1kZW1vIiwiaWF0IjoxNzAwMDAwMDAwMDAwLCJraW5kIjoic2VnbWVudCIsIm5vbmNlIjoi"
          "Zml4ZWQtbm9uY2UteHl6IiwicHNldWRvbnltIjoiYWJjZDEyMzRhYmNkMTIzNCIsInYiOjF9."
          "FvmBQZcGIJIPdPYSoKRi-Mbp50FayqBq_1DwDCICheI")
IAT = 1700000000000
EXP = 1700000600000


def test_python_verifies_a_js_signed_token():
    payload = L.verify_lease(GOLDEN, key=KEY, now_ms=IAT)
    assert payload is not None
    assert payload["chunk_id"] == "seg__a__s1__0"
    assert payload["pseudonym"] == "abcd1234abcd1234"
    assert payload["kind"] == "segment"


def test_expired_token_rejected():
    assert L.verify_lease(GOLDEN, key=KEY, now_ms=EXP + 1) is None
    assert L.verify_lease(GOLDEN, key=KEY, now_ms=EXP - 1) is not None


def test_wrong_key_rejected():
    assert L.verify_lease(GOLDEN, key=b"not-the-shared-key", now_ms=IAT) is None


def test_tampered_body_rejected():
    body_b64, sig_b64 = GOLDEN.split(".")
    # Flip a character in the body; signature no longer matches.
    bad_body = ("X" + body_b64[1:]) if body_b64[0] != "X" else ("Y" + body_b64[1:])
    assert L.verify_lease(f"{bad_body}.{sig_b64}", key=KEY, now_ms=IAT) is None


def test_garbage_rejected():
    for t in ["", "nodot", "a.b.c", None, 42]:
        assert L.verify_lease(t, key=KEY, now_ms=IAT) is None


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")


if __name__ == "__main__":
    _run_all()
