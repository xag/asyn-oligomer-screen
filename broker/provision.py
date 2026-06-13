"""Create/update the broker HF Space and set its secrets.

    python broker/provision.py <path-to-lease-key-file>

The lease key is read from a file (not the command line) so it never lands in
shell history. Uses the cached `hf auth login` token for both the API calls and
the Space's HF_TOKEN secret.
"""
from __future__ import annotations

import pathlib
import sys

from huggingface_hub import HfApi, get_token

SPACE = "xagg/asyn-dwell-broker"
DATASET = "xagg/asyn-dwell-results"
APP_URL = "https://xagg-asyn-dwell-broker.hf.space"


def main() -> None:
    lease = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8").strip()
    here = pathlib.Path(__file__).resolve().parent
    api = HfApi()

    api.create_repo(SPACE, repo_type="space", space_sdk="docker", exist_ok=True)

    files = {
        "app.py": here / "app.py",
        "Dockerfile": here / "Dockerfile",
        "requirements.txt": here / "requirements.txt",
        "README.md": here / "README.md",
        "lease.py": here.parent / "screen" / "lease.py",
    }
    for dest, src in files.items():
        api.upload_file(path_or_fileobj=str(src), path_in_repo=dest,
                        repo_id=SPACE, repo_type="space")
        print(f"  uploaded {dest}", flush=True)

    api.add_space_secret(SPACE, "HF_TOKEN", get_token())
    api.add_space_secret(SPACE, "SCREEN_LEASE_KEY", lease)
    api.add_space_secret(SPACE, "SCREEN_DATASET_REPO", DATASET)
    print("  secrets set: HF_TOKEN, SCREEN_LEASE_KEY, SCREEN_DATASET_REPO", flush=True)

    print(f"\nspace:   https://huggingface.co/spaces/{SPACE}")
    print(f"app url: {APP_URL}")


if __name__ == "__main__":
    main()
