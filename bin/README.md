# bin

Third-party binaries the pipeline shells out to.

<details>
<summary><b>Plain English</b></summary>

This folder is for the docking program (AutoDock Vina). Vina is the
standard open-source tool that, given a protein and a small molecule,
figures out where on the protein the molecule sits and how tightly it
binds. We don't ship it with this repository because of licensing — you
download the binary for your operating system from the AutoDock website
and drop it here before running the screen. If you're not running the
pipeline, you don't need this folder.
</details>

| file | source | used by |
| --- | --- | --- |
| `vina.exe` | AutoDock Vina 1.2.5, [github.com/ccsb-scripps/AutoDock-Vina/releases](https://github.com/ccsb-scripps/AutoDock-Vina/releases) | [screen/stage3.py](../screen/stage3.py) |

The Vina binary is not vendored — download the release for your platform
and drop it here (`vina.exe` on Windows, `vina` on Linux/macOS — adapt
the path in [screen/stage3.py](../screen/stage3.py)). Tested against 1.2.5; earlier 1.2.x is
likely compatible.
