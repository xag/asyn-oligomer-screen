# bin

Third-party binaries the pipeline shells out to.

| file | source | used by |
| --- | --- | --- |
| `vina.exe` | AutoDock Vina 1.2.5, [github.com/ccsb-scripps/AutoDock-Vina/releases](https://github.com/ccsb-scripps/AutoDock-Vina/releases) | [screen/stage3.py](../screen/stage3.py) |

The Vina binary is not vendored — download the release for your platform
and drop it here (`vina.exe` on Windows, `vina` on Linux/macOS — adapt
the path in [screen/stage3.py](../screen/stage3.py)). Tested against 1.2.5; earlier 1.2.x is
likely compatible.
