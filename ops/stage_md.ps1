#requires -Version 5.1
<#
.SYNOPSIS
  Stage the OpenMM-MD runtime on one GPU node over password SSH (Posh-SSH):
  find a writable scratch dir, clone the repo, create a uv venv with the `md`
  group (pip-only openmm + pdbfixer), verify the CUDA platform, run the geometry
  selftest. Idempotent — re-run is a fast `git pull` + `uv sync`.

.EXAMPLE
  . .\ops\local.ps1
  .\stage_md.ps1 -TargetHost host01                       # uses $env:CLUSTER_SCRATCH
  .\stage_md.ps1 -TargetHost host01 -Scratch /data/scratch/you
#>
[CmdletBinding()]
param(
  [Parameter(Mandatory)][string]$TargetHost,
  [string]$User = $(if ($env:CLUSTER_USER) { $env:CLUSTER_USER } else { $env:USERNAME }),
  [string]$Domain = $env:CLUSTER_DOMAIN,
  [string]$Scratch = $env:CLUSTER_SCRATCH,                      # blank = auto-discover
  [string]$RepoUrl = 'https://github.com/xag/asyn-oligomer-screen.git',
  [string]$PyVersion = '3.12',
  [int]$ConnectTimeout = 10,
  [int]$CommandTimeout = 900                                # clone + uv sync can take minutes
)

if (-not (Get-Module -ListAvailable -Name Posh-SSH)) {
  Write-Error 'Posh-SSH not installed. Run:  Install-Module Posh-SSH -Scope CurrentUser'; return
}
Import-Module Posh-SSH -ErrorAction Stop

if (-not $Domain) { Write-Error 'Cluster config not set. Run:  . .\ops\local.ps1   (copy from ops/local.ps1.example)'; return }

$h = "$TargetHost.$Domain"
$cred = Get-Credential -UserName $User -Message "SSH password for $User@$h"
if (-not $cred) { Write-Error 'no credentials entered'; return }

# Staging block, run by the remote shell. Placeholders are substituted below.
$stage = @'
set -u
echo "host: $(hostname -f 2>/dev/null || hostname)"

echo "-- pick scratch --"
SCRATCH=""
for c in "__SCRATCH__" /data/work/$USER /data/$USER /data/scratch/$USER /data/tmp/$USER "$HOME/asyn-dwell"; do
  [ -z "$c" ] && continue
  if mkdir -p "$c" 2>/dev/null && touch "$c/.w" 2>/dev/null; then rm -f "$c/.w"; SCRATCH="$c"; break; fi
done
if [ -z "$SCRATCH" ]; then
  echo "  !! no writable scratch found. Ask for /data/work/$USER (chmod u+rwx), or pass -Scratch."; exit 3
fi
echo "  scratch: $SCRATCH"
df -h "$SCRATCH" 2>/dev/null | awk 'NR==2{print "  fs free: "$4" of "$2" ("$5" used)"}'

echo "-- repo --"
REPO="$SCRATCH/asyn-oligomer-screen"
if [ -d "$REPO/.git" ]; then
  echo "  exists -> git pull"; git -C "$REPO" pull --ff-only 2>&1 | sed 's/^/    /'
else
  echo "  cloning (shallow)"; git clone --depth 1 __REPOURL__ "$REPO" 2>&1 | sed 's/^/    /'
fi
cd "$REPO" || { echo "  !! cannot cd $REPO"; exit 4; }

echo "-- uv venv (managed Python __PYVER__) + md group (pip-only openmm/pdbfixer) --"
uv venv --python __PYVER__ .venv 2>&1 | sed 's/^/    /'
uv sync --group md 2>&1 | tail -6 | sed 's/^/    /'

echo "-- verify openmm + CUDA platform --"
OK=0
out=$(.venv/bin/python - <<'PYEOF' 2>&1
import openmm
from openmm import Platform
print("openmm", openmm.version.version)
print("platforms", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
PYEOF
)
echo "$out" | sed 's/^/    /'
echo "$out" | grep -q "CUDA" && OK=1

echo "-- dwell geometry selftest (no GPU) --"
.venv/bin/python screen/dwell_time.py selftest 2>&1 | sed 's/^/    /'

echo "-- micromamba (for one-time ligand parametrisation) --"
command -v micromamba >/dev/null 2>&1 && echo "  present: $(command -v micromamba)" || echo "  absent (will bootstrap in scratch at the prep step)"

echo "-- DONE --"
if [ "$OK" = 1 ]; then
  echo "  => STAGED. repo=$REPO  venv=$REPO/.venv  (CUDA platform present)"
else
  echo "  => STAGED but NO CUDA platform — check openmm wheel vs driver"
fi
'@

$stage = $stage.Replace('__SCRATCH__', $Scratch).Replace('__REPOURL__', $RepoUrl).Replace('__PYVER__', $PyVersion)

try {
  $s = New-SSHSession -ComputerName $h -Credential $cred -AcceptKey -ConnectionTimeout $ConnectTimeout -ErrorAction Stop
  try {
    $r = Invoke-SSHCommand -SessionId $s.SessionId -Command $stage -TimeOut $CommandTimeout
    if ($r.Output) { $r.Output }
    if ($r.Error)  { $r.Error | ForEach-Object { "  [stderr] $_" } }
  }
  finally { Remove-SSHSession -SessionId $s.SessionId | Out-Null }
}
catch { "  !! ${h}: $($_.Exception.Message)" }
