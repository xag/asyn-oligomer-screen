#requires -Version 5.1
<#
.SYNOPSIS
  Bootstrap the `asyn-md` conda env on a GPU node via micromamba (no admin),
  over password SSH. This env's conda-forge `openmm` ships the CUDA platform the
  pip wheel lacked, AND carries the OpenFF stack for ligand parametrisation — so
  it covers every MD step the dwell-time pilot shells out to.

.EXAMPLE
  . .\ops\local.ps1
  .\bootstrap_md_env.ps1 -TargetHost host01
#>
[CmdletBinding()]
param(
  [Parameter(Mandatory)][string]$TargetHost,
  [string]$User = $(if ($env:CLUSTER_USER) { $env:CLUSTER_USER } else { $env:USERNAME }),
  [string]$Domain = $env:CLUSTER_DOMAIN,
  [string]$Scratch = $env:CLUSTER_SCRATCH,             # blank = same auto-discovery as staging
  [int]$ConnectTimeout = 10,
  [int]$CommandTimeout = 1800                       # micromamba solve+download can take a while
)

if (-not (Get-Module -ListAvailable -Name Posh-SSH)) {
  Write-Error 'Posh-SSH not installed. Run:  Install-Module Posh-SSH -Scope CurrentUser'; return
}
Import-Module Posh-SSH -ErrorAction Stop

if (-not $Domain) { Write-Error 'Cluster config not set. Run:  . .\ops\local.ps1   (copy from ops/local.ps1.example)'; return }

$h = "$TargetHost.$Domain"
$cred = Get-Credential -UserName $User -Message "SSH password for $User@$h"
if (-not $cred) { Write-Error 'no credentials entered'; return }

$boot = @'
set -u
echo "host: $(hostname -f 2>/dev/null || hostname)"

# Re-discover the same scratch the staging step picked.
SCRATCH=""
for c in "__SCRATCH__" /data/work/$USER /data/$USER /data/scratch/$USER /data/tmp/$USER "$HOME/asyn-dwell"; do
  [ -z "$c" ] && continue
  if [ -d "$c" ] || mkdir -p "$c" 2>/dev/null; then
    if touch "$c/.w" 2>/dev/null; then rm -f "$c/.w"; SCRATCH="$c"; break; fi
  fi
done
[ -n "$SCRATCH" ] || { echo "  !! no writable scratch"; exit 3; }
REPO="$SCRATCH/asyn-oligomer-screen"
echo "  scratch: $SCRATCH"

echo "-- micromamba --"
MAMBA="$SCRATCH/bin/micromamba"
if [ -x "$MAMBA" ]; then
  echo "  present: $($MAMBA --version)"
else
  echo "  downloading micromamba"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$SCRATCH" bin/micromamba 2>&1 | sed 's/^/    /'
  [ -x "$MAMBA" ] && echo "  installed: $($MAMBA --version)" || { echo "  !! micromamba download failed"; exit 4; }
fi

echo "-- create asyn-md env (conda-forge openmm[CUDA] + OpenFF) --"
ENVDIR="$SCRATCH/mm/envs/asyn-md"
if [ -x "$ENVDIR/bin/python" ]; then
  echo "  env exists: $ENVDIR"
else
  "$MAMBA" create -y -p "$ENVDIR" -c conda-forge \
    python=3.11 openmm pdbfixer openff-toolkit openff-nagl openff-nagl-models openmmforcefields rdkit \
    2>&1 | tail -8 | sed 's/^/    /'
fi
[ -x "$ENVDIR/bin/python" ] || { echo "  !! env creation failed"; exit 5; }

echo "-- verify CUDA platform in asyn-md --"
"$ENVDIR/bin/python" - <<'PYEOF' 2>&1 | sed 's/^/    /'
import openmm
from openmm import Platform
print("openmm", openmm.version.version)
plats = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
print("platforms", plats)
fails = openmm.Platform.getPluginLoadFailures()
if fails:
    print("plugin load failures:")
    for f in fails: print("  ", f)
print("CUDA_OK" if "CUDA" in plats else "CUDA_MISSING")
PYEOF

echo "-- space --"
df -h "$SCRATCH" 2>/dev/null | awk 'NR==2{print "  fs free: "$4" of "$2" ("$5" used)"}'

echo "-- DONE --"
echo "  ASYN_MD_PYTHON=$ENVDIR/bin/python"
echo "  (export that before running the pilot so the orchestrator finds this env)"
'@

$boot = $boot.Replace('__SCRATCH__', $Scratch)

try {
  $s = New-SSHSession -ComputerName $h -Credential $cred -AcceptKey -ConnectionTimeout $ConnectTimeout -ErrorAction Stop
  try {
    $r = Invoke-SSHCommand -SessionId $s.SessionId -Command $boot -TimeOut $CommandTimeout
    if ($r.Output) { $r.Output }
    if ($r.Error)  { $r.Error | ForEach-Object { "  [stderr] $_" } }
  }
  finally { Remove-SSHSession -SessionId $s.SessionId | Out-Null }
}
catch { "  !! ${h}: $($_.Exception.Message)" }
