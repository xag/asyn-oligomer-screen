#requires -Version 5.1
<#
.SYNOPSIS
  Probe the idle GPU nodes for OpenMM-MD readiness over password SSH (Posh-SSH).
  Read-only except for a single write-test under /data/work/$USER. Decides whether
  the dwell-time MD sweep can actually launch on these hosts tonight.

.DESCRIPTION
  Checks, per host: system python + whether `import openmm` works and which
  Platforms it exposes (CUDA needed), conda/mamba + env list (for the one-time
  OpenFF ligand parametrisation), which GPUs are idle *right now*, a writable
  scratch dir under /data, core tooling (git/uv/pip/nvcc), and outbound net to
  GitHub. Prints a READY / missing-bits verdict per host.

  Setup (once):  Install-Module Posh-SSH -Scope CurrentUser

.EXAMPLE
  . .\ops\local.ps1
  .\probe_md.ps1 -Hosts host01,host02
#>
[CmdletBinding()]
param(
  [string]$User = $(if ($env:CLUSTER_USER) { $env:CLUSTER_USER } else { $env:USERNAME }),
  [string]$Domain = $env:CLUSTER_DOMAIN,
  [string[]]$Hosts = @(),
  [string]$Out = "md_probe_$(Get-Date -Format yyyyMMdd_HHmmss).txt",
  [int]$ConnectTimeout = 10,
  [int]$CommandTimeout = 60
)

if (-not (Get-Module -ListAvailable -Name Posh-SSH)) {
  Write-Error 'Posh-SSH not installed. Run:  Install-Module Posh-SSH -Scope CurrentUser'
  return
}
Import-Module Posh-SSH -ErrorAction Stop

if (-not $Domain) { Write-Error 'Cluster config not set. Run:  . .\ops\local.ps1   (copy from ops/local.ps1.example)'; return }
if (-not $Hosts -or $Hosts.Count -eq 0) { Write-Error 'Pass -Hosts (e.g. -Hosts host01,host02).'; return }

$cred = Get-Credential -UserName $User `
  -Message "SSH password for $User@{$($Hosts -join ',')}.$Domain"
if (-not $cred) { Write-Error 'no credentials entered'; return }

# Remote readiness probe. Single-quoted here-string: evaluated by the remote shell.
$probe = @'
echo "host: $(hostname -f 2>/dev/null || hostname)"
OK_OMM=0; OK_CUDA=0; OK_SCRATCH=0

echo "-- python + openmm --"
PY=""
for cand in python3 python; do command -v $cand >/dev/null 2>&1 && { PY=$cand; break; }; done
if [ -n "$PY" ]; then
  echo "  $PY -> $($PY --version 2>&1) @ $(command -v $PY)"
  omm=$($PY - <<'PYEOF' 2>&1
try:
    import openmm
    from openmm import Platform
    plats = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    print("openmm", openmm.version.version)
    print("platforms", plats)
except Exception as e:
    print("openmm: NOT available (%r)" % (e,))
try:
    import openff.toolkit
    print("openff.toolkit", openff.toolkit.__version__)
except Exception as e:
    print("openff.toolkit: NOT available")
PYEOF
)
  echo "$omm" | sed 's/^/    /'
  echo "$omm" | grep -q "^openmm [0-9]" && OK_OMM=1
  echo "$omm" | grep -q "CUDA" && OK_CUDA=1
else
  echo "  no python found"
fi

echo "-- conda / mamba (for OpenFF parametrisation) --"
CFOUND=0
for c in conda mamba micromamba; do
  if command -v $c >/dev/null 2>&1; then
    echo "  $c -> $(command -v $c)"; CFOUND=1
    $c env list 2>/dev/null | sed 's/^/    /'
    break
  fi
done
for d in "$HOME/miniconda3" "$HOME/miniforge3" "$HOME/mambaforge" /opt/conda; do
  [ -d "$d" ] && { echo "  dir: $d"; CFOUND=1; }
done
[ "$CFOUND" = 0 ] && echo "  none (OpenFF/--prepare-only would need an install; apo+segment chunks still run pip-only)"

echo "-- idle GPUs right now (<1 GiB used) --"
nvidia-smi --query-gpu=index,name,memory.used,memory.free,utilization.gpu --format=csv,noheader 2>/dev/null \
  | awk -F',' '{u=$3+0; if(u<1024) print "    GPU"$1" FREE (used"$3", free"$4", util"$5")"}'
nfree=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader 2>/dev/null | awk '{u=$1+0; if(u<1024) n++} END{print n+0}')
echo "  idle count: $nfree"

echo "-- scratch under /data/work/$USER --"
SD="/data/work/$USER/asyn-dwell"
if mkdir -p "$SD" 2>/dev/null && touch "$SD/.probe_write_test" 2>/dev/null; then
  rm -f "$SD/.probe_write_test"; OK_SCRATCH=1
  echo "  writable: $SD"
else
  echo "  NOT writable: /data/work/$USER"
fi
df -h /data 2>/dev/null | awk 'NR==2{print "  /data free: "$4" of "$2" ("$5" used)"}'

echo "-- tooling --"
for t in git uv pip pip3 gcc nvcc; do command -v $t >/dev/null 2>&1 && echo "    $t -> $(command -v $t)"; done

echo "-- outbound net --"
if command -v curl >/dev/null 2>&1; then
  curl -m 5 -sSI https://github.com >/dev/null 2>&1 && echo "    github.com: reachable" || echo "    github.com: NOT reachable (limited egress)"
else
  echo "    curl absent"
fi

echo "-- READINESS --"
miss=""
[ "$OK_OMM" = 1 ]     || miss="$miss openmm"
[ "$OK_CUDA" = 1 ]    || miss="$miss cuda-platform"
[ "$OK_SCRATCH" = 1 ] || miss="$miss writable-scratch"
if [ -z "$miss" ]; then echo "  => READY to run MD chunks"; else echo "  => missing:$miss"; fi
'@

$report = foreach ($n in $Hosts) {
  $h = "$n.$Domain"
  '============================================================'
  "# $h"
  '============================================================'
  try {
    $s = New-SSHSession -ComputerName $h -Credential $cred -AcceptKey `
           -ConnectionTimeout $ConnectTimeout -ErrorAction Stop
    try {
      $r = Invoke-SSHCommand -SessionId $s.SessionId -Command $probe -TimeOut $CommandTimeout
      if ($r.Output) { $r.Output }
      if ($r.Error)  { $r.Error | ForEach-Object { "  [stderr] $_" } }
    }
    finally { Remove-SSHSession -SessionId $s.SessionId | Out-Null }
  }
  catch {
    "  !! ${h}: connection/auth failed -- $($_.Exception.Message)"
  }
  ''
}

@(
  'MD-readiness probe'
  "date: $(Get-Date)"
  "user: $User   hosts: $(( $Hosts | ForEach-Object { "$_.$Domain" } ) -join ', ')"
  ''
) + $report | Tee-Object -FilePath $Out

Write-Host "report written to $Out"
