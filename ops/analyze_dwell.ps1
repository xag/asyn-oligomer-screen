#requires -Version 5.1
<#
.SYNOPSIS
  Ship ops/asyn_diag.py to the node and run it (CPU-only, pip venv) to diagnose
  why the dwell-time gate came back null — occupancy, basin-threshold fit, and
  apo-vs-control time-courses on the trajectories already on disk. No GPU/MD.

.EXAMPLE
  . .\ops\local.ps1
  .\analyze_dwell.ps1 -TargetHost host01
#>
[CmdletBinding()]
param(
  [Parameter(Mandatory)][string]$TargetHost,
  [string]$User = $(if ($env:CLUSTER_USER) { $env:CLUSTER_USER } else { $env:USERNAME }),
  [string]$Domain = $env:CLUSTER_DOMAIN,
  [string]$Scratch = $env:CLUSTER_SCRATCH,
  [string]$Shape = 'fusco_parallel_3mer_core70-88_relaxed',
  [string]$Ligands = 'silibinin dhea caffeine',
  [switch]$Rescore = $true,
  [int]$ConnectTimeout = 10,
  [int]$CommandTimeout = 600
)

if (-not (Get-Module -ListAvailable -Name Posh-SSH)) {
  Write-Error 'Posh-SSH not installed. Run:  Install-Module Posh-SSH -Scope CurrentUser'; return
}
Import-Module Posh-SSH -ErrorAction Stop
if (-not $Domain)  { Write-Error 'Cluster config not set. Run:  . .\ops\local.ps1'; return }
if (-not $Scratch) { Write-Error 'CLUSTER_SCRATCH not set — see ops/local.ps1.'; return }

$h = "$TargetHost.$Domain"
$cred = Get-Credential -UserName $User -Message "SSH password for $User@$h"
if (-not $cred) { Write-Error 'no credentials entered'; return }

$pyLocal = Join-Path $PSScriptRoot 'asyn_diag.py'
if (-not (Test-Path $pyLocal)) { Write-Error "missing $pyLocal"; return }
$b64 = [Convert]::ToBase64String([IO.File]::ReadAllBytes($pyLocal))

$rescoreFlag = if ($Rescore) { ' --rescore' } else { '' }
$cmd = @'
set -e
SCRATCH="__SCRATCH__"
REPO="$SCRATCH/asyn-oligomer-screen"
echo '__B64__' | base64 -d > "$SCRATCH/asyn_diag.py"
cd "$REPO"
env ASYN_REPO="$REPO" SHAPE="__SHAPE__" LIGANDS="__LIGANDS__" \
  "$REPO/.venv/bin/python" "$SCRATCH/asyn_diag.py"__RESCORE__
'@
$cmd = $cmd.Replace('__B64__', $b64).Replace('__SCRATCH__', $Scratch).
            Replace('__SHAPE__', $Shape).Replace('__LIGANDS__', $Ligands).
            Replace('__RESCORE__', $rescoreFlag)

try {
  $s = New-SSHSession -ComputerName $h -Credential $cred -AcceptKey -ConnectionTimeout $ConnectTimeout -ErrorAction Stop
  try {
    $r = Invoke-SSHCommand -SessionId $s.SessionId -Command $cmd -TimeOut $CommandTimeout
    if ($r.Output) { $r.Output }
    if ($r.Error)  { $r.Error | ForEach-Object { "  [stderr] $_" } }
  }
  finally { Remove-SSHSession -SessionId $s.SessionId | Out-Null }
}
catch { "  !! ${h}: $($_.Exception.Message)" }
