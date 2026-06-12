#requires -Version 5.1
<#
.SYNOPSIS
  Low-cost blocked dwell-time validation. `run` launches one block (one
  conformer; matched apo/test/decoy seeds) detached on idle GPUs — sized for a
  short opportunistic session; `pool` re-reads every accumulated block and prints
  the pre-registered sequential decision. `status`/`stop` manage a running block.

.DESCRIPTION
  Each night: pick a conformer and `run` (skips replicas already on disk, so
  re-running tops a block up); change -Conformer across nights to add blocks.
  Each morning: `pool`. Blocks accumulate under results/blocks/<conformer>/.

.EXAMPLE
  . .\ops\local.ps1
  .\block.ps1 -TargetHost host01 -Conformer fusco_parallel_3mer_core70-88_relaxed
  .\block.ps1 -TargetHost host01 -Action status
  .\block.ps1 -TargetHost host01 -Action pool
#>
[CmdletBinding()]
param(
  [Parameter(Mandatory)][string]$TargetHost,
  [ValidateSet('run','pool','status','stop')][string]$Action = 'run',
  [string]$User = $(if ($env:CLUSTER_USER) { $env:CLUSTER_USER } else { $env:USERNAME }),
  [string]$Domain = $env:CLUSTER_DOMAIN,
  [string]$Scratch = $env:CLUSTER_SCRATCH,
  [string]$Conformer = 'fusco_parallel_3mer_core70-88_relaxed',
  [string]$Ligands = 'silibinin caffeine',
  [int]$Replicas = 5,
  [int]$SeedBase = 7000,
  [string]$Gpus = '4 5 6 7',
  [int]$PerGpu = 2,
  [double]$ProdNs = 2.0,
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

function B64([string]$name) {
  $p = Join-Path $PSScriptRoot $name
  if (-not (Test-Path $p)) { throw "missing $p" }
  [Convert]::ToBase64String([IO.File]::ReadAllBytes($p))
}

if ($Action -eq 'run') {
  $b64Block = B64 'asyn_block.py'
  $cmd = @'
set -e
SCRATCH="__SCRATCH__"
REPO="$SCRATCH/asyn-oligomer-screen"
MDPY="$SCRATCH/mm/envs/asyn-md/bin/python"
echo '__B64BLOCK__' | base64 -d > "$SCRATCH/asyn_block.py"
[ -x "$MDPY" ] || { echo "!! asyn-md python missing"; exit 2; }
[ -x "$REPO/.venv/bin/python" ] || { echo "!! pip venv missing"; exit 2; }
if [ ! -x "$REPO/.venv/bin/vina" ]; then
  curl -fLs -o "$REPO/.venv/bin/vina" \
    https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 \
    && chmod +x "$REPO/.venv/bin/vina" || { echo "  !! vina fetch failed"; exit 3; }
fi
MKR=$(ls "$REPO/.venv/bin/"mk_prepare_receptor* 2>/dev/null | head -1)
cd "$REPO"
nohup env ASYN_REPO="$REPO" ASYN_MD_PYTHON="$MDPY" MK_PREPARE_RECEPTOR_BIN="$MKR" \
  CONFORMER="__CONFORMER__" LIGANDS="__LIGANDS__" REPLICAS=__REPLICAS__ SEED_BASE=__SEEDBASE__ \
  GPUS="__GPUS__" PER_GPU=__PERGPU__ PROD_NS=__PRODNS__ \
  "$REPO/.venv/bin/python" "$SCRATCH/asyn_block.py" < /dev/null > "$SCRATCH/block.log" 2>&1 &
echo "BLOCK_PID $!"
sleep 4
echo "-- first log lines --"; head -n 18 "$SCRATCH/block.log" 2>/dev/null | sed 's/^/  /'
'@
  $cmd = $cmd.Replace('__B64BLOCK__', $b64Block).Replace('__SCRATCH__', $Scratch).
              Replace('__CONFORMER__', $Conformer).Replace('__LIGANDS__', $Ligands).
              Replace('__REPLICAS__', "$Replicas").Replace('__SEEDBASE__', "$SeedBase").
              Replace('__GPUS__', $Gpus).Replace('__PERGPU__', "$PerGpu").Replace('__PRODNS__', "$ProdNs")
}
elseif ($Action -eq 'pool') {
  $b64Pool = B64 'asyn_pool.py'
  $b64Pre = B64 'prereg.py'
  $cmd = @'
set -e
SCRATCH="__SCRATCH__"
REPO="$SCRATCH/asyn-oligomer-screen"
echo '__B64POOL__' | base64 -d > "$SCRATCH/asyn_pool.py"
echo '__B64PRE__'  | base64 -d > "$SCRATCH/prereg.py"
cd "$REPO"
env ASYN_REPO="$REPO" "$REPO/.venv/bin/python" "$SCRATCH/asyn_pool.py"
'@
  $cmd = $cmd.Replace('__B64POOL__', $b64Pool).Replace('__B64PRE__', $b64Pre).Replace('__SCRATCH__', $Scratch)
}
elseif ($Action -eq 'status') {
  $cmd = @'
SCRATCH="__SCRATCH__"
echo "-- block running? --"
pgrep -af "$SCRATCH/asyn_block.py" | sed 's/^/  /' || echo "  (not running)"
echo "-- finished replicas per block --"
for d in "$SCRATCH/asyn-oligomer-screen/results/blocks"/*/; do
  [ -d "$d" ] || continue
  n=$(ls -1 "$d"*_final.pdb 2>/dev/null | wc -l)
  echo "  $(basename "$d"): $n *_final.pdb"
done
echo "-- GPU 4-7 --"
nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader | awk -F',' '($1+0)>=4{print "  GPU"$0}'
echo "-- block.log tail --"
tail -n 18 "$SCRATCH/block.log" 2>/dev/null | sed 's/^/  /'
'@
  $cmd = $cmd.Replace('__SCRATCH__', $Scratch)
}
else {  # stop
  $cmd = @'
SCRATCH="__SCRATCH__"
pkill -f "$SCRATCH/asyn_block.py" && echo "killed block launcher" || echo "no block launcher"
pkill -f "$SCRATCH/asyn-oligomer-screen/screen/md_relax.py" && echo "killed md_relax workers" || echo "no md_relax workers"
'@
  $cmd = $cmd.Replace('__SCRATCH__', $Scratch)
}

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
