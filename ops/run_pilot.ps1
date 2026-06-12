#requires -Version 5.1
<#
.SYNOPSIS
  Ship ops/asyn_launch.py to a GPU node and run it detached (nohup), or check
  status / stop it. The launcher pins replicas to idle GPUs and packs PER_GPU
  per card. Default is apo-only (no Vina/OpenFF); pass -Ligands for the gate.

.EXAMPLE
  . .\ops\local.ps1                                          # set CLUSTER_* (gitignored)
  .\run_pilot.ps1 -TargetHost host01                          # launch apo
  .\run_pilot.ps1 -TargetHost host01 -Action status
  .\run_pilot.ps1 -TargetHost host01 -Action stop
  .\run_pilot.ps1 -TargetHost host01 -Ligands "silibinin dhea trehalose caffeine"
#>
[CmdletBinding()]
param(
  [Parameter(Mandatory)][string]$TargetHost,
  [ValidateSet('launch','status','stop','diag')][string]$Action = 'launch',
  [string]$User = $(if ($env:CLUSTER_USER) { $env:CLUSTER_USER } else { $env:USERNAME }),
  [string]$Domain = $env:CLUSTER_DOMAIN,
  [string]$Scratch = $env:CLUSTER_SCRATCH,
  [string]$Shape = 'fusco_parallel_3mer_core70-88_relaxed',
  [string]$Gpus = '4 5 6 7',
  [int]$Replicas = 10,
  [int]$PerGpu = 2,
  [double]$ProdNs = 2.0,
  [string]$Ligands = '',
  [int]$ConnectTimeout = 10,
  [int]$CommandTimeout = 120
)

if (-not (Get-Module -ListAvailable -Name Posh-SSH)) {
  Write-Error 'Posh-SSH not installed. Run:  Install-Module Posh-SSH -Scope CurrentUser'; return
}
Import-Module Posh-SSH -ErrorAction Stop

if (-not $Domain)  { Write-Error 'Cluster config not set. Run:  . .\ops\local.ps1   (copy from ops/local.ps1.example)'; return }
if (-not $Scratch) { Write-Error 'CLUSTER_SCRATCH not set — see ops/local.ps1.'; return }

$h = "$TargetHost.$Domain"
$cred = Get-Credential -UserName $User -Message "SSH password for $User@$h"
if (-not $cred) { Write-Error 'no credentials entered'; return }

if ($Action -eq 'launch') {
  $pyLocal = Join-Path $PSScriptRoot 'asyn_launch.py'
  if (-not (Test-Path $pyLocal)) { Write-Error "missing $pyLocal"; return }
  $b64 = [Convert]::ToBase64String([IO.File]::ReadAllBytes($pyLocal))
  $cmd = @'
set -e
SCRATCH="__SCRATCH__"
REPO="$SCRATCH/asyn-oligomer-screen"
MDPY="$SCRATCH/mm/envs/asyn-md/bin/python"
echo '__B64__' | base64 -d > "$SCRATCH/asyn_launch.py"
echo "shipped asyn_launch.py ($(wc -c < "$SCRATCH/asyn_launch.py") bytes)"
[ -x "$MDPY" ] || { echo "!! asyn-md python missing at $MDPY"; exit 2; }
[ -x "$REPO/.venv/bin/python" ] || { echo "!! pip venv missing"; exit 2; }
# Docking (complex gate) needs the Vina binary — gitignored, so fetch once when ligands are requested.
if [ -n "__LIGANDS__" ] && [ ! -x "$REPO/.venv/bin/vina" ]; then
  echo "fetching AutoDock Vina 1.2.5 (linux) -> .venv/bin/vina"
  curl -fLs -o "$REPO/.venv/bin/vina" \
    https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 \
    && chmod +x "$REPO/.venv/bin/vina" \
    && "$REPO/.venv/bin/vina" --version 2>&1 | sed 's/^/  /' \
    || { echo "  !! vina fetch/exec failed"; exit 3; }
fi
# meeko's receptor-prep CLI name varies by build (mk_prepare_receptor or .py);
# stage3 checks $MK_PREPARE_RECEPTOR_BIN first, so point it at whatever installed.
MKR=$(ls "$REPO/.venv/bin/"mk_prepare_receptor* 2>/dev/null | head -1)
if [ -n "__LIGANDS__" ]; then
  if [ -n "$MKR" ]; then echo "receptor tool: $MKR"
  else echo "  !! no mk_prepare_receptor* in .venv/bin — run: -Action diag"; fi
fi
cd "$REPO"
nohup env ASYN_REPO="$REPO" ASYN_MD_PYTHON="$MDPY" MK_PREPARE_RECEPTOR_BIN="$MKR" \
  SHAPE="__SHAPE__" GPUS="__GPUS__" \
  REPLICAS=__REPLICAS__ PER_GPU=__PERGPU__ PROD_NS=__PRODNS__ LIGANDS="__LIGANDS__" \
  "$REPO/.venv/bin/python" "$SCRATCH/asyn_launch.py" < /dev/null > "$SCRATCH/launch.log" 2>&1 &
echo "LAUNCHER_PID $!"
sleep 4
echo "-- first launcher log lines --"
head -n 20 "$SCRATCH/launch.log" 2>/dev/null | sed 's/^/  /'
'@
  $cmd = $cmd.Replace('__B64__', $b64).Replace('__SCRATCH__', $Scratch).Replace('__SHAPE__', $Shape).
              Replace('__GPUS__', $Gpus).Replace('__REPLICAS__', "$Replicas").Replace('__PERGPU__', "$PerGpu").
              Replace('__PRODNS__', "$ProdNs").Replace('__LIGANDS__', $Ligands)
}
elseif ($Action -eq 'status') {
  $cmd = @'
SCRATCH="__SCRATCH__"; SHAPE="__SHAPE__"
SD="$SCRATCH/asyn-oligomer-screen/results/dwell/$SHAPE"
echo "-- launcher running? --"
pgrep -af "$SCRATCH/asyn_launch.py" | sed 's/^/  /' || echo "  (not running)"
echo "-- finished replicas (*_final.pdb) --"
n=$(ls -1 "$SD"/*_final.pdb 2>/dev/null | wc -l); echo "  count: $n"
ls -1 "$SD"/*_final.pdb 2>/dev/null | sed 's#.*/##;s/^/    /'
echo "-- GPU 4-7 right now --"
nvidia-smi --query-gpu=index,utilization.gpu,memory.used,memory.free --format=csv,noheader | awk -F',' '($1+0)>=4{print "  GPU"$0}'
echo "-- launcher log tail --"
tail -n 22 "$SCRATCH/launch.log" 2>/dev/null | sed 's/^/  /'
'@
  $cmd = $cmd.Replace('__SCRATCH__', $Scratch).Replace('__SHAPE__', $Shape)
}
elseif ($Action -eq 'diag') {
  $cmd = @'
SCRATCH="__SCRATCH__"
REPO="$SCRATCH/asyn-oligomer-screen"
echo "-- meeko CLIs in .venv/bin --"
ls -1 "$REPO/.venv/bin/" | grep -i 'mk_' | sed 's/^/  /' || echo "  (none)"
echo "-- meeko version --"
"$REPO/.venv/bin/python" -c "import meeko, sys; print(meeko.__version__)" 2>&1 | sed 's/^/  /'
echo "-- vina --"
{ [ -x "$REPO/.venv/bin/vina" ] && "$REPO/.venv/bin/vina" --version 2>&1; } | sed 's/^/  /' || echo "  vina not present"
echo "-- mk_prepare_receptor* --"
MKR=$(ls "$REPO/.venv/bin/"mk_prepare_receptor* 2>/dev/null | head -1)
if [ -n "$MKR" ]; then echo "  found: $MKR"; "$MKR" --help 2>&1 | head -6 | sed 's/^/    /'; else echo "  none found in .venv/bin"; fi
'@
  $cmd = $cmd.Replace('__SCRATCH__', $Scratch)
}
else {  # stop
  $cmd = @'
SCRATCH="__SCRATCH__"
pkill -f "$SCRATCH/asyn_launch.py" && echo "killed launcher" || echo "no launcher process"
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
