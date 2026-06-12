#requires -Version 5.1
<#
.SYNOPSIS
  Windows/PowerShell twin of check_cluster.sh. Surveys CPU/GPU hardware and
  current utilization across <prefix>{1..8}.<domain> (set via env; see
  ops/local.ps1.example) over PASSWORD SSH, using the Posh-SSH module (no
  sshpass / WSL / PuTTY needed). Read-only: it runs only query commands on each
  host and writes a timestamped report locally.

.DESCRIPTION
  Why: sizing the dwell-time MD sweep (screen/dwell_time.py) needs to know how
  many GPUs each host has and how loaded they are right now — the run packs many
  small (~55k-atom) replicas per GPU, so GPU count, free VRAM, and the compute
  processes already running are the numbers that decide replicas-per-card.

  Setup (once, needs PSGallery/internet):
      Install-Module Posh-SSH -Scope CurrentUser

.EXAMPLE
  . .\ops\local.ps1
  .\check_cluster.ps1
  .\check_cluster.ps1 -Range 1,2,3 -Out report.txt
#>
[CmdletBinding()]
param(
  [string]$User = $(if ($env:CLUSTER_USER) { $env:CLUSTER_USER } else { $env:USERNAME }),
  [string]$Domain = $env:CLUSTER_DOMAIN,
  [string]$Prefix = $env:CLUSTER_PREFIX,
  [int[]]$Range = (1..8),
  [string]$Out = "cluster_survey_$(Get-Date -Format yyyyMMdd_HHmmss).txt",
  [int]$ConnectTimeout = 10,
  [int]$CommandTimeout = 60
)

if (-not (Get-Module -ListAvailable -Name Posh-SSH)) {
  Write-Error 'Posh-SSH not installed. Run:  Install-Module Posh-SSH -Scope CurrentUser'
  return
}
Import-Module Posh-SSH -ErrorAction Stop

if (-not $Domain -or -not $Prefix) { Write-Error 'Cluster config not set. Run:  . .\ops\local.ps1   (copy from ops/local.ps1.example)'; return }

# Prompt once; reuse the same credential for every host. The password lives only
# in the in-memory PSCredential (SecureString) — never on the command line.
$cred = Get-Credential -UserName $User `
  -Message "SSH password for $User@$Prefix[$($Range[0])-$($Range[-1])].$Domain"
if (-not $cred) { Write-Error 'no credentials entered'; return }

# Remote survey. Single-quoted here-string: PowerShell does NOT touch it, so every
# $(...) and quote is evaluated by the remote POSIX shell. Same checks as the .sh.
$survey = @'
echo "host: $(hostname -f 2>/dev/null || hostname)"
echo "-- os / kernel --"
( . /etc/os-release 2>/dev/null && echo "  ${PRETTY_NAME:-unknown}" ); echo "  kernel $(uname -r)"
echo "-- uptime / load --"
uptime | sed 's/^/  /'
echo "-- cpu --"
if command -v lscpu >/dev/null 2>&1; then
  lscpu | grep -E '^(Model name|Socket\(s\)|Core\(s\) per socket|Thread\(s\) per core|CPU\(s\)):' | sed 's/^/  /'
else
  echo "  $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //')"
  echo "  logical CPUs: $(nproc)"
fi
echo "-- memory --"
free -h | awk 'NR==1 || /^Mem:/ {print "  " $0}'
echo "-- busiest processes (top 5 by %CPU) --"
ps -eo pcpu,pmem,user,comm --sort=-pcpu | head -6 | sed 's/^/  /'
echo "-- disk --"
df -h 2>/dev/null | awk 'NR==1 || $6 ~ /^(\/|\/tmp|\/scratch|\/data|\/home)$/ {print "  " $0}'
echo "-- gpu --"
if command -v nvidia-smi >/dev/null 2>&1; then
  echo "  driver/CUDA: $(nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -1) / $(nvidia-smi | sed -n 's/.*CUDA Version: \([0-9.]*\).*/\1/p' | head -1)"
  echo "  idx, name, mem.total, mem.used, mem.free, util.gpu, util.mem, temp:"
  nvidia-smi --query-gpu=index,name,memory.total,memory.used,memory.free,utilization.gpu,utilization.memory,temperature.gpu --format=csv,noheader | sed 's/^/    /'
  echo "  compute processes now:"
  apps=$(nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name,used_memory --format=csv,noheader 2>/dev/null)
  if [ -n "$apps" ]; then echo "$apps" | sed 's/^/    /'; else echo "    (none -- GPUs idle)"; fi
else
  echo "  no nvidia-smi (no NVIDIA GPU, or driver not installed)"
  command -v lspci >/dev/null 2>&1 && lspci | grep -iE 'vga|3d|display' | sed 's/^/    pci: /'
fi
'@

$report = foreach ($x in $Range) {
  $h = "$Prefix$x.$Domain"
  '============================================================'
  "# $h"
  '============================================================'
  try {
    $s = New-SSHSession -ComputerName $h -Credential $cred -AcceptKey `
           -ConnectionTimeout $ConnectTimeout -ErrorAction Stop
    try {
      $r = Invoke-SSHCommand -SessionId $s.SessionId -Command $survey -TimeOut $CommandTimeout
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
  'cluster hardware / utilization survey'
  "date: $(Get-Date)"
  "user: $User   hosts: $Prefix[$($Range[0])..$($Range[-1])].$Domain"
  ''
) + $report | Tee-Object -FilePath $Out

Write-Host "report written to $Out"
