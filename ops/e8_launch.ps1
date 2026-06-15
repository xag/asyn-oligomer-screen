#requires -Version 5.1
<#
.SYNOPSIS
  Launch the E8 MSM occupancy swarm: many short, velocity-seeded MD replicas
  from a library of distinct 3-mer basin seeds (both registers) across the
  live-idle half of the GPUs on px103/104/105. Writes only under
  /data/infint/xgrehant/e8. Needs Posh-SSH and $env:PX_PASSWORD.

  Idempotent: the node runner skips replicas that already carry a DONE marker,
  so a re-run resumes rather than restarts.

.PARAMETER Reps        velocity replicas per seed conformation (default 10)
.PARAMETER EquilPs     equilibration ps per replica (default 200)
.PARAMETER SegmentPs   production ps per segment (default 5000 = 5 ns)
.PARAMETER NSegments   segments per replica (default 3 -> 15 ns production)
.PARAMETER MaxGpus     cap on selected GPUs per host (default 4; politeness)
.PARAMETER SmokeName   if set, write a tiny 1-seed/1-rep manifest to ONE host
                       (px104) with short times and launch it foreground-style
                       for verification, then stop.
#>
[CmdletBinding()]
param(
  [int]$Reps = 10,
  [double]$EquilPs = 200,
  [double]$SegmentPs = 5000,
  [int]$NSegments = 3,
  [double]$TrajIntervalPs = 50,
  [int]$MaxGpus = 4,
  [int[]]$HostNums = @(4,3,5),                 # <prefix><num>; half-idle GPUs taken per host
  [string]$Prefix = $env:CLUSTER_PREFIX,
  [string]$Domain = $env:CLUSTER_DOMAIN,
  [string]$User = $env:CLUSTER_USER,
  [string]$Base = $env:CLUSTER_SCRATCH,
  [switch]$Smoke
)

Import-Module Posh-SSH -ErrorAction Stop
if (-not $env:PX_PASSWORD) { Write-Error 'set $env:PX_PASSWORD'; return }
if (-not $Domain -or -not $Prefix -or -not $Base) { Write-Error 'Cluster config not set. Run:  . .\ops\local.ps1   (copy from ops/local.ps1.example)'; return }
if (-not $User) { $User = $env:USERNAME }
$Hosts = @($HostNums | ForEach-Object { "$Prefix$_" })
$sec = ConvertTo-SecureString $env:PX_PASSWORD -AsPlainText -Force
$cred = New-Object System.Management.Automation.PSCredential($User, $sec)
$repo = Split-Path -Parent $PSScriptRoot
$work = "$Base/e8"
$mdpy = "$Base/miniforge3/envs/md/bin/python"

# --- seed library: distinct relaxed-basin conformations (core58-102 chunks) ---
# name carries the register (par_/anti_) for per-register MSM grouping later.
$seeds = @(
  @{ n='par_c7088';      p="$repo\results\funnel_states\inputs\fusco_parallel_3mer_core70-88_core58-102.pdb" }
  @{ n='par_c7088_s123'; p="$repo\results\funnel_states\inputs\fusco_parallel_3mer_core70-88_s123_core58-102.pdb" }
  @{ n='par_c7088_s777'; p="$repo\results\funnel_states\inputs\fusco_parallel_3mer_core70-88_s777_core58-102.pdb" }
  @{ n='par_c6583';      p="$repo\results\funnel_states\inputs\fusco_parallel_3mer_core65-83_core58-102.pdb" }
  @{ n='par_c7391';      p="$repo\results\funnel_states\inputs\fusco_parallel_3mer_core73-91_core58-102.pdb" }
  @{ n='par_basin1';     p="$repo\results\stable_states\analysis\fusco_parallel_3mer_core70-88_basin1_medoid.pdb" }
  @{ n='anti_c7088';      p="$repo\results\funnel_states\inputs\fusco_antiparallel_3mer_core70-88_core58-102.pdb" }
  @{ n='anti_c7088_s123'; p="$repo\results\funnel_states\inputs\fusco_antiparallel_3mer_core70-88_s123_core58-102.pdb" }
  @{ n='anti_basin1';     p="$repo\results\stable_states\analysis\fusco_antiparallel_3mer_core70-88_basin1_medoid.pdb" }
  @{ n='anti_basin2';     p="$repo\results\stable_states\analysis\fusco_antiparallel_3mer_core70-88_basin2_medoid.pdb" }
)
foreach ($s in $seeds) { if (-not (Test-Path $s.p)) { Write-Error "missing seed $($s.n): $($s.p)"; return } }

function New-Session($h) {
  New-SSHSession -ComputerName "$h.$Domain" -Credential $cred -AcceptKey -ConnectionTimeout 12 -ErrorAction Stop
}
function Invoke-Remote($sid, $cmd, $timeout = 120) {
  (Invoke-SSHCommand -SessionId $sid -Command $cmd -TimeOut $timeout)
}
function Get-IdleGpus($sid) {
  # idle := <2 GiB used and <10% util right now
  $r = Invoke-Remote $sid 'nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader,nounits'
  $idle = @()
  foreach ($line in ($r.Output -split "`n")) {
    if ($line -match '^\s*(\d+),\s*(\d+),\s*(\d+)') {
      if ([int]$Matches[2] -lt 2000 -and [int]$Matches[3] -lt 10) { $idle += [int]$Matches[1] }
    }
  }
  return ,$idle
}

# ---------------- SMOKE TEST: tiny 1-rep job on px104, verify end-to-end ------
if ($Smoke) {
  $h = $Hosts[0]; $s = New-Session $h; $sid = $s.SessionId
  try {
    $idle = Get-IdleGpus $sid
    if (-not $idle) { Write-Error "no idle GPU on $h for smoke test"; return }
    $g = $idle[0]
    Write-Host "[smoke] $h GPU $g; staging tiny job..."
    Invoke-Remote $sid "mkdir -p $work/inputs $work/logs $work/traj $work/prep; cp -n $Base/e7/md_relax.py $work/md_relax.py" | Out-Null
    Set-SCPItem -ComputerName "$h.$Domain" -Credential $cred -AcceptKey -Path "$PSScriptRoot\e8_swarm_runner.py" -Destination $work -Force
    $sd = $seeds[0]
    Set-SCPItem -ComputerName "$h.$Domain" -Credential $cred -AcceptKey -Path $sd.p -Destination "$work/inputs" -Force
    $remotePdb = "inputs/" + (Split-Path $sd.p -Leaf)
    $man = @{ work=$work; py=$mdpy; md_relax="$work/md_relax.py"; gpus=@($g); per_gpu=1;
              equil_ps=20.0; segment_ps=100.0; n_segments=1; traj_interval_ps=50.0;
              jobs=@(@{ seed=$sd.n; pdb=$remotePdb; rep=0 }) } | ConvertTo-Json -Depth 5
    $manB64 = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($man))
    Invoke-Remote $sid "echo $manB64 | base64 -d > $work/manifest.json" | Out-Null
    Write-Host "[smoke] running runner foreground (short)..."
    $r = Invoke-Remote $sid "cd $work && timeout 1200 $mdpy e8_swarm_runner.py 2>&1 | tail -25" 1300
    $r.Output
    $chk = Invoke-Remote $sid "ls -la $work/traj/$($sd.n)_rep00/ 2>&1; echo '--- platform line ---'; grep -i 'platform\|CUDA' $work/logs/$($sd.n)_rep00.log 2>/dev/null | head -3"
    $chk.Output
  } finally { Remove-SSHSession -SessionId $sid | Out-Null }
  return
}

# ---------------- FULL SWARM ------------------------------------------------
# Probe hosts, select half of live-idle GPUs (capped), assign whole seeds to
# hosts weighted by GPU count (so each host prepares only its seeds).
$plan = @()
foreach ($h in $Hosts) {
  $s = New-Session $h; $sid = $s.SessionId
  try {
    $idle = Get-IdleGpus $sid
    $take = [Math]::Min([Math]::Max([int][Math]::Floor($idle.Count / 2), 1), $MaxGpus)
    $sel = @($idle | Select-Object -First $take)
    $plan += [pscustomobject]@{ host=$h; sid=$sid; session=$s; gpus=$sel; seeds=@() }
    Write-Host "$h : idle=[$($idle -join ',')] -> using [$($sel -join ',')]"
  } catch { Write-Host "$h : probe failed - $($_.Exception.Message)"; if ($s) { Remove-SSHSession -SessionId $s.SessionId | Out-Null } }
}
if (-not $plan) { Write-Error 'no usable hosts'; return }

# weighted round-robin seed assignment, preserving register interleave
$totalGpu = ($plan | Measure-Object -Property { $_.gpus.Count } -Sum).Sum
$idx = 0
foreach ($sd in $seeds) {
  # pick host with the highest (gpu_share - assigned_share) deficit
  $best = $plan | Sort-Object @{ E = { ($_.seeds.Count + 1) / [double]$_.gpus.Count } } | Select-Object -First 1
  $best.seeds += $sd
  $idx++
}

$launched = @()
foreach ($p in $plan) {
  if (-not $p.seeds) { Write-Host "$($p.host): no seeds assigned, skipping"; continue }
  $sid = $p.sid; $h = $p.host
  Write-Host "==== staging $h : gpus[$($p.gpus -join ',')] seeds[$(( $p.seeds | ForEach-Object { $_.n }) -join ',')] ===="
  Invoke-Remote $sid "mkdir -p $work/inputs $work/logs $work/traj $work/prep $work/tmp $work/nvcache; cp -n $Base/e7/md_relax.py $work/md_relax.py" | Out-Null
  Set-SCPItem -ComputerName "$h.$Domain" -Credential $cred -AcceptKey -Path "$PSScriptRoot\e8_swarm_runner.py" -Destination $work -Force
  $jobs = @()
  foreach ($sd in $p.seeds) {
    Set-SCPItem -ComputerName "$h.$Domain" -Credential $cred -AcceptKey -Path $sd.p -Destination "$work/inputs" -Force
    $remotePdb = "inputs/" + (Split-Path $sd.p -Leaf)
    for ($r = 0; $r -lt $Reps; $r++) { $jobs += @{ seed=$sd.n; pdb=$remotePdb; rep=$r } }
  }
  $man = @{ work=$work; py=$mdpy; md_relax="$work/md_relax.py"; gpus=@($p.gpus); per_gpu=1;
            equil_ps=$EquilPs; segment_ps=$SegmentPs; n_segments=$NSegments;
            traj_interval_ps=$TrajIntervalPs; prep_workers=4; jobs=$jobs } | ConvertTo-Json -Depth 5
  $manB64 = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($man))
  Invoke-Remote $sid "echo $manB64 | base64 -d > $work/manifest.json" | Out-Null
  # launch detached; survives the SSH session
  $launch = "cd $work && nohup $mdpy e8_swarm_runner.py > runner.log 2>&1 & echo PID=`$!"
  $r = Invoke-Remote $sid $launch
  $pidLine = ($r.Output -split "`n" | Where-Object { $_ -match 'PID=' }) -join ''
  Write-Host "  launched: $pidLine ($($jobs.Count) replicas across $($p.gpus.Count) GPUs)"
  $launched += [pscustomobject]@{ host=$h; gpus=$p.gpus; replicas=$jobs.Count; pid=$pidLine }
  Remove-SSHSession -SessionId $sid | Out-Null
}

Write-Host ""
Write-Host "=== E8 swarm launched ==="
$launched | Format-Table -AutoSize
Write-Host "monitor:  .\ops\e8_status.ps1     (or: ssh <host> 'tail -f $work/runner.log')"
