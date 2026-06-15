#requires -Version 5.1
<#
.SYNOPSIS
  Progress of the E8 MSM swarm across px103/104/105: replicas DONE / running /
  queued per host, runner alive, GPUs in use. Needs Posh-SSH + $env:PX_PASSWORD.
#>
[CmdletBinding()]
param(
  [int[]]$HostNums = @(4,3,5),
  [string]$Prefix = $env:CLUSTER_PREFIX,
  [string]$Domain = $env:CLUSTER_DOMAIN,
  [string]$User = $env:CLUSTER_USER,
  [string]$Base = $env:CLUSTER_SCRATCH
)
Import-Module Posh-SSH -ErrorAction Stop
if (-not $Domain -or -not $Prefix -or -not $Base) { Write-Error 'Cluster config not set. Run:  . .\ops\local.ps1'; return }
if (-not $User) { $User = $env:USERNAME }
$Hosts = @($HostNums | ForEach-Object { "$Prefix$_" })
$Work = "$Base/e8"
$sec = ConvertTo-SecureString $env:PX_PASSWORD -AsPlainText -Force
$cred = New-Object System.Management.Automation.PSCredential($User, $sec)
$probe = @"
done=`$(find $Work/traj -name DONE 2>/dev/null | wc -l)
started=`$(ls -d $Work/traj/*/ 2>/dev/null | wc -l)
total=`$(grep -c '"rep"' $Work/manifest.json 2>/dev/null)
alive=`$(pgrep -fc e8_swarm_runner.py)
md=`$(pgrep -fc md_relax.py)
echo "  runner_alive=`$alive  md_relax=`$md  replicas: `$done DONE / `$started started / `$total total"
tail -1 $Work/runner.log 2>/dev/null | sed 's/^/  last: /'
"@
$allDone = $true
foreach ($h in $Hosts) {
  "######## $h ########"
  try {
    $s = New-SSHSession -ComputerName "$h.$Domain" -Credential $cred -AcceptKey -ConnectionTimeout 12 -ErrorAction Stop
    try {
      $r = (Invoke-SSHCommand -SessionId $s.SessionId -Command $probe -TimeOut 40).Output
      $r
      $rt = ($r -join "`n")
      if ($rt -match 'runner_alive=([1-9])') { $allDone = $false }
      if ($rt -match 'replicas: (\d+) DONE / \d+ started / (\d+) total' -and $Matches[1] -ne $Matches[2]) { $allDone = $false }
    } finally { Remove-SSHSession -SessionId $s.SessionId | Out-Null }
  } catch { "  !! $h : $($_.Exception.Message)"; $allDone = $false }
}
if ($allDone) { "`nALL HOSTS COMPLETE" } else { "`nstill running" }
