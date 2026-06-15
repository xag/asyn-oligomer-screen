#requires -Version 5.1
<#
.SYNOPSIS
  Pull the E8 swarm's trajectory frames (the small solute-only seg PDBs, not the
  9 MB state files) from px103/104/105 into results/msm_states/traj locally, so
  e8_msm.py can featurize + build the MSM. Tars seg_*.pdb on the remote and SCPs
  the tarball. Needs Posh-SSH + $env:PX_PASSWORD.
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
$repo = Split-Path -Parent $PSScriptRoot
$dest = Join-Path $repo 'results\msm_states'
$tdir = Join-Path $dest 'traj'
New-Item -ItemType Directory -Force -Path $tdir | Out-Null

foreach ($h in $Hosts) {
  Write-Host "==== collecting $h ===="
  try {
    $s = New-SSHSession -ComputerName "$h.$Domain" -Credential $cred -AcceptKey -ConnectionTimeout 12 -ErrorAction Stop
    try {
      # tar only finished replicas' frames (DONE present) to keep it clean
      $cmd = "cd $Work && rm -f e8_segs.tgz && find traj -name DONE | sed 's#/DONE##' | tar czf e8_segs.tgz -T - --transform 's#^#$h/#' 2>/dev/null; ls -la e8_segs.tgz | awk '{print `$5}'"
      $sz = (Invoke-SSHCommand -SessionId $s.SessionId -Command $cmd -TimeOut 300).Output
      Write-Host "  tarball bytes: $sz"
    } finally { Remove-SSHSession -SessionId $s.SessionId | Out-Null }
    Get-SCPItem -ComputerName "$h.$Domain" -Credential $cred -AcceptKey -PathType File -Path "$Work/e8_segs.tgz" -Destination $tdir -Force
    $tgz = Join-Path $tdir 'e8_segs.tgz'
    if (Test-Path $tgz) {
      & tar -xzf $tgz -C $tdir
      Remove-Item $tgz -Force
      Write-Host "  extracted under $tdir\$h"
    } else { Write-Host "  !! no tarball pulled from $h" }
  } catch { Write-Host "  !! $h : $($_.Exception.Message)" }
}
$n = (Get-ChildItem -Path $tdir -Recurse -Filter 'seg_*.pdb' -ErrorAction SilentlyContinue).Count
Write-Host "`ncollected $n seg PDBs into $tdir"
