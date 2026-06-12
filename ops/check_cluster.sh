#!/usr/bin/env bash
#
# check_cluster.sh — survey CPU/GPU hardware and current utilization across
# ${CLUSTER_PREFIX}{LO..HI}.${CLUSTER_DOMAIN} over password SSH. Read-only: it runs only
# query commands on each host and writes a timestamped report locally.
#
# Why: sizing the dwell-time MD sweep (screen/dwell_time.py) needs to know how
# many GPUs each host actually has and how loaded they are *right now* — the run
# packs many small (~55k-atom) replicas per GPU, so GPU count, free VRAM, and
# the compute processes already running are the numbers that decide how many
# replicas can go in flight.
#
# Requires: sshpass  (apt install sshpass | dnf install sshpass |
#                      brew install hudochenkov/sshpass/sshpass)
#           Without it, OpenSSH cannot take a password non-interactively; the
#           script falls back to plain ssh, which prompts once per host.
#
# Config via env (set your cluster here, e.g. export CLUSTER_DOMAIN=... CLUSTER_PREFIX=...):
#   CLUSTER_DOMAIN  (required)  e.g. cluster.example.com
#   CLUSTER_PREFIX  (required)  host-name prefix, e.g. node  -> node{1..8}
#   CLUSTER_USER    (optional)  defaults to $USER
#
# Usage:  CLUSTER_DOMAIN=cluster.example.com CLUSTER_PREFIX=node ./check_cluster.sh [-u user] [-o report.txt] [-r LO-HI]
#
set -uo pipefail

USER_NAME="${SSH_USER:-${CLUSTER_USER:-$USER}}"
DOMAIN="${CLUSTER_DOMAIN:?set CLUSTER_DOMAIN (e.g. cluster.example.com)}"
PREFIX="${CLUSTER_PREFIX:?set CLUSTER_PREFIX (host-name prefix, e.g. node)}"
OUT=""
LO=1
HI=8

while getopts "u:o:r:h" opt; do
  case "$opt" in
    u) USER_NAME="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    r) LO="${OPTARG%%-*}"; HI="${OPTARG##*-}" ;;
    h) sed -n '2,/^set /{/^set /d;s/^# \{0,1\}//p}' "$0"; exit 0 ;;
    *) echo "usage: $0 [-u user] [-o report.txt] [-r LO-HI]" >&2; exit 2 ;;
  esac
done

# Prompt once and reuse for every host. sshpass -e reads the password from the
# SSHPASS env var so it never shows up in 'ps' output (unlike sshpass -p).
read -r -s -p "SSH password for ${USER_NAME}@${PREFIX}[${LO}-${HI}].${DOMAIN}: " SSHPASS
echo
export SSHPASS

HAVE_SSHPASS=1
command -v sshpass >/dev/null 2>&1 || {
  HAVE_SSHPASS=0
  echo "warning: sshpass not found — falling back to interactive ssh (one password prompt per host)." >&2
}

[ -n "$OUT" ] || OUT="cluster_survey_$(date +%Y%m%d_%H%M%S).txt"

# Force password auth (the hosts only accept it), keep first-connect non-blocking
# with accept-new, and time out fast on dead hosts so one outage can't stall the
# whole sweep.
SSH_OPTS=(
  -o ConnectTimeout=10
  -o StrictHostKeyChecking=accept-new
  -o PreferredAuthentications=password
  -o PubkeyAuthentication=no
  -o NumberOfPasswordPrompts=1
)

# ---- remote survey: emitted to stdout, piped into `bash -s` on each host. ----
# Single-quoted heredoc so every $ and $(...) is evaluated ON THE REMOTE HOST.
remote_survey() {
cat <<'REMOTE'
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
  nvidia-smi --query-gpu=index,name,memory.total,memory.used,memory.free,utilization.gpu,utilization.memory,temperature.gpu \
             --format=csv,noheader | sed 's/^/    /'
  echo "  compute processes now:"
  apps=$(nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name,used_memory --format=csv,noheader 2>/dev/null)
  if [ -n "$apps" ]; then echo "$apps" | sed 's/^/    /'; else echo "    (none — GPUs idle)"; fi
else
  echo "  no nvidia-smi (no NVIDIA GPU, or driver not installed)"
  command -v lspci >/dev/null 2>&1 && lspci | grep -iE 'vga|3d|display' | sed 's/^/    pci: /'
fi
echo
REMOTE
}

run_one() {
  local host="$1"
  if [ "$HAVE_SSHPASS" -eq 1 ]; then
    remote_survey | sshpass -e ssh "${SSH_OPTS[@]}" "${USER_NAME}@${host}" 'bash -s'
  else
    remote_survey | ssh "${SSH_OPTS[@]}" "${USER_NAME}@${host}" 'bash -s'
  fi
}

{
  echo "cluster hardware / utilization survey"
  echo "date: $(date)"
  echo "user: ${USER_NAME}   hosts: ${PREFIX}{${LO}..${HI}}.${DOMAIN}"
  echo
  for x in $(seq "$LO" "$HI"); do
    host="${PREFIX}${x}.${DOMAIN}"
    echo "============================================================"
    echo "# ${host}"
    echo "============================================================"
    run_one "$host"
    rc=$?
    [ "$rc" -eq 0 ] || echo "  !! ${host}: ssh failed (exit ${rc}) — unreachable, auth rejected, or timed out"
    echo
  done
} 2>&1 | tee "$OUT"

echo "report written to ${OUT}" >&2
