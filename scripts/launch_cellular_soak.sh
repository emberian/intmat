#!/usr/bin/env bash
# Run on hbox after the source and pinned Julia dependencies are provisioned.
# The source tree and any referenced checkpoint must remain immutable.
set -euo pipefail
if (( $# < 4 || $# > 6 )); then
  echo 'usage: launch_cellular_soak.sh SOURCE PROFILE OUTPUT CPU_PAIR [SECONDS=43200] [CHECKPOINT]' >&2
  exit 2
fi
source_root=$(realpath "$1")
profile=$(realpath "$2")
output=$(realpath -m "$3")
cpus=$4
seconds=${5:-43200}
[[ $cpus =~ ^[0-9]+,[0-9]+$ && $seconds =~ ^[0-9]+$ ]] || { echo 'CPU_PAIR needs two CPU indices and seconds must be integral' >&2; exit 2; }
: "${SOAK_DEPOT:?Set SOAK_DEPOT to the isolated writable Julia depot}"
julia=${SOAK_JULIA:-/othersys/tmp/intmat-probe-20261014/julia-1.12.4/bin/julia}
export JULIA_DEPOT_PATH="$SOAK_DEPOT:/othersys/tmp/intmat-probe-20261014/depot:"
export JULIA_NUM_THREADS=2 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_PKG_PRECOMPILE_AUTO=0 OPENBLAS_NUM_THREADS=1
mkdir -p "$(dirname "$output")"
stamp=$(date -u +%Y%m%dT%H%M%SZ)-$$
log="$output.launch-$stamp.log"
args=("$julia" --startup-file=no --compiled-modules=yes --threads=2 --gcthreads=1 --project="$source_root"
      "$source_root/scripts/cellular_soak.jl" "$profile" "$output" --seconds "$seconds"
      --budget-bytes "${SOAK_ARTIFACT_BUDGET:-4000000000}" --metrics-seconds 20 --archive-seconds 600)
if (( $# == 6 )); then args+=(--resume "$(realpath "$6")"); fi
[[ ${SOAK_THROUGHPUT:-0} =~ ^[01]$ && ${SOAK_ALLOW_SOURCE_CHANGE:-0} =~ ^[01]$ ]] || { echo 'Soak mode environment values must be 0 or 1' >&2; exit 2; }
if [[ ${SOAK_THROUGHPUT:-0} == 1 ]]; then args+=(--throughput); fi
if [[ ${SOAK_ALLOW_SOURCE_CHANGE:-0} == 1 ]]; then
  (( $# == 6 )) || { echo 'Source migration requires an explicit checkpoint' >&2; exit 2; }
  args+=(--allow-source-change)
fi
# timeout is a secondary guard; the Julia driver exits and saves at its own wall
# deadline. INT gives its finally block a chance to finish before the kill grace.
nohup env SWARM_MEM_MAX="${SOAK_MEM_MAX:-8G}" /home/hbox/bin/swarm-build \
  taskset -c "$cpus" timeout --signal=INT --kill-after=60s "$((seconds+900))" \
  "${args[@]}" </dev/null >"$log" 2>&1 &
pid=$!
printf 'launcher_pid=%s\nlog=%s\noutput=%s\ncpus=%s\nseconds=%s\n' "$pid" "$log" "$output" "$cpus" "$seconds" | tee "$output.launch-$stamp.txt"
