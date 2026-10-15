#!/usr/bin/env bash
set -euo pipefail
# Run only after acquiring the shared AMD window. This launcher never downloads,
# updates, or patches the copied environment or installed compiler packages.
: "${DEVICE_SOURCE_ROOT:?Set the immutable source root containing src/ and scripts/}"
: "${DEVICE_OUTPUT:?Set a fresh output directory}"
PROBE_ROOT=/othersys/tmp/intmat-probe-20261014
export JULIA_DEPOT_PATH="$PROBE_ROOT/depot:"
export JULIA_NUM_THREADS=2 JULIA_NUM_GC_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1
export ROCM_PATH="$PROBE_ROOT/rocm-debs/opt/rocm-6.4.3"
export LD_LIBRARY_PATH="$ROCM_PATH/lib" HSA_OVERRIDE_GFX_VERSION=10.3.0
export DEVICE_BACKEND=amd DEVICE_AMD_OPTLEVEL=1
: "${DEVICE_PROJECT:?Set the copied pinned environment matching DEVICE_SOURCE_ROOT}"
: "${DEVICE_SECONDS:=600}"
SWARM_MEM_MAX=8G /home/hbox/bin/swarm-build timeout "$DEVICE_SECONDS" \
  taskset -c 10,11 "$PROBE_ROOT/julia-1.12.4/bin/julia" \
  --startup-file=no --compiled-modules=yes --project="$DEVICE_PROJECT" \
  "$DEVICE_SOURCE_ROOT/scripts/benchmark_device_support.jl" "$DEVICE_OUTPUT"
