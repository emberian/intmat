#!/usr/bin/env bash
# Run on hbox after provisioning an isolated Julia + ROCm directory.
# All Julia writes stay in PROBE_ROOT; source is a copied Interior checkout.
set -euo pipefail
: "${PROBE_ROOT:?Set the isolated probe directory}"
export JULIA_DEPOT_PATH="$PROBE_ROOT/depot"
export JULIA_NUM_THREADS=2
export JULIA_NUM_PRECOMPILE_TASKS=1
export JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1
export ROCM_PATH="$PROBE_ROOT/rocm-debs/opt/rocm-6.4.3"
export LD_LIBRARY_PATH="$ROCM_PATH/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PATH="$ROCM_PATH/bin:$PATH"
julia="$PROBE_ROOT/julia-1.11.7/bin/julia"
run() {
    SWARM_MEM_MAX=8G /home/hbox/bin/swarm-build taskset -c 0,1 \
        timeout --signal=TERM --kill-after=15s 900 "$julia" \
        --startup-file=no --compiled-modules=existing \
        --project="$PROBE_ROOT/environment" "$@"
}
case "${1:-probe}" in
    install)
        run -e 'using Pkg; Pkg.develop(path=ENV["PROBE_ROOT"] * "/source"); Pkg.add("AMDGPU")'
        ;;
    diagnose)
        run -e 'using AMDGPU; AMDGPU.versioninfo(); println("functional=", AMDGPU.functional())'
        ;;
    probe)
        run "$PROBE_ROOT/source/scripts/benchmark_accelerator.jl" amd 32 10
        ;;
    verify)
        run -e 'using AMDGPU; AMDGPU.versioninfo(); append!(ARGS, ["amd", "32", "10"]); include(ENV["PROBE_ROOT"] * "/source/scripts/benchmark_accelerator.jl")'
        ;;
    *) echo 'usage: bootstrap_hbox_probe.sh install|diagnose|probe|verify' >&2; exit 2 ;;
esac
