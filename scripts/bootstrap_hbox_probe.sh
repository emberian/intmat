#!/usr/bin/env bash
# Run on hbox after provisioning an isolated Julia + ROCm directory.
# All Julia writes stay in PROBE_ROOT; source is a copied Interior checkout.
set -euo pipefail
: "${PROBE_ROOT:?Set the isolated probe directory}"
# The trailing colon retains bundled system depots, excluding the normal user
# depot. Writes remain in the isolated first depot.
export JULIA_DEPOT_PATH="$PROBE_ROOT/depot:"
export JULIA_NUM_THREADS=2
export JULIA_NUM_PRECOMPILE_TASKS=1
export JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1
export ROCM_PATH="$PROBE_ROOT/rocm-debs/opt/rocm-6.4.3"
export LD_LIBRARY_PATH="$ROCM_PATH/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PATH="$ROCM_PATH/bin:$PATH"
julia="$PROBE_ROOT/julia-${PROBE_JULIA_VERSION:-1.11.7}/bin/julia"
project="${PROBE_PROJECT:-$PROBE_ROOT/environment}"
command_timeout="${PROBE_TIMEOUT_SECONDS:-900}"
run() {
    SWARM_MEM_MAX=8G /home/hbox/bin/swarm-build taskset -c 0,1 \
        timeout --signal=TERM --kill-after=15s "$command_timeout" "$julia" \
        --startup-file=no --compiled-modules="${PROBE_COMPILED_MODULES:-yes}" \
        --project="$project" "$@"
}
case "${1:-probe}" in
    install)
        # Reproduce the released diagnostic baseline; this stack did not pass
        # Interior GPU compilation. Other retained manifests are selected via
        # PROBE_PROJECT and must not be silently replaced by this mode.
        run -e 'using Pkg; Pkg.develop(path=ENV["PROBE_ROOT"] * "/source"); Pkg.add([Pkg.PackageSpec(name="AMDGPU",version="2.8.0"), Pkg.PackageSpec(name="GPUCompiler",version="2.6.0"), Pkg.PackageSpec(name="LLVM",version="9.13.2")]); Pkg.pin(["AMDGPU","GPUCompiler","LLVM"])'
        ;;
    prefetch)
        # Install the declared lazy artifact without loading the affected JLL.
        run -e 'using Pkg, LazyArtifacts; info=only(filter(p->p.name=="LLVMDowngrader_jll",collect(values(Pkg.dependencies())))); Pkg.Artifacts.ensure_artifact_installed("LLVMDowngrader",joinpath(info.source,"Artifacts.toml"))'
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
    *) echo 'usage: bootstrap_hbox_probe.sh install|prefetch|diagnose|probe|verify' >&2; exit 2 ;;
esac
