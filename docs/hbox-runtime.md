# hbox Julia AMD runtime (verified 2026-10-15)

Julia AMDGPU now runs the original Interior dense workload on the RX 6750 XT.
The 8³ CPU/device comparison passed five steps with maximum absolute field error
**3.73e-9**, below the 1e-4 gate. Ten warmed 64³ steps took **10.759 ms**, or
**243.64 million cell updates/s**. This measures the dense fields workload,
excluding continuous arena ecology and the implicit FFT phase. No simulation
kernel or reduction rewrite was needed. All simulation and GPU kernel code is
Julia; the default runtime uses standalone ROCm libraries.

The current three-array `DoseThermalWriter` also passed the Julia 12³/25-step
pulsed/shutdown fixture on this backend: maximum temperature error 7.45e-9,
dose error 4.66e-10, nutrient error 5.96e-8 and deposit error 1.75e-10. A bounded
384³ run then took 92.67 ms for 25 timed advances (3.707 ms/step), after an
excluded warm advance. `scripts/probe_hbox_julia_dose.jl` is only a generic host
launcher: it constructs Interior's existing writer and calls `step!` and
`writer_material!`; it contains no thermal algorithm implementation. Evidence
is in `gpu/hbox-julia-dose-verification.json`.

The same Julia implementation also completed a bounded 768³ run: 25 timed
advances in 728.63 ms (29.145 ms/step), with the same fixture checks passing.
Eight times the voxels took 7.86 times the measured advance time. This is a
short, stable-timestep writer benchmark, not a completed five-unit specimen,
grid-convergence study or whole-arena result. Its three-field base state is
5.44 GB; peak GPU memory was not instrumented in this Julia probe.

This fixture run uses a separate current-source snapshot at
`source-dose-20261015/` and copied pinned environment
`environment-dose-20261015/`, leaving the older dense benchmark source and live
producers unchanged. Its fixture is
`/tank/chreatures/probes/intmat-20261015/dose-fixture/fixture.json`.

For the retained dose snapshot, the following host-only launcher calls the
existing Julia implementation; select 384 or 768 for the bounded measurement:

```sh
export PROBE_ROOT=/othersys/tmp/intmat-probe-20261014
export JULIA_DEPOT_PATH="$PROBE_ROOT/depot:"
export JULIA_NUM_THREADS=2 JULIA_NUM_PRECOMPILE_TASKS=1
export ROCM_PATH="$PROBE_ROOT/rocm-debs/opt/rocm-6.4.3"
export LD_LIBRARY_PATH="$ROCM_PATH/lib" HSA_OVERRIDE_GFX_VERSION=10.3.0
SWARM_MEM_MAX=8G /home/hbox/bin/swarm-build timeout 180s taskset -c 0,1 \
  "$PROBE_ROOT/julia-1.12.4/bin/julia" --startup-file=no --compiled-modules=yes \
  --project="$PROBE_ROOT/environment-dose-20261015" \
  "$PROBE_ROOT/source-dose-20261015/scripts/probe_hbox_julia_dose.jl" \
  /tank/chreatures/probes/intmat-20261015/dose-fixture/fixture.json \
  "$PROBE_ROOT/julia-dose-results.json" 384 25
```

## Working configuration and replay

The isolated directory is `/othersys/tmp/intmat-probe-20261014` on hbox:

| Component | Verified version/location |
| --- | --- |
| Julia | 1.12.4 |
| AMDGPU | upstream commit `c77460e860b35f21277d8f1253e37304bd8b8993` (reports 2.8.0) |
| GPUCompiler | 2.14.1 |
| LLVM.jl | 10.2.1 |
| AMDGPU LLVM backend JLL | 23.1.3+0 |
| LLVMDowngrader JLL | 0.12.0+0 |
| KernelAbstractions | 0.9.45 |
| Project / manifest | `environment112/Project.toml`, `environment112/Manifest.toml` |
| Interior source | `source/` (preserved dense benchmark snapshot) |
| Default ROCm | extracted standalone 6.4.3; HIP 6.4.43484 |
| Device target | `gfx1030`, with `HSA_OVERRIDE_GFX_VERSION=10.3.0` |

The physical device is `gfx1031`, AMD Radeon RX 6750 XT, 11.984 GiB usable VRAM.
AMDGPU also detects its `navy_flounder` topology and applies the architecture
compatibility override. This is an experimentally verified configuration, not a
claim of official support for this GPU/Ubuntu pairing.

The copied `bootstrap_hbox_probe.sh` matches the repository helper:

```sh
export PROBE_ROOT=/othersys/tmp/intmat-probe-20261014
bash "$PROBE_ROOT/bootstrap_hbox_probe.sh" verify 64 10
```

`install` instantiates the retained pinned manifest, prefetches its declared
LLVMDowngrader lazy artifact, and precompiles AMDGPU/Interior. It does not select
new package versions. It requires the saved source and Julia/ROCm installation.
`diagnose` prints runtime information; `probe N STEPS` runs parity and timing;
`verify N STEPS` combines runtime information and the benchmark in one process.
`PROBE_PROJECT`, `PROBE_JULIA_VERSION`, `PROBE_ROCM_PATH`, and
`PROBE_TIMEOUT_SECONDS` can select another already prepared environment.

Normal compiled modules are enabled. **Keep the trailing colon** in
`JULIA_DEPOT_PATH=$PROBE_ROOT/depot:`: writes stay in the isolated first depot,
while bundled system depots remain available and the normal user depot is
excluded. The recovered run created AMDGPU and dependency package images before
GPU execution. The earlier cold `--compiled-modules=existing` path prevented
creating these images and failed GPU compilation.

The [structured results](hbox-accelerator-results.json) preserve manifest SHA256
hashes, exact package/source hashes, all measurements, and historical failures.
The manifests themselves remain in the isolated directory. The default minimal
ROCm runtime has no rocFFT/BLAS/SOLVER libraries; these are not needed for this
dense benchmark.

## Recovery evidence

The successful Julia sequence was:

1. A no-op `@roc` kernel launched and synchronized.
2. A Float32 device pointer store compiled, executed, and passed readback;
   `@device_code_warntype` captured its compiler trace.
3. A KernelAbstractions arithmetic kernel passed exact readback.
4. `minimum(vec(array))`, 3D `minimum`, and 3D `extrema` all passed.
5. Original Interior CPU/device parity and dense timing passed.

The first recovery used read-only ROCm 6.3 libraries already bundled at
`/tank/chreatures/envs/rocm-dev/lib/python3.12/site-packages/torch/lib`, exposed
through the fresh `rocm-bundled/` symlink root. No Torch module was imported or
used as a simulation backend, and that environment was not modified. A second
control switched back to the original standalone ROCm 6.4.3 extraction and passed
all the same Julia gates. Thus those bundled libraries were not required to cure
the compiler failure; the standalone configuration is the replay default.

| Runtime | Grid | Warmed steps | Seconds | Maximum parity error |
| --- | --- | --- | --- | --- |
| Bundled ROCm 6.3 / HIP 6.3.42134 | 32³ | 10 | 0.005935074 | 3.73e-9 |
| Standalone ROCm 6.4.3 / HIP 6.4.43484 | 32³ | 10 | 0.005933782 | 3.73e-9 |
| Standalone ROCm 6.4.3 / HIP 6.4.43484 | 64³ | 10 | 0.010759461 | 3.73e-9 |

These are bounded samples on a shared host, not a universal GPU ranking.
Logs: `recover-julia.log`, `recover-julia-rocm64.log`, and
`recover-julia-64.log`.

The earlier InvalidIR contained dynamic LLVM pointer-store/load calls; boxing
then pulled hostcall allocation into GPU code. Direct host generated-function
expansion of LLVM 10.2.1 `_pointerset` and `_pointerref` emitted valid LLVM IR with
no hidden generator exception (`trace-llvm-generator.log`). The same pinned
Julia/compiler stack now passes with normal package images on both ROCm
runtimes. This identifies the earlier cold source/module-image route as the
practical failure boundary. The exact lower-level generator/context defect has
not been fully isolated, so a more specific causal claim would be premature.

## Historical failures (2026-10-14)

These failures are retained as diagnostics, superseded by the successful replay:

- Julia 1.11.7, released AMDGPU 2.8.0/GPUCompiler 2.6/LLVM 9 failed at Interior's
  `Excitable` extrema reduction. Separate scalar reductions failed too. The
  temporary scalar-reduction edit was removed, preserving original physics.
- Julia 1.12.4 with that released compiler stack failed 1D `minimum(vec(array))`.
- Latest upstream AMDGPU required KA 0.10 and could not resolve against Interior's
  KA 0.9 constraint. The pinned compatible commit above avoided this migration.
- The compatible LLVM 10/compiler 2.14 stack, loaded without cached package
  images, failed a tiny KA arithmetic kernel before Interior. Its missing lazy
  LLVMDowngrader artifact was separately prefetched and SHA256-verified; no JLL
  source was patched.
- A released AMDGPU 2.8.0/GPUCompiler 2.10/LLVM 9 normal-cache control was rejected
  by dependency resolution: compiler 2.10 requires backend JLL 23.1.1+, while the
  AMDGPU release requires 22.1.8. That control never executed GPU code.

An earlier native HIPRTC arithmetic diagnostic also passed exact readback on
native `gfx1031`; the current successful Interior measurements use Julia.
Historical logs and compiler/package hashes are listed in the results JSON.

## Host and resource controls

Host: Ubuntu 24.10, i9-12900 (24 logical CPUs), about 123 GiB RAM. Existing
`render`/`video` membership permits `/dev/kfd` and DRM render access. Julia,
package/depot writes, extracted runtime, and source copies stay in the isolated
probe directory. No apt installation, driver changes, existing-cache deletion,
or changes to the working bundled runtime were made.

Every replay uses `/home/hbox/bin/swarm-build` under the existing shared
`swarm.slice`, 8 GiB maximum RAM, no swap, nice 15, CPUs 0 and 1, two Julia threads,
one precompile task, and a 900-second timeout. Setup/cache creation used a
1200-second cap. The previous busy-host sample showed no cgroup CPU throttling;
all existing jobs and the live arena were left alone.
