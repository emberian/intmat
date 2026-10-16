# Implicit phase runtime and zero-diffusion coupling

The CPU process keeps its original physical struct and serialized schema. Its
Float64 real/complex work arrays, Float32 physical composition and forcing,
finite-difference spectrum, stabilization and arithmetic remain unchanged.
`Interior.ImplicitPhaseRuntime(p)` owns the original three FFTW ESTIMATE plans
outside that state. `Interior.advance!(w,p,dt,runtime)` lets a persistent worker
reuse them. The ordinary advance still creates fresh plans. A runtime belongs
to its worker and compatible buffers; concurrent mutation of one process is
unsupported, as it was before.

Plan reuse did not justify an automatic cache. Interleaved complete brine steps
at 96³ gave median 84.59 ms with reuse and 84.75 ms with replanning. At 192³,
682 versus 656 ms and substantial host variance likewise showed no reliable
benefit. The earlier sequential pilot is preserved alongside the paired report:
`benchmarks/2026-10-16-implicit-phase-runtime-cpu-{pilot,paired}.json`.
No automatic global cache was adopted. Physical copying/serialization never
captures native FFT resources; explicit-runtime continuation matches fresh plans.

The useful optimization is in the existing phase coupling kernel. When nutrient
diffusion is exactly zero, skip the neighbor/permeability loop while retaining
its zero flux and original final nutrient/activity assignments. Brine uses this
regime, so the previous kernel evaluated permeability and tanh repeatedly for
no nutrient transport. Nonzero diffusion uses the same loop as before.

The complete coupled brine benchmark uses the original initialization, forcing,
phase advance, harmonic-face diffusion twice, and bounded mineral reaction.
An automatically source-derived oracle removes only the zero-D guard, retaining
all production equations and the complete phase schedule. Interleaved local
2-thread measurements gave 226→122 ms at 96³ and 1370→566 ms at 192³. Absolute
host timings varied; these are bounded workload measurements, not isolated FFT
or single-kernel claims. The original report is
`benchmarks/2026-10-16-implicit-phase-zero-diffusion-cpu.json`.
At 192³ dt=.01, compared with .04 at 96³, preserves the explicit solute positivity
bound under the same extent and diffusion constants. Owned arrays were about
69 MB and 553 MB respectively, excluding native FFTW allocations.

`test/implicit_phase_runtime.jl` adds 30 assertions covering explicit-runtime
continuation, physical deepcopy/serialization, D=0/nonzero, relaxation=0/nonzero,
all World fields/scratch, and complete coupled brine equality. Existing implicit
phase tests pass 187 assertions, including signed conservation, temporal
refinement, odd/singleton FFT layouts and checkpoint alias continuity.
`studies/implicit_phase_runtime.jl` contains the bounded timing interfaces.

Native Float64 FFT capability and a full coupled GPU process are separate gates;
no CPU timing here claims a GPU speedup or changes precision.

## Genuine Julia Float64 native FFT capability

`benchmarks/2026-10-16-rocfft-capability.json` records successful AMDGPU/rocFFT
real transforms at 9×7×5, 16³ and 96³. Maximum forward absolute discrepancy
against FFTW was 7.50e-11 at 96³; inverse roundtrip error was 1.75e-15. The
96³ warmed transform pair took approximately .44 ms. These transform-only
numbers establish capability and are not a complete GPU process speedup.

The standalone runtime had no rocFFT. The matching official
[rocFFT 1.0.32 ROCm 6.4.3 package](https://repo.radeon.com/rocm/apt/6.4.3/pool/main/r/rocfft/rocfft_1.0.32.60403-128~22.04_amd64.deb)
was unpacked under the hbox probe root's `rocfft-standalone-6.4.3` directory.
Its `overlay` directory symlinks the established runtime's components and adds
only this private library. Setting ROCM_PATH to that overlay before Julia starts
lets AMDGPU discover it. Established runtime files, system packages and drivers
remain unchanged. Set LD_LIBRARY_PATH to the overlay's lib directory and retain
the documented HSA override, optlevel 1 and resource controls.

Normalized AbstractFFTs inverse plans call rocBLAS for scaling, which is absent
from this standalone runtime. The working probe uses the existing unnormalized
`plan_brfft`, followed by ordinary Julia broadcast normalization. It introduces
no replacement FFT algorithm and needs no BLAS installation. The replay driver
is `scripts/probe_rocfft.jl`; its output path is its first argument. Versions 1–4
of the remote probe/log preserve launcher, load and normalized-scaling failures;
version 5 completed under a 300-second/two-thread/8-GiB bound. The GPU was released.

## Shared native process and complete workload gate

The current implementation supports a backend with Float64 arithmetic and real
FFT plans. Physical arrays remain Float32 and FFT work arrays remain
Float64/ComplexF64. Metal reports no Float64 support and receives an explicit
error; no reduced-precision mode is selected. The physical struct field layout
is unchanged. CPU retains its original split normalized inverse. The device
runtime uses the backend's existing real forward and unnormalized inverse plans,
then a Julia kernel applies normalization. Inline nonlinear-potential and
spectral-update functions are shared by the CPU loops and device kernels. One
advance schedule owns their ordering. Runtime resources bind to their backend,
buffer identities, element types, sizes and strides; incompatible reuse errors.
Invalid dt is rejected before default plan construction.

`medium_diffuse!` now has a genuine device path using the same harmonic-face
flux and scalar stencil as the generic CPU path. The optimized CPU array method
and its face order remain intact. GPU scratch must not alias either input;
such calls error before mutation. Historical CPU alias behavior is preserved.
Recipe scalar conversions in the reusable brine execution helper are performed
outside broadcasts, avoiding a GPU type-object argument without changing math.

The immutable hbox gate `source-implicit-phase-backend-20261016-v6` completed
with its environment sibling, isolated FFT overlay, optlevel 1, CPU 10–11,
two Julia threads, 8 GiB/no-swap and nice 15. Results are in
`benchmarks/2026-10-16-implicit-phase-backend-gpu.json`; replay with
`scripts/probe_implicit_phase_backend.jl`. The unchanged CPU equations passed
523 assertions across phase, runtime/checkpoint and compartment tests. Native
boundary checks add 21 phase and 120 transport assertions covering odd and
singleton layouts, zero mode, constant forcing, anisotropic periodic transport,
closed/unequal permeability, zero diffusion, positivity, mass and rejected aliases.
Earlier frozen versions/logs preserve composition, alias and accuracy-gate
failures; no historical source or producer was changed.

The paired timings below are medians of three warm blocks, each five complete
steps followed by the same eight-field physical download on both backends.
They include validation, clocks, observation collection, phase coupling, and for
brine the existing forcing/permeability updates, two conservative transports and
bounded reaction. They exclude frame-file encoding/writes and checkpoint I/O.

| Grid | Complete phase CPU / GPU | Complete brine CPU / GPU | Brine ratio |
|---|---:|---:|---:|
| 96³ | 162 / 21.1 ms | 334 / 16.5 ms | 20.2× |
| 192³ | 1402 / 88.9 ms | 2640 / 87.6 ms | 30.1× |

Brine CPU/GPU setup took 1.92/2.01 s at 96³ and 10.36/1.87 s at 192³. Those
components include initialization, transfers, runtime planning and first-use
specializations, so they are not clean incremental compilation measurements.
Cold first brine steps are recorded separately. The complete gated job took
64.97 s, including both workloads and resource-rebuilding restarts. It is a
short numerical/performance gate, not a morphological long-run conclusion.

At 192³ the largest primary composition/A/B discrepancy was 1.19e-7; mineral
was 4.47e-8 and interface 3.58e-7. Derived chemical potential differed by
2.54e-5 because its Laplacian amplifies composition ULPs at h=.125. Its gate uses
the measured composition difference times
`bulk*(3maxabs(phi)^2+1) + 4kappa*sum(h^-2)`, plus Float32 evaluation rounding;
the observed discrepancy is below that bound. Each field's raw error and bound
is reported. Signed phase mean error stayed below 4.5e-10 and relative chemical
mass error below 2.8e-9. Clocks matched exactly. dt=.04 at 96³ and .01 at 192³
preserve the explicit solute bound; an implicit phase solver does not remove
that solute constraint.

## Physical checkpoints and alias preservation

Raw GPUArrays serialization loses repeated device-array aliases. Use
`Interior.implicit_phase_checkpoint(w,p)`, serialize that host payload, then
`Interior.implicit_phase_restore(payload; backend=...)`. The reusable transfer
helper `implicit_phase_transfer(w,p,backend)` uses one identity memo across all
World and process arrays and reconstructs one canonical process. The checkpoint
has no FFT plans or workareas; callers rebuild `ImplicitPhaseRuntime(restored)`.
CPU and native tests verify actual writes through restored aliases, new buffer
identities, exact immediate physical state, and continued coupled evolution.

The strict payload pins Julia version and the implicit-phase, phase and core
source files, which contain this transfer, grid/World definitions, launch and
phase equations. It does not claim to checkpoint arbitrary additional processes
or their source laws. A coupled producer must checkpoint its private A/B arrays
and identify its recipe/source separately. The gate does so explicitly.

## Memory ceiling before larger allocations

`benchmarks/2026-10-16-implicit-phase-memory-plan.json` gives byte counts, not
allocated larger runs. Let N=n³ and H=(n/2+1)n². Carried fields cost 28N bytes,
derived phase observables 8N, other scratch 16N, FFT arrays 8N+36H, and the two
native plan buffers 32H. Their observed additional rocFFT workareas were zero
at both measured sizes. Higher-size FFT workareas remain unmeasured.

| Grid | Carried fields | Derived / other scratch | FFT arrays / plan buffers | Active native total |
|---|---:|---:|---:|---:|
| 384³ | 1.48 GiB | .42 / .84 GiB | 1.38 / .85 GiB | 4.97 GiB |
| 576³ | 4.98 GiB | 1.42 / 2.85 GiB | 4.64 / 2.86 GiB | 16.75 GiB |
| 1152³ | 39.87 GiB | 11.39 / 22.78 GiB | 37.06 / 22.82 GiB | 133.92 GiB |

The device reports 11.98 GiB total. Measured 192³ active owned arrays/plans were
.62 GiB; total driver/allocator use after restarts was 2.55 GiB, including a
.92-GiB baseline and retained allocations. The 192³ whole gate peaked at
4.47 GiB host cgroup memory. These measurements include CPU reference states,
full observations, serialization buffers and repeated native reconstruction.

A single 384³ producer's 4.97-GiB active device state leaves plausible headroom
on this GPU, but full CPU initialization/reference, host staging and checkpoint
serialization cannot simply be scaled into the 8-GiB host limit. One full
384³ observation is 1.69 GiB; one field is .211 GiB. Its physical checkpoint
arrays cost 4.12 GiB before serialization's output buffer. A bounded producer
needs native/shared initialization or released host temporaries, streamed field
capture, and a checkpoint plan avoiding simultaneous full staging and a second
serialized copy. Device restart temporarily owning two states also needs a
separate headroom check. No 384³ or larger allocation has been attempted.

576³ and 1152³ exceed this device's capacity in the current representation even
before driver overhead. Reaching them requires representation/workspace changes
or a larger/distributed/out-of-core execution plan; changing FFT precision alone
would not make this Float64-reference workload fit. Higher resolution also
reduces the explicit solute timestep, so time-to-equal-physical-time must be
budgeted separately from raw step throughput.
