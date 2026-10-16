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
