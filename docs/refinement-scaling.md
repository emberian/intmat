# Scaling native material simulations

The next useful redesign is **localized thermal computation with dense output
only when needed**, backed by an implicit or exact linear thermal propagator.
The writing equation has considerably more exploitable structure than a generic
coupled simulation: temperature obeys a linear diffusion/cooling equation with
a prescribed moving Gaussian source. Deposit and nutrient do not feed back into
temperature in this writer-only experiment.

## What has actually run

The same physical extent (2,2,2) and five-unit duration were simulated at 384³
and 768³. The 768³ Metal run used 1600 steps of .003125, took 651.77 seconds of
evolution, and exported finite native deposit, temperature and dose fields.
Its compact carrier uses five Float32 volumes,9.06 GB; the former nineteen-array
layout would require 34.43 GB. Default versus compact 96³ runs matched exactly
for 120 steps. This allocation improvement is implemented, not prospective.

| Native grid | Five-array base | Three-field export | Explicit runtime estimate |
|---|---:|---:|---:|
| 768³ | 9.06 GB | 5.44 GB | **10.86 min measured** |
| 1152³ | 30.58 GB | 18.35 GB | ~1.4 h |
| 1536³ | 72.48 GB | 43.49 GB | ~5.8 h |
| 2304³ | 244.61 GB | 146.77 GB | ~44 h |

Estimates scale the measured 768³ time by `(n/768)^5`: cubic cell count times
quadratic explicit-diffusion step count. They exclude compilation, exports and
rendering, and assume comparable effective GPU bandwidth/utilization. They are
planning estimates, not measured limits. Thermal stability bounds at those four
grids are approximately .0033405, .0014970, .0008445 and .0003761. The 768 run
used a slightly smaller step than its bound.

Another 2× or 3× per axis therefore costs roughly 32× or 243× with this explicit
thermal algorithm. Acoustic CFL instead gives approximately 16× or 81× for the
same physical duration. Larger memory on another accelerator, distributed
execution, streaming, or a better representation can change feasibility;
this 96 GB host is not a universal resolution ceiling. Base state also is not
peak memory: staging, solver scratch, runtime and rendering need headroom.

## Three improvements, in order

**1. Remove the diffusion timestep restriction while retaining the physical
field.** For constant diffusivity and periodic boundaries, the discrete thermal
operator separates by axis. Directional implicit diffusion solves can retain
real Float32 arrays and use GPU batched periodic tridiagonal solves. Alternatively,
a Fourier propagator can apply the exact exponential of the discrete linear
operator. This is the same diffusion/cooling model with a changed time integrator;
source and threshold-dose quadrature still need accuracy checks. Fourier state
and work buffers can erase memory savings, especially full complex storage.
A directional implementation should benchmark scratch and bandwidth before
choosing FFTs by default.

This can replace the quadratic step-count growth with an accuracy-controlled
step count. At fixed acceptable temporal error, another 2×/3× per axis approaches
8×/27× work, rather than 32×/243×. It does **not** permit arbitrary large steps:
beam motion, pulsing, threshold crossings and dose integration remain. The
current path's speed is bounded by about 3.18 distance units/time, and its Gaussian
width is .018. At dt=.003125 it can move .0099 units, over half a source width.
Thus the existing 768 step is already a temporal-accuracy question, not free
headroom for a much larger implicit step. Source quadrature along each step
and adaptive threshold integration are more valuable than merely disabling CFL.

A practical first implementation is a positivity-preserving directional
backward-Euler thermal step with explicit source quadrature, followed by a
higher-order variant if error/time comparisons justify it. Require homogeneous
cooling/decay checks, Fourier-mode diffusion, a stationary Gaussian source,
moving-path timestep refinement, and nonnegative temperature/dose. Compare
physical volume integrals and local profiles, not only finite output.

**2. Exploit localized support and eliminate redundant evolving arrays.**
Measured 768³ dose exceeds .01 in only .323% of the domain; deposit exceeds .18
in .312%. Those fractions are evidence of locality, not a guarantee that only
those cells need thermal computation. Allocate active physical tiles around the
source and thermal wake, including diffusion halos; retain written dose tiles
when the hot region moves away. Bin activation must look ahead along the swept
source path. Halo exchange precedes each diffusion stage; shared tile faces
must use the same stencil and boundary values. Dense output can be streamed
from tiles, filling certified inactive regions with their physical defaults.

A Gaussian and diffusion have nonzero mathematical tails. Skipping tiles is
therefore a controlled approximation unless a bound proves their temperature
can never cross the writing threshold. Choose an absolute temperature/dose
error budget; bound omitted source plus propagated heat and cooling, and expand
halos until below it. Test the bound against dense 768 data. A halo chosen only
from observed deposit occupancy is not sufficient. Preserve the original
periodic domain; cropping the thermal domain and imposing zero edges changes
its boundary conditions.

For this isolated writer, nutrient obeys `dn/dt=-conversion*n*max(T-threshold,0)`
and dose integrates the same nonnegative exposure. With initial nutrient 1 and
no other material process, `nutrient=exp(-conversion*dose)` and
`deposit=1-nutrient`. Thus dose can be the sole persistent material variable;
reconstruct nutrient/deposit during output. This removes two evolving arrays
and repeated material writes. It is an analytic invariant of this writer-only
model, but will differ slightly from the reference's repeated Float32 rounding.
Validate the difference and conserve the reconstructed pair by construction.
It does not apply unchanged to mixed metabolic, erosive or externally fed worlds.

**3. Use the linear PDE directly where that beats stepping a full grid.**
The continuous infinite-space temperature from a past Gaussian source at ageτ
is another Gaussian with variance `width²+2 Dτ`, amplitude factor
`(width²/(width²+2 Dτ))^(3/2)*exp(-cooling*τ)`, centered at the past beam position.
Integrating these contributions over the prescribed path gives temperature;
adaptive integration above threshold gives dose. Periodic images restore the
periodic model when needed. This is a concrete route to evaluating only the
written support, without a global diffusion timestep or evolving dense thermal
array. It changes spatial discretization from the current finite-difference
solver, so assess continuum convergence rather than expecting bitwise parity.

Naively summing every past source sample at every voxel can be worse than the
current solver. It needs temporal clustering, Gaussian splatting, tile bounds
or another fast evaluation strategy. Cooling limits useful history at a chosen
error tolerance; source motion and threshold exposure guide adaptive quadrature.
A bounded stationary/moving-source prototype would determine the crossover
before committing to a full implementation. This option is promising because
of the actual equation's structure; it need not be dismissed merely because a
short unoptimized prototype is slow.

## Execution choice

First implement the tiled writer with dose-only material state and a controlled
thermal propagator; keep a dense reference path for comparison. Independently
prototype analytic Gaussian-history evaluation on the same physical paths.
These target empty-space work and the timestep bottleneck together. Multi-GPU
slabs or out-of-core GPU tiles remain useful when dense states are required:
they need halo exchange, source ownership, overlap and synchronized output, and
do not by themselves remove the explicit quadratic step-count cost.

Native 384/768 outputs and measured timing records remain under `outputs/refinement-20261015/`; the comparison includes smaller
timesteps at 768 and is not an isolated spatial convergence-order measurement.

## Implemented optimization and current run

`DoseThermalWriter` now implements the three-array alternative: temperature,
thermal scratch and dose. It shares the exact thermal/exposure helper with
`ThermalWriter`; its private fixed uniform material inventory is reconstructed
only for observations. World nutrient/deposit fields are untouched. An ordinary
material-writing process can coexist on its separate World reservoir; this does
not exchange material with the private writer. Scheduling the same embedded
thermal state twice is rejected through existing ownership preflight.

Thirty-four CPU assertions pass, including exact temperature/dose parity, pulsing,
shutter stop, fractional conversion, private inventory accounting and slab
readout. A 2000-step saturation test reached zero nutrient; maximum absolute
nutrient difference from incremental exchange was 3.70e-12 and deposit difference
1.02e-6. These are measured rounding differences, not an all-parameter error bound.
A warm 96³ CPU comparison measured 49.13 ms original versus 35.50 ms dose-only
(1.38×). A bounded 12³/25-step Metal comparison also passed temperature, dose and
material readout parity. Warm Metal throughput measurements appear below.

`studies/refinement_dose.jl` and `experiments/refinement/dose-writing1152.toml` drive the active run:
three-array 18.35 GB base state,22 GB recipe ceiling, dt=5/3600 below .001497,
and four-plane exports without allocating whole nutrient/deposit readout arrays.
A five-step 16³ CPU smoke run actually exported the expected native fields.

The bounded Green-function probe uses the existing dense writer at 32³/64³/96³
and a stationary Gaussian analytic reference. At dt=.001 and t=.2, temperature
relative differences at sampled points fell to about .14%, while dose remained
about .61% different. Thus finite time quadrature visibly matters even when
spatial resolution improves; this prototype supports the linear-PDE route but
is not a validated moving-source replacement. Exact data are in
`outputs/refinement-20261015/thermal-green-probe.json`.

The Julia dose writer now also swaps two preconstructed thermal handles instead
of copying the temperature field after every step. No extra grid arrays are
allocated. Its temperature/scratch handles are borrowed for one step: reacquire
`p.writer.temperature` after advancing, or copy to retain a snapshot. Dose keeps
its array identity. Both handles participate in ownership checks; readout rejects
aliases of temperature, scratch or dose. The original ThermalWriter keeps stable
buffer identities. Thirty-four CPU tests cover these contracts.

This removes eight bytes/cell/step of copy traffic: 12.23 GB per 1152³ step,
or about 44 TB over 3600 steps, before cache effects. All future implementation
paths described here are Julia work; diagnostic measurements from other stacks
do not authorize a second algorithm implementation or production backend.

The paired warm Julia Metal benchmark now measured 96³ original 1.887 ms versus
three-array ping-pong 1.423 ms (1.326×), and 384³ 36.225 versus 29.926 ms
(1.210×), with exactly equal temperatures after matched steps. These are
shared-host component measurements. The 1152³ four-step allocation pilot also
passed: 18.346 GB state, 907 ms median, all finite. The authorized 3600-step
production is now running at native 1152³ with renderer overlap explicitly
allowed below 40 GB combined working memory; no completed 1152 specimen is
claimed until final fields and diagnostics exist.

The native 1152³ run is recorded in `outputs/refinement-20261015/writing-1152/`;
its loaded source and recipe are frozen there. It passed step 1440/3600 with a
recent median of .966 s/step. This is progress, not a final runtime. The host has
96 GiB RAM. At 1536³, three evolving arrays alone require 43.49 GB; feasibility
depends on scratch, output staging and competing jobs, not a categorical size cap.

Original Julia kernels now also execute on hbox through AMDGPU. The dose writer
measured 3.707 ms/step at 384³ and 29.145 ms/step at 768³ in bounded warmed probes.
These are component timings, not complete bake or cellular-arena throughput.
[Runtime and pinned environment](hbox-runtime.md)
