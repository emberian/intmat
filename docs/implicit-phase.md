# Conserved semi-implicit phase evolution

`ImplicitPhaseSeparation(w; kwargs...)` preserves the signed composition,
interface observations, bulk double-well energy, periodic finite-difference
spatial stencil, and nutrient/activity coupling of `PhaseSeparation`. It accepts
the reference constructor's parameters, plus `stabilization=2` and an optional
`chemical_forcing` array. It currently runs on CPU; the explicit reference
continues to support KernelAbstractions backends.

For the nonnegative eigenvalue q of minus the discrete Laplacian, each step is

```
phi_hat_new = ((1 + dt*M*q*S)*phi_hat_old
               - dt*M*q*FFT(bulk*(phi_old^3-phi_old) + chemical_forcing))
              / (1 + dt*M*q*(S + kappa*q))
q = sum_d 4*sin(pi*k_d/n_d)^2/h_d^2
```

This is a stabilized first-order semi-implicit step: the stiff biharmonic term
and linear stabilizer are implicit, while the cubic chemical potential and
external chemical forcing are explicit. The q=0 mode is unchanged, so forcing
remodels composition through conservative flux rather than adding phase mass.
Float32 state introduces small rounding drift; real-to-complex transforms use
Float64 work space and a ComplexF64 Hermitian half-spectrum. The conjugate modes
are redundant for this real material, so the representation changes neither the
discrete operator nor its time integration. Native FFT plans stay local to each
advance, preserving ordinary serialization and deterministic checkpoint restart.
The FFT diagonalizes the *finite-difference* operator, so small-step comparison
with the explicit reference tests temporal error without changing spatial
physics. This follows the stabilized semi-implicit approach analyzed by
[Shen and Yang (2010)](https://www.aimsciences.org/article/id/b4cba61a-377b-449a-b226-e85e2c00cc6a).

The solver requires constant face mobility. Nonzero `mineral_pinning` is
supported exactly for spatially uniform mineral deposit, with effective mobility
`mobility/(1 + mineral_pinning*deposit)`. It rejects spatially varying deposit
when pinning is nonzero. Use the explicit reference for that constitutive law;
choosing `mineral_pinning=0` explicitly defines an unpinned constant-mobility
material suitable for active remodeling. The solver never silently ignores
mineral pinning.

`stability_limit` retains bounds for the explicitly updated nutrient transport
and activity relaxation. With those disabled it returns infinity: this removes
the stiff *linear* CFL restriction, not temporal error or nonlinear restrictions.
There is no unconditional energy-stability guarantee for the unbounded quartic
potential with a fixed stabilizer. Large amplitudes, large timesteps or strong
forcing can produce inaccurate or unstable trajectories. Choose dt by refinement,
monitor `phase_energy`, and increase stabilization or reduce dt when necessary.
External forcing can legitimately increase passive free energy. Composition is
not clipped, which would violate the conserved model.

The persistent public fields `phi`, `interface`, `chemical_potential` and the
writable `chemical_forcing` compose with active processes. Write forcing before
advancing phase. `phase_optical_fields`, `phase_energy` and
`record_process_fields` support the same observation/recording paths as the
reference. `chemical_potential` records the passive material chemical potential;
external forcing is recorded separately. Recordings omit FFT scratch buffers. FFT plans are ephemeral and never enter
the process object graph; ordinary Julia Serialization preserves aliases and
checkpoint restart reproduces the next step exactly. For portable recordings,
construct a fresh process and restore recorded arrays and parameters.
Binary checkpoints made with the earlier full-complex process layout require
their original pinned source. The new real-FFT layout has no implicit binary
migration; its own serialization roundtrip and alias-preserving continuation are
tested.

Validation in `test/implicit_phase.jl` covers 2D and 3D unequal spacing,
conservation, monotone passive energy at dt=0.2, uniform stationary states,
constant/nonuniform forcing, exact uniform pinning and rejected variable
pinning. Halving dt in a same-stencil comparison with explicit Euler halves
the discrepancy, as expected for first-order time integration. The 187 focused
assertions also cover full-complex representation parity on odd, anisotropic,
singleton-axis and true 2D grids with forcing and uniform mineral pinning, plus
current-source checkpoint continuity.

A 192³ same-process comparison on the local Apple M2 Max used two Julia/FFTW
threads and one BLAS thread. Three interleaved warm samples gave median times of
0.643 seconds for the real-FFT core and 1.187 seconds for the previous full-complex
implementation (1.85× throughput), with 14,544 Julia allocation bytes per core
step. The cached production prototype measured 0.671 seconds in that comparison.
The core splits its inverse into complex axes 2–3 and real axis 1 to avoid
FFTW.jl's full-input copy during multidimensional real-inverse planning; no global
plan cache is required. Evidence is retained in
`outputs/refinement-20261015/phase-backend-comparison-192.json`.

An earlier full-complex CPU measurement at 64³ with extent=(32,32,32), bulk=1,
kappa=0.7, mobility=1, pinning=0 and explicit coupling disabled used the seeded
noise plus a smooth cosine mode:

| Solver | dt | Seconds per step |
| --- | ---: | ---: |
| Explicit reference | 0.000317412 | 0.05048 |
| Semi-implicit | 0.025 | 0.04746 |

Four hundred implicit steps traversed time 10 in 18.98 seconds. One hundred
explicit steps took 5.05 seconds, giving about 84 times greater simulated-time
throughput (equal-time explicit cost extrapolated from the measured steps).
Passive energy decreased from 8189.05 to 6812.47; composition remained finite
in [-0.595,0.609]. Signed mass drift was 5.08e-5 summed over 262144 cells,
or a mean drift of 1.94e-10. FFTW planning is included in these warmed timings.

A separate 16², spacing-one trajectory to time 10 compared dt=0.05 and 0.025
with dt=0.0125: maximum field discrepancies were 0.06483 and 0.02276,
respectively. Thus longer timesteps buy throughput but change the transient;
refine dt for quantitative conclusions. This benchmark is not a claim that
implicit and explicit fields agree at those different step sizes or a
hardware-independent guarantee. For arena grids with cell spacing one,
dt=0.1–0.2 is a useful exploratory starting point; preserve physical
extent/spacing explicitly when comparing resolutions.
