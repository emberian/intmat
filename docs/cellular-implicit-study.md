# Cellular implicit mechanics study

The graph-local midpoint prototype has **no demonstrated throughput advantage at
matched accuracy** over the production adaptive explicit method. It remains an
unloaded study, not an arena process or a checkpoint migration.

The actual fixture is the constructive-volume pilot at step 1500, time 150:
2,393 cells, 6,242 bonds, fixed radii/controller inputs/funds. Each measurement
advances the same physical interval 0.1. Refinement caps internal mechanical steps
while preserving the outer motor-funding horizon. Repeated public calls with a
smaller `dt` would change the drag reservation horizon near energy exhaustion and
would therefore confound a numerical comparison.

The solver performs symmetric funded motor half-drifts around a passive,
forward/reverse body sweep. Each body uses a three-dimensional implicit-midpoint
local solve over all incident contact/bond potentials. A conservative periodic
contact cache bounds each accepted local move; objective/passivity gates and the
existing positive motor-work debit remain active. Failure publishes no partial
state. Coordinate splitting does not conserve the center exactly; its drift is
measured rather than projected away.

## Controlled result

The source-pinned follow-up ran on hbox CPU 14–15, Julia 1.12.4, 2 Julia threads,
1 GC/BLAS thread, 8 GiB/no-swap/nice 15 controls and a 180-second limit. It completed
in about 30 seconds. AMD art production was running concurrently; no isolated GPU
timing gate overlapped. Timings are medians of three warmed mechanics-only
measurements, not complete ecology steps. No broad speed claim follows from them.
The 512→1024-step explicit reference difference is RMS 2.27e-6, maximum 1.89e-5.

| Method | Time(ms) | Position RMS | Position maximum | Paid motor work | Allocated MB |
|---|---:|---:|---:|---:|---:|
| Adaptive explicit,62 substeps |244.1|3.52e-5|2.92e-4|0.939925|16.25|
| Midpoint,1 sweep pair |84.0|8.70e-4|6.60e-3|1.160521|27.49|
| Midpoint,2 sweep pairs |125.3|2.88e-4|2.17e-3|1.031086|49.96|
| Midpoint,4 sweep pairs |202.7|7.99e-5|6.01e-4|0.969479|89.48|
| Midpoint,8 sweep pairs |364.5|2.28e-5|3.77e-4|0.951510|159.51|

Every case retained the same bond set and had zero reported work-ledger residual.
Soft overlap remains present; this is not a hard nonpenetration method. The
four-sweep case is faster but less accurate in position and paid work; the
eight-sweep case improves RMS but is slower, has larger maximum error, and spends
about 1.23% more motor work than adaptive explicit. Allocation is also materially
higher. Smooth two-body tests show second-order position and chemical-work
convergence, which does not establish the same order across crowded topology
changes.

An earlier normal-Hessian/global Newton route was substantially more expensive;
an actual 2393-cell warm solve exhausted its bounded 180-second study window.
Repeated global contact discovery inside line search dominated. Radial pair
splitting also failed to earn a useful accuracy/cost advantage. These unsuccessful
prototypes and outputs are retained.

## Sources and reproduction

Study implementations live in `studies/cellular_implicit_mechanics.jl`,
`studies/cellular_pair_proximal.jl`,
`studies/cellular_coordinate_proximal.jl`, and
`studies/cellular_capped_reference.jl`. They are not loaded by Interior.
The corresponding four scoped test files pass 318 assertions, including exact
uncapped-reference parity, analytic convergence, finite-funding cases, periodic
cache checks and transactional failure behavior.

Run the study with a fresh empty output directory:

```sh
julia --project --threads=1 scripts/study_cellular_implicit.jl \
  --fixture outputs/constructive-volume-20261015/mechanics-fixture-middle.jls \
  --output NEW_EMPTY_DIRECTORY --coordinate --midpoint \
  --reference-splits 512 --verify-reference
```

The final captured source/results are under
`outputs/cellular-implicit-study-20261015/host-tailfix-037aa/` locally and
`/othersys/tmp/intmat-coordinate-midpoint-20261015/runtime-tailfix` plus
`probe-002` remotely. The manifest verifies 136 files. The parent exact pin and
first failed 4/8-sweep results are retained separately. The only numerical change
in this fork absorbs a floating-point-sized final time remainder into the final
actual step; it neither drops physical time nor skips work. A two-body regression
requires exactly 4/8 accepted steps over 0.1.

The driver captures loaded default-source hashes, exact study/generated-reference
bytes, fixture hash, runtime, raw measurements and final mechanical states.
These fixed-body outputs are not resumable ecological checkpoints. The original
arenas, source pins and long-running pilots were not modified.
