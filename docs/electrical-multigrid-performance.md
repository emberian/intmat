# Electrical multigrid execution study

The native 513³ amber archive runs from immutable source under
`outputs/art/runtime-amber-archive-20261015`. This study does not modify, restart,
or migrate that run. All kernel and simulation code is Julia.

## Cost map

The active harmonic solver is `studies/branch_metal_multigrid.jl`, not the core
lexicographic SOR reference in `src/electrical.jl`. A 513³ hierarchy contains
513, 257, 129, 65, 33, 17, 9 and 5-cell levels. Every noncoarse level has three
red-black pre-sweeps and three post-sweeps, a residual pass, a child clear, a
27-point full-weighting restriction, and an eight-point trilinear prolongation.
That is 16 kernel dispatches at each of seven levels, one serial coarse solve,
and one final residual dispatch: approximately 114 dispatches per V-cycle,
before the stopping reduction. About 87% of hierarchy cell volume is finest-level
volume. Fine smoothing performs six full-grid-equivalent updates per V-cycle.

The only host readback in `device_solve!` is the scalar maximum residual; it
synchronizes once per V-cycle. There is no host coarse-grid solve or full-field
transfer in its loop. The residual needed for restriction must be materialized;
the separate residual used solely for stopping offers a different memory-saving
opportunity. Finest-level scalar fields are each approximately 515 MiB at 513³.

## Layout change

The shared solver now accepts `workgroup=(32,4,2)`. It moves more neighboring
threads along Julia's contiguous first array dimension than `(8,8,4)`. Arithmetic,
color order, sweep counts, restriction/prolongation, conductor masks, boundaries,
and stopping tolerance are identical. The historical `(8,8,4)` default remains.
There is one implementation of the solver rather than a duplicate optimized law.

The reserved AMD window used the verified RX 6750 XT Julia environment from
`docs/hbox-runtime.md`, CPUs 10–11, two Julia threads, an 8 GiB memory cap and no
swap under `swarm-build`. No other coordinated GPU lane ran concurrently. Both
paths were warmed for the exact grid before timing; measured cold solves reset
the same initial field after compilation. Three or nine warm solves each add one
identical, face-connected conductor site chosen from the baseline harmonic field.
Timed regions contain the complete converged solve and device synchronization;
setup, frontier selection, readback for parity, and independent residual checks
are excluded. Solve order alternates between the two paths.

On native 513³, nine warm paired solves gave a median paired speedup of **1.197×**.
The baseline median was **0.13517 seconds per V-cycle**, versus **0.11272 seconds**
with `(32,4,2)`. Both used the same four or five cycles for every changed geometry.
Potential arrays and reported residuals matched exactly; the independent final
normalized harmonic residual was **6.93e-6**, below the **1.6e-5** solve tolerance.
The 257³ comparison gave approximately 1.28× over three changed geometries.
`(64,2,2)` gave roughly 1.20× at 513³ and did not justify a different default.

A separate 33³, 64-event growth fixture used tolerance 3e-7. CPU and AMD layout
comparisons retained identical growth choices and exactly equal final potentials;
the independent residual was 1.59e-7. The focused CPU tests also compare against
the independent core SOR solve and check exact conductor/grounded values.

Individual timing outliers occurred even within reserved windows. Reported
speedup is the median of *paired* baseline/candidate ratios, not the ratio of
unpaired medians, which can be misleading when four- and five-cycle samples mix.
Earlier cold timings that accidentally included a new reduction specialization
were not used to claim speedup. AMD measurements do not establish a Metal speedup
or a complete electrothermal artwork speedup. Metal remains occupied by the
immutable production; its layout comparison awaits a separate window.

## Fused stopping norm

An optional `DeviceMGNorm(first(levels))` workspace computes the exact stopping
residual and its workgroup maxima together. At 513³ its maxima array is about
2 MiB, versus a 515 MiB residual field. The 256-thread kernel pads its launch,
initializes all reduction lanes, and preserves the original residual arithmetic.
The residual used before restriction is unchanged. Reuse one workspace across
all event solves; its shape/backend must match the fine level.

With `(32,4,2)` already selected, the fused norm provided another **1.086×** median
paired full-solve speedup over nine changed native 513³ geometries. The directly
measured combination against the original layout and full residual was **1.302×**:
**0.13494 → 0.10365 seconds per V-cycle**. The cold physical solve, with compilation
excluded, took **1.6172 → 1.2432 seconds**, with the same twelve cycles. All nine
warm solves used the same four or five cycles and had bitwise-equal potentials
and reported residuals. A separate 64-event AMD growth fixture retained all
choices exactly at tolerance 3e-7. Eighteen focused CPU tests passed, including
unrelaxed residual fields and final workgroups that are not full.

For a **new** producer, the explicit opt-in is:

```julia
workspace = DeviceMGNorm(first(levels))
stats = device_solve!(levels; tol=1.6e-5, workgroup=(32,4,2),
                     norm_workspace=workspace)
```

The returned norm is current. When this workspace is selected, `level.residual`
remains V-cycle scratch rather than a published final residual field; call
`device_residual!(first(levels))` if that full field is needed. Existing calls
without a workspace retain the old behavior. The live artwork is unchanged.

Exact source hashes, raw paired samples, independently checked residuals, growth
trajectories, and derived medians are tracked in
`benchmarks/2026-10-15-electrical-multigrid.json`. Baseline/candidate source snapshots
are retained locally under `outputs/mg-performance-20261015` and in the isolated
hbox probe directory recorded in that report.

## Rejected candidate

A coarse execution prototype batched 33→17→9→5 into one workgroup with barriers,
keeping the original red-black operators and converged solution. It matched
potentials exactly but was slower: native 513³ warm speedup approximately 0.986×,
and much worse at small sizes. Removing dispatches did not compensate for
serializing coarse work into one workgroup. Its frozen source and raw results
remain under `outputs/mg-performance-20261015`; the duplicate implementation is
not retained in the shared solver.

## Reproduce

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. --threads=1 test/electrical_multigrid.jl
MG_EVENTS=9 julia --project=YOUR_AMDGPU_ENV scripts/benchmark_electrical_multigrid.jl amd results.json 513
```

Use the documented hbox environment and its resource controls for the AMD command;
coordinate an exclusive GPU window first. The benchmark also accepts `cpu`.
`MG_RBGROUP` and `MG_DENSEGROUP` select the shared shape and must agree;
`MG_BASE_GROUP` selects the baseline shape; `MG_FUSED_NORM=1` selects the preallocated stopping-norm workspace. These are execution experiments on the
same harmonic equation, not a license to migrate an existing artwork.
