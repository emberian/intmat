# Exact material-feedback sensing fusion

This study removes a redundant cell-support traversal and launch from enabled
native material feedback. Production first launches `dc_sense!`, then immediately
launches `dh_material_sense!`; both traverse the same support and read the same
pre-controller state. The second pass averages raw food, bound food and catalyst
concentrations, then writes the three optional material receptor slots.

`studies/feedback_factorization.jl` derives a fused kernel from current production
source. Raw, bound and catalyst sums retain independent accumulators and original
per-support arithmetic order. Saturation still follows averaging. The original
chemical, recurrent, bond, contact and internal sensors retain their expressions.
Material receptors are written after native sense finishes, so its initialization
cannot erase them. No inventories or physical laws change. Both materialized and
axis descriptors remain supported. Inactive cells and prelatched failures retain
the native outer guard. New failure latching during a nonfinite sensor evaluation
can suppress other cells differently under parallel scheduling; exact failure
snapshot equivalence is not claimed for invalid inputs.

The adapter requires explicit enabled material receptors. Ordinary disabled
execution is left on the production path. The copied enqueue schedule is a
source-derived temporary benchmark, with hard source markers and recorded SHA256
hashes. It is not a second maintained simulator. Production adoption should use a
narrow shared sensing hook or generated schedule rather than carrying a fork.
The generated methods have `ff_` names and are installed only on explicit study
inclusion; production files are never modified.

The test checks complete snapshot equality over twelve consecutive steps for
materialized/axis support and two channel sets covering totals, signed balances
and individual species. The full fixture includes inactive capacity, births and
paid constructive chemistry. The benchmark compares full ordinary and fused
steps, including status synchronization, for four and sixty-four founders on
16³ with capacity256/1024 and production mechanics. It alternates measurement order,
excludes compilation, compares full snapshots after every measured block, and
records all samples. This measures single-trajectory complete-step throughput.
No native AMD or Metal claim follows from CPU measurements.

Run with a direct installed Julia binary and two threads:

```sh
JULIA_NUM_THREADS=2 julia --project=. --startup-file=no test/feedback_factorization.jl
JULIA_NUM_THREADS=2 julia --project=. --startup-file=no -e 'include("studies/feedback_factorization.jl"); benchmark_feedback_factorization("benchmarks/2026-10-16-feedback-factorization.json")'
```

The completed CPU fusion gate passed 52/52 assertions: 48 exact complete snapshots
and four periodic-boundary/disabled-path checks. The first gate took 3m26 including
compilation; a later warm-source repeat took 1m47 plus 54s for boundary/F32 checks.

Two-thread local CPU whole-step brackets gave 1.077× for four founders and 0.929×
for 64 founders, with wide 93–240 ms/step brackets under shared load. Every compared
full snapshot matched exactly. These noisy brackets establish no CPU speedup.
Compact evidence: `benchmarks/2026-10-16-feedback-factorization.json`.

The native AMD 64³ gate completed normally on the RX 6750 XT with documented Julia
1.12.4/opt-level1 runtime, CPUs 10–11, 8 GiB/no-swap/nice15, two Julia threads.
It compared five alternating 20-step brackets, compilation excluded:

| Founders / capacity | Ordinary median | Fused median | Ratio |
| --- | ---: | ---: | ---: |
|256 /16384|8.608ms|8.597ms|1.0013×|
|2048 /32768|11.387ms|11.369ms|1.0016×|

All twelve receptor channel codes passed exact same-state input comparisons
under both materialized and axis support (eight three-channel checks). Sparse
complete snapshots were exact. Dense trajectories retained exact identities and
bond arrays; maximum checked physical-state/field difference was 5.96e−8 over
101 steps. Scheduling can alter native atomic accumulation order, so this is not
an assertion of bitwise dense trajectory identity. The physical comparison checks
positions, masses, hidden/adaptation/action states, metabolites, batteries,
resource/feedstock fields, bound food and catalyst; it is not all scratch arrays.

The benefit is negligible. Keep this as a negative optimization study; do not
adopt the copied enqueue schedule into production for this result. The useful
next question is work avoidance for zero diffusion coefficients, which currently
still traverse every voxel's six neighbors in feedstock transport.

GPU evidence: `benchmarks/2026-10-16-feedback-factorization-gpu.json`. Frozen
hboxsource/environment suffix is `feedback-factorization-20261016-v3`; its
`feedback-factorization-gpu.jl` contains the exact native gate driver. The v1
startup failure (JSON3 imported directly rather than through Interior) and v2
comparison-helper typo (`positions` versus native `position`) remain preserved
with their logs. Neither failure was a numerical fusion failure. The completed
v3 source was never mutated while running. GPU was released directly to the
learning lane after completion; no additional GPU probe was launched.


## Zero-diffusion follow-up candidate

`install_zero_diffusion_factorization!()` explicitly installs a separate source-
derived transport kernel and schedule. It retains every output expression,
scratch write and commit launch. A species with declared D=0 skips coordinate,
neighbor, face-permeability and stencil work, then writes the original
`value + dt*zero(value)` expression. Positive coefficients follow the original
stencil unchanged. Both base resource and habitat feedstock passes use it.
This differs from omitting the whole transport pair: scratch/checkpoint arrays
can remain exactly equal, including the sign behavior of zero.

The all-zero feedstock law in the constructed garden currently performs eighteen
neighbor visits per voxel per step (three storage columns ×six neighbors) despite
zero physical transport. The third column is a padded feedstock role. The new
branch removes these neighbor visits without claiming an engine speedup from an
isolated transport benchmark. Finite valid fields/permeabilities are the gate;
invalid NaN propagation details can differ as with other work-avoidance paths.
CPU tests passed 14 unique assertions: ten complete Float32/Float64 coupled
snapshots and four bitwise stencil oracles with all-zero/mixed coefficients and
signed zero. Four native AMD bitwise stencil oracles also passed. Native
complete-step brackets show a useful gain:

| Grid / founders / capacity | Ordinary median | Zero-D median | Ratio |
| --- | ---: | ---: | ---: |
|64³ /256 /16384|8.602ms|8.287ms|1.0380×|
|64³ /2048 /32768|11.412ms|11.081ms|1.0299×|
|128³ /256 /16384|33.879ms|31.527ms|1.0746×|

Each case uses five alternating 20-step brackets after an excluded warm step,
with full status synchronization in every complete step. 64³ and 128³ sparse
full snapshots matched exactly. Dense 64³ identities/bonds matched exactly and
maximum checked physical-state/field difference was 1.49e−8; all scratch arrays
are not claimed bitwise equal under native atomic scheduling. This is a single-
trajectory complete-step gain, not aggregate replicas or isolated controller
throughput. Changing grid resolution is a scaling comparison, not a claim that
64³ and 128³ trajectories agree.

Evidence: `benchmarks/2026-10-16-feedback-zero-diffusion-gpu.json`. Frozen hbox
source/environment suffix is `feedback-zero-diffusion-20261016`; its driver is
`feedback-zero-diffusion-gpu.jl`. The run completed normally under the same
CPU 10–11/8 GiB/no-swap/nice15/opt-level1 controls, and the GPU was released.
The narrow D=0 branch is now integrated into the shared native transport
kernel, preserving the single production schedule. Its whitespace-independent
kernel token SHA256 is
`3634c451d8c6331adda4173f6f6d034edc4aaf06cd7d3fcc2c2e2231e666f8ca`,
identical to the source-derived kernel in the completed native benchmark.
The study installer verifies that fingerprint before creating its explicit old-
stencil reference (`if true` instead of the D=0 guard). Production CPU tests
compare the new native schedule with that reference, including all scratch and
physical arrays. The test is registered in an isolated module in `test/runtests.jl`.
Post-integration CPU validation passed all 14 unique reference/oracle assertions
(ten complete snapshots and four bitwise stencils). The native benchmark used
the identical arithmetic body under its own frozen source.
