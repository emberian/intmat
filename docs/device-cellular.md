# Experimental device-native cellular process

`DeviceCellular` is a separate experimental Julia/KernelAbstractions process,
exported by Interior. Construct it from a CPU `CellularEcology` reference and a
World on the selected backend. The World may have no fields or scratch buffers:
all chemical/material arrays are process-owned. Existing nonzero World nutrient
is rejected to prevent confusing the two inventories. Generic `step!`,
`Experiment`, process ports and `record_state` are supported; `Arena` deliberately
rejects it until an adapter supplies the full statistics/recording contract.
This is not an interchangeable backend switch for CPU `CellularEcology`.

The complete step keeps resource fields, capacity-major recurrent states and
genomes, physical exchange supports, spatial contact bins, bonds, division
admission, stable identities and inheritance on the device. Capacity is fixed;
a failed radius bound, mechanical resolution, identity range, nonfinite write,
or unrepresentable accounting debit latches a status. Later queued submissions
leave all device arrays and bindings unchanged. The failing step can be partial;
there is no transactional rollback. `attempted_steps` is a host submission count,
not a completed physical step count.

The controller, adaptive hidden state, local plastic marks, complementary
precursor metabolism and work/heat accounting follow the current cellular laws.
Changes are explicit: junctions use simultaneous donor-limited network exchange;
contacts use mutual bounded-degree proposals; mechanical evolution uses checked
fixed substeps and conservatively funded motors instead of CPU backtracking.
Narrow spatial bins use the actual maximum body radius rebuilt on device. The
query stencil covers every eligible periodic pair, including thin grids;
persistent stretched bonds are enumerated separately. The projected execution
recipe caches narrow typed stage arguments. `execution=:staged` retains a wider
argument reference recipe. A fusion candidate provided no useful additional
whole-loop gain and was removed.

Genotype mutation and developmental partition use a counter-based seeded device
random stream. This is inherited stochastic variation, not evidence of evolved
fitness. Optional typed anabolic callbacks provide pre-respiration planning, a shared
growth/export allocator, post-division inheritance, and habitat completion before
the step is marked complete. The base `nothing` hook preserves its physiology;
paid catalyst/sorption composition is developed separately. It shares finite
work/precursor requests rather than requiring growth to exhaust precursors
before construction.

Physical arrays have fixed canonical roles, with guarded device commits from
scratch. References remain valid across steps. The process owns a private
ResourcePool representation; If World nutrient is present, it must explicitly be zero. Replacing
owned arrays, their shapes, process parameters or grid geometry fails preflight.
Direct array-content edits are the caller's responsibility. Both finite resource
budgets and cumulative per-cell ledgers have high/low Float32 compensation.
Waste smaller than a representable grid increment is retained in a compensated
inert residue. Grid/scatter arithmetic remains Float32 and is audited with
explicit tolerances, not claimed mathematically exact.

`device_cellular_step!(w,p,dt;check_status=true,execution=:projected)` checks a four-byte status each
step by default. Batching with `check_status=false` requires a bounded explicit
`device_cellular_status!` or `device_cellular_statistics` boundary. Statistics
reduce on device and download 244 bytes per observation. Their Float32 reduction
roundoff makes them monitoring data; `device_cellular_audit(
 device_cellular_snapshot(p))` downloads all arrays and sums in Float64 for budget
checks. Neither full downloads nor controller-only work are hidden inside the
reported full-step timing.

`device_cellular_checkpoint(w,p)` synchronizes and emits a complete host-array
payload, including scratch, bin/free-list buffers, compensated budgets, IDs and
RNG/submission counters. `device_cellular_restore(w,payload)` allocates independent
arrays on the chosen backend and checks the exact source pin, Julia version,
array roles/shapes/types and complete-state boundary. A failed partial state
cannot be checkpointed. The payload may be serialized by the caller using Julia
Serialization within the pinned execution environment; never deserialize an
untrusted file. Exact replay across different GPU hardware is not promised.

## Initial native AMD evidence

The first populated 3D run used RX 6750 XT, two launch-host threads pinned to
CPUs 10,11, an 8 GiB process budget, 16³ resource grid, capacity 1024 and 16 mechanical
substeps. Over 100 steps of dt .01, population 512→960, with 480 births, 32 deaths and
1,998 final bonds. Complete-step throughput was 39.32/s, including 244-byte compact
observations/failure checks every 8 steps. Process-owned resident arrays occupied
8,417,656 bytes; measured host allocation was229,789,552 bytes over 100 steps.
Construction, compilation and full audit downloads were outside that timer.
Strict Float64 material residual was 5.70e-6 over 4362.88 initial material;
potential/work/heat residual was 4.89e-5 over 14110.72 initial potential.

Evidence is `outputs/device-cellular-native-v8/benchmark.json`; immutable remote
source and log are `/othersys/tmp/intmat-probe-20261014/device-cellular-v8`.
Model hash 831b0ac1b2784657e59c9d6ecf370c26633e42db2c16eb2d5efcdbf8827265b6;
kernel hash 825d6b5bdf2eeec57e771eecfdd7c55f2bea1f6b19c51a137facdf37b596fbf7.
This is a first working whole-loop sample, not a speedup relative to the existing
CPU adaptive law.

The pinned AMD runtime hits an LLVM 18 EarlyCSE/MemorySSA compiler segmentation
fault in some junction specializations at device optimization level 2. The probe
and benchmark expose `DEVICE_AMD_OPTLEVEL=1`: it changes only the current
process's cached AMDGPU compiler configuration, avoiding that pass. It neither
patches installed packages nor changes Julia's host `-O` flag. This uses private
AMDGPU internals and is a pinned-runtime workaround, not a portable public API.
The normal-optimization 2D 20-step fixture also passed with two births, one death
and exact CPU-KA/device primary array agreement. Broader 3D parity and sustained
checkpoint/restart trajectories remain required before production integration.

Run focused tests with `julia --project test/device_cellular.jl`; execute the
portable probe with `scripts/probe_device_cellular.jl cpu|amd|metal NEW_OUTPUT`.
The full benchmark uses `DEVICE_BACKEND`, `DEVICE_CELLS_SIDE`, `DEVICE_STEPS`,
`DEVICE_BATCH`, and `DEVICE_MECHANICAL_SUBSTEPS`; output destinations must be fresh.

## Sustained closed system and robust bin measurements

The earlier capacity-1024 trajectory saturated and recorded 130,053 blocked
birth attempts; it is a throughput fixture, not an unconstrained ecological
result. A separate capacity-20,000 run covered a conservative material-based
live-count bound of 17,451.6. Over 2,000 steps, 512 founders became 1,418 living
cells with 938 births, 32 deaths, 3,651 bonds and **zero blocked births**.
Whole trajectory rate was 29.53 steps/s, including 81.3 MB of selected-field
downloads and 6.42 MB of compact active-row frame files. This run predates the
narrow-bin optimization. Strict material residual was 1.10e-5 and the
potential/work/heat residual was −0.001283 over initial potential 14,110.72.
It is a finite, closed food supply, not a claim of long-term ecological balance.

Every restored checkpoint array matched exactly. A subsequent 20-step AMD
continuation differed by at most 9.54e-7 in position and 1.19e-7 in resources;
only linked-bin list order differed among integer arrays. Future GPU trajectories
are not bitwise reproducible under atomic floating-point accumulation. The
runner declares short-replay tolerances of 1e-4 for primary physical arrays and
1e-3 for all floating arrays, with identity/bond arrays unchanged. Checkpoints
pin source files, controller layout and Julia version; run metadata records
backend/device, compiler optimization and project/manifest hashes.

Five alternating pairs of 200 complete steps measured narrow versus wide bins
with the same projected recipe and 16 mechanical substeps. Median speedups were
**2.722×** at 512 founders (113.20 versus 41.61 steps/s), and **3.003×** at 4,096
founders (88.07 versus 29.28 steps/s). All samples exercised births, deaths and
contacts, had zero blocked births, matched active identities and had identical
final bond sets. These are domain-enlargement comparisons at fixed spacing,
not spatial refinement. They accelerate this declared experimental device law,
not a controller-only surrogate or a claimed equivalent CPU adaptive arena.

Exact candidate coverage is separate from Float32 trajectory agreement.
Independent brute-force checks found no missing pairs among 9,190 eligible pairs,
including 600 threshold/boundary cases. In the paired trajectories, maximum
position differences were 3.15e-5 (512 founders) and 0.001135 (4,096 founders);
maximum position RMS values were 1.63e-6 and 1.16e-5. The larger fixture's maximum
hidden-state/resource differences were 0.00330/0.000653. Different bin traversal
and atomic accumulation order can affect threshold histories.

Compact tracked measurements are in `docs/device-cellular-results.json`.
Full source-copied trajectories are under `outputs/device-cellular-native-v12/`
and paired-bin evidence under `outputs/device-cellular-native-v15/`.
`scripts/run_device_cellular.jl` writes initial/final complete checkpoints and
bounded actual frames. `DEVICE_CAPACITY` is independent of founder count;
`DEVICE_REQUIRE_HEADROOM=1` checks the conservative material forecast. Frames
store active rows, original slot indices, IDs/parents/lineages/generations,
positions/radii, neural states, metabolites and bonds; bond slots refer to the
original slot indices. `DEVICE_CAPTURE` bounds capture frequency.

The base cellular implementation has been exercised on CPU-KA and native AMD.
The script's Metal selector is available for a future compatibility gate; it is
not evidence that this full cellular pipeline has already passed on Metal.

The final default 50-step random-controller 3D CPU-KA/AMD probe passed with
position error 1.31e-6, hidden-state error 7.45e-8 and resource error 9.69e-8,
including two births, one death and nonzero recurrent activity. A separate
500-step capacity-20,000 default run sustained 82.26 steps/s including actual
captures; it restored every checkpoint array exactly and its 20-step continuation
had maximum floating error 1.19e-7. Evidence is source-pinned in
`outputs/device-cellular-native-v16/`. Focused base tests pass 25 assertions;
protocol and independent coverage tests are maintained separately.
