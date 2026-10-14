# Arena throughput and accelerator placement

## Current complete-step comparison

Matched 2,400-agent CPU runs after workspace integration:

| Path | Step | Allocation/step |
|---|---:|---:|
| Neural, default | 19.13 ms | 28.42 MB |
| Neural, reusable workspace | 11.20 ms | 4.32 MB |
| Composed metabolic arena | 12.74 ms | 5.10 MB |

The workspace shares one synchronous sensing/controller frame between metabolism
and ecology. Birth arrays grow once per step; dead contribution entries are
pruned. These are bounded warm measurements, not sustained live-arena throughput.
The metabolic row adds mechanisms and is not a matched speedup comparison.

## Original allocation baseline

Measured on the local Mac, Julia 1.12.4, one OpenBLAS thread, 64×64×1,
seed 11, fixed founder populations, warmed functions, median of eight calls.
Run `OPENBLAS_NUM_THREADS=1 julia --project scripts/profile_arena.jl results.json`.
The JSON includes allocated bytes and per-process timings. Independent stage
medians are diagnostic and cannot be added: GC and changing ecology state affect
short samples. Births are suppressed using zero initial storage and ample energy.

| Population | Full step | Allocation/step | Sensors+neighbors | Controller | Phase FFT | Snapshot construction |
|---:|---:|---:|---:|---:|---:|---:|
| 256 | 6.45 ms | 9.44 MB | 2.36 ms | 0.98 ms | 0.40 ms | 10.46 ms |
| 1000 | 17.06 ms | 37.19 MB | 15.49 ms | 2.14 ms | 0.40 ms | 31.88 ms |
| 2400 | 49.77 ms | 90.40 MB | 22.87 ms | 7.01 ms | 0.40 ms | 84.00 ms |

This baseline identified population-dependent CPU ecology and allocation, not
64² FFT execution, as the bottleneck. Its sensors created temporary vectors and
sets of nearby cells; its controllers sliced genome arrays and added copies. A GPU FFT port alone therefore has a very small ceiling at this size.
`src/neural_batch.jl` supplies allocation-reduced CPU controller evaluation and a
KernelAbstractions backend-native controller kernel. The backend entry point
requires all arrays already on the same backend; it does not hide transfers.
CPU tests compare both paths with the original matrix-vector expression.
The existence of a kernel is not an end-to-end accelerator speed claim.

The original measured baseline implies about 20 steps/s for 2400 cells, or roughly
1.74 million steps/day before output and ecology population changes. At dt=.05,
a million steps represent 50,000 model-time units and take about 14 hours.
Use a bounded warm run on the actual chosen settings before committing a day:
birth/death array resizing, densities and sensor radius change throughput.

Snapshot construction is synchronous. At 2400 cells, one snapshot each two
seconds consumes about 4% of the measured step budget before JSON encoding and
filesystem cost. Every .1 seconds would consume most of the budget. Prefer
2–5 seconds for live observation, 30–60 seconds for archived samples, and
checkpoint every 5–15 minutes depending on acceptable lost work. Profile
serialization on the full state before shortening checkpoint cadence.
The writer atomically replaces latest state/checkpoints, but append-only history
can grow indefinitely: measure serialized snapshot bytes, multiply by 86400 / archive
interval, and set an explicit retention budget. Controller state alone scales as
2×(hidden+9)×(hidden+sensors)×population×4 bytes for genome plus plastic marks;
checkpoints also retain positions, fields, IDs and arena state.

The local CPU phase implementation explicitly uses FFTW complex arrays.
GPU/TPU execution needs backend-native field and ecology state plus compatible
FFT or a deliberately changed numerical implementation. Host/device transfers
of every agent every step are particularly unattractive for these small controllers.
Batch many independent environments or larger worlds before expecting accelerator
utilization; measure compilation, transfers, output and complete wall time separately.

A matched controller-only follow-up (20 warm calls, random Float32 inputs,
identical genomes) reduced 2400-agent median time from 8.28 to 3.54 ms and
allocation from 23.33 MB to .328 MB; maximum absolute difference was 3.80e-7.
At 256 and 1000 agents it reduced .968→.366 ms and 3.32→1.16 ms.
These are controller results, not full-step speedups or GPU timings.

## Reusable observation architecture

`NeuralWorkspace` now keeps contiguous position/input/output arrays and a dense
CSR cell index (counts, prefix offsets, agent indices). Neighborhood traversal
visits each periodic cell once without a Set, and computes displacement scalars
without allocating a vector per neighbor. Environmental arrays are captured once
behind a function barrier, avoiding boxed scalar access through `medium::Any`.
Buffers grow geometrically and survive population shrink; active frame views
restrict all consumers to the current population. Resource competition has a CSR
implementation preserving the reference proportional Float64 calculation and
Float32 downward rounding. Sensor tests cover tiny wrapped grids, 2D, 3D,
empty populations, controller reference parity, and resource allocations.

Complete sensing plus controller preparation measured **3.39 ms / 1280 bytes**
at 2400 agents and **22.65 ms / 1280 bytes** at 10,000 agents (20 warm calls,
64² grid). The earlier 2400-agent separate sensor/controller costs were about
30 ms and 53.58 MB. These are complete frame costs, not complete simulation
steps; metabolism, plasticity, division, field updates and outputs remain.
Use `ARENA_PROFILE_POPULATIONS=2400,10000` with the profiler to inspect larger
populations. Fixed contiguous arrays and CSR bins make a future resident GPU
implementation possible; they already remove the largest CPU allocation source.

Actual Metal controller-only measurements by `scripts/benchmark_controller.jl`
were 3.32 ms CPU versus .890 ms resident GPU at 2400 agents, but 11.77 ms when
uploading all controller arrays and downloading output. At 10,000 agents those
were 9.51 / 2.17 / 23.09 ms. Maximum Float32 reference error was 1.19e-7.
Keeping controller state on-device is required for this path to be useful.

## Contract for a complete resident GPU ecology

Implemented now: reusable CPU observation buffers, contiguous controller inputs
and outputs, dense CSR bins, synchronous complete-frame sensing, and fairshare
CSR allocation. Genome/plasticity arrays remain exact-population CPU arrays;
birth/death resizing is not yet a fixed-capacity GPU carrier.

The remaining resident implementation needs these explicit stages:

1. Allocate persistent genotype, plastic marks, agent state and active-slot arrays
   with spare capacity; retain stable IDs separately from physical slots.
2. Rebuild cell counts, prefix offsets and stable CSR occupancy on-device.
   Complete this barrier before sensors read positions or occupancy.
3. Gather sensors and infer all controllers from the same pre-exchange world and
   agent frame. Generate nutrient/erosion requests without changing field stock.
4. Reduce requests per cell and compute proportional allocations. Preserve the
   Float64 summation and downward Float32 rounding contract, or explicitly test
   and document any changed GPU reduction precision/order. Complete allocation
   before agents commit actions; otherwise execution order changes competition.
5. Commit funded actions, plastic marks and material/work ledgers. Reduce shared
   signal, remodeling and deposit writes per cell, avoiding racing atomic updates
   whose order silently changes the physical result.
6. Scan birth/death flags, reserve child slots, and compact or maintain active
   slots. Carry parent/lineage/ID relationships through the permutation. Seed
   mutations from the established seed/event/parent-ID contract rather than GPU
   thread scheduling. Finish this barrier before the next frame observes state.

Phase/medium fields and compatible FFT operations must also stay resident;
controller residence alone is insufficient. Snapshots should copy a bounded,
consistent post-step frame asynchronously to a host staging buffer, with at most
one pending transfer and explicit ownership until writing finishes. Checkpoints
must serialize a completed step including RNG/event counters and active IDs.
These stages describe remaining implementation work, not an existing full GPU
arena. The CPU workspace is the implemented array/interface foundation.

Runtime guards are also part of the execution contract. Current `step!` and individual `advance!` preflights can repeat validation; GPU reductions that return extrema to the host synchronize execution. A resident backend must retain the numerical checks while eliminating redundant scans and unnecessary host round trips—for example, device error flags at explicit synchronization points. Disabling checks without a replacement would change the contract rather than optimize it.
