# Arena throughput and accelerator placement

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

The current bottleneck is population-dependent CPU ecology and allocation, not
64² FFT execution. Each sensor creates many temporary vectors and a set of nearby
cells; each original controller slices two genome arrays and adds them into a
third. A GPU FFT port alone therefore has a very small ceiling at this size.
`src/neural_batch.jl` supplies allocation-reduced CPU controller evaluation and a
KernelAbstractions backend-native controller kernel. The backend entry point
requires all arrays already on the same backend; it does not hide transfers.
CPU tests compare both paths with the original matrix-vector expression.
The existence of a kernel is not an end-to-end accelerator speed claim.

At the measured baseline, 2400 cells sustain about 20 steps/s, or roughly
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
