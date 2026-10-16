# Living-step representation experiments

The next support representation should be chosen from complete-step costs. In
current 64³ living scenes, compact support neighborhoods alone have a small
throughput ceiling. A CPU grants-fusion prototype saves retained scratch memory
and preserves successful trajectories exactly, but has no measured throughput
win. Neither result is a production representation change.

## Current CPU cost and support geometry

`studies/living_support_worklists.jl` profiles the unchanged complete native step,
including materialized support, finite feed, native contacts, division and paid
private learning. The fixtures warm for 220 steps (time 4.4), beyond the paid
actor's doubled reward horizon. Sparse and denser founder populations are 24 and
46, capacity 512, extent 2, physical radius .18 and 64³ grid. Both have births,
actual geometric contacts and paid actor updates.

The captured states have 32/58 residents, 25,590/46,464 support edges and
24,651/44,638 distinct support sites (9.40%/17.03% of the grid). Materialized
support retains 11,239,424 bytes. Support construction costs about 2–3% of the
whole step; the five grid passes that might exploit a support-union worklist
cost about 3–4%. A union's construction/clear/scatter overhead must be subtracted
from this ceiling. Reconstructing axis support is already known to lose complete
GPU throughput (`device-support-performance.md`), so no new union implementation
is justified by these CPU measurements.

Resource grants and their compensated reductions account for about 23–26% of
CPU step time. Reaction/sorption, diffusion and whole-state validation are also
substantial. The source profile compares 23,129,418 array elements exactly per
case, including private learned state, status, scratch, resource levels and
lifecycle watermark. These are intrusive local attribution measurements, not
isolated throughput claims; timing variability was visible.

## Implicit grants prototype

`studies/living_resource_fusion.jl` invokes the authoritative current
`device_resource_voxel!` once per voxel/species and directly accumulates the first
resource reduction level. Supplied and removed channels retain their individual
strided traversal and compensated addition order. Later levels and budget
updates are unchanged. Unlike the older `device_resource_fused.jl` study, this
prototype does not repair endpoint arithmetic or change numerical semantics.

The experimental carrier replaces only the reconstructible V×6 grants backing
with a zero-row array. It retains every physical field, compensated level and
resource/lifecycle contract. This saves 12,582,912 retained bytes at 64³ with
Float64 inventories; core construction still temporarily allocates the dense
workspace, so this is not peak-initialization evidence. The omitted grants are
transient and are deliberately excluded from parity; all other arrays and clocks
match bitwise. The explicit carrier is required: the compact process must never
be passed to an ordinary unfused resource phase. No historical checkpoint is
migrated and no core constructor is changed.

Three alternating 12-step comparisons for each founder population compare
64,669,662 array elements per case exactly, including paid updates and births.
A separate direct gate passes 240 checks for Float32/Float64, volumes
1/7/127/128/129/16,641, and zero/source/sink/mixed exchange. Endpoint fields and
every high/low reduction level match exactly. This establishes successful-state
parity; partial arrays after an invalid-input halt are not claimed identical
because fused work can latch the failure earlier.

Local median speed ratios are .995× and .960× with substantial paired variance.
There is no supported CPU speedup. Keep the prototype for its concrete scratch
saving and as a basis for an isolated device experiment, rather than adopting it
as a throughput default. Raw local reports are
`benchmarks/2026-10-16-living-resource-fusion-cpu-local.json` and
`benchmarks/2026-10-16-living-support-profile-cpu-local.json`.

Run the focused gate with Julia2/OpenBLAS1 and the project environment:

```
julia --project=. --compiled-modules=existing studies/probe_living_resource_fusion_edges.jl
julia --project=. --compiled-modules=existing studies/probe_living_resource_fusion.jl output.json
```

## GPU complete-step attribution

The isolated AMD profile reverses the CPU priorities. At 64³ native median
steps take 21.202 ms (24 founders) and 23.344 ms (46 founders). The warmed
states have 32/58 residents, 8/12 births and 72/139 paid actor updates. Three
alternating native/profile blocks of eight steps each compare every physical,
private, event and workspace array bitwise (23,129,418 elements per pair).

The profiler brackets both core and policy kernels. It records 1,065 launches
per eight steps, including the decision cadence: about 133 per step. Native
host enqueue occupies .70–2.60 ms per step; final synchronization/status read
occupies 18.83–23.13 ms per step. There are no observer frames or scientific
field downloads inside this boundary. This workload is predominantly waiting
for device execution, rather than exhausting host submission capacity.

HIP events add roughly 9–11% wall overhead and include driver/submission gaps;
they are attribution evidence, not isolated instruction timings. The middle
paired samples give the following inclusive times per complete step:

| Kernel group | 24 founders | 46 founders |
| --- | ---: | ---: |
| support construction (two calls) | 5.556 ms | 6.037 ms |
| contact forces (16 mechanical substeps) | 4.308 ms | 5.025 ms |
| support sensing | 3.900 ms | 3.931 ms |
| uptake regulation/scatter | 1.737 ms | 1.902 ms |
| metabolism | 1.586 ms | 1.577 ms |
| resource grants + all reductions | .585 ms | .536 ms |

The resource fusion is therefore a CPU/memory experiment, not the main GPU
throughput opportunity. Conversely, GPU support construction warrants a different
factorization despite its small CPU share. The existing kernel maps one lane per
resident; only a few dozen lanes execute thousands of candidate geometry tests.
The next study parallelizes candidate geometry, then compresses and sums in the
original order and normalizes edges in parallel. It keeps materialized support
consumers unchanged and initially uses extra scratch, rather than claiming an
unproven compact representation.

Hardware/resource pin: RX 6750 XT, verified AMDGPU2.8.0/Julia1.12.4 opt-level1,
CPU10–11 are SMT siblings of **one physical P-core**, Julia2/OpenBLAS1,
nice15, 8GiB MemoryMax and zero allowed swap. The successful profile job costs
60.92 seconds wall, 66.59 seconds user +9.54 system CPU, maxRSS1,608,560KiB,
zero swaps. Two startup/warm-wrapper repairs are preserved separately in the
private frozen stage; no failed timings enter the report. Source/environment:
`/othersys/tmp/intmat-probe-20261014/{source,environment}-living-gpu-profile-20261016-v1`.
Raw report: `benchmarks/2026-10-16-living-gpu-profile.json`.

## Parallel materialized-support construction: useful GPU result

`studies/living_support_factorization.jl` factors only construction into three
kernels: candidate-parallel geometry and unnormalized weights; resident-wise
stable compaction and sequential total in the original x-fastest order; then
edge-parallel normalization. Consumers still read the original compressed
capacity-major support arrays. Inactive counts, fallback support and untouched
stale tails retain their original values. This is execution factorization, not
axis reconstruction or a new physical kernel. It needs 11,241,472 additional
scratch bytes at 64³/capacity 512, including the candidate sites/weights and
per-resident total. It does not currently reduce memory.

CPU source proof compares every existing array in six fixtures: 4³/16³/64³,
each with 24/46 founders. This covers inactive slots and small-grid fallback.
The subsequent isolated AMD gate warms 220 native steps, warms both exact timed
call paths on disposable states, then performs five alternating 20-step pairs
per population. All timed blocks have zero compilation and zero GC time.

| Founders | Final residents | Births | Paid actor updates | Native step | Factored step | Throughput ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
|24|32|8|125|20.953 ms|17.409 ms|1.2036×|
|46|59|13|241|23.233 ms|19.511 ms|1.1908×|

Sparse paired ratios range 1.1999–1.2088×; denser ratios 1.1866–1.2008×.
Every pair matches 23,129,418 existing array elements bitwise: all physical,
private learned, identity/contact, event, status, resource-level/highwater and
transient state. Ten pairs compare 231,294,180 elements. World/private clocks and
device completion metadata agree; attempted host clocks use the same increments
in the shared loop. Initial cloning/strict checkpoint restore, scratch allocation,
compilation, all-state downloads and observer serialization are outside timings.
This is single-trajectory complete-step throughput, not aggregate replica rate.

The successful gate costs 74.64 seconds wall, 83.38 user+14.80 system CPU,
maxRSS 1,658,156KiB (1.58GiB), zero swaps. The resource pin is the same exclusive
RX 6750 XT and one physical P-core/SMT 10–11 described above. A standalone-entry
parse error occurred before simulation and remains in its separate private log;
it contributes no timing samples. The immutable executed source/environment is
`/othersys/tmp/intmat-probe-20261014/{source,environment}-living-support-factorization-20261016-v1`.
The raw report is `benchmarks/2026-10-16-living-support-factorization-amd.json`;
it records the exact study hash and core pins. The local harness is
`studies/living_support_factorization_gate.jl`.

The useful next integration is an explicit parallel materialized-support
execution option with candidate scratch in the process's validated stage
ownership. Keep the serial CPU default: this gate proves GPU throughput, not a
CPU improvement. Checkpoint scientific schema can stay unchanged because the
extra arrays are reconstructible construction scratch; restored processes must
allocate fresh scratch. Full production integration should validate support mode,
backend, dimensions and no aliasing, retain strict current-source checkpoints,
and run ordinary complete-step composition/parity for the selected execution
option. No production source or historical run has been changed by this study.
The GPU/CPU slot was explicitly released to the queued developmental-learning
gate before writing this report.

## Reusable native implementation

The tested factorization is now integrated behind an explicit choice:

```julia
cells = DeviceCellular(world, reference;
    support_strategy=:materialized, support_execution=:parallel)
# Direct founder initialization accepts the same two keywords.
```

`support_execution=:serial` remains the default for every existing constructor
and native recipe. It allocates no candidate scratch. `:parallel` is supported
with `:materialized` only; `:axis` plus `:parallel` is rejected. This selects
execution of the same materialized support law, not a new support representation.
The choice is currently intended for native single-world CPU/AMD execution;
experimental replica-batch schedulers retain their own serial support schedule.
The AMD speed claim remains the frozen prototype evidence above; the production
routing is checked separately without repeating that timing sweep.

The three kernels now live in `src/device_cellular_support.jl`, and all ordinary
`device_cellular_step!`/habitat/living callers reach them through the common
native launch route. Candidate sites, weights and resident totals are owned by
`ownership.stages.support`, with backend, identity, shape, element type and
construction-time nonalias checks against physical arrays, resource workspace
and lifecycle watermark. The ordinary step preflight checks the scratch contract.
The two support rebuilds per living step use the same reusable buffers.

Checkpoint payloads record `support_execution` outside the scientific snapshot.
Restore validates the choice and allocates fresh reconstruction scratch. Existing
strict source and Julia-version checks remain; historical snapshots/runs are not
migrated. The source-pinned checkpoint format version remains 4: its scientific
snapshot/geometry layout is unchanged, and the source pin distinguishes this
implementation metadata addition.

`experiments/native/living-patch-parallel-support-64.toml` selects the production
option for the existing 64³ recruitment law. The native recipe adapter forwards
`[native] support_execution = "parallel"`; older TOMLs implicitly select serial.
Run it through the ordinary living study/workbench, with the usual resource and
immutable-run controls. This recipe is a concrete new-run option, not a resumed
historical trajectory.

`test/device_support_parallel.jl` checks full-state serial/parallel composition
and real division for 2D/3D coarse/fallback and finer grids, Float32/Float64
inventories, stale tails, checkpoint restore/continuation, invalid choices,
workspace alias rejection and preflight shape rejection. It is included in the
normal test entry point. The original DeviceCellular suite also remains a focused
compatibility check for default serial and axis execution.

The production routing/restart gate passes for both founder populations: eight
serial/parallel complete steps compare 23,129,418 existing elements per case,
then three parallel/restored continuation steps compare another 23,129,418 per
case, all bitwise exact (92,517,672 elements total). Final birth counts are 8/13
and paid actor updates 98/190. The checkpoint-selected parallel route survives
restore and uses fresh owned scratch. Its 62.43s bounded job costs 65.28 user+8.85
system CPU seconds, maxRSS 1,816,200KiB (1.73GiB), zero swaps, on the same
exclusive AMD/one-P-core resource pin. No additional throughput sweep was run.
Report: `benchmarks/2026-10-16-living-support-production-amd.json`; harness:
`studies/living_support_production_gate.jl`; frozen executed stage:
`/othersys/tmp/intmat-probe-20261014/{source,environment}-living-support-production-20261016-v1`.
The 91 focused CPU checks and 69 existing DeviceCellular checks all pass.
