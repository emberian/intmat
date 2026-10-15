# Cellular step performance

The cellular implementation remains Julia on CPU. These measurements cover complete
`arena_step!` calls, including resource transport, validation, neural development,
metabolism, contacts and lifecycle. They are distinct from controller-only GPU
throughput and from runner throughput including observations, checkpoints and sleeps.

## Matched checkpoint workloads

On 2026-10-15, immutable copies of three actual broad-landscape checkpoints were
profiled on hbox (Intel Core i9-12900) with Julia 1.12.4. The two existing overnight soaks and their source
snapshot were left running unchanged. The benchmark used two allowed CPU cores,
one GC thread, one BLAS thread, nice 15, the existing `swarm-build` memory controls
and an 8 GiB memory ceiling. Other host work remained active; reported medians are
five warmed samples of eight complete steps, with state cloning outside timing.

The sparse mature case starts at step 206665 with 207 cells and no bonds. The
bonded mature case starts at step 55965 with 148 cells and 71 bonds. The clustered
early case replays the original initial checkpoint for 475 steps, yielding 187
cells and 345 bonds. The exact source snapshot and checkpoint digests belong in
the accompanying machine-readable benchmark evidence.

## Complete-step results

| Starting state | Baseline ms/step | Optimized ms/step | Speedup |
| --- | ---: | ---: | ---: |
| Sparse mature: 207 cells, 0 bonds | 29.598 | 2.585 | 11.45× |
| Bonded mature: 148 cells, 71 bonds | 21.228 | 1.658 | 12.80× |
| Clustered early: 187 cells, 345 bonds | 29.862 | 3.185 | 9.38× |
| Larger initial: 1,024 cells, 0 bonds | 154.185 | 12.641 | 12.20× |

The larger case uses the same initial founder density and grid spacing: 256 × 192
sites over a 128 × 96 domain, compared with 64 × 48 sites over 32 × 24. It retains
the broad profile's rates and one complementary landscape wave. This is a larger
physical domain, not an overlapped population placed into the original area.

All four cases produce **bitwise identical full serialized Arena states after
100 steps**, including lifecycle, bonds, ledgers and array aliases. The clustered
case improves substantially but remains below the 10× target. Exact metadata,
source/checkpoint hashes and sample timings are in
[cellular-performance-results.json](cellular-performance-results.json).

These byte comparisons use the same restore-and-`deepcopy` harness on each side.
Julia serialization is not a canonical value fingerprint: an Arena and its
unchanged deep copy can serialize differently because of internal container
layout. A separately integrated runtime must compare persisted values and alias
relationships, or reproduce the exact harness, before treating a byte mismatch
as a trajectory difference.
The integrated `361d0a9` runtime passed that check: 2,206 persisted values and
89 mutable alias relationships match the baseline, including bitwise array
contents. Its direct-restored state differed only in configuration dictionary
iteration order; reproducing the benchmark's deep-copy fixture also reproduced
the exact 100-step serialization digest.


A separate three-second continuous loop, without cloning or artificial sleeps,
measured approximately 362, 569, 405 and 45 complete steps/s respectively. These
are evolving workloads, not matched-state speedup ratios: the large case grew
from 1,024 to 1,745 cells and formed 1,734 bonds during 136 steps. CPU time was
approximately wall time. This establishes a busy simulation thread, not saturation
of the i9-12900 or its accelerator. Production runner rates also include snapshots
and checkpoint I/O; `scripts/cellular_soak.jl --throughput` removes its configured
batch sleep and records that choice in the session metadata.

## Removing dynamic numerical execution

The initial sparse profile placed about 73% of active samples in the plasticity
loop. The checkpoint-compatible `parameters::NamedTuple` field hid its concrete
type from inference, boxing scalar accesses throughout physiology and development.
Function barriers now specialize the main numerical loop, sensor construction,
porosity and division on the actual parameter tuple. Medium access is specialized
at the same boundary; density is passed explicitly to the radius calculation.
The serialized struct layout is unchanged.

Mechanics had a second inference problem. Its Boolean force-evaluation option
produced a union of energy-only and force-bearing states, widening loop-carried
positions, timesteps and energy terms. Internal `Val(true)` / `Val(false)` dispatch
now gives each path a concrete result. Let-bound snapshots also prevent reassigned
loop variables from being captured dynamically in reduction closures. All force,
energy, resource and lifecycle equations retain their previous operation order.

Contact broadphase now uses integer linear bins with linked heads. A bounded
array replaces hash lookup for ordinary bin counts; a sparse integer dictionary
handles tiny interaction radii without a huge allocation. Pair output is still
sorted canonically and uses the original physical distance predicate. Direct
all-pairs oracles cover both storage paths and periodic neighborhood aliasing.

Mechanics also reuses each conservative skin candidate list for the bond refresh
after an accepted move. Each body's motor and passive displacements are separately
bounded by `travel_limit`, so a pair can close by at most `4travel_limit` during
that split step. The skin includes the configured bond range, and refresh applies
the original exact distance cutoff before sorting candidates by stable IDs. This
avoids a fresh spatial bin rebuild for every adaptive substep; the clustered
checkpoint requires 13 substeps per frame. No integration timestep or energy
acceptance criterion was relaxed.

The periodic displacement calculation retains its original Float64 arithmetic
order but uses subtraction/addition for a single adjacent period, avoiding a
general floating-point remainder on wrapped coordinates. Inputs outside that
range still use `Base.mod`. Signed-zero and seam cases plus 600,060 random and
boundary comparisons match the general modulo bitwise. The focused transport,
contact and cellular suite passes 348 assertions, including direct contact oracles
and cached versus fresh bond formation after motion.

The final combined implementation passes the four 100-step exact checkpoint
comparisons above. Historical benchmark snapshots and raw profiles remain in the
scoped hbox probe directory; live overnight runs were not migrated or restarted.

## Further hardware work

A full device-resident cellular backend is not implemented. It must keep resource
allocation, periodic contacts and funded mechanics, developmental state, stable-ID
birth/death compaction and the associated ledgers resident together. At these small
populations, transferring only controller inputs and outputs would address the
wrong bottleneck. Batched independent replicas could improve aggregate accelerator
occupancy, but their throughput is distinct from the speed of one trajectory.

## Larger-workload follow-up

The [second performance pass](cellular-performance-pass2.md) adds measured
1.08–1.27× whole-step throughput over the already optimized `dc30e40` baseline
on 1,024-, 1,489- and 4,096-cell workloads, with 23–43% fewer allocated bytes.
It also fixes uninitialized proposal scratch and provides a replayable exact
physical-state and alias comparator. These new measurements use their own
matched baseline and should not be multiplied into the earlier ratios.
