# Constructive independent-world batches

This experimental carrier packs multiple existing `DeviceConstructiveEcology`
worlds into one kernel launch schedule. It reuses the production kernel bodies,
shared cellular enqueue, catalytic hooks and sorption law. Raw, free and bound
food, enzymes, matrix, recurrent controllers, reservoir budgets, clocks, random
seeds and failure latches remain independent. Common physical parameters and
geometry are required; each world can have its own cellular and habitat seed.
This is an aggregate ensemble throughput experiment, not a faster single world.

`studies/constructive_world_batches.jl` extends the base study adapter. Private
kernel wrappers change indexing and parameter views, not numerical bodies.
Source hashes pin the adapter and production habitat implementation. CPU and AMD
are supported; other backends are rejected rather than changing habitat types.
No default runtime changes are made.

## CPU evidence

The composition tests passed 140 exact complete-state/budget assertions across
axis and materialized support and closed and finite-flux reservoirs. Another 15
checkpoint/serialization checks and eight failed-world isolation checks passed.
The base adapter regression suite passed 61 assertions. Complete snapshots include
all cellular, raw-food and habitat arrays, with canonical alias ownership checked
on restore.

A four-world, 160-step CPU pilot paired seeds 11 and 29 under closed and finite
food supplies. All sampled complete snapshots matched separate ordinary worlds
exactly. The maximum material residual was 8.26e-6; potential/work/heat residual
was 7.54e-5. Both founders divided in every arm, leaving four cells and no deaths
or capacity hits. Supplies changed chemistry and work slightly, without a measured
population effect. This short fixture establishes composition and independent
budgets; it does not establish ecological adaptation or sustained turnover.

A separate Julia process loaded the midpoint checkpoint from disk and advanced
80 more steps. Final populations, work, heat, supplies, binding, enzymes and
construction counters exactly matched the original 160-step run in all worlds.
The pilot and original reader are retained under
`outputs/world-batches-20261015/constructive-cpu-pilot-1/source`.

## Checkpoint ownership

The experimental envelope contains ordinary per-world production habitat
checkpoints, explicit member indices, Julia version and source hashes. Those
members retain an ordinary reader path independent of the batch adapter.
Restoration validates on CPU before packing/uploading, rebuilds aliases and
requires distinct empty compatible World objects. All-member checkpointing
rejects a failed or partial member; explicitly selecting healthy members permits
recovery without resuming the failed member. Successful-world clocks advance
independently, including when another world has latched a failure.

Run the CPU study with `scripts/study_constructive_world_batches.jl`; resume a
trusted archived checkpoint using `scripts/resume_constructive_world_batch.jl`
with its frozen project. The bounded native driver is
`scripts/probe_constructive_world_batches.jl`. Timing excludes construction,
compilation, bulk snapshots and checkpoint transfers, and includes the entire
coupled enqueue plus endpoint status observation. Native constructive throughput
must be measured independently of the earlier base-cellular batching result.

## Native AMD composition gate

The frozen `constructive-world-batches-20261015-v1` gate passed on RX 6750 XT
with Julia 1.12.4/AMDGPU optlevel 1, two CPU threads on hbox CPU10–11,
8 GiB/no swap and a 300-second process bound. It compared two small worlds and
four worlds with 64 founders each, capacity 256 per world, 16³ voxels over extent
8³, exchange radius 1, hidden width 8 and 16 mechanical substeps. Cellular seeds
were 101–104; habitat seeds were 1101–1104.

For the moderate fixture, median time for four complete steps in all four worlds
was 0.276569 s serial versus 0.119375 s packed: **2.3168× aggregate throughput**
(57.85 versus 134.03 world-steps/s). Three samples were measured after warming.
This includes paid construction, raw/free/bound chemistry, supply accounting,
controller/physiology, mechanics and endpoint status. Bulk state download and
checkpoint round trips are outside timing. This small-grid ensemble result does
not predict a single larger world's speed.

All compared floating arrays agreed within 5.96e-8; discrete identities and
counters matched exactly. Linked-bin scratch ordering was checked for valid
per-world coverage and acyclic in-world links rather than byte equality. The
encoded radius metadata was compared as Float32. Native checkpoint upload and
one post-restore complete step passed. A separately latched member neither froze
nor contaminated its healthy neighbor, which could be exported independently.

Every moderate world reached 128 cells with no capacity pressure and had positive
paid secretion, solubilization and bound food. Maximum absolute material residual
was 1.48e-4 over roughly 2706 material units; potential/work/heat residual was
7.66e-4 over roughly 15291 potential units, below the two-parts-per-million gate.
This is active construction and replication, not evidence of long-term evolution.
Raw samples, per-world audits, hashes and backend details are recorded in
`benchmarks/2026-10-15-constructive-world-batches.json`; the frozen source and
native log are retained under `outputs/world-batches-20261015/constructive-gpu-v1`.

## Shared resource-reduction adaptation

The current experimental adapter also packs the production typed
`DeviceResourceWorkspace`: grants and every compensated reduction level have a
separate final world axis. The shared `device_resource_phase!` still owns the
schedule and arithmetic. Private kernel wrappers project each world's grants,
partials and status; one writer per species per world updates its ledger.
Scratch is reconstructed on restore and never serialized as physical inventory.
Backing arrays are checked for aliasing against each other and all packed process
state, rather than trusting wrapper identity. Slice tests write distinct markers
through each world view and check isolation and nonoverlapping offsets.

This adaptation changes the authoritative resource source hash. The earlier AMD
throughput numbers above refer to its retained v1 pin; they are not a measurement
of the new resource-reduction composition.

The shared-resource adapter passed 249 CPU assertions: 229 whole-composition and
protocol regressions, 11 backing-storage/slice checks, six multi-level 3D
checks and three habitat/raw backing-parent alias checks. These include closed and fed inventories, failed-world isolation and
healthy-member checkpoint recovery. The native composition gate below verifies this source version.


### Shared-resource native gate (`d0e41d1`)

The bounded v2 AMD gate passed on the same 4-world × 64-founder, 16³ fixture,
including the actual-parent ownership guards. Three complete-workload samples
(16 world-steps each) took serial `[0.173556, 0.201348, 0.203551]` seconds versus
packed `[0.045500, 0.051979, 0.052393]` seconds. Median aggregate throughput was
**79.46 versus 307.81 world-steps/s**, a **3.8736×** ratio. The smaller two-world
fixture measured 1.9704×. These are the v2 serial-versus-packed comparisons;
changes from the historical v1 timings are not an isolated resource-kernel study.

All compared floating state differed by at most 5.96e-8 (the small fixture was
exact); identities, counters and world ownership passed. Every moderate world
again reached 128 cells with positive paid construction, raw conversion and
binding. Maximum material residual was 8.40e-5 and potential/work/heat residual
5.11e-4, both below the two-ppm gates. Independent failure, selected healthy-member
checkpoint export, native restore and one complete post-restart step passed.
No scratch is checkpointed. The source is archived from exact commit `d0e41d1` at
`outputs/world-batches-20261015/constructive-gpu-v2-d0e41d1`, including the native
log and raw JSON. The process exited successfully and released the GPU; no long
ensemble was launched.

## Lifecycle watermark adaptation

The current adapter supplies the production typed `DeviceLifecycleWorkspace`
with a separate Int32 high-water marker for every world. Kernel wrappers project
that marker alongside state and resource views. The shared enqueue and shared
prefix/dense lifecycle bodies remain authoritative. The batch constructor copies
owned per-world markers; backend re-homing preserves them. Strict checkpoint
restore reconstructs ordinary per-world markers from restored IDs before packing,
so the marker adds no physical checkpoint payload. Backing-storage guards cover
base state, resource scratch, habitat and raw-food arrays.

The carrier preserves each array's element type, including optional Float64
battery storage; mixing battery representations within one packed batch is
rejected by the common geometry/type contract. Historical native throughput
above remains attached to its exact source pins and does not measure this newer
lifecycle version.

Run the related CPU gates in one process with
`JULIA_NUM_THREADS=1 julia --project test/cellular_world_batches_all.jl` to share
study-kernel compilation across the base, constructive and lifecycle checks.

The lifecycle adaptation passed 278 CPU assertions: the existing 249 regression
checks, 21 marker ownership/birth/death/failure/healthy-checkpoint checks, and
eight Float64 battery checks comparing batched prefix lifecycle with standalone
dense lifecycle. Invalid marker ranges and a marker aliased to a process role are
rejected. Logs are retained under
`outputs/world-batches-20261015/lifecycle-cpu-1`. No native GPU measurement was
made for this lifecycle adaptation.

A further 21 full-composition alias assertions passed, bringing the cumulative
CPU count to 299. Lifecycle markers backed by strided/reinterpreted views into
enzyme, bound food or raw food are rejected before either sibling's complete
physical state, clock or attempted-step count changes. Selected-member checkpoint
export also rejects the invalid owner. These preserve the earlier 278 results;
no additional numerical implementation or GPU run was needed.
