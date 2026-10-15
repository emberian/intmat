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

The shared-resource adapter passed 246 CPU assertions: 229 whole-composition and
protocol regressions, 11 backing-storage/slice checks and six multi-level 3D
checks. These include closed and fed inventories, failed-world isolation and
healthy-member checkpoint recovery. A fresh native gate remains to be scheduled
for this source version.
