# Independent paid-learning world batches

`studies/learning_world_batches.jl` combines the existing constructive-world batch
adapter with the existing paid lifetime actor/critic. It does not copy a simulator
schedule, chemistry, controller law or learning equation. Original production
kernels and original `NativeLifetimeCredit.choose!`/`reward!` kernels gain only
world-axis projections through parsed Julia syntax.

Each member owns its physical fields, cells, material receptors, random seed,
learned readout, eligibility trace, held action, predictor, reward baseline and
paid-work account. Actor rate/sign, perturbation scale, learning flags, critic
flags, payment coefficient and noise seed may differ between members. Physical
parameters, array geometry, physics timestep, decision/reward strides, action
rows and submitted step clocks must agree. Distinct initial fields, inherited
controllers and identities remain independent; different physical laws currently
need separate cohorts.

```julia
include("studies/learning_world_batches.jl")
learners = [NativeLifetimeCredit.sidecar(p; dt=.02, seed=6201,
    rate=(0., .003, -.003)[i], learning=i>1) for (i,p) in enumerate(processes)]
batch = learning_world_batch(processes, learners;
    identities=["control", "credit", "reversed"])
learning_world_batch_advance!(worlds, batch, 100, .02)
observations = learning_world_batch_observation(batch)
checkpoint = learning_world_batch_checkpoint(worlds, batch)
restored = learning_world_batch_restore(fresh_worlds, checkpoint)
```

The single-step adapter validates ownership at each step. The block adapter
validates at its boundary, enqueues the existing production schedule repeatedly,
and synchronizes once at the end. Every original kernel retains its own member's
failure latch. Each `World.time` advances only by actually completed member
steps. A failed world's private learner freezes even though the cohort's shared
host submission counter continues; failed worlds cannot become restart payloads.
Explicitly selecting healthy members preserves their clocks and experiment names.

Strict checkpoints include ordinary production physical checkpoints, all private
learning arrays and options, source fingerprints, and a per-member journal of
experiment identity, completed/attempted steps, learning step, physical time,
seed and failure status. Restore validates these before changing caller clocks.
A checkpoint between reward windows retains the partially accumulated trace,
held drive, prior inventory value, predictor features and paid account. Restoring
only learned weights would change the next return and is not supported here.

Snapshot/observation collection downloads each packed parent once, then slices
on the host. This avoids downloading all members repeatedly for each member.
Production's ordinary checkpoint writer and reader remain the physical schema
and validation oracle. Sidecar restart adds its own array/domain/source checks.

The CPU gate passed 125 assertions, including exact independent physical/private
state over ten coupled steps, different actor signs and critic modes, nonzero
paid updates, partial reward-window restart, block-versus-step continuation,
source/journal rejection, failure isolation and selected healthy restart. It
includes births through the constructive fixture. The native AMD gate passed
23 checks, with zero observed serial/batch discrepancy, exact immediate strict
restore and complete failed-sidecar isolation. Each member paid nonzero learning
work. Material residuals were at most 3.55e−15 and potential-plus-work/heat
residuals at most 2.68e−14.

The throughput driver reports initialization, complete stepping, observation,
end-to-end time and logical owned bytes separately. Its serial baseline enqueues
all independent worlds before one backend synchronization, then reads each
ordinary owner's status. It does not synchronize on every simulation step.
Warm-up uses each exact capacity/world partition because KernelAbstractions can
specialize launch geometry; an initial CPU pilot exposed substantial first-sample
compilation and is retained separately rather than called a steady-state result.

Measured native cohorts have 8 or 32 worlds and a fixed total capacity of 4096.
Sparse/dense cases retain respectively 32 or 256 total founders. Every member
has an 8³ field, so total voxel count grows with member count and is reported.
These partitions are not ecologically equivalent trajectories. The assay uses
closed, nondividing production worlds with material receptors and explicit
.02 physics/.04 decision/.08 reward clocks; faster learning windows make paid
updates occur in a bounded gate. Performance is aggregate independent-world
throughput, not a faster single trajectory or evidence of useful learning.


## Native throughput and end-to-end cost

The RX 6750 XT gate used the documented Julia 1.12.4/AMDGPU environment, device
opt-level 1, CPUs 10–11, two Julia threads, 8 GiB/no-swap/nice 15. Production sources
were frozen before the parallel lane's subsequent request-helper extraction.
Three alternating serial/batch samples follow each exact cohort's warm-up.
The serial baseline queues all independent owners before one backend barrier;
its speed difference does not come from synchronizing each simulation step.

| Worlds | Total founders | Aggregate step gain | Setup +16steps + observation gain | Batch shared-step wall |
| ---: | ---: | ---: | ---: | ---: |
|8|32|7.19×|3.28×|4.66ms|
|32|32|26.67×|6.90×|4.84ms|
|8|256|7.32×|3.57×|6.08ms|
|32|256|28.14×|7.77×|5.28ms|

All twelve measured serial/batch comparisons had zero discrepancy across the
checked physical floating arrays, private learner arrays and identity/bond
arrays. Paid work matched exactly. These are independent-world throughput gains;
no isolated single-trajectory speedup is claimed. All worlds advance their own
model time, and shared-step wall time is explicitly different from dividing that
wall time by the number of worlds.

At 32 worlds, the sparse median costs are serial setup 203 ms + stepping 2067 ms +
observation 124 ms, versus batch setup 251 ms + stepping 77.5 ms + observation 18.2 ms.
The dense median total drops 2.708 s→.348 s. Each 8-world cohort owns 39,645,368 bytes
(37.81 MiB), each 32-world cohort 44,401,376 bytes (42.34 MiB); this includes private
learning and resource/lifecycle workspaces, not compiler/runtime/driver memory.
Cohorts retain 4096 total cell slots; total field voxels are 4096 or 16384 respectively.

The additional warm batch initialization cost repays after 1–2 shared steps in
these samples. Cold process loading took 6.99 s; the strict native gate plus its
compilation took 67.38 s. The first benchmark's combined serial/batch initial warm-up
took 13.58 s, and exact-case warm-ups took .65–2.50 s. A subsequent benchmark in the
same loaded process needed only .127 s for its common warm-up. These components
include both execution modes and validation; they do not identify an incremental
batch-only compile cost, so no compilation break-even number is invented. A
persistent loaded Julia worker amortizes those costs across interventions.

Compact evidence is in `benchmarks/2026-10-16-learning-world-batches-gpu.json`.
The completed frozen source/environment suffix is `learning-world-batches-20261016-v2`.
Its `learning-world-batches-gpu.jl` is the exact driver. The earlier source/log,
whose throughput stage lacked the workbench recipe after passing the native
gate, remains preserved. The completed run released its GPU slot directly to the
construction lane; no existing history was resumed or reset.

The local two-thread warmed CPU pilot gives 1.243× aggregate/1.200× end-to-end for
four worlds, and .864×/.908× for two. All compared arrays and payments match. This
small shared-host pilot is not a general CPU acceleration claim. Its full samples
are in `benchmarks/2026-10-16-learning-world-batches-cpu-warm.json`; the initial
case-warming failure is retained in the separate CPU pilot report.

The main remaining short-assay cost is setup, now larger than the native batched
stepping. Reusing validated resident cohort storage and retaining the Julia
process are natural follow-ups; they must preserve explicit experiment identity,
source pins and independent branch state rather than silently resetting runs.

## Original host hooks and policy registration

The current adapter executes the sidecar's original begin/sense/prepare/finish
hooks. `NativeLifetimeCredit.policy_launch!` is an execution seam: its native
launcher calls the same kernel/backend/workgroup as before, while a batch
launcher appends the independent world axis. Kernel equations remain in their
original study modules. Registration lists original kernel functions, member
options, member properties, cohort-common properties, private array caches,
constructor/restart functions and source fingerprints. It does not implement a
second host schedule or a pretend device backend.

`register_construction_world_batches!(C)` handles the two-actuator construction
policy; `register_joint_construction_world_batches!(J)` handles the investment,
allocation and reserve policy. Memory ablation is private to each world. Original
source masks are private cached state because the current masks may already have
been permuted. Strict checkpoint schema version 2 records the policy identity,
source, caches and properties, and rejects the prior schema rather than silently
changing the continuation.

The separate WIP extension `studies/construction_lag_world_batches.jl` supplies
`register_construction_lag_world_batches!(P)`, which composes the joint policy with a
held-command intervention. Axis, amplitude, producer and memory can differ per
member. Pulse origin and duration must agree across the cohort. Saved commands
and original source masks survive restart. The pulse's original prepare hook
restores the exact sampled drive after submitting the perturbed requests/reserve;
the batch adapter does not add/subtract commands itself.

The refactored lifetime CPU gate retains all 125 checks. The joint cohort adds
30 checks through 60 steps at dt=.1, crossing three reward windows and a source
permutation, comparing original serial hooks against original batched hooks and
strictly restarting between blocks. Construction's launch seam separately passes
16 actuator/forcing boundary checks and 8 full coupled equivalence/count checks.

The original-hook joint native gate is recorded in
`benchmarks/2026-10-16-joint-world-batches-gpu.json`; its replay driver is
`scripts/probe_joint_world_batches.jl`. Eight actual independent construction
worlds (capacity 256 each, memory enabled/disabled, noise 6201/6277) ran 300 steps
at dt=.02, crossing three paid rewards and the four-unit source exchange.
Every compared physical/private floating array matched exactly (maximum error 0).
Strict restart matched every member, and material/potential accounting residuals
were below 1e-8 with positive paid updates in every member.

Warm 100-step blocks including the same full snapshot/observation schedule took
3.492/3.481 seconds for serial owners versus .462/.455 seconds for the cohort
(7.56–7.65× aggregate throughput). The serial baseline enqueued every owner's
trajectory before a common synchronization. Owned packed arrays totaled about
19.1 MiB. Loading took 22.18 seconds; the first one-step serial and batch calls
took 15.37 and 19.98 seconds respectively, including separate compilation. That
cold difference is not a measured incremental compilation cost: specializations
and the order of loading differ. Initialization was not separately timed in
this gate; use the earlier lifetime end-to-end evidence for that narrower claim.

The immutable native source/environment are
`source-joint-world-batches-20261016-v2` and its environment sibling under the
standard hbox probe root. It used optlevel 1, CPU 10–11, two Julia threads,
8 GiB/no-swap and nice 15. Version 1 stopped before numerical execution because a
script dependency was absent; its source/log remain preserved. The subsequent
cohort-common pulse-property schema extension was validated on CPU, separate
from this frozen native joint gate; it does not change joint numerical kernels.

The pulse cohort CPU gate adds 25 assertions: distinct axes, signs, producers and
memory flags match serial hooks exactly; saved commands match, and checkpoints
restart exactly halfway through the pulse and after it. Mixed pulse origins are
rejected. Native pulse cohorts still require their own gate before scientific use.

## Decision-aligned return queues

The separate extension `studies/cohort_construction_world_batches.jl` registers
`register_cohort_construction_world_batches!(Q)`. It reuses joint construction's
launches and adds the original cohort choose/mature/diagnostic kernels. Queue
length, horizon and action clocks remain common; memory, learning/critic flags,
noise, actor rate and critic gain may differ per member. All queue arrays live
in the existing private state and survive the ordinary strict batch checkpoint.

The policy projector now world-projects chemistry arguments too, because they
may contain private scalars such as critic gain. The physical batch engine's
common-parameter behavior remains unchanged. CPU lifetime 125 and joint 30 checks
pass after this projection change. The cohort adds 28 exact physical/private
parity, partial-queue restart, original-mask, maturity and positive paid-update
checks through H=2, action time 4, drain time 6, with distinct memory and actor-only
versus critic-only members.

The independent native batch gate passed 11 checks on four worlds through 60
steps/time 6, including positive/reversed actors and critic-only gains .01/.03.
Maximum CPU/AMD difference was 7.63e-6 (gate 2e-5); partial-queue restart at step
19 and subsequent native continuations were exact. All members paid and matured
eligible records with strict material/potential accounting. This is a parity gate,
not a new throughput claim. Compact evidence is
`benchmarks/2026-10-16-cohort-world-batches-native.json`.
