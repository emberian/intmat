# Independent-world batching and biological clocks

The experimental batch adapter preserves the tested law, but does not justify
a CPU default change: the reserved moderate fixture ran at 0.586× serial
throughput. GPU throughput remains unverified.

The existing cellular implementation already runs sensing, controller, plasticity,
resource allocation, metabolism, junction exchange and lifecycle once per outer
step. CPU contact mechanics takes adaptive internal steps; DeviceCellular takes
16 fixed internal steps by default. Moving biology to a separate slow clock is
therefore already implemented, rather than a new performance opportunity.

Spatially disconnected contact components are not currently independent CPU
mechanical systems: positive motor-induced potential change is funded using a
global sum of motor squared speeds, and passive trial acceptance uses global
potential. Splitting components changes that allocation and acceptance law.
Calling the existing mechanics function repeatedly on dt/k also changes its
remaining-outer-interval motor reservation horizon. Local clocks need a stated
new approximation and full trajectory/work convergence against a reference
that preserves the original funding horizon. The separate implicit-mechanics
studies investigate that problem; this study does not duplicate them.

## Implemented decomposition

`studies/cellular_world_batches.jl` packs independent, identically configured
worlds into a final array dimension. A typed world view supplies each kernel
with its original array shape. A compact per-world seed vector preserves each
standalone world’s mutation and keyed lifecycle random streams. The NDRange gains a world coordinate. Each world
retains its own resources, contact bins, status latch, identities, lifecycle
random streams and material/work ledgers. No resource or force crosses worlds.
The existing generic mechanical host loop is reused through launch dispatch.

The physiology is embedded in kernel definitions, so this prototype parses the
source AST and creates private kernel names. Its only kernel edits are the index
binding, insertion of per-world state/environment views, and substitution of
the matching world’s seed in the parameter view. The checked batch wrapper calls the shared generic
`device_cellular_enqueue!` schedule directly, with an early rejection of
unsupported anabolic extensions. No host schedule is generated. It never replaces a production method.
Both source SHA256 values are recorded with benchmark output. This is a study
adapter, not the proposed permanent abstraction: extraction of scalar kernel
bodies would make independent-world execution maintainable without generation.

The constructor requires matching physical parameters, geometry and
attempted-step clocks; RNG seeds may differ. Worlds can have different initial states. Construction
copies the existing worlds and is excluded from warmed execution timings. The
first prototype is CPU-KA verified; it does not establish GPU compilation,
throughput, or compatibility of GPU atomic operations with the world views.
The batch helper must currently be included from Main. The ordinary CPU
CellularEcology representation and checkpoint ABI are unchanged.

The proposed benefit is aggregate ensemble throughput. At 16 substeps, mechanics
alone takes 82 launches per world: five per substep plus final bin clearing and
geometry. A batch submits that schedule once for all worlds. A single world's
trajectory does not become physically faster; memory use still grows with the
number of independent worlds. HIP graph replay is a separate study that can
potentially compose with batching after both are validated.

## Verification and reproduction

Run with one Julia thread:

```
julia --project test/cellular_world_batches.jl
julia --project scripts/benchmark_cellular_world_batches.jl
```

Tests compare every state and environment array after complete steps, including
resource and compensated accounting arrays, rather than comparing only positions.
An independent failed-world latch is tested without aborting its healthy peer.
The benchmark warms both complete schedules, includes status synchronization,
uses the default 16 mechanical substeps, reports all-array parity, and excludes
packing and compilation. `BATCH_WORLDS`, `BATCH_CELLS`, `BATCH_SIDE`, and
`BATCH_STEPS` bound its workload. No default is changed.

## Representation result and outstanding gate

An initial eager projection of every state/environment field into a SubArray at
all kernel threads was unsuitable: approximately 104 MB allocated per four-world
macrostep. Lazy typed field projection reduced this to 0.40 MB in the local
single-thread probe, versus 1.02 MB for four serial worlds, with exact complete
state equality. On the reserved hbox CPU 14/15 pair, the lazy SubArray version
still took 0.1153 s versus 0.03576 s serial per ten complete macrosteps (0.31×
aggregate throughput). It allocated 0.656 MB versus 1.056 MB. The source-pinned
probe passed 50 assertions, including real birth/extinction and failure-latch
isolation. Maximum closed material balance error was 5.88e-7 and maximum
potential-plus-work-plus-heat balance error 1.43e-6. These earliest measurements
used a common RNG seed, so their worlds were physically uncoupled but their
per-ID random streams were correlated. They are not independent ecological
ensemble samples.

The hbox probe used two Julia threads, one GC thread, one BLAS thread, cores 14/15,
nice 15, 8 GiB memory cap and no swap. Artist AMD production ran concurrently;
no GPU was used by this study. Logs and the frozen source are retained under
`outputs/world-batches-20261015/` locally and
`/othersys/tmp/intmat-probe-20261014/world-batches-20261015-v1/` on hbox.

The final candidate replaces general SubArray indexing with a concrete
contiguous offset slice, and supplies distinct seeds per world. Results are
reported separately from the common-seed reference. GPU adaptation, atomic
access compatibility, complete-step parity and end-to-end aggregate throughput
remain explicit gates before GPU ecological use. No GPU performance claim is
supported by the CPU measurements.

The distinct-seed contiguous-offset candidate passed 51 complete-state and
accounting assertions locally. With one Julia thread, four worlds of eight
cells on 8² grids, default 16 mechanical substeps and three samples of ten full
steps, median batch time was 0.02663 s versus 0.04492 s serial (1.69× aggregate).
Allocations were 0.418 MB versus 1.022 MB per macrostep. Every state/environment
array matched; maximum material error was 2.09e-6 and potential/work/heat error
2.39e-6. These local timings ran alongside other work and require a reserved
host confirmation, reported below. Explicit logical-slice bounds checks were added after this
measurement and checked separately; they preserve ordinary array bounds
semantics without exposing adjacent worlds.

## Reserved host confirmation

The final source-pinned hbox follow-up passed 56 assertions with seeds 12, 13, 14, 15,
including bounds checks. Both complete workloads matched every serial state and
environment array. Timings below are median milliseconds for one four-world
macrostep, including all 16 mechanical substeps and status synchronization.

| Representation | Cells/world | Grid; physical extent | Serial ms | Batch ms | Aggregate ratio |
|---|---:|---|---:|---:|---:|
| Lazy SubArray, common seed |8|8²;4×4|3.576|11.530|0.310×|
| Contiguous offset, distinct seeds |8|8²;4×4|3.987|3.837|1.039×|
| Contiguous offset, distinct seeds |64|16²;8×8|10.001|17.065|0.586×|

The candidate allocated 0.677 MB per ensemble macrostep, compared with 1.056 MB
serial for the tiny fixture and 1.092 MB for the moderate fixture. The maximum
material balance error across the final fixtures was 2.22e-6; the maximum
potential/work/heat error was 1.44e-5. These errors were also present in their
exactly matching standalone runs. This is not a refinement comparison: cell
count, physical extent and grid all differ between fixtures.

The concrete slice repaired the initial representation cost, but the moderate
CPU workload still loses. Keep the adapter experimental and retain serial CPU
execution. The remaining question is whether fewer host GPU submissions offset
view/indexing and adaptation costs on launch-bound independent GPU worlds.
That requires a separate GPU gate; the current evidence cannot answer it.
Structured measurements are in
`benchmarks/2026-10-15-cellular-world-batches.json`, with final source/logs under
`outputs/world-batches-20261015/final/` and the untouched hbox v2 source pin.

## Frozen GPU probe

`probe_device_world_batches.jl` uses the common setup/runtime pin checks and
an optional `cellular_world_batches_amdgpu.jl` adapter. It adapts the packed parent
arrays to device descriptors, then constructs borrowed `ROCDeviceArray` world
slices inside kernels. Using the native array type preserves AMDGPU’s floating
atomic dispatch and avoids routing GPU scalar accesses through the CPU slice
wrapper. This is an unverified GPU implementation until the coordinated probe
runs; importing or compiling it is not evidence of GPU parity.

The bounded probe first tests two worlds with distinct seeds, a failed-world
latch and a healthy peer. Its native 3D fixture includes contact, resource flux,
birth and death. It requires exact identity/discrete-state agreement and bounded
Float32 differences across all floating arrays against corresponding standalone
worlds. Atomic linked-list order is scratch and is checked for ownership and
coverage rather than byte equality; newly born cells enter bins at the next
ordinary geometry stage. The later moderate fixture defaults to four 64-cell
worlds on 16³ grids, extent 8 and fixed physical exchange radius 0.8, with capacity
only twice the founder count. The probe records full sample times, resource and
work/heat balances, seeds, capacity, and runtime/source pins. It does not shrink
the exchange radius or construct a capacity-sized CPU spatial stencil.

The queued command, only after GPU handoff, is:

```
DEVICE_BACKEND=amd DEVICE_AMD_OPTLEVEL=1 \
  BATCH_PROBE_RESULT=outputs/device-world-batches-probe.json \
  julia --project=PINNED_ENV scripts/probe_device_world_batches.jl
```

Status checks occur at each timed sample endpoint for both execution recipes;
packing, JIT and full audit copies are excluded. The adapter currently rejects
anabolic hooks. Future composition with constructive helpers or HIP graphs needs
its own interface/parity gate. The GPU question remains open until actual native
execution is measured.

The shared-enqueue adapter passed all 56 CPU assertions after narrowing the host
interface. A tiny CPU dry-run of the queued native driver (two four-cell worlds,
8³ and 16³ grids) also passed: all floating fields matched exactly, identity and
failure isolation held, maximum material error 2.28e-7 and work/heat/potential
error 2.66e-6. This does not compile or validate the AMD-specific borrowed-view
method. The complete frozen GPU probe source is staged at
`/othersys/tmp/intmat-probe-20261014/world-batches-20261015-v3/source/`, with its
copied pinned environment beside it; source, dry-run JSON and test log are also
retained in `outputs/world-batches-20261015/gpu-probe-source/`.
