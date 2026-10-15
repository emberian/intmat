# Exact storage specialization for coupled device ecology

A CPU1 prototype reduces the coupled grid's dense data from **268 to 168 bytes
per voxel**, retaining the existing DeviceCellular/DeviceConstructiveEcology
carriers, numerical kernels and launch schedule. It introduces storage views,
not a parallel simulation algorithm or a new process wrapper. Production source
and source-pinned runs are unchanged. Native GPU adaptation/measurement remains
a gate; the prototype deliberately rejects GPU execution.

Run `JULIA_NUM_THREADS=1 julia --project studies/run_device_memory_specialization.jl result.json`.
Code: `studies/device_memory_specialization.jl`; gates:
`test/device_memory_specialization.jl`. Retained results are in
`docs/device-memory-specialization-results.json`.

## Earned dense savings

| Storage change | Saved bytes/voxel | Exactness condition |
|---|---:|---|
| Four uniform source/sink mask matrices | 48 | Each column checked bitwise uniform; nonuniform columns remain dense. Masks are immutable. |
| Feedstock inactive third species field | 4 | Core feed contract proves third species zero; nonzero or negative-zero tail is rejected. |
| Feed diffusion scratch shares base arena | 12 | Feed diffuse/commit finishes before base diffusion starts. |
| Base scratch/demand/scale/returns/occupancy share arena | 24 | Diffusion commits before clear_grid; scale overwrites demand sitewise; occupancy clears after waste commit. |
| Uniform fixed feed permeability | 4 | Checked uniform; this role is read-only in the current schedule. |
| Habitat ledger's unused twelfth column, high/low | 8 | Initial zero verified; any nonzero write fails. Columns 1–11 remain fully stored. |
| **Total for uniform-mask fixture** | **100** | Fallback retains all nonuniform masks/permeability. |

A single V×3 Float32 backing arena supplies the five existing base transient
roles and feed diffusion. Views keep their original logical shapes; the two raw
species remain a logical three-species matrix with a checked zero tail. Every
persistent array, compensated ledger entry, field, cell and counter is retained.

This depends on the current sequential schedule. Do not run food/feed diffusion
concurrently or retain scratch values as observations. Caller must transfer
ownership to the specialized carrier rather than step both source and result:
they borrow the same persistent state. A production typed workspace should own
the backing arena and provide stage leases; overlapping active leases must fail.

## Actual complete coupled gates

All 46 assertions pass. They cover persistent-state `isequal` after every step,
checked inactive writes including signed zero, nonuniform mask fallback, page
storage, and source-pinned core checkpoint materialization/restore followed by
five further matching steps. Sixteen additional complete steps were measured
for each 16³ closed/globally-fed recipe at 64-cell capacity:

| Representation | Actual retained array backing bytes | Host graph including metadata | Exact persistent trajectory |
|---|---:|---:|---|
| Original | 1,520,372 | not sampled | reference |
| Dense specialized | 1,110,772 | 1,119,364 | yes |
| Optional paged ledgers, after 16 steps | 725,748 | 736,612 | yes |

The grid-only reduction is 37.3%; the full small carrier saves 26.9% because its
fixed-capacity cell arrays are a large share. These are actual reachable backing
arrays counted once through views and aliasing, not a formula applied to dense
logical shapes. The reference and candidate coexist during comparison; this is
not process-RSS or constructor-peak evidence. Existing fixture constructors still
allocate densely before specialization, so they cannot bootstrap an oversized
GPU run. Production needs an allocator hook/direct specialized initialization.

Scratch contents intentionally differ. Excluded base roles are scratch, demand,
scale, returns and occupancy; feed scratch is also transient. Each is overwritten
before the next read by the declared schedule. Core checkpoints still receive
ordinary materialized arrays in the existing schema; the tested restore has an
identical future persistent trajectory. This is not byte-identical serialized
scratch. A direct specialized GPU restore must allocate the new representation
from host payloads without first allocating the dense GPU state.

## Optional exact sparse ledgers: separate proof, not GPU readiness

The habitat's 11 live compensated per-voxel ledger columns plus two inert arrays
consume 96 bytes/voxel even when most are zero. An optional CPU PagedMatrix stores
only nonzero pages, returning exact zero elsewhere and preserving every history
bit. After16 steps the sampled high-ledger page coverage was 2.734% for both
closed and globally fed fixtures. Global external feed changes raw inventory at
all sites but does not by itself create catalyst/sorption ledger history there.
Long tissue spread, imported scaffold or globally distributed catalysts can
activate many or all pages; the short fixture does not bound future occupancy.

The Dict implementation requires one CPU thread. A GPU implementation needs a
bounded atomic page map/pool, safe first-writer initialization, and a latched
failure before pool exhaustion loses any history. No silent dropping, approximate
aggregation or reset is allowed. Account both the page-map and pool capacity.
Worst-case data returns to 168 bytes/voxel, plus map/allocator overhead; sparse
storage is not a worst-case memory reduction. CPU metadata can be substantial
for tiny pages and is reported separately. Future GPU page/tile size needs an
occupancy-versus-map-overhead study on real habitats.

## Current 12GB hardware budget and staged resource reduction

The documented hbox device has 11.984 GiB usable VRAM, approximately 12.87 GB.
At 512³, before cell state, head maps, allocator metadata and runtime:

| Layout | Dense grid data |
|---|---:|
| Original 268 B/voxel | 35.97 GB |
| Exact dense 168 B/voxel | 22.55 GB |
| Paged-ledger floor 72 B/voxel | 9.66 GB + occupied ledger/inert pages |
| Floor with staged-grant arena 84 B/voxel | 11.27 GB + pages/partials |

Thus the earned dense reduction is useful at smaller sizes but does not fit
512³ on this GPU. The sparse floor alone leaves limited room; there is no
512³ fit or throughput claim. A real total-residency gate must include chosen
capacity/support storage, actual page occupancy/pool reservation, context and
reduction workspace. More memory, distributed storage or a separate explicit
coarse-chemistry approximation remain legitimate alternatives.

Recovery's staged actual-grant reduction needs V×6 Float32 values, live from
exchange through its finishing reduction. Food/feed phases are sequential. Its
six-column arena can later serve diffusion and the existing transient views,
so the combined maximum is six columns, not separate 6+3 buffers. Relative to
this prototype it adds 12 bytes/voxel plus small compensated-tree partials.
Neither prototype has yet demonstrated their combined native execution. The
integration target is one typed reusable workspace shared by ordinary and
constructive-world-batch schedules, with source/checkpoint and failure gates,
not a growing family of alternate carrier wrappers.

Promotion gates also include backend adaptation of the dense storage views,
device-safe invariant failure reporting (the CPU prototype throws), bulk snapshot
materialization without forbidden scalar GPU reads, and alias-aware validation.
Constant-column specialization may increase compiler type/dispatch complexity;
benchmark compilation and complete queued steps, not only allocation. No speedup
is inferred from these memory savings. Preserve the existing constructive-batch
readout/layout semantics and independently latched world failures when sharing
workspace capacity across batched worlds.
The measured 16³ pilot uses an empty World carrier, so its world contributes no
additional dense grid arrays. Existing reference-based habitat construction can
also retain large CPU initialization fields; a high-resolution initializer must
stream or generate those roles directly within its host-memory budget. Neither
retained-state accounting nor a successful small checkpoint proves that startup
and restore peaks fit the target device/host limits.
Concurrent batched worlds require disjoint arena slices for their simultaneous
live intervals; sharing one world's storage across parallel worlds is invalid.
A single max-size arena can be reused between worlds only when the launch
schedule deliberately serializes those epochs. This is a capacity/throughput
choice, not an implicit property of the storage views.
Matched warm CPU1 timings are retained in the JSON, with the first three steps
excluded. Dense closed was 5.41→5.58ms and fed 3.98→3.90ms; paged closed was
7.22→12.37ms and fed 3.91→3.72ms. Shared-host variation and tiny fixtures prevent
a throughput conclusion; the paged case demonstrates that memory savings can
carry substantial lookup overhead. Native whole-step timing is required.
