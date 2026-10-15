# Physical exchange support at high resolution

Fixed physical exchange radius must cover more grid sites as resolution grows.
The normalized kernel is `(1 - distance²/radius²)²` inside its periodic sphere,
zero outside, with the existing nearest-center fallback for unresolved supports.
No radius was reduced with resolution in this study. All implementation and
measurements use Julia, on one CPU thread; no GPU benchmark or full-arena speedup
is claimed.

## Work and memory

For voxel spacing `h`, a resolved sphere has approximately
`4π/3 * (R/h)^3` live sites per cell. The materialized DeviceCellular strategy reserves the bounding box
`S = product(min(dim, 2ceil(R/h)+2))` for every capacity slot, including inactive
ones. `sites::Int32` and `support::Float32` cost `8 * capacity * S` bytes.
The CPU CellStencil stores sites, Float64 weights, owners, linked edges and
proposal scratch: approximately 40 bytes per live edge, plus a site dictionary
and site-order vector. Vector capacity and dictionary slack add overhead.

At physical extent 4 and radius 0.18:

| Grid | Bounding-box slots/cell | Expected live sites/cell | Device support at capacity 20k | Float32 axis descriptor at capacity 20k |
| --- | ---: | ---: | ---: | ---: |
| 64³ | 512 | 100 | 81.92 MB | 4.08 MB |
| 128³ | 2744 | 800 | 439.04 MB | 6.96 MB |
| 256³ | 17576 | 6404 | 2812.16 MB | 12.72 MB |

The current device grid arrays add approximately 84 bytes per voxel: 1409 MB at
256³, before constructive fields, genomes, host initialization, or captures. With
extent 16 and radius 1, 256³ support has 39304 reserved slots per capacity cell:
6289 MB at capacity 20k, versus about 16.56 MB for Float32 axis geometry. Capacity
40k doubles capacity-dependent storage even when only 20k cells are alive.

There is also an initialization trap. The pre-`a3db6c2` DeviceCellular constructor accepted a complete
CPU CellularEcology, which builds its CSR support, and then copied every device
state array to host with `map(Array, state)` before initialization. A high-resolution
launch can exceed the 8 GiB host budget before useful GPU execution. Optional
axis support needs a fresh initialization payload that omits unused CPU CSR and
uploads only the small physiological/genome/bond arrays, rather than merely
removing the device support arrays.

## Isolated representation prototype

`studies/physical_support_prototype.jl` caches the one-dimensional periodic site
coordinates and squared offsets for each axis, plus a normalization and fallback
per cell. It regenerates the same sphere and normalized weights while iterating.
Storage is `O(cells * sum(axis candidates))`; computation still visits the
physical support. The spherical cutoff prevents a simple globally separable
three-dimensional kernel. This is not a whole-process replacement or a new
checkpoint format.

The helper supplies a scalar support visitor, fused concentration/gradient
sampling, and a synchronous proportional allocation prototype. The allocator
retains the original overflow-resistant largest/scaled proposal calculation,
rounds Float32 withdrawals toward each request cap, and credits actual field
losses to grants. It uses three Float64 voxel buffers plus a unique active-site
list, avoiding per-edge proposal/owner/link storage. Its workspace can approach
`28 * voxel_count` logical bytes before vector-capacity slack, so geometry savings must not be described as savings for
the entire allocator.

## Measurements

The fixed extent-4, radius-0.18 fixture uses nonuniform concentration fields and
seeded positions. At 2000 cells:

| Grid | Existing stencil logical bytes | Axis geometry bytes | Existing / axis build time | Existing / axis sample time |
| --- | ---: | ---: | ---: | ---: |
| 64³ | 16.91 MB | 0.608 MB | 0.0408 / 0.00475 s | 0.0254 / 0.0614 s |
| 128³ | 139.83 MB | 1.040 MB | 0.559 / 0.0246 s | 0.220 / 0.163 s |
| 256³ | 1126.03 MB | 1.904 MB | 3.967 / 0.142 s | 1.891 / 1.029 s |

Sampling includes concentrations and all central-difference gradient channels
for three species. The prototype fuses the two sampling passes and regenerates
geometry; the comparison does not isolate compression from loop fusion. Sampling
allocated zero bytes after warmup. Maximum concentration and gradient differences
were 2.94e-15 and 5.14e-16. These are shared-host CPU timings with recorded sample
variation; the 64³ result shows that the new representation is not universally
faster.

Complete fair allocation at 2000 cells gave median paired speedups of 0.81×,
0.93× and 1.31× at 64³, 128³ and 256³. The largest grant difference was 2.52e-17,
field difference 6.78e-21, and grant/withdrawal ledger discrepancy 2.49e-14. At
256³ the CPU prototype needs up to approximately 470 MB of allocation workspace
in addition to 1.9 MB of geometry, rather than the 1.13 GB reference stencil.

A real 20k-cell 256³ descriptor held 128.08 million live support edges implicitly
in 19.04 MB of Float64 geometry. It built in 2.484 s and sampled in a median
16.12 s. The 20k reference was not materialized. This confirms bounded storage,
while exposing the remaining physical-support work; it is not a GPU throughput
measurement.

## Validation and integration boundary

355 tests passed: CPU sampling/gradient/allocation comparisons across periodic,
large-radius and unresolved cases, plus bitwise sites/order/weights comparisons
against the existing Float32 KA `dc_support!` on 99 cells, including 64³, 128³ and
256³ grids at unchanged physical radius, plus four input-failure checks before exchange. Float64 geometry is mathematically
identical to the reference's sqrt-based expression, with tiny roundoff changes;
no identical-bit claim is made for CPU Float64 allocation.

This isolated prototype leaves the CPU CSR unchanged. The subsequent optional
device strategy, integrated at `a3db6c2`, uses capacity-major axis arrays `(capacity, axis_candidate, axis)`
and a typed scalar iterator shared by the base and constructive consumers. It
must preserve the normalizer, fallback, periodic images, edge ordering, and
pre-/post-movement support schedule. Reuse concentration data only while resource
fields and support are unchanged. Current sensing and regulation use different
Float32 arithmetic groupings; an exact cache can compute both accumulators during
one field read rather than silently changing arithmetic.

CPU CellularEcology has a concrete `stencil::CellStencil` field and several direct
CSR consumers (porosity, observation, habitat secretion and anabolic sensing).
A type-parametric replacement would change its Serialization ABI. That broader
change is deferred: preserve old readers/source pins and checkpoint claims.
The integrated device strategy uses an explicit option, source-pinned checkpoint
format 2, and direct founder initialization without the CPU CSR. Its 69 CPU-KA
checks pass; AMD measurements are a separate pending gate. See
[device-cellular.md](device-cellular.md).

Full settings, forecasts, sample distributions and source hashes are in
`benchmarks/2026-10-15-physical-support.json`. Reproduce the bounded CPU work with:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. --threads=1 test/physical_support_prototype.jl
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. --threads=1 scripts/benchmark_physical_support.jl FRESH_OUTPUT
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. --threads=1 scripts/benchmark_support_allocation.jl FRESH_OUTPUT
```
