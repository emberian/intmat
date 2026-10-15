# Zero resource atomics experiment, 2026-10-15

Skipping zero-valued read-modify-write atomics reduces the measured closed coupled step from38.70ms to29.59ms (1.308× throughput,23.5% less wall time). All compared physical arrays, identities, contacts and counters match exactly. This is an opt-in study; no production kernel changed.

The source-derived `dc_resource_exchange!` variant retains every field write and same-kernel failure check. When supplied/removed is exactly zero, it atomically loads the ledger value instead of adding zero. Nonzero and NaN amounts retain the original atomic increment. The load is race-free alongside positive writers. Quota resets accumulators to positive zero before exchange, so reachable signed-zero behavior is preserved. No resource stage is skipped.

The complete coupled fixture uses default materialized support,128³ voxels over extent4³,512 active founders, capacity8192, radius.18,dt=.001 and16 mechanics substeps. Five settling steps precede three measured steps. The native wall sample excludes per-kernel profiling; separate event samples bracket the unchanged shared schedule.

| Actual source law | Baseline complete step | Guarded complete step | Resource-pass event baseline / guard |
|---|---:|---:|---:|
|Closed pools|38.701ms|29.595ms|10.472 /1.385ms|
|Finite positive raw sources|5053.257ms|4903.894ms|4892.034 /4751.572ms|

The fed variant supplies .5 total mass/time of each raw precursor, distributed over every voxel, with finite100-unit budgets. It exposes severe positive atomic contention: millions of credits target the same ledger scalar. Zero guards do not solve this. Its observed3% wall difference is not a robust speedup claim: only one short pair was collected and profiled/native wall ratios vary by about2.5%. A staged block reduction of actual grants is the next architectural candidate; it must preserve finite budgets, actual representable source/sink amounts, failures and conservation, while explicitly allowing changed Float32 summation order.

Closed profile/native wall ratios were1.035 and1.045; fed ratios were .974 and .975. Event samples are intrusive and collected separately, so these ratios can be below1 from timing variation. Closed comparisons are bitwise exact; fed primary-array differences are at most2.33e-10 across implementations and6.41e-10 between baseline profile/native samples. IDs, contacts and counters are exact.

CPU tests pass63 assertions covering closed/source/sink/mixed masks, deliberate NaN accumulators, prior failure latches and24 complete coupled closed/fed steps. A separate native AMD4096-voxel mixed-mask gate checks all four stream modes, nonnegative fields/effective compensated remaining budgets and status. GPU fields, masks and quotas match exactly; atomic-order ledger differences are at most1.94e-7. Positive atomic cases deliberately use a documented tolerance rather than promising order-independent bit parity. Independent static review found no race or reachable signed-zero change.

The bounded AMD process completed all four cases and mixed-stream gates with exit0 within300s, then released the GPU. It used CPUs10,11, two Julia threads, one GC/BLAS thread,8GiB/no swap/nice15, Julia1.12.4, AMDGPU2.8.0, process-local optlevel1, ROCm6.4.3 and gfx1030 override on RX6750XT. Source/runtime pins and every event are in [machine-readable evidence](../benchmarks/2026-10-15-device-resource-zero.json). The immutable stage is `/othersys/tmp/intmat-probe-20261014/device-resource-zero-575-carry`; prior archives and live runs were untouched.

Reproduce with `scripts/profile_device_resource_zero.jl FRESH_OUTPUT` in that pinned AMD environment under an exclusive GPU slot. The new study/carrier changes dispatch only for resource exchange and leaves production methods intact.
