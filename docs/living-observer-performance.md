# Living frame observation

`NativeLivingConstruction.frame` preserves its existing record format, field
precision, native ID-validated bonds, and complete physical closure audit. It
borrows CPU physical arrays only while extracting a frame, after the cellular
status check synchronizes the backend. Every returned material column, cell
position/store vector, bond endpoint, and cumulative request vector owns storage.
No source state or recorded history is changed.

Stepping and observation must be serialized by the caller. The existing batch
capture loop advances and then observes synchronously. The workbench worker
processes one command at a time; its `serve!` loop does not launch concurrent
stepping tasks. Backend synchronization completes outstanding kernels; it is not
a lock against an independently launched future simulation task.

GPU observation retains the full physical host snapshot required by the audit.
Only the learner fields used by the public frame (`learned`, `paid`, `requested`)
are transferred, rather than copying private eligibility queues and all other
learner arrays. CPU observation also reads only these learner fields.

Run the bounded CPU gate from the repository root:

```sh
JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=. --compiled-modules=existing studies/living_observer_gate.jl 16 32 64
```

The reference function is recovered from commit
`1ff4acea0673b9697c0dd9e543f0c8638c0a570b`. The gate compares complete records
and exact serialized JSON bytes, then overwrites every physical and learner
source array and checks the previously captured JSON hash. It includes nonzero
learner values and an ID-valid native bond. Timing alternates warmed reference
and optimized calls, with the origin audit computed once outside timing. It
reports both complete frame extraction and frame extraction plus JSON encoding,
including allocations and GC samples. These are static observation fixtures,
not simulation throughput or ecological results. GPU transfer behavior is covered by the separate device execution gate below.

The earlier 167–172 ms hbox observation number included an extra origin audit
inside the timed argument expression and must not be described as production
frame latency. The corrected performance lane isolates that boundary separately.

The local two-thread macOS arm64 gate passed 16³, 32³, and 64³, including the
native bond and overwrite checks. At 64³, median frame allocation fell from
171.40 MB to 67.27 MB, and frame plus JSON allocation from 213.38 MB to
109.25 MB. Median times were 14.83 → 10.14 ms for frame extraction and
43.64 → 39.18 ms including JSON encoding. These nine alternating samples ran
at nice 15 and are local observations, not isolated hbox performance claims.
Full samples are in `benchmarks/2026-10-16-living-observer-cpu.json`.

The gate hashes explicit byte vectors. Hashing a large JSON String directly in
the local Julia 1.12 SHA path repeatedly hashed its immutable CodeUnits owner
during alias checking; that probe-only overhead is excluded from all capture
and serialization measurements.

The bounded AMD 16³ gate also passed exact record types, complete record and JSON
bytes, unchanged source arrays, and post-overwrite return independence for 24
residents and one native bond. It transfers nonzero learned/paid/requested values.
Evidence: `benchmarks/2026-10-16-living-observer-amd.json`. This is a correctness
gate, with no GPU latency or simulation-performance claim. It ran on the RX 6750
XT using the pinned Julia 1.12.4 / AMDGPU 2.8.0 runtime at device opt level 1,
CPU10–11, two Julia threads, and the existing 8 GB / 300 s host controls.
The preserved separate stages are
`/othersys/tmp/intmat-probe-20261014/source-living-observer-20261016-v1/` and
`environment-living-observer-20261016-v1/`. To reproduce, pass the baseline
`native_living_construction.jl` source file and an output JSON path to
`studies/living_observer_device_gate.jl` under the normal native device launcher.
