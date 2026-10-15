# Coupled AMD launch profiles, 2026-10-15

The matched default materialized-support ecology takes **38.67ms per complete step**, versus127.83ms with optional carry-axis support: materialized is3.305× faster on the same128³ physical fixture. All physical arrays, identities, contacts and counters match exactly across descriptors and the exact-zero variants. The earlier carry-only profile below must not be generalized to the default engine.

|128³,512 active cells,capacity8192|Native ms/step|Profile/native wall|
|---|---:|---:|
|Materialized|38.672|1.034|
|Carry-axis|127.827|1.013|
|Materialized + exact-zero support guards|38.663|1.035|
|Carry-axis + exact-zero support guards|125.023|1.011|

The zero guards skip empty waste/export traversals. They show no demonstrated whole-step gain on materialized support; the axis sample improves about2.2%. These single short samples do not establish a small speedup robustly, so the guards are not nominated as a default optimization here.

For **materialized support**, per-step event intervals are: dense resource passes10.49ms, support construction5.53ms, mechanics4.69ms, sensing3.97ms, diffusion3.07ms, habitat reaction/sorption2.10ms, regulation1.94ms, metabolism/paid export1.77ms, lifecycle1.37ms, waste/porosity1.16ms and uptake gather.82ms. Both `dc_resource_exchange!` launches together dominate resource passes: their implementation atomically credits zero supplied/removed amounts to shared ledger scalars even in this closed recipe. Exact-zero atomic guards are a concrete next hypothesis requiring a matched conservation/trajectory gate. Repeated geometry recomputation explains much of the optional axis cost; it is not evidence that default support work irreducibly dominates the engine.

The matched runner is `scripts/profile_device_support_paired.jl`. It pins the same575 enqueue/carry numerical source in a fresh sibling archive, adds the source-derived exact-zero experiment, and compares physical state while excluding representation-specific support buffers. It requires exact equality, not a floating tolerance. All four cases completed within the300s controlled window; no production numerical kernels changed.

## Earlier carry-axis scaling diagnostic

Repeated per-cell physical-support traversals dominate the measured 128³ coupled ecology **using optional carry-axis support**. This is not a profile of default materialized support; axis repeatedly regenerates geometry and can be substantially slower. A matched materialized comparison is required before generalizing. Dense diffusion and Julia launch overhead are smaller targets on this fixture. This is an intrusive stage profile, not a new performance optimization.

The unchanged Julia schedule runs 122 launches per step, including 16 mechanics substeps, paid catalyst secretion, feedstock conversion and sorption. Five settling steps precede three measured steps at dt=.001. The closed physical recipe has extent4³, exchange radius.18, density128,512 founders, and capacity8192; all512 cells remain active. Both measured paths start from identical seeded states. Every compared floating array, identity/contact array and counter agrees exactly afterward.

| Carry-axis grid / capacity | Native complete step | Profile/native wall | Lifecycle event time |
|---|---:|---:|---:|
|64³ /8192|12.53ms|1.146|1.39ms|
|128³ /8192|128.28ms|1.010|1.37ms|
|128³ /32768|132.22ms|1.011|5.35ms|

At128³/C8192, event intervals per step are:

| Stage | ms | Share of profiled stream |
|---|---:|---:|
|Metabolism, including funded catalyst export|27.13|20.9%|
|Regulation and support uptake scatter|26.72|20.6%|
|Support uptake gather|14.78|11.4%|
|Waste and porosity support scatter|14.32|11.1%|
|Sensing, gradient gather and physical cues|13.76|10.6%|
|Dense resource passes|10.47|8.1%|
|Habitat support readout|7.87|6.1%|
|Mechanics|4.70|3.6%|
|Dense diffusion, both pools|3.06|2.4%|
|Dense habitat reaction/sorption|2.10|1.6%|
|Support construction, both builds|.73|.6%|

The six mixed support-consuming groups including metabolism account for about81% of the128³ stream. Their code serially traverses each cell's support, often repeatedly; axis support also regenerates radial candidates. This supports investigating exact-zero export/waste guards, concentration reuse with preserved arithmetic, then parallel per-cell support reductions/scatters. A parallel reduction changes Float32 ordering and needs explicit trajectory/conservation gates. These timings do not separate the support instructions from plasticity, metabolism or physical-cue work inside mixed kernels.

At64³, mechanics instead occupies34% of the profiled stream. Quadrupling unused capacity at128³ adds about4ms chiefly in lifecycle; it barely changes the dominant support stages. These are short settled founder states, not mature dense-contact scaling results.256³ detached parity was omitted to respect the8GiB host cap.

The observer allocates its HIP event pool before timing, records pairs around existing launches without per-launch waits, then waits once at the end. Event intervals include recording/driver gaps; the64³ overhead is material. Native and instrumented wall measurements are retained separately. These single three-step samples are diagnostic, not uncertainty estimates or a claimed speedup.

Reproduction uses `scripts/profile_device_cellular_launches.jl OUTPUT` in the existing controlled AMD environment, with `DEVICE_BACKEND=amd DEVICE_AMD_OPTLEVEL=1`, two host threads, one GC/BLAS thread, CPUs10,11,8GiB/no swap and an exclusive GPU slot. The immutable run is `/othersys/tmp/intmat-probe-20261014/device-launch-profile-575-carry`: carry source7d3caeb plus the575c14b enqueue-only model, new observer and runner. Original source archives and live runs were untouched. Julia1.12.4, AMDGPU2.8.0, ROCm6.4.3 and gfx1030 override ran on RX6750XT. Full hashes, runtime and all366 event intervals per case are preserved in [machine-readable results](../benchmarks/2026-10-15-device-launch-profile.json).

The portable observer's34 CPU assertions cover shared-schedule tracing, exact coupled trajectories, fake-event ordering/pool bounds and failed-step clock behavior. No numerical production kernel changed.
