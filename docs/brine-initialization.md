# Exact brine initialization with concrete slab work

`brine_initialize!` resolves recipe options, conversions and Fourier modes once,
then fills independent z slabs through `Threads.@threads`. Each voxel retains
the original arithmetic and mode order; `phase_observe!` runs after all slabs
finish. The signature, initial values and chemical storage selection are
unchanged. No evolving or historical state is migrated.

The source-derived oracle loads the original initializer from commit `9c9f1a5`,
or from `BRINE_INITIALIZER_REFERENCE_FILE` for a frozen source without Git
history. It is a benchmark/test dependency, not a simulation runtime dependency.
Eight CPU cases cover Float32/Float64 chemistry, spinodal/lamellar phase and
lens/sheet reservoirs. All physical fields are bitwise identical, including
phase interface and chemical potential.

The isolated hbox gate used CPU 10–11, which are SMT siblings of one physical
core, two Julia/FFTW threads, OpenBLAS 1, 8 GiB/no swap and nice 15. Two
interleaved pairs at each resolution included the complete initializer and
`phase_observe!`, excluding allocation and FFT planning:

| Resolution | Original seconds | Specialized slab seconds | Median ratio |
|---|---|---|---|
| 96³ | 3.931, 2.148 | 0.1824, 0.1508 | 18.25× |
| 192³ | 19.142, 18.801 | 0.7373, 0.7365 | 25.75× |

Both large cases retained exact bitwise physical parity. This is a complete
initialization improvement, not a two-core scaling claim. Moving recipe work
out of the loop and into concrete threaded closures changes specialization as
well as scheduling; this gate does not separate those contributions. An earlier
serial-hoist experiment under concurrent local load was stopped and is not
performance evidence. Preserved compact results, source fingerprints and timing
boundaries are in `benchmarks/2026-10-16-brine-initialization.json`.

The first specialized factory call at each size, including allocation and FFT
planning, took 0.432 s and 1.142 s respectively, after the small warmup and
boundary tests. These are not cold process startup measurements.

The completed immutable 384³ producer used the original initializer: 174.06 s
of its 1036.94 s total wall time. Even free initialization could improve that
whole workload by at most about 1.20×. No new 384³ producer was run for this
change, and sustained coupled-step throughput is unchanged.
