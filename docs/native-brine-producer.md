# A separate 384³ porous brine experiment

`studies/native_brine_producer.jl` reuses the original brine initializer, reaction
and coupled execution helper with the shared Float64 FFT backend. It creates a
new state from `experiments/art/brine-porous-384.toml`; the 128³ producer and its
artifacts remain independent. The selected porous recipe retains its fixed
physical Fourier modes, thick finite lens, material parameters and feedback.
Only resolution, timestep, step count and execution/output budgets change.

The planned 5760 steps at Float32 dt=.0025 span 14.399999678 physical time. The
explicit solute positivity number is .6144, below .9. Phase stabilization does
not remove that transport constraint. A 1100-second internal deadline reserves
room for final outputs before the enforced 1200-second launcher limit. A budget
stop records its actual completed steps/time rather than claiming the target.

The small-only preflight first exercises the actual shared 16³ phase/transport/
reaction path, full-field Float64 scalar reductions and replay sampler. The
successful report is `benchmarks/2026-10-16-native-brine-384-preflight.json`:
11.66 GB device free versus 5.33 GB active arrays/plans plus a 2.15-GB margin;
6.41 GB estimated host peak versus the enforced 8-GiB limit; 962 MB planned
artifacts versus 2 GiB. Only after this gate may the large allocation begin.

CPU initialization is transferred once, then its entire fixture is released.
The worker retains native state and explicit runtime resources. Scalar histories
reduce the full native fields; a 64³ replay is a clearly identified point sample,
not a mass estimate. Final phi, A, B and mineral remain full-resolution Float32
files in Julia column-major order. Each final field transfers/writes separately.
There are no duplicated full checkpoints or eight-field host observations.

The immutable source uses committed core `aae69fe` plus hashes of the actual
brine study, execution helper, producer, setup and recipe. Import-time phase
fingerprints are checked against the frozen files. Physical fields remain
Float32 and FFT arrays Float64/ComplexF64. This is a new resolution experiment,
not continuation or migration of an earlier state. Conservation checks concern
finite chemical inventory and signed phase composition; no thermodynamic energy
closure or biological interpretation is added.

The launcher uses the verified isolated rocFFT overlay, optlevel 1, CPU 10–11,
two Julia threads, 8 GiB/no-swap and nice 15. It enforces fresh output directories
and the artifact estimate before allocation, then checks the final file total.

The first Float32-chemical epoch (`native-brine-384-20261016-v1`) failed its
unchanged 2e-6 chemical conservation guard after time 7.8 (the next check was
8.1). Its last passing relative error was -1.9847095648e-6. Phase composition
remained stable. Its source, log and replay artifacts are preserved; the original
exception escaped before a full final capture, so this epoch has no resumable
final state. It is not a completed 14.4-time experiment.

A separate explicit recipe, `brine-porous-384-f64.toml`, selects Float64 storage
for A, B, mineral and chemical transport scratch through
`phase_runtime_fixture(...; chemical_type=Float64)`. The default remains
Float32. Initial values, coefficient values, reaction and conservative transport
formulas are shared with the original Julia study. Chemical feedback and
permeability expressions evaluate with the Float64 mineral and store into their
Float32 fields; phase and its interface remain Float32, FFT Float64/ComplexF64.
This changes arithmetic precision explicitly, not the physical equations.

For that variant, the preflight adds 16N bytes to physical/native state; native
active storage is 76N+68H, where N=n³ and H=(n/2+1)n². Full final storage is
28N bytes (phi .f32; A/B/mineral .f64), about 1.585 GB at 384³ before sampled
replays and metadata. Each frame declares its output dtypes. Sampled replay
remains an explicit Float32 observation of the original storage.

The revised producer catches numerical-guard, wall-budget and execution stops,
then attempts a streamed final diagnostic capture. `status.json` records failure
and export errors even when export fails. Only a valid planned-step completion
is labeled successful. Captures are observations, explicitly nonresumable, and
must never be treated as partial checkpoints. The original 2e-6 conservation
guard is retained for both variants.

The explicit-storage CPU suite passed 54 checks (47 existing runtime checks and
7 chemical-storage checks). The paired native 16³/32³ × 200-step gate had exact
phase parity, maximum chemical field difference 4.76e-10 and relative inventory
error at most 1.20e-16, with exact Float64 mineral checkpoint restoration and
matching clocks. Evidence: `benchmarks/2026-10-16-brine-f64-native-gate.json`.

Before the fresh large epoch, native producer tests exercised normal completion,
a wall-budget stop, and an injected failed final export. They verified successful
completion only for the complete case, nonresumable metadata in all cases,
status preservation on export failure and exact Float64/Float32 output lengths.
The final small-only F64 preflight measured 7,305,797,632 bytes estimated host
peak, 6,238,568,448 active native bytes plus a 2-GiB safety margin against
11,651,776,512 free device bytes, and 1,641,021,440 estimated artifact bytes.
Evidence: `benchmarks/2026-10-16-native-brine-384-f64-preflight.json`.

The F64-chemical memory forecast is separate from the earlier F32-chemical
backend table. Active arrays plus native FFT plans require approximately 5.81,
19.60 and 156.71 GiB at 384³, 576³ and 1152³ respectively; only 384³ fits the
present 12-GB device. Full four-field files alone require 1.59, 5.35 and 42.81 GB.
These larger forecasts involve no allocation:
`benchmarks/2026-10-16-brine-f64-memory-ceiling.json`.

The separate Float64-chemical epoch completed all 5760 steps to physical time
14.399999678 in 1036.94 seconds. Warm complete coupled steps had median 144.45 ms;
CPU initialization took 174.06 s and transfer/native runtime setup 4.22 s.
Evolution including observation/capture took 839.42 s. These are single-run
execution measurements, not a paired CPU speedup. The measured host cgroup peak
was 6.317 GB (below 8 GiB), final device free memory 5.276 GB, and reported artifact
bytes 1,639,983,706 (before its compact run report). Final chemical relative error
was -9.90e-16, signed phase mean drift -9.36e-11, and all fields passed the finite,
positivity and unchanged conservation guards. Source-pinned complete evidence:
`benchmarks/2026-10-16-native-brine-384-f64-complete.json`.

All 64 original output files (1,640,032,730 bytes including the completed run
report) were copied to private local art storage and verified against remote
SHA256 hashes. The actual eight executed source files were separately verified
against the run's fingerprints. Compact inventory:
`benchmarks/2026-10-16-native-brine-384-f64-inventory.json`.

A later exact initializer specialization is measured separately in
`docs/brine-initialization.md`. It does not alter the frozen source or fields of
the completed 384³ epoch.
