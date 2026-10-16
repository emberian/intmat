# Experimental compact brine execution

`studies/compact_brine_runtime.jl` is a separate explicit adapter for the brine
recipe's zero nutrient diffusion and zero activity relaxation. The ordinary
`World`, `ImplicitPhaseSeparation`, constructors, recordings and checkpoint
interfaces retain their meanings. This adapter exposes only carried phase,
interface, chemical forcing and FFT resources; it does not pretend to carry an
up-to-date chemical-potential diagnostic.

The complete step still uses the existing Float64 real FFT stages and scalar
spectral update, conservative harmonic-face solute transport, bounded reaction,
Float64 A/B/mineral, Float32 phase/forcing/permeability and actual Float32 dt
accumulated in the Float64 clock. The adapter keeps finite composition/forcing
validation and the same explicit solute positivity bound. The producer's
unchanged 2e-6 signed-phase and chemical guards remain requirements for any
future production integration.

The dependency audit removes or shares only fields whose complete lifetimes
permit it:

| Storage | Dependents in the original brine step | Compact representation |
| --- | --- | --- |
| World nutrient, activity, nutrient scratch (12N bytes) | Zero nutrient/activity coupling only | Absent; deposit remains F64 |
| Phase chemical potential and explicit-step scratch (8N) | Diagnostic mu and unused explicit solver | Absent |
| Solute transport scratch (8N) | A then B transport, after phase inverse | Shares F64 FFT real workspace |
| Full-grid eigenvalues (4H) | Spectral update only | Three F64 axis arrays, summed in original order then rounded to F32 |
| Interface (4N) | Actual reaction gate | Retained, using extracted original centered-gradient function |

The reference phase observer still computes the same chemical potential and
interface. Its interface formula is factored into `phase_interface_value` so the
compact observer does not duplicate the gradient law. Runtime creation delegates
to a generic internal factory; the existing typed public constructor is unchanged.

The real workspace's alias is explicit: phase forward/nonlinear/inverse fully
finish before transport overwrites it; the next phase step overwrites it before
reading. Every step rejects a different scratch binding. A single identity memo
spans world/phase/solute transfer. Checkpoints serialize only host state, pin
Julia/source/recipe, validate this alias on restore and reconstruct native plans.
No native resource is serialized. No existing long run is migrated.

For N=n³, H=(n/2+1)n², owned arrays change from 76N+36H to 48N+32H+O(n).
Adding the previously measured 32H plan buffers estimates active native storage
at 48N+64H+O(n), versus 76N+68H. This saves approximately 27.2%; it is a storage
forecast, not a GPU peak measurement. At 384³ it is 4.23 GiB, 480³ 8.25 GiB,
512³ 10.02 GiB and 576³ 14.26 GiB. 512³ remains tight against driver overhead and
unmeasured size-specific FFT workareas.576³ still exceeds the present device.
Full final four-field outputs remain 28N bytes; 512³ alone would exceed the 2-GiB
producer artifact budget, requiring an explicit changed output policy. No large
allocation or new producer is authorized by these forecasts.

The focused local two-thread CPU suite passes 264 assertions across 9³/16³/24³:
12 complete matched steps give bitwise equal phi/interface/forcing/A/B/mineral/
permeability and identical clocks. Serialized restore preserves the workspace
alias, independent buffers and exact three-step future. Owned array byte counts
match the formula; incompatible source/schema/recipe is rejected. Native gate
and whole-step performance are separate pending evidence. Replay:
`scripts/probe_compact_brine.jl` (small 16³/32³/96³ only).

The same bounded gate ran locally on CPU at 16³/32³/96³. All seven physical
field comparisons were exact before and after interleaved timing; restart futures
were exact, and relative chemical error was at most 2.39e-16. At 96³ the median
complete step changed from 66.71 ms to 59.69 ms (1.12× in this local two-thread
measurement). It excludes capture/download/file IO and makes no GPU speed claim.
Evidence: `benchmarks/2026-10-16-compact-brine-cpu.json`.

After the timing snapshot, the adapter additionally pins its included
`implicit_phase_runtime.jl` helper and enforces its supported recipe contract at
fixture creation, step and restore. Nonzero nutrient diffusion, activity
relaxation/interface feedback, mineral pinning, nonoriginal bulk, and explicit
non-F64 chemicals are rejected. Added negative source/alias/contract tests bring
the focused suite to 264 assertions. These checks do not change physical
calculations; CPU timing evidence retains the original source hashes. The queued native v2
snapshot includes these additions; the unlaunched v1 snapshot is preserved.

The isolated AMD native v2 gate passed the original early recipe at 16³/32³/96³.
All seven retained physical arrays and clocks were bitwise equal before/after
interleaved timings; immediate serialized restore and three-step futures were
exact. Chemical relative error was at most 1.42e-16 and signed phase mean drift
at most 1.58e-10. Actual forward/inverse plan buffers totaled the forecast 32H
bytes, with zero additional rocFFT workarea at these tested sizes. No large-size
workspace measurement follows from this small gate.

At 96³ its median complete physical step changed from 2.460 ms to 2.181 ms
(1.128×), excluding field download/capture/file IO. This recipe is still before
mineral production at this short time; this timing does not measure the later
384³ specimen. Evidence: `benchmarks/2026-10-16-compact-brine-native-v2.json`.
A separate funded-mineral fixture explicitly debits A/B, activates positive
reaction and tests the retained coupling under nonzero forcing/permeability.

The separate native v3 correctness fixture starts with synthetic funded mineral:
mineral=.05*min(A,B), debited equally from A and B, and explicitly selects
basal_reaction=.15 so the original reaction is active. Initial total chemical
inventory is 30485.2987 voxel-sum; mineral rises from 83.1201 to 179.3971 through
the full bounded fixture, with relative chemical error −1.19e-16. Nonzero
interface, signed forcing and heterogeneous permeability are observed. All seven
physical fields, clocks and both immediate/future restart comparisons remain
bitwise equal. This verifies active coupling and workspace reuse; it does not
change the ordinary recipe or establish later-state 384³ performance. Evidence:
`benchmarks/2026-10-16-compact-brine-native-v3.json`.

The complete-fixture CPU constructor allocates full FFT resources before upload.
At 512³ its owned arrays alone are 8.598 GB, exceeding the 8-GiB host limit before
baseline/initialization/planning temporaries. Thus the device forecast alone is
insufficient: a data-only host initializer followed by native resource creation,
or native initialization, is also needed before that resolution is executable.

## Data-only host seed and target resource construction

The additional explicit constructor pair `compact_brine_seed(c,n)` and
`compact_brine_seed_upload(seed,c,backend)` now provides that host-side path.
The seed contains only phi, interface, A, B and mineral (32N bytes), plus its
recipe and small grid/clock metadata. It runs the shared physical initializer;
it contains no Fourier spectra, real workspace or FFT plan. Upload allocates
all FFT arrays/plans directly on the requested target and creates the shared
real-workspace/transport-scratch alias there. Recipe mismatches, changed chemical/phase precision and nonzero
seed clocks are rejected; this is initialization, not long-run migration.

The new CPU suite passes 43 checks at 9³/16³/24³: exact initialization and
five-step future parity with the existing compact fixture, independent target
arrays, scratch alias, serialized seed reconstruction, recipe/time rejection
and exact array-byte accounting. The seed is 4 GiB at 512³; target GPU
headroom, large-size FFT workspace and output budgets still require a measured
preflight. Native validation of this constructor is pending. Existing compact
native v2/v3 evidence concerns the preceding complete-fixture upload path; no
new large allocation or producer has been launched.
