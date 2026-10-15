# Local correction before harmonic multigrid

A bounded native 257³ AMD study reduced 36 global V-cycles to 30 across nine
recorded growth additions, with a 1.328× median paired harmonic-solve speedup.
This is a candidate for a fresh producer, not a change to the amber artwork.
Coupled histories differ at finite tolerance; tighter common-state references
support the patch choice at the first divergence. Running art is unchanged.

The active amber-confluence solver already keeps its potential and multigrid
hierarchy between events. Each event adds one occupied Dirichlet cell and updates
its coarse masks. Its Laplace operator has constant coefficients: thermal dose
and nutrient affect site selection, not conductivity. A future variable-
conductivity law would need a different invalidation rule; it is not a rank-one
update merely because one conductor site changed.

The new isolated helper, `studies/electrical_local_correction.jl`, performs a
bounded red-black patch correction before the ordinary global V-cycle solve.
The patch's exterior retains the previous warm potential. It does not impose
new physical boundaries or accept a patch residual as a global certificate.

A prior certified global residual bound remains valid for every stencil outside
the patch plus its one-cell halo: no input to those stencils changed. Checking
the affected stencils and taking their maximum with the previous bound therefore
provides a global bound without an extra full-volume pass. This contract requires
that the new occupied site and patch were the **only ** fine potential/mask edits
since the previous constant-coefficient homogeneous harmonic solve. Coefficient,
RHS, boundary or unknown field edits invalidate the certificate. Without a supplied
certificate the helper checks the full residual. If the bound exceeds tolerance,
the unchanged global multigrid solver supplies the correction and convergence
check. An independent full-grid residual check is also performed after every
study event, outside the work comparison.

## Initial evidence

Two one-thread CPU fixtures each followed 32 deterministic growth events on a
forked conductor specimen. Each pair applied the same sites and separately
checked what site the candidate would have selected. The public residual target
was 2e-5; both arms solved at 0.8 times that target, matching the artwork's existing
roundoff margin. All candidate choices matched their control and all independent
residuals were below 1.6e-5. Potential arrays are not bitwise equal: their largest
observed difference was 9.60e-5, consistent with separately converged finite-
tolerance solutions. This is not evidence of indefinite trajectory equivalence
or a coupled thermal-growth run.

| Grid | Patch radius/sweeps | Baseline cycles | Candidate cycles | Added patch cell updates |
|---|---|---:|---:|---:|
|33³|4/8|132|122|183,320|
|33³|8/16|132|113|2,457,776|
|33³|12/24|132|104|10,586,712|
|65³|4/8|159|142|183,560|
|65³|8/16|159|134|2,501,088|
|65³|12/24|159|131|11,953,728|

A 65³ V-cycle performs roughly 1.7 million red-black cell updates alone, plus
restriction/prolongation/residual work. Radius 8 removes 25 such cycles while
adding about 2.5 million local updates across the entire 32-event experiment.
That supports a device experiment; it does not price dispatches, synchronizations,
whole-event thermal/frontier costs, or a native 257³ workload. No local wall-time
claim is made because the workstation was under unrelated contention.

The focused test covers untouched-stencil identity, the residual-bound argument,
global convergence, conductor/ground boundary values and fallback after an edit
without a prior certificate. The study helper assumes homogeneous harmonic RHS,
matching the current finest level; it is not a general Poisson patch smoother.

## Reproduce and next gate

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project --threads=1 test/electrical_local_correction.jl
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project --threads=1 studies/electrical_local_correction.jl 65
```

Source hashes, captured study source and raw 32-event rows are retained under
`outputs/electrical-local-correction-20261015`. The bounded native-size device experiment below includes patch dispatch cost,
full-solve synchronization and independent residual checks. The native coupled trajectory comparison and common-state accuracy follow-up
are recorded below.

## Single-workgroup device candidate

`studies/electrical_patch_kernel.jl` assigns the radius 8 patch and its halo to one
256-thread workgroup. Its 19³ potential and mask arrays plus 256 reduction floats
use 35,319 bytes of explicit local memory before compiler alignment. Threads loop
through cells, use physical-grid parity and synchronize between each red/black
color. The halo stays fixed. The final affected-stencil residual uses shared
neighbors, and unchanged global values only for neighbors outside the halo.
It is reduced in the same dispatch before writing the patch back.

The CPU KA backend matches the scalar patch bitwise in 38 assertions, covering
17³/33³ grids, centers near faces/corners and 0/1/16 sweeps. Its local residual also
matches within Float32 rounding. This checks the workgroup implementation; it does
not establish AMD performance or hardware resource availability.

`scripts/study_electrical_patch_device.jl` prepares a bounded native fixture from
the artist's exact recorded integer positions and birth events. The chosen prefix
is reconstructed with a cold converged harmonic field; this is explicitly not
the production's saved finite-tolerance iterate (no such potential capture was
made). Both arms start from that identical field and receive nine identical
recorded additions. The candidate's timed region includes its single patch
kernel, scalar readback and any ordinary global fallback solve. Both paths
synchronize at timing boundaries. Full residual, conductor and grounded boundary
checks follow every event outside timing. These recorded-site comparisons do not
replace a coupled thermal-growth trajectory check.

The future AMD command first queries the actual selected HIP device's shared-
memory and thread limits, requiring at least 36 KiB and 256 threads per workgroup.
It then checks scalar CPU versus GPU patch fields bitwise before measuring.
The native gate below verified these limits on the selected device before
compilation. A two-event 17³ CPU harness smoke also passed.


## Native AMD gate

The fixed prefix was chosen before any native timing: actual amber-confluence
births through event 4096, followed by recorded additions 4097–4105. No prefix
search or parameter tuning selected this result. The exact occupied geometry was
cold-solved once, and both arms started from that same captured potential.
The baseline was the existing `(32,4,2)` multigrid layout with fused stopping norm;
the candidate added the radius 8, 16-sweep patch and retained that solver.
Solve tolerance was 1.6e-5; independent checks required 2e-5.

On the actual selected RX 6750 XT, HIP reported 65,536 bytes shared memory and
1,024 threads per block. The kernel's 35,319 explicit bytes and 256 threads fit.
Its scalar CPU versus AMD patch field comparison passed bitwise before the native
fixture. The bounded job exited zero and released the GPU to the next queued lane.
It used CPU 10–11, two Julia threads, one GC/BLAS thread, 8 GiB/no-swap/nice 15 and a
300-second limit; there was no other coordinated GPU job during this gate.

| Measure | Baseline | Patch + multigrid |
|---|---:|---:|
| Total V-cycles, nine additions |36|30|
| Sum of complete timed solves |0.455883s|0.373053s|
| Median complete solve |0.051581s|0.037661s|

The median of paired speedups was **1.328×**; the ratio of summed solve times was
**1.222×**. Six additions saved one cycle each. Three saved none and were slower
(ratios 0.921, 0.977, 0.926). Every independent normalized residual was below 1.43e-5;
conductor and grounded values remained exact. The largest baseline/candidate
potential difference was 1.033e-4. This is convergence to the same harmonic
operator/tolerance, not bitwise trajectory preservation. Food/dose-dependent
selection was not rerun, so these results do not establish coupled growth parity
or a complete electrothermal-event speedup. No production state was migrated.

The source pin is commit `a6ec1a8` plus the explicitly hashed driver-only Julia
world-age loader fix (SHA256
`6c8cf63c722e5e7228d62d9e9797b886e9c3ac3cde756cd128d67b80ca39fd9a`).
Numerical kernel sources were unchanged. All 209 source/environment/geometry
hashes were verified before launch. The fresh remote root is
`/othersys/tmp/intmat-electrical-patch-a6ec1a8-worldage-20261015`;
its original a6 pin and an unused partial full-tree transfer remain preserved
separately. Local logs, environment, immutable fixture and source copies are under
`outputs/electrical-local-correction-20261015/native-a6ec1a8`.
Tracked raw timings/provenance are in
`benchmarks/2026-10-15-electrical-local-correction.json`.

## Optional tracked wrapper and coupled gate

`studies/electrical_temporal_solver.jl` provides `ElectricalTemporalSolver` for a
fresh producer. `electrical_temporal_mark!` updates every hierarchy mask and tracks
one known Dirichlet addition. `electrical_temporal_solve!` can patch only when the
previous global certificate is valid and that addition was the only edit. Multiple
additions force ordinary MG. A call with no intervening edits may return its
cached global bound.

Every external geometry, fine/coarse field, RHS, boundary or operator edit must
call `electrical_temporal_invalidate!(solver, reason)`. This clears the certificate
and disables local reuse. Direct writes to exposed Julia arrays cannot be detected
without a scan; the wrapper makes this ownership contract explicit. To re-enable,
`electrical_temporal_recertify!(solver; homogeneous_operator=true)` requires an
explicit assertion that the supported homogeneous constant-coefficient grounded
operator and all coarse masks have been restored. It rejects nonzero fine RHS and
performs an ordinary global solve before enabling patches. It does not implement
variable conductivity or make an incompatible operator correct by invalidation.
Eighteen focused assertions cover this state machine and its fallback paths.

A fresh pair of independent coupled histories used the actual amber-confluence
recipe at overridden 33³ resolution for 256 events. Both arms performed thermal
writing, material depletion, dose feedback and their own site selection. All 256
choices matched; occupied Jaccard was 1; final solid totals were 68.524230659 in both,
and relative material ledger errors were 3.95e-10. Independent harmonic residuals
were below 1.50e-5. V-cycles decreased 1082→908. No local CPU timing claim follows;
the workstation was contended. The narrow native heat width is underresolved on
this coarse grid, so this is a composition check, not a native physical result or
a guarantee of indefinite trajectory agreement.

The new `scripts/study_electrical_temporal_coupling.jl` retains both independent
histories after any changed choice, records the first choice divergence, overlap,
field differences, material budgets, residual checks and complete-event timings.
The first 16 events are explicitly excluded from steady timing to avoid kernel-
compilation bias; raw timings remain. No rescue or forced choice-matching occurs.
The native 257³/1024-event bounded AMD gate below has completed. Producer recipe
wiring remains an opt-in engineering choice; the active artwork is unchanged.

## Native coupled outcome: finite-tolerance history changed

The fresh `d2ae74f` gate completed all 1024 events on native 257³ in about 141 seconds,
within its 300-second/8 GiB/no-swap/nice 15/CPU 10–11 bounds. The artist had released
its pilot; recovery explicitly yielded the profiler slot; root granted the run.
The GPU was released directly to recovery when it exited zero. No running artwork
was modified.

The first changed growth choice was **event 85**. Only 109/1024 subsequent ordinal
choices matched overall, and final occupied Jaccard was **0.1184**. Both histories
continued independently, with no forced choices or rescue. Independent harmonic
residuals stayed below 1.582e-5. Relative material ledger errors were −6.70e-13 and
−1.87e-12. Final solid totals were 4019.739 and 4029.889 (about 0.25% apart), but spatial
thermal/material fields differed substantially. Their final differences compare
different occupied geometries, not errors against a common exact solution.

After the planned 16-event warmup exclusion, events 17–84 retained the same choices,
geometry and thermal/material history. Their complete-event times summed to
4.0153 seconds baseline versus 3.5008 seconds candidate (ratio 1.147× across 68 events).
The full remaining steady totals were 54.411 versus 47.836 seconds, but after event 85
these compare **different geometries** and cannot establish a trajectory-preserving
complete-event speedup. The separate fixed-recorded-geometry solve result remains
valid at its stated tolerance and scope.

This rejects a claim that local correction preserves the existing deterministic
coupled history at the current stopping tolerance. It does **not** establish that
either arm is closer to a tighter harmonic reference; both passed the same
residual criterion. The common-state follow-up below replays the common prefix
through 84 and inspects event 85 against tighter solves. Occupied overlap alone is
not an accuracy oracle. The optional wrapper remains unloaded; no production
recipe selects it.

Tracked raw rows, quality metrics and provenance are in
`benchmarks/2026-10-15-electrical-temporal-coupling-native.json`. Logs, source copies,
environment and manifest are under
`outputs/electrical-local-correction-20261015/coupled-native-d2ae74f`.
The 683,228,034-byte final electrical/material array records remain remotely at
`/othersys/tmp/intmat-electrical-coupled-d2ae74f-20261015/results-001/final-states.jls`,
SHA256 `894408677ec98f4b0524861fe786630ec83cc1691db28e3e320bb2df17d3be30`.
They are observation arrays, not an executable archive checkpoint. No new GPU
campaign is running from this lane.


## Common-state accuracy: the patch choice agrees with tighter solves

The targeted follow-up reloaded the unchanged `d2ae74f` numerical sources, with a
single guarded driver injection to capture the state immediately before event 85.
All 84 preceding baseline/candidate choices and both original event 85 choices were
reproduced exactly. Occupied geometry and temperature, dose, food and solid fields
were verified bitwise equal between the arms before either next choice. Both
potential iterates and the exact common material inputs were saved before any
reference refinement. This bounded 90-second job exited zero and handed the GPU
directly to recovery's queued profiler.

The diagnostic extracts the propensity calculation from the exact archived
`choose_archive` source, retaining its candidate order, arithmetic, nutrient gate,
dose multiplier and susceptibility envelope. It changes only the return value to
inspect weights; there is no separately handwritten growth law. Ten focused CPU
assertions cover the distribution/CDF audit. The stateless target at event 85 was
**0.7641862771655998**.

| Iterate | Selected site | Distance to nearest normalized CDF boundary |
|---|---|---:|
| Ordinary MG, original tolerance |(126,123,148)|2.936e-5|
| Patch + MG, original tolerance |(127,123,148)|2.316e-6|
| Refine ordinary iterate to 1e-6 |(127,123,148)|4.113e-5|
| Then refine it to 3e-7 |(127,123,148)|4.329e-5|
| Refine patch iterate to 1e-6 |(127,123,148)|3.849e-5|
| Then refine it to 3e-7 |(127,123,148)|4.314e-5|

The original distributions differed by total variation **8.723e-5** (about 0.009%).
That small probability change moved the target across a neighboring site's CDF
boundary. Both tighter starting paths selected the patch's original site, at both
reference tolerances. The final tight CDFs differed by only **9.81e-7**, well below
their roughly 4.3e-5 target margins. Independent residuals of those final references
were 1.589e-7 and 2.384e-7. Against them, ordinary's original potential maximum
difference was about 1.28e-4; patch's was about 6.1e-5.

At this common state, the patch iterate is closer to the tighter references and
selects their site. The ordinary tolerance-level iterate chooses differently.
This directly refutes treating deterministic history divergence as a scientific
disqualification of the patch. It is still one common-state experiment using
finite-precision references, not exact PDE truth, universal superiority, or a
proof that later histories should match. The method targets the same harmonic
operator with the same public residual tolerance; a fresh producer may reasonably
select it as a numerical recipe with explicitly different realizations. Existing
source pins must retain their own solver semantics.

Raw CDFs, intervals, refinement residuals and provenance are tracked in
`benchmarks/2026-10-15-electrical-divergence-cdf.json`. Captured source/logs are under
`outputs/electrical-local-correction-20261015/divergence-d2ae74f`. The 409,538,219-byte
common-state input remains at
`/othersys/tmp/intmat-electrical-cdf-d2ae74f-20261015/results-001/capture/divergence-input.jls`,
SHA256 `9a46a7ea204c689b1df406f64404f6ce7713bad90ae9a1b825d1f7c5ef984ebe`.
