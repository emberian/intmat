# Local correction before harmonic multigrid

A bounded native 257³ AMD study reduced 36 global V-cycles to 30 across nine
recorded growth additions, with a 1.328× median paired harmonic-solve speedup.
This is a candidate for a fresh producer, not a change to the amber artwork.
Coupled growth parity and complete-event speedup remain unmeasured.

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
that the new occupied site and patch were the **only** fine potential/mask edits
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
full-solve synchronization and independent residual checks. A coupled archive
trajectory comparison remains the next gate before a producer opt-in.

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
