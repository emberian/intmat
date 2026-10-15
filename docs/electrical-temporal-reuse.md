# Local correction before harmonic multigrid

A small CPU study shows that correcting the newly occupied site's neighborhood
can reduce the number of subsequent global multigrid cycles. This is a candidate
for a fresh producer, not a change to the active amber artwork. No GPU timings or
complete-event speedup have been measured.

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
`outputs/electrical-local-correction-20261015`. The next useful check is a separate
native-size, paired device experiment including patch dispatch cost, full solve
synchronization and independent residual checks, followed by a coupled archive
trajectory comparison. That needs the coordinated GPU window after the artist,
base-ecology and HIP-graph gates; none has been launched by this lane.

## Device candidate prepared (not yet measured)

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
This query and GPU launch remain pending the coordinated window; no hardware
limit is assumed from the GPU model name. A two-event 17³ CPU harness smoke passed.
