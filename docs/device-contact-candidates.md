# Contact-search reuse: isolated device study

`studies/device_contact_candidates.jl` is an unloaded Julia/KernelAbstractions
experiment. It changes no production method or existing source pin. The current
mechanical schedule clears bins and rebuilds geometry before each of 16 force
passes and once afterward: 82 launches. Cell identities, masses, metabolites,
radii and bonds do not change inside that mechanical interval. Positions and
funds do change; biology and habitat retain their ordinary outer schedule.

The study builds a bounded per-cell candidate list once, in the initial bin-walk
order. A conservative **axis-aligned cube** includes pairs within each cell's
radius sum plus a physical skin and Float32 rounding padding. Coverage-mode bin
search is explicitly enlarged to `ri + maximum_radius + skin + padding`, with
unique periodic bins; it does not borrow the old bond-capture margin. Each cell
must remain within half the skin on every periodic axis, less rounding padding.
Validation is fused into the original position-commit kernel. No force can use
an invalid list on the next substep: invalidation is global and monotone until
the next mechanical call. Radius/identity/activity changes also invalidate it.

Two modes expose the arithmetic boundary:

- **Strict reference:** retain the original geometry rebuilds; require unchanged
  bin membership and linked-list `next` pointers as well as the skin certificate.
  The cache is then a subsequence of the exact current force walk, omitting only
  guaranteed zero-overlap pairs. Original contribution order is preserved. GPU atomic insertion order can change
  even without bin crossings, so this mode can invalidate on every rebuild; it
  is an exact-order reference, not the nominated GPU optimization.
- **Coverage mode:** reuse initial geometry through all 16 substeps, preserve the
  captured candidate order, then rebuild final geometry for ordinary bond
  proposals and lifecycle. Mathematical pair coverage and force formulas are
  preserved, but floating summation can differ from rebuilt atomic bin lists.
  This is an explicit arithmetic-order change, not a timestep change.

A list that overflows its fixed per-row capacity uses the original force helper
in strict mode, or complete deterministic all-pairs traversal in coverage mode.
An expired skin uses the same fallback. No candidates are silently truncated,
no overflowed row is treated as complete, and bonds are always evaluated in the
original bond-slot order. The all-pairs fallback can be expensive, especially
with a large unused capacity; this is a correctness fallback, not a dense-case
performance guarantee. Fallback/hit/overflow counts are recorded. Public Arena
checkpointing rejects the unloaded carrier; an experiment cannot silently label
its numerical policy as the ordinary source law.

Coverage mode has **53 mechanics launches versus 82**, removing 30 geometry
launches and adding one list build. Strict mode has 83. A capacity 8192 / 64-neighbour
cache occupies about 2.38 MiB. The source-derived force arithmetic and funding
remain unchanged: 16 substeps, original checked conservative motor cost, original
move helper, same biochemical clock. This avoids the remaining-outer-horizon
pitfall of invoking a complete mechanical funding schedule repeatedly.

## Checks and bounded CPU result

The scoped CPU suite passes **298 assertions**: strict exact snapshots, complete
contact coverage, overflow and expired-skin fallback, unchanged full-cellular
composition, schedule counts, and adversarial approaching pairs outside the old
bin stencil, including across a periodic seam. Those pairs enter contact while
remaining inside the skin; cached forces/rows match a freshly rebuilt oracle.

One approved <=180s host job used CPUs 20,21, two Julia threads, one GC/BLAS thread,
8 GiB/no swap/nice15 via `swarm-build`. Core source is commit 14f54e1 plus hashed
study-only overlays. It measured the **complete coupled cellular, catalytic and
sorptive loop**, materialized support, Float64 battery, 12³ grid and capacity 512,
16 mechanical substeps, five alternating pairs of six dt=.001 steps. Sparse has
64 lattice founders; crowded has 128, with the same initial coordinates
contracted to one quarter their distance from the box centre. Every arm starts
from the same fresh physical state; no old run is migrated.

| Case | Median paired speedup | Default median (6 steps) | Candidate median | Max cached neighbours |
|---|---:|---:|---:|---:|
| Sparse64 | 1.148× | 223.27 ms | 188.94 ms | 0 |
| Crowded128 | 0.958× | 167.69 ms | 167.59 ms | 26 |

Ratios vary substantially: sparse `[1.148,1.032,1.590,2.284,1.009]`, crowded
`[.958,1.054,.882,.766,1.082]`. This is **not evidence of a robust crowded speedup**;
median ratios and ratios of medians differ. It is not a GPU throughput result.
All ten paired terminal snapshots match bitwise, including physical/scratch
arrays, identities, bonds, counters and work. In the recorded final macrostep every active force call hits its
cache; none overflows or falls back. Earlier macrostep hit/fallback counts were
not accumulated by this driver, so this is not an all-history cache-hit claim. Closed mass error is at most 1.20e-7 and
energy error 4.78e-7. Six short steps do not establish ecological equivalence.

Full samples/errors are in
[device-contact-candidates-results.json](device-contact-candidates-results.json).
The immutable host root is
`/othersys/tmp/intmat-contact-candidates-14f54e1-20261015`; local logs and state
observations are under `outputs/device-contact-candidates-20261015`. Mac tar
included harmless, unexecuted AppleDouble sidecars; actual source bytes passed
the manifest check and loaded fingerprints record the extra sidecars. Nothing
in that pin was removed or rewritten. Future packaging should use Python/git
tar archives to avoid those metadata files.

Native GPU profiling is still needed: at the earlier materialized 128³ fixture,
mechanics was about 4.7 ms while resource/support stages dominated; even removing
all mechanics cannot give a large whole-loop multiplier there. At 64³ mechanics
was a much larger share. A GPU gate must measure the complete coupled loop,
cache/fallback counts, compiler/register cost and contact/work/chemical errors
under the declared summation-order boundary. Existing live garden/artwork pins
remain unchanged. No GPU slot was used in this study.

The CPU probe completed with exit0 and CPUs20,21 were released. All 298 scoped
assertions passed before that run; no further GPU work is authorized or running
from this lane. The existing garden remains on its immutable implementation.
