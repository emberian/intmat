# Block-reservoir constructive ecology prototype

This experimental factorization keeps physical extent, cell radii, exchange
radius, neural/physiological rates, identity, inherited controllers, fair uptake,
paid construction and mechanical laws. Free and raw pools retain a fine transport
lattice. Cells and local catalyst/sorption reactions operate on a coarser nested
lattice. It deliberately changes spatial resolution of biochemical interaction;
it is not an exact optimization of the existing fine model.

The fine pools own mass. Each step first transports on the fine lattice with the
existing diffusion kernel. Fine body porosity is rebuilt once per step using
the existing fine axis-support kernels; the previous coarse matrix permeability
is multiplied blockwise. Block sums are rounded down into coarse representative pools. Existing
complete DeviceConstructiveEcology kernels execute with coarse diffusion disabled.
The actual change of each coarse pool is reconciled to its fine block: withdrawals
scale donors proportionally; additions distribute uniformly by volume. Fine
residual mass below a coarse Float32 ULP remains owned and is not erased. Coarse
representatives are excluded from total inventory, avoiding double counting.

This preserves nonnegativity and mass to transfer-rounding precision, while the
existing paid chemistry and physiology retain their conservative ledgers. It
retains fine within-block heterogeneity under transport and proportional removal,
but suppresses sub-block selectivity: cells in one block share its reservoir,
sensed gradients are coarse, export shapes are blockwise, and catalyst/reaction/
sorption kinetics use block-average concentrations. Nonlinear reaction of an
average is generally different from averaging fine reactions. Matrix permeability remains piecewise constant within blocks; body porosity
uses the original fine spatial kernel. These errors must be
measured against resolution and physical support radius.

The initial CPU-only prototype uses Float64 fine transfer pools and the existing
Float32 device cellular/chemistry law. That separates transfer roundoff from
spatial approximation; it does not claim GPU parity or a pure Float32 speedup.
Nested blocks use an O(1) descriptor and an index-carry iterator, with no
fine-grid index map. Fine body geometry uses O(capacity × axis candidates)
storage and one traversal of physical support per step, rather than repeated
fine support scans for every biochemical cell consumer. There are no per-cell host transfers because this prototype is CPU-KA,
not an alleged device-resident implementation.

`studies/block_reservoir_ecology.jl` reuses the shared pilot fixture and production
kernels, including the complete controller, competition, metabolism, mechanics,
lifecycle, paid secretion, catalysis and sorption stages. It constructs a fresh
private carrier with zero coarse diffusion instead of mutating an existing
process's ownership contract. No production source is changed. This experimental
carrier has no checkpoint/restore protocol yet.

## Initial complete-path measurements

Seventeen CPU assertions pass: conservative transfers (including retained
sub-Float32 residuals), three complete closed steps, implicit nested-grid cover,
and initial fine-body porosity agreement within1e-6 with the fine reference.
Short same-seed comparisons use128 founders, capacity512, unchanged physical
extent4/radius.18/rates and all complete cellular/habitat stages. These sparse, short initial-state probes do not assess lifecycle/contact-rich
ecological performance. They exercise coupled stages rather than establish
an ecological trajectory.

| Fine / interaction grid | Support storage, baseline → prototype* | Complete-step speed ratio | Food A/B relative L2 | Catalyst relative L2 |
| --- | ---: | ---: | ---: | ---: |
| 64³ / 32³ | 2.10 → .991 MB | 1.90–2.05× warm | 5.84% / 5.81% | 55.2% |
| 128³ / 64³ | 11.24 → 2.277 MB | 1.34–2.17× | 2.60% / 2.55% | 28.5% |

*Prototype support includes both coarse materialized interaction support and fine
axis support used once for body porosity. The block descriptor costs24 bytes.
Fine128³ pools use218.1 MB (Float64), fine body occupancy8.39 MB and transfer
workspace6.29 MB; support savings are not whole-process memory reductions.

The64³ comparison advances12 steps;128³ advances8 steps at dt=.001. The first64³
sample was slower than baseline due to compilation and is retained in raw data,
not presented as a warm result. The128³ target received a separate two-step
warmup, but timing still varies materially; no robust universal speed factor is
claimed from two samples. All timing includes fine transport and transfers,
not only the controller. No GPU run has occurred.

The128³ prototype material residual is-9.27e-8 and potential/work/heat residual
-2.41e-5. Food approximation improves at the finer interaction resolution, while
waste-field L2 error remains28.4%, bound-field error31.3% and maximum combined
permeability error.112. These spatial differences matter even when total mass is
nearly identical. The short runs do not establish long-term morphology or
selection equivalence. [All samples and errors](block-reservoir-results.json)
retain the cold sample, timing variation and approximation limits.

## Next useful factorization

This establishes a working conservative block-budget path with a small support
footprint. It is a candidate *new interaction law* with explicit communication
resolution, not a drop-in performance improvement. A narrower follow-up should
keep catalyst, bound pools and reaction/sorption on the fine lattice while
restricting only cell sensing and shared uptake budgets. Existing typed anabolic
hooks permit that split: prepare/readout samples coarse aggregates, paid export
reconciles to fine catalyst receivers, and the finish stage invokes the existing
fine reaction/sorption kernels. The body's fine geometry already follows this
once-per-step pattern. This avoids crediting coarse local chemistry as an exact
optimization and tests whether reduced support scans still pay after chemistry
stays fine.

A native GPU prototype would launch one independent block/species owner for
restriction and reconciliation, using the same implicit nested descriptor.
No inter-block atomics are necessary. It must retain explicit transfer residue
accounting if fine pools switch toFloat32, and must add ownership, failure-latch
and strict checkpoint contracts before generic composition. The present CPU
prototype is intentionally outside those production protocols.
