# Cellular performance: larger workloads

This follow-up measures additional improvements against the already optimized
`dc30e40` implementation, rather than reusing the original slower baseline.
Complete `arena_step!` calls include resource transport, runtime checks, recurrent
control, metabolism, contacts and lifecycle. Raw samples, exact source maps and
replay commands are in [cellular-performance-pass2-results.json](cellular-performance-pass2-results.json).

The bounded hbox probes used CPU affinity 8,9 (one physical i9-12900 core), two
Julia threads, one GC thread and one BLAS thread, with the existing low-priority
`swarm-build` controls, 8 GiB memory ceiling and no swap. Existing soaks on CPUs
2–7 were left pinned and unchanged. These are CPU results, not GPU throughput.

Three identical-input comparisons cover 1,024 initial cells, a real 1,489-cell
birth-heavy continuation with 1,199 bonds, and 4,096 initial cells at the same
physical density and spacing in a larger domain. The 4,096-cell recipe inherits
the 4,096 capacity guard: it measures larger transport/contact work, not unlimited
lifecycle growth. The separate birth-heavy case exercises allocations during
actual reproduction. Each timing is the median of five warmed eight-step samples;
state cloning and explicit GC occur outside timing.

## Complete-step results

| Workload | Baseline ms/step | Final ms/step | Throughput ratio | Fewer allocated bytes |
| --- | ---: | ---: | ---: | ---: |
| 1,024 initial cells | 13.240 | 11.107 | 1.19× | 38.9% |
| 1,489 birth-heavy cells | 27.561 | 25.606 | 1.08× | 23.4% |
| 4,096 initial cells | 59.989 | 47.064 | 1.27× | 43.1% |

Final per-step allocation is 2.61, 15.91 and 10.41 MB respectively. The larger
case benefits most from reusing conservative contact candidates. The birth-heavy
case remains more expensive and improves less; these results do not establish a
uniform speedup across population structures or host loads.

All three final 100-step continuations pass the exact physical-state and alias
comparison described below. No live simulation was restarted or source-migrated.

## Changes

Mechanics now allocates force, velocity, motor, trial, weight and spending scratch
once per call, rather than on every adaptive substep or rejected budget attempt.
The buffers are transient locals, so checkpoint layouts and process ownership do
not change. This alone reduces allocations substantially but produces little
wall-time improvement; it is not presented as a speedup by itself.

CPU Float32 diffusion uses direct periodic neighbor indices with the same ordered
face contributions and arithmetic as the generic axis/side implementation.
Bitwise reference tests cover 2D, 3D, singleton axes and aliasing supported by the
generic method. Other element/storage types retain the generic path.

Resource exchange caches the four ledger scalars during a serial field sweep.
Every source/sink rounding and ledger update retains its original order. If any
ledger vectors share memory, the original interleaved-storage loop remains active.
The public field layout, ledgers and reservoir API are unchanged.

Contact candidates may now survive multiple adaptive substeps. The existing
`4travel_limit + bond_range` skin is retained. Before a current, proposed motor or
accepted trial position is evaluated, the cache is reanchored when any body's
measured displacement from its reference exceeds `travel_limit`. This conservative
threshold is half the per-body skin allowance. If a trial triggers a reanchor at
the current positions, the established motor/passive bounds still place it within
`2travel_limit` per body. The original exact distance filter and canonical pair
order remain unchanged. `reuse_candidates=false` supplies the always-rebuilt
reference mode; 160 exact mechanics assertions exercise seam crossings, 2D/3D,
funded/unfunded motors, bond identities and energy diagnostics.

## State comparison

This pass exposed a pre-existing uninitialized scratch tail:
`rebuild_cell_stencil!` resized `stencil.proposals` without initializing newly grown
slots. Resource allocation overwrites every proposal before reading it, so these
incidental bytes were not physical state but could enter checkpoints. Rebuild now
zeros that scratch explicitly; tests poison it before allocation and verify
unchanged grants, fields and ledgers.

The replayable [comparator](../scripts/compare_cellular_states.jl) deserializes
copies and zeros only this documented scratch vector. Its type, length and alias
relationships are still checked. Every other persisted field, Float32/Float64
bit pattern, dictionary value and two-way mutable alias relationship must match
after 100 steps. There is no numerical tolerance or broad scratch exclusion.
