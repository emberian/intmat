# Empty-material chemistry in complete living steps

The exact empty-material shortcut raises measured complete-step CPU throughput
by 7.4–7.6% in the two tested native 64³ living scenes. AMD throughput is unchanged
within the measurement variation. The corrected all-voxel catalyst controls stay
at parity. This is a bounded same-resource result, not a controller-only timing
or a claim about all habitats.

## Law and representation

If both catalyst channels and both bound inventories are zero, reaction, decay,
binding and collapse are the identity. `dh_react_sorb!` now takes that branch
explicitly. It still resets both release channels and performs the same
compensated-zero normalizations in grid-ledger channels 1–4 and 9–11 and the
inert pair. It avoids sorption arithmetic and soluble-field writes; the safety
guard still reads all five relevant soluble inventories. Diffusion, finite
sources, contacts, division, learning and validation remain in the timed step.

The shortcut requires finite soluble inventories with nonnegative sign bits.
Nonfinite values and negative zero retain the original endpoint path, including
its diagnostics and zero canonicalization. Public cellular steps reject
nonfinite/nonpositive timesteps, and habitat parameter validation requires a
positive half-saturation constant. No state schema or physical parameter changed.

## Isolated complete-step gate

`studies/living_step_native_gate.jl` compares only the reaction kernel against
its original implementation at Git revision
`9c9f1a51f277f66949da00fdb63c7a20c2fcaed5`; every other kernel and schedule is
shared. `LIVING_BASELINE_KERNEL` supplies an immutable extracted baseline on
non-Git frozen stages. The underlying general harness is
`studies/living_step_performance.jl`.

All cases use physical extent 2³, capacity 512, materialized exchange support,
dt .02, finite patch-recruitment feeding, native motion/contact mechanics and
paid inheritance/learning. The 46-founder fixture increases contact population
while respecting the same conservative funded-capacity bound. Warmup reaches
4.4 model time, beyond twice the two-time credit horizon, and the native gate
asserts actual paid actor updates before timing. There are five alternating
baseline/candidate pairs of 30 complete steps from matched independent states.
Clone-triggered collection occurs before timing. Timed blocks contain no
compilation and include native status synchronization. Transfers, observation,
cloning, warmup and startup are separate.

The assigned hbox CPUs 10–11 are **two SMT threads of one physical i9-12900
P-core**, not two physical cores. Both sides use that same pin, two Julia
threads, OpenBLAS 1, enforced 8 GiB/no swap and actual nice 15. AMD uses the
verified standalone ROCm 6.4.3 environment, RX 6750 XT, gfx1030 compatibility
configuration and the established opt-level-1 compiler configuration. The
completed brine producer and its output transfer released before these timings.

| Backend | Fixture | Residents after block | Baseline ms/step | Candidate ms/step | Speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| CPU, one core/two SMT threads | 24 founders, native material | 33 | 93.944 | 87.448 | 1.0743× |
| CPU, one core/two SMT threads | 46 founders, native material | 59 | 99.520 | 92.485 | 1.0761× |
| CPU, one core/two SMT threads | 24 founders, full catalyst coverage | 33 | 111.509 | 111.416 | 1.0008× |
| AMD | 24 founders, native material | 33 | 21.091 | 21.039 | 1.0025× |
| AMD | 46 founders, native material | 59 | 23.312 | 23.299 | 1.0005× |
| AMD | 24 founders, full catalyst coverage | 33 | 21.192 | 21.224 | 0.9985× |

Every CPU native paired block improves: ratios are 1.071–1.090 for the sparse
scene and 1.053–1.093 for the denser scene. The corrected CPU control ranges
0.991–1.007. Native CPU cases execute 152 and 293 paid actor updates and nine
and thirteen births, respectively. Thus the CPU gain resolves above this
control's variation. The GPU differences are too small to claim a useful gain.

All mutable scientific, private-learning, event, scratch, status, identity and
contact arrays, resource-reduction workspace, lifecycle watermark and clocks
match exactly: **416,329,524 compared array elements per backend**, or
832,659,048 across the six valid cases. The AMD comparator permits small
floating discrepancies for device scheduling, but none occurred in these
comparisons; all reported error maps are empty.

Reports are `benchmarks/2026-10-16-living-step-{cpu,amd}-isolated.json` for native
cases and `benchmarks/2026-10-16-living-step-{cpu,amd}-control.json` for the
corrected controls. Their runtime/project/manifest hashes, all paired samples,
activity counters and actual resource metadata are retained.

## Boundary checks and measurement corrections

`studies/living_chemistry_edge_gate.jl` compares the original and candidate
kernels directly on CPU and AMD: 14 cases and 140 exact array/status checks on
each backend. Cases cover NaN, Inf and negative zero in raw, food and waste;
negative-zero bound inventory and catalyst bits; positive bound/catalyst fallback;
noncanonical nonzero compensated low ledger/inert values; and nonempty release
buffers. The empty case explicitly proves that normalization occurs and release
7 resets to zero. Reports are
`benchmarks/2026-10-16-living-chemistry-edge-{cpu,amd}.json`.

The first artificial full-coverage fixture overwrote existing catalyst with
1e-7. Strict AMD checkpoint creation correctly rejected the resulting bound
inventory above its reduced capacity. The repaired fixture uses
`max(existing_catalyst,1e-7)`, preserving existing material and capacity. Only the
corrected control was repeated on each backend; native cases were not rerun.
The original CPU control is explicitly marked superseded, and the failed AMD
phase and its cost remain recorded.

The v1 observer timing also included an extra origin-audit argument that
production computes once outside capture. Its 167–172 ms values are not the
production frame cost. With origin precomputed, the corrected control's warmed
frame costs 34.73 ms on CPU and 44.27 ms on AMD; complete-array copy probes cost
14.91 and 44.60 ms, respectively. These are separate single observer samples,
not a serialization benchmark or part of the simulation speedup.

All bounded hbox phases together consumed 441.88 wall seconds and 661.03 process
CPU seconds, including startup, compilation, the superseded fixture, corrected
controls, parity, transfers and the AMD edge gate. Peak RSS was 1,760,064 KiB
(1.68 GiB), all phases had zero swaps, and each used the same 8 GiB/no-swap
resource limit. GPU VRAM peak was not instrumented; the peak above is host RSS.
Phase costs and exits are stored in their reports. No worker
from this lane remains active; CPU10–11 and GPU were explicitly released to
the initializer lane.

## Source history and earlier local evidence

The exact pre-comment-change performance/64³-scene kernel has SHA256
`856670a8279e22c456012c8984036a3c6f0c981c95ca686be9d29653f0ecf39a`.
Its private local executed copy is
`outputs/provenance/living-step-gate-20261016-v1/src/device_habitat_kernels.jl`.
The immutable hbox source/environment locators are
`/othersys/tmp/intmat-probe-20261014/source-living-step-gate-20261016-v1/` and
`environment-living-step-gate-20261016-v1/`; repaired control uses the separate
`-v2` siblings. Existing scene histories keep their executed source.

The final kernel hash is
`b9774fc4b57b8d6dd844dc996979b007cf4f3e8160660a2961525514dad52c9c`.
Only the traffic comment changed after the complete-step gates; CPU and AMD
edge gates cover this final source. Its isolated AMD edge stage is
`/othersys/tmp/intmat-probe-20261014/source-living-chemistry-edge-20261016-v1/`.

Earlier contended local timings are explicitly inconclusive in
`benchmarks/2026-10-16-living-step-cpu-local.json`. Those initial twelve
16³/32³/64³ array gates passed but warmed only 2.2 model time, before the first
possible paid actor update at twice the horizon. They do not establish paid-actor
coverage. The isolated gate above corrects that boundary. The original local
profile identified chemistry as the largest named kernel (3255 sampled stack
frames versus diffusion 1389 and habitat validation 956); these counts guided
the candidate and are not additive wall-time percentages.
