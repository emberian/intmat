# Selective battery precision

The new AMD constructed garden explicitly selects `energy_type=Float64`. Geometry, material stores, controller inputs/outputs, physiological requests, and motion coordinates remain Float32. The default device-cell constructor still selects Float32 energy. Float64 energy is supported on CPU/AMD and rejected on Metal before allocation.

The motivation is a physical accounting error in the previous Float32 battery: a motor charge smaller than half an energy ULP could leave the battery unchanged while increasing the work ledger. A direct free-motion reproduction with the garden's timestep and substeps retained energy 1 while reporting positive motor work. Maintenance can round in the opposite direction; these errors need not cancel.

The wider battery records actual endpoint debits as work. Float64 work/heat amounts are split into two Float32 increments before entering the existing compensated ledgers. Respiration preserves the Float32 material/credit request, records the actual battery credit, and sends the remaining released energy to heat. Default growth uses a separate Float64-budget specialization while retaining Float32 material arithmetic; catalytic growth retains its existing representable-material and paid-export law. Motor affordability retains a Float32 motion factor with a bounded downward correction.

`geometry.energy_bits` declares the battery representation. Device checkpoint version 3 validates this declaration against the energy array dtype and reconstructs the matching allocator. Physical energy arrays remain full precision in checkpoints; scalar/visual observations intentionally round to Float32. Source-pinned historical runs are unchanged and require their original runtime.

This repairs the battery-resolution failure, not every possible finite-precision endpoint effect. The resource grant subtraction limitation documented in [device-resource-reduction.md](device-resource-reduction.md) remains independent. The native accounting gate below supplies the bounded evidence for this repair; longer runs must continue monitoring the same closure limits.

CPU gates passed: 50 independent assertions covering actual motor microsteps, maintenance, legacy default parity, respiration, material types, division/death and Float64 checkpoint continuation; 69 existing device assertions; 72 existing coupled habitat/audit assertions. Separate fresh-factory whole-step and 20-step coupled Float64 restart smokes passed. The native comparison and longer ledger trend use `scripts/gate_device_energy_precision.jl`, with results below.

## Native AMD gate

The source-pinned gate completed on the RX 6750 XT: core `e44b20c` with gate-only overlay `2b25091` (SHA256 `5f84771b1aafd7b80fe3b394e31fbbd1b1f5252c47d698ec2c41b12f7db94a76`). It used the actual 64³ garden, 256 founders and capacity 16,384, under the usual 8 GiB/no-swap, CPUs 10–11, two Julia threads, one GC/BLAS thread, nice 15 and AMD compiler optlevel 1 controls. Total gate elapsed time was 169.88 seconds.

With the frozen e44 dense-lifecycle baseline, after 20 warm steps three alternating 100-step samples measured **10.2802 ms/step for Float32 energy and 10.2881 ms/step for Float64 energy**, a **0.076% overhead**. These compare battery representations with the same initial recipe; they do not assert equivalent trajectories.

The Float64 trajectory then completed **10,000 additional steps**, reaching step 10,320 / time 10.32. Potential closure error was **−0.0009176** (5.46×10⁻⁷ of initial potential), and material error was **−0.0002838**. All intermediate checks retained the existing absolute/relative guard. The prior Float32 gate had reached +0.0384 potential error by time 1.1. The remaining small drift is measurable; this is not an arbitrary-precision conservation claim.

The prefix trajectory sustained **111.02 steps/s** across its 10,000 timed steps (90.075 seconds, excluding audits/checkpoints). This is trajectory throughput, not a matched prefix speedup: the short 20-step paired prefix timings were contaminated by JIT/GC and are excluded from performance claims. Earlier separately frozen prefix evidence remains independent.

Float64 native/prefix arrays matched exactly, including forced simultaneous death and slot reuse. Strict restore, continuation, reconstructed high-water mark and failure-latch checks passed; the final long-trend checkpoint restored bitwise exactly. The battery consumes 131,072 bytes, **65,536 bytes more** than Float32 at this capacity. No long garden was launched by this gate.

Compact evidence: [benchmark JSON](../benchmarks/2026-10-15-device-energy-precision.json). Raw JSON is retained under `outputs/device-energy-precision-20261015/`; the approximately 244 MB strict checkpoint remains on hbox at `/othersys/tmp/intmat-probe-20261014/device-energy-e44b20c/results/checkpoint.jls`.

After production lifecycle integration, the replay driver explicitly requests `lifecycle=:prefix` for future battery-overhead measurements. That harness change has not been used to reinterpret the frozen dense-lifecycle timing above.
