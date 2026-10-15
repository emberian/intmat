# Selective battery precision

The new AMD constructed garden explicitly selects `energy_type=Float64`. Geometry, material stores, controller inputs/outputs, physiological requests, and motion coordinates remain Float32. The default device-cell constructor still selects Float32 energy. Float64 energy is supported on CPU/AMD and rejected on Metal before allocation.

The motivation is a physical accounting error in the previous Float32 battery: a motor charge smaller than half an energy ULP could leave the battery unchanged while increasing the work ledger. A direct free-motion reproduction with the garden's timestep and substeps retained energy 1 while reporting positive motor work. Maintenance can round in the opposite direction; these errors need not cancel.

The wider battery records actual endpoint debits as work. Float64 work/heat amounts are split into two Float32 increments before entering the existing compensated ledgers. Respiration preserves the Float32 material/credit request, records the actual battery credit, and sends the remaining released energy to heat. Default growth uses a separate Float64-budget specialization while retaining Float32 material arithmetic; catalytic growth retains its existing representable-material and paid-export law. Motor affordability retains a Float32 motion factor with a bounded downward correction.

`geometry.energy_bits` declares the battery representation. Device checkpoint version 3 validates this declaration against the energy array dtype and reconstructs the matching allocator. Physical energy arrays remain full precision in checkpoints; scalar/visual observations intentionally round to Float32. Source-pinned historical runs are unchanged and require their original runtime.

This repairs the battery-resolution failure, not every possible finite-precision endpoint effect. The resource grant subtraction limitation documented in [device-resource-reduction.md](device-resource-reduction.md) remains independent. A fresh native sustained accounting gate is required before launching the longer garden.

CPU gates passed: 50 independent assertions covering actual motor microsteps, maintenance, legacy default parity, respiration, material types, division/death and Float64 checkpoint continuation; 69 existing device assertions; 72 existing coupled habitat/audit assertions. Separate fresh-factory whole-step and 20-step coupled Float64 restart smokes passed. The native comparison and longer ledger trend are prepared in `scripts/gate_device_energy_precision.jl`; they are not yet measured.
