# Native inventory precision

New CPU/AMD experiments can request `inventory_type=Float64` together with
`energy_type=Float64`. Native TOML accepts `native.inventory_type = "Float64"`.
The default remains Float32. Historical runs and their strict source-pinned
checkpoints are unchanged; this is a new numerical implementation, not a resume
or a correction applied to a saved trajectory.

The separate [audited64 recipe](../experiments/native/constructed-garden-audited-64.toml)
selects both wide inventories and batteries. It bounds a fresh observation window
at64 model time or900 wall seconds, with .2-time coarse captures and16-time native
captures. It has not been launched. The runner accepts this relative cap without
an obsolete calendar deadline; an explicitly provided UTC deadline still applies
as an additional earlier bound. Use the documented hbox resource controls and
a fresh source stage/output directory.

The old garden's accounting stop exposed material creation by independently
rounded Float32 diffusion endpoints. Widening only a field would leave narrower
uptake, scatter and budget transfers. The option therefore widens cell biomass,
metabolites, soluble food/waste/feedstock, demand, allocation, grants, junction
exchange scales, waste returns, source/sink grants, reduction pairs and resource
budgets. Controllers, geometry, catalyst and bound-food storage remain Float32.
The habitat's existing conservative endpoint rules still account catalyst and
bound-food conversion residues as inert material and their potential as heat.
No guard tolerance was increased and no balancing charge was added to an audit.

A separate error was exposed in waste scattering: Float32 support weights need
not sum to one. Wide waste scatter caps each transfer by its remaining inventory
and puts the remainder on the final support edge. This closes the actual donor
transfer, including after a cell dies; its support was built before metabolism.
It does not renormalize the controller's observations. Generic source/sink
reduction now keeps both accumulator components at the selected inventory
precision. A mixed high/low prototype was rejected during review and is retained
only as diagnostic evidence.

Checkpoint version 4 records `inventory_bits`, verifies owned array types, and
restores the same precision. Existing checkpoints require their own original
source, as before. Native captures can still down-convert display fields; those
images are observations, not restart inventory.

## Bounded evidence

RX 6750 XT, Julia 1.12.4, pinned AMDGPU environment, CPUs 10–11, 8 GiB host limit,
no swap, nice 15, two Julia threads and one GC thread. The established process-local
`DEVICE_AMD_OPTLEVEL=1` workaround was needed: the first optlevel-2 attempt hit the
previously diagnosed LLVM18 EarlyCSE/MemorySSA compiler crash. No package or driver
files were patched.

The same fresh 64³ garden recipe ran 100 warm steps plus three 200-step timed
blocks in each precision. Both used Float64 batteries. Measurements include the
whole native coupled step and a final synchronization, excluding checkpoint and
image exports.

| Inventory | Median ms/step | Material residual, step 700 | Potential closure residual | Owned arrays + resource workspace |
|---|---:|---:|---:|---:|
| Float32 | 7.764 | −6.405e−6 | −1.766e−5 | 250.14 MB |
| Float64 | 8.350 | −1.137e−12 | +3.723e−12 | 275.14 MB |

Here wider accounting costs **7.56% time and 10.0% owned memory**, not an assumed
hardware FP64 peak-rate penalty. This is one short workload on a shared host.
It is not a long-run ecological result or a proof of indefinite conservation.
More accurately represented source delivery also changes the numerical trajectory;
this is not a bit-identical speed comparison.

A 1,000-step closed, reproducing/construction/sorption fixture on native AMD had
material residual zero and potential closure 5.37e−14. CPU/device identities and
numeric parity passed, checkpoint restoration was exact, and a seven-step
continuation passed a 1e−10 tolerance on biomass, metabolites and energy. That
continuation check does not assert equality of every future GPU array. CPU
seven-step continuation compared complete snapshots exactly. CPU precision regressions passed 29 assertions, including
species-wise heterogeneous diffusion, sub-Float32-ULP source delivery, three-level
budget reductions and strict precision-aware restart. A broader habitat,
sorption, resource and battery regression batch passed 14,283 assertions before
the final tree-type correction; the final precision tests and native tests cover
that correction. Float32 production defaults remain covered by that batch. A final 800-assertion
CPU batch/lifecycle gate passed, including 304 batch assertions and a new
positive-source wide-inventory serial/batch parity and strict-restart case.
Batch reduction storage and backend rehoming now follow physical inventory
precision. Wide batched execution has not yet been checked on native AMD.

## The actual stopped garden fields

The read-only diffusion fixture hash-checked the preserved final checkpoint
`171707c2a7a7e0a614e7e8d509f88a04264027591b33ca2869729d9c6a385716`,
deserialized its named arrays and copied only soluble fields and permeability.
It did not call restore, metabolism, source delivery or any ecological step, and
wrote nothing into the historical run. Both precisions began with the identical
Float32 voxel values, promoted exactly for the wide case. Four thousand native
AMD diffusion/commit steps used the recorded coefficients and timestep 0.001.

| Inventory | Precursor A residual | Precursor B residual | Waste residual | Total residual | Time |
|---|---:|---:|---:|---:|---:|
| Float32 | −7.124e−4 | −7.675e−4 | −3.933e−3 | −5.413e−3 | 1.051 s |
| Float64 | 0 | +1.421e−14 | 0 | +1.421e−14 | 1.563 s |

All fields stayed finite and nonnegative. Wide potential closure was +5.684e−14.
The Float32 loss here differs in sign from the prior single-step CPU diagnostic
and the original integrated positive residual. Rounding depends on state and
compilation, and the field changes during these 4,000 steps. This isolates a
transport defect and validates its wider representation; it does not attribute
every term of the original long trajectory's error. Transport alone cost 1.49×,
while the complete coupled step cost 1.076× in the separate native benchmark.

[Exact-field numerical evidence](../benchmarks/2026-10-15-native-final-field-diffusion.json).
The initial standalone driver missed its `KernelAbstractions` import; that failed
attempt is preserved. The successful invocation supplied the import before
including the same frozen driver. The repository driver now imports it directly.

Structured evidence: [native whole-step measurements](../benchmarks/2026-10-15-native-inventory-precision.json).
The immutable validation stage and earlier diagnostic failures are retained under
`/othersys/tmp/intmat-probe-20261014/inventory-precision{-v2,-v3,-v4}-20261015`.
Source hashes in the measurements identify the actual executed numerical code.
Raw local logs remain in `outputs/inventory-precision-20261015/`.

A longer integrated trajectory remains necessary to validate closure through the
recovery regime, dense contacts, many generations and depletion. That run has
not been launched. The new option provides a trustworthy numerical baseline for
that decision; a faster compensated representation can later be compared against
it without weakening the accounting contract.
