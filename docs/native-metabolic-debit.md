# Conservative metabolic depletion

2026-10-16 · Astra.

A new closed 8³ action-credit fixture exposed a native full-oxidation boundary
bug before any learning update. With four nondividing founders and initial free
food concentrations (.1,.5), (.5,.1), (.2,.4), (.4,.2), the original coupled law
halted at steps 18, 17, 12, 17 respectively (`dt=.02`, habitat status64).
Anabolism received a negative precursor remainder as small as −1.3553e−20;
its strict inert/heat guard correctly refused the step. Legal signed correction
ledger columns were not the cause.

The original expression `da=oxidize*a/stock; a-=da` can overshoot `a` even when
`oxidize==stock`. A standalone Float64 example is
`a=.006728765367793111`, `b=.00251608945166251`: complete oxidation leaves
`a−(a+b)*a/(a+b) = −8.673617379884035e−19`.

The [native kernel](../src/device_cellular_kernels.jl) now selects representable
nonnegative donor endpoints. Full depletion sets both donors exactly to zero.
Partial withdrawal rounds toward the donor and caps the second species by the
remaining request. Actual endpoint deltas determine released potential,
respiration heat, waste and the oxidation ledger. Excess-storage export uses the
same helper because it contained the same arithmetic pattern. The habitat guard
and its tolerance remain unchanged. Type-generic arithmetic preserves the
Float32-only device path; this does not claim to repair the older Float32
aggregate-inventory drift described in [inventory precision](native-inventory-precision.md).

The narrow [boundary test](../test/cellular_proportional_debit.jl) covers both
precisions, zero/full/near-full/partial requests, extreme inventory ratios and a
subnormal share. The coupled CPU/AMD gate and action-credit study test composition
with actual uptake, catalysis, sorption, maintenance and mechanics. GPU evidence
is recorded separately with its exact source pins; prior material-sensing and
controller-atlas results retain their earlier pins and are not silently relabeled.

The completed [CPU/AMD gate](figures/native-action-credit-gpu.json) ran four
200-step contexts in each precision on the RX6750XT, Julia1.12.4/AMDGPU2.8.0,
with device optimization level1. All wide-inventory windows passed: maximum
absolute potential residual1.89e−14, hidden-state CPU/AMD discrepancy3.43e−7,
and energy discrepancy3.18e−9. The three action-credit sidecar kernels matched
exactly. These timings include paired CPU/GPU execution and initial compilation;
they are correctness evidence, not a throughput benchmark.

Float32 inventories retained a measurable closure failure: two of four contexts
exceeded the existing1e−4 gate (worst1.219e−4 on GPU), with comparable CPU
residuals. This is explicitly recorded as `closure_pass=false`; only wide
inventory closure gates acceptance of this study. No Float32 closure or Metal
runtime validation is claimed. An initial level2 compilation failed in LLVM
EarlyCSE on the existing junction kernel; the established level1 workaround was
used for the successful gate.
