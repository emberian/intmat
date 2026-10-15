# Split finite-genotype bank: memory win, no measured throughput win

This unloaded CPU study factors the lattice model's 12 additive binary loci into
six low and six high bits. The low table accumulates the common Float32 base and
bits 0–5; the high table accumulates bits 6–11 from zero. Each coefficient is one
Float32 addition before the original dot-product FMA. Both tables use contiguous
`(input,output,halfcode)` storage. The full 4096×14×11 bank occupies 2,523,136 bytes;
two 64×14×11 tables occupy **78,848 bytes**, exactly 32× less.

Splitting changes Float32 addition order. It is an explicit arithmetic variant,
not an exact-law or trajectory-replay claim. The original frozen lattice source,
existing specimens and GPU runs are unchanged. The study's ordinary checkpoint
entry points reject saving/loading rather than pretending its custom bank is a
compatible baseline checkpoint. Source-captured observations remain inspectable.

The source-derived generator preserves the original RNG sequence and quantized
base/basis. An independent exhaustive Float64 additive target covers all 630,784
coefficients: original-vs-split maximum difference 4.768e-7; original error against
target 7.451e-7; split error 4.098e-7. The conservative pairwise forward-error bound
is 2.899e-6, with zero violations. It uses unit roundoff u=2^-24 and gamma(n)=nu/(1−nu),
covering both accumulation paths. An independent bounded-input controller sample
covers 2,816 neurons, including extreme codes/input signs: maximum dot difference
4.768e-7, hidden-target difference 2.384e-7 and logistic-action difference 1.192e-7.
The FMA dot bound is 1.119e-5 with zero violations. Real tanh/logistic Lipschitz
bounds exclude their implementation rounding. Inputs are explicitly bounded by
one here; this is not a proof for arbitrary damaged or out-of-domain states.
The scoped suite passes **27 assertions**, including ownership, zero-strength
identity, index bounds, closed ledgers and checkpoint rejection.

## Complete-loop controlled comparison

Each of three arms starts from the **same baseline 64-step warmed physical state**:
original strided bank, existing contiguous full-bank control, and split bank.
The loop retains diffusion, neural metabolism, paid construction, contact
matching, lifecycle events, validation and bookkeeping. Every measurement copies
that warm state outside timing; compilation is warmed separately. Three rounds
alternate arm ordering. Allocation is measured by `@timed` in each complete window.
This compares against the useful packed layout, not only the weakest baseline.

| Grid | Window steps | Baseline median | Packed median | Split median | Split speed vs baseline / packed |
|---|---:|---:|---:|---:|---:|
|32²|100|164.69 ms|148.33 ms|216.14 ms|0.762× / 0.686×|
|64²|100|338.38 ms|339.09 ms|348.77 ms|0.970× / 0.972×|
|32³|20|354.23 ms|322.84 ms|341.54 ms|1.037× / 0.945×|

Split allocates 62,800 bytes per complete step versus 55,600 for both controls.
The custom AbstractArray adapter and validation contribute representation costs;
these measurements do not isolate CPU hardware cache misses. A warmed specimen
uses a small subset of the 4,096 codes, so a contiguous full bank can already have
a small working set. **No measured size beats the packed control.** The weak 3D
baseline improvement does not establish a useful throughput optimization.

After each window, packed-control dynamic arrays match baseline bitwise. Split
integer arrays—including genotypes, lineages, generations, event counts, stable
IDs and parents—match exactly on all three cases. Split floating maxima are
8.35e-7, 2.39e-6 and 9.54e-7 respectively; the largest is compensated heat. Work,
resource, cell, hidden and action discrepancies are retained separately. The 3D
baseline final state has 32,602 inhabitants, cumulative 16,140 births, 9,150 hops
and 655 mutations: this is an active full loop, not an inert controller benchmark.
Closed material residuals stay below 2.4e-5 and accounted-energy residuals below
1.01e-4 on total inventories up to 565,601 mass / 1.90 million energy. Neither
short matched events nor small floating errors establish long-term equivalence.

The job used hbox CPUs20,21, 2 Julia threads, 1 GC/BLAS thread, 8 GiB/no swap/nice15
through `swarm-build`, with a 300-second cap. It completed exit0; the lane is
released. Core commit is df33fa743c435f5355e430494955490b7cfe017b plus hashed
study/driver overlays. Python tar preserves the captured source without macOS
sidecars. Every captured source hash was verified after completion. No GPU was
used. The immutable root is
`/othersys/tmp/intmat-lattice-split-bank-df33fa7-20261015`; local source, logs and
results are under `outputs/lattice-split-bank-20261015`.

[All measurements and errors](lattice-split-bank-results.json).
Reproduce using `scripts/study_lattice_split_bank.jl FRESH_OUTPUT` with the pinned
study and ordinary project environment. The factorization remains useful as a
memory/representation experiment; these results do not justify replacing the
packed full bank for CPU throughput. No additional variants were launched.
