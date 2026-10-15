# Lattice complete-step layout experiment, 2026-10-15

A contiguous shared policy bank is a plausible CPU optimization with exact tested trajectories. The earlier almost-full 64² throughput collapse did **not** reproduce, so the old result cannot establish an unavoidable density or cache scaling law. The ecological law, its original study, and their source pins remain unchanged.

The source-pinned diagnostic uses one Julia thread on Apple M2 Max, Julia 1.12.4. Each grid starts at occupancy 0.5 with the same seed and warms for 64 complete physics steps. All measured windows restart from copies of that exact warmed state. Three ABBA rounds each execute 100 full steps per arm; every round includes diffusion, neural chemistry, contact events, normal validation, and bookkeeping. Compilation and copying are outside timing. Allocation is measured in separate cloned windows. The retained profiler runs another 200 steps per arm. No GPU launch was performed.

The only accepted variant changes bank storage from `(genotype,input,output)` to `(input,output,genotype)` and changes its constructor, shape check and indexed read correspondingly. Multiplication/addition order, weights, positions, hidden state, events, funding and physics are unchanged. The 14 weights for a neuron become contiguous rather than separated by 4,096 Float32 entries. This is a layout hypothesis; these samples do not isolate cache hardware counters.

| Grid | Arm | Median voxel steps/s | Range across six windows | Allocated bytes/step |
|---|---|---:|---|---:|
| 32² | baseline | 1,112,346 | 1,042,964–1,227,602 | 50,704 |
| 32² | packed | 1,424,605 | 1,320,653–1,577,965 | 50,704 |
| 64² | baseline | 1,195,438 | 967,716–1,995,578 | 50,704 |
| 64² | packed | 1,410,887 | 967,511–1,507,319 | 50,704 |

At 32² all packed windows exceed all baseline windows; the ratio of medians is 1.28. At 64² the ratio of medians is 1.18, but ranges overlap substantially and the fastest baseline window exceeds every packed window. These are descriptive measurements under shared-host scheduling, not confidence intervals or a demonstrated universal speedup. In particular, current 64² baseline median 1.20 million voxel steps/s is over eight times the historical single 0.145-million measurement under the same nominal fixture; without controlled host conditions, that historical slowdown cannot be assigned to representation alone.

After 100 matched future steps at each size, **every dynamic array is bit-identical**, including environment, diffusion scratch, body/energy/age, hidden state, actions, compensated heat, work, genotype, founder lineage, generation, event counts, IDs and parents; ticks and world time also match. The transformed banks map back exactly. Equality is much stronger evidence here than matching final population alone. It does not prove equivalence for all possible inputs or GPU backends.

Flat sampling profiles locate the main cost in the metabolism kernel: 211/259 baseline snapshots at 32² and 718/859 at 64² include that kernel. At 64², 197 snapshots include the weight dot-product helper, versus 142/747 under the packed layout. Counts include nested/inlined frames and must not be added as exclusive time shares. Normal validation is a small share. The result points to policy evaluation/data access rather than recipient matching as the useful next optimization target. There is no evidence here that the KernelAbstractions backend failed to specialize or introduced per-voxel allocation: both layouts allocate the same approximately 50.7 KB per complete step, independent of grid size.

An additional experimental variant placed `@inbounds` directly on the dot-product loop. It failed the exact-trajectory assertion and is excluded from accepted performance results. Removing a check can alter compiler lowering; a failed equivalence check is not a performance win. A follow-up diagnosis found only Float32 differences through 100 steps (largest absolute difference 2.39×10⁻⁶, in compensated heat); identity arrays and event counts remained exact. This is small roundoff on these fixtures, but a tolerance criterion and longer trajectory check were not established before accepting that variant, and its aborted timing batch is not reported. Its source, per-array discrepancy observations and failure are retained separately under `outputs/alife-20261015/lattice-layout-bounds` for diagnosis.

A practical next change would make shared-bank layout explicit at construction and dispatch weight reads on that static layout, retaining the current bank layout as a replay option. The contiguous form should be evaluated as a CPU choice; GPU neighborhood/genotype locality may favor a different access pattern. Keep the shared finite bank: this experiment gives no reason to allocate a full genome per voxel or factor the network with a different floating-point contraction order. Any production integration should retain exact-state regressions and benchmark longer windows on an isolated designated host before advertising a throughput factor.

Reproduce the accepted diagnostic with:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. studies/lattice_performance_experiments.jl outputs/my-layout-probe
```

Raw data, archived driver, both law variants, exact final states and flat/raw profiles are retained in `outputs/alife-20261015/lattice-layout-profile-v2`. The machine result pins all loaded shared core sources. The earlier `lattice-layout-profile` directory records a path-resolution launch error and has no measured results. The optional `--bounds` diagnostic intentionally retains a strict equality assertion; its observed failure is not suppressed.

Baseline SHA256: `e6d6cab2b0624965a5e7ec66b7f9f9ce44a01f9ded5653db660de61368cde9b5`. Packed variant: `20cbb7054fe0531c4860e6f205ec078c6817a0b03603daf1f11ed412652e1528`. Accepted-run driver: `e79daddf9ff3a0e0e6fdfae706d622772c200ec3b21adc0a6fdb3fbb4a698375`. The current driver additionally archives itself, resolves its dynamic module binding through `invokelatest`, and supports the rejected bounds experiment; those driver changes do not modify either accepted law variant.
