# BondNetwork scalar-loop performance

The central-force `BondNetwork` now avoids allocating endpoint slices, displacement vectors, force vectors and per-node update vectors during mechanical steps. Three-coordinate tuples retain Julia's robust `norm` calculation; coordinate loops preserve each force and Verlet operation's order. The external-work dot product still uses the original strided BLAS reduction with one reusable three-element scratch vector: a tuple dot changed rounding in a numerical probe. State layout, persistent IDs, failure law, timestep bound and energy ledgers are unchanged.

On an Apple M2 Max, Julia 1.12.4, two Julia threads, one GC thread and one BLAS thread, five interleaved samples of 20 complete mechanical steps gave:

| Actual weave snapshot | Original ms/step | Optimized ms/step | Speedup | Allocated bytes/step |
| --- | ---: | ---: | ---: | ---: |
| Initial: 1,440 nodes, 13,490 bonds | 20.97 | 3.66 | 5.73× | 40.10 MB → 1.27 MB |
| Step 6,000: 14,859 bonds, 1,369 repair bonds | 18.16 | 2.23 | 8.14× | 41.62 MB → 1.50 MB |

The timed operation includes validation, tangent bound, grip motion, two force evaluations, fracture, damping and work accounting. Cloning, warm-up and explicit garbage collection are outside timing. Existing repair bonds remain in the second snapshot; the short benchmark does not insert additional repair bonds. The machine had other active work, so raw samples are retained in [the evidence](bond-network-performance-results.json).

Both snapshots passed exact serialized-state comparisons after 100 driven steps, using copies from the same construction. This is a local equality check, not a portable serialization fingerprint. The focused tests add a 250-step bitwise trajectory with nonzero loads, fractures and newly inserted stress-free repair bonds; all 271 existing and new mechanics assertions pass.

The full, unchanged 12,000-step repair-weave study was also replayed from a new immutable source pin. **All 61 geometry files were byte-identical**, and every saved history, energy, repair-stock and fracture metric matched bit for bit. The final state has 1,632 new repair bonds and 2,201 broken bonds, including 360 broken repair bonds. This run completed in 49.42 seconds; the historical reference took 308.77 seconds under different host conditions, so the matched speedups above are the performance comparison.

Run the matched benchmark with:

```sh
JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 julia --gcthreads=1 --project \
  scripts/benchmark_bond_network.jl outputs/art/repair-weave-pilot-v2 RESULTS.json
```

The retained full replay and source pin are under `outputs/bond-network-performance-20261015/`. The original artist outputs and pins were not modified. The allocation-heavy mechanics oracle is preserved in `test/fixtures/bond_network_reference.jl`; `test/mechanics_allocations.jl` exercises it against the production implementation.
