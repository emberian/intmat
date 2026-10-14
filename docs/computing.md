# Computing the material

The dense fields use Float32 three-dimensional arrays and KernelAbstractions kernels. A grid with `nz=1` has two-dimensional transport. Boundaries are periodic. Explicit diffusion rejects steps beyond its conservative CFL limit; excitable dynamics also cap the step at 0.05 nondimensional units. Tuple order is the coupling schedule, and every step preflights its process limits before mutating fields.

The CPU environment is the root Julia project. The optional Metal environment is separate, so ordinary CPU use does not load GPU dependencies:

```sh
julia --project=gpu -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=gpu scripts/benchmark_gpu.jl
```

The benchmark verifies all eight fields against CPU on a 12³ grid, then times 50 already-compiled steps on a 64³ grid. It writes measured errors and timings to `gpu/verification.json`. Kernel synchronization is included. Compilation and host initialization are excluded from the timings. This is an end-to-end dense-field comparison with four processes, including array copies; it is not a peak hardware throughput claim.

Construct a GPU world with `World(Grid(64,64,64); backend=MetalBackend())` after `using Metal`. Dense excitable activity, nutrient diffusion, deposition, and imposed wave forcing share their CPU and GPU implementations. Particle exchange currently requires CPU Arrays: it does not execute on Metal. Exporting fields or rendering on the host also requires copying with `Array(field)`.

The imposed pressure and velocity are analytic external forcing, with zero vertical velocity; they do not solve a fluid or elastic wave equation. Excitable activity is an invented active rule. Nutrient and deposit are conserved inventories per cell; particle storage joins that budget on CPU. Activity itself is not a conserved concentration and may be negative during recovery.

## Measured on this host

On 14 October 2026, Julia 1.12.4 and Metal 1.11.1 successfully executed all four dense-field processes on the Apple M2 Max (38 GPU cores). After 21 steps on 12³, the largest CPU/Metal field discrepancy was 1.49×10⁻⁸; nutrient and deposit matched exactly. Fifty warmed steps at 64³ took 2.231 seconds on CPU and 0.909 seconds on Metal, a 2.45× speedup for this workload. These timings include the state-dependent excitation stability guard (GPU extrema reductions), synchronization, and scratch-array copies and are a single bounded comparison. Initial dependency precompilation took approximately three minutes and is excluded. Full machine-readable measurements are in `gpu/verification.json`.

A separate recurrent neural-medium CPU/Metal comparison also passed, with maximum activity discrepancy 7.45×10⁻⁹ in its tested short run. The reported benchmark above contains the four dense-field processes, rather than the neural medium or particles.
