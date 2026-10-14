# Bigruns: resolve material structure before enlarging its photograph

The repository now has executable 128³ and 256³ profiles for the dense field study and phase separation. They retain physical extents and material coefficients across resolution. The default run exports state and selected scalar fields; it skips optical construction and rendering. No full 256³ study has been launched.

```sh
# Inspect timestep, step count, memory and raw storage without allocating the big grid.
julia --project=gpu examples/bigrun.jl examples/bigruns/dense-128.toml outputs/bigruns/dense-128-study --dry-run

# A bounded pilot, with an explicit fresh output directory.
julia --project=gpu examples/bigrun.jl examples/bigruns/dense-128.toml outputs/bigruns/dense-128-pilot --steps 100

# The complete requested physical interval, when ready.
julia --project=gpu examples/bigrun.jl examples/bigruns/dense-128.toml outputs/bigruns/dense-128-study

# CPU override works for either profile.
julia --project=. examples/bigrun.jl examples/bigruns/phase-128.toml outputs/bigruns/phase-cpu-pilot --backend cpu --steps 20
```

Output directories must be new. The runner prints progress at fifteen-second wall intervals even when the next physical-time capture is distant. It streams frame metadata to JSONL and selected raw Float32 fields at the first completed step crossing each capture time. Initial, periodic, and final `record_state` observations retain effective rules and process state. These observations are inspection exports, not restart checkpoints. The runner never accumulates all captured volumes in memory. `--steps` caps work explicitly; omitted, it advances through the configured `t_end`, with at most one timestep of overshoot.

The timestep is computed as 0.8 times the minimum current process stability limit, then fixed and checked at every step. State-dependent bounds can reject it later; such a failure is reported rather than silently changing the schedule or coefficients. The dry run builds an 8³ prototype at the requested cell spacing, so it avoids allocating 256³ arrays. It reports the initial bound, not a promise about future nonlinear state. World time is Float64 while kernels and timesteps remain Float32; the recorder also reports the exact scheduled time `step × Float64(dt)`. Periodic forcing phases are reduced on the host before transfer to GPU Float32 kernels.

Run `julia --project=gpu scripts/benchmark_scaling.jl` to repeat the bounded measurements. It performs twenty warmed steps per workload/resolution, two warm-up steps, CPU comparisons through 64³, and Metal comparisons through 128³. The measured results are in `outputs/bigruns/scaling.json`; a checked-in snapshot is in [the benchmark record](../benchmarks/2026-10-14-m2max-scaling.json).

| Workload | Grid | Metal seconds/step | Initial dt | Estimated interval duration |
|---|---:|---:|---:|---:|
| Dense fields | 32³ | 0.00635 | 0.04 | t=10: 1.6 s |
| Dense fields | 64³ | 0.00789 | 0.01172 | t=10: 6.7 s |
| Dense fields | 128³ | 0.02091 | 0.002930 | t=10: 71 s |
| Phase separation | 32³ | 0.00268 | 0.002890 | t=30: 28 s |
| Phase separation | 64³ | 0.01001 | 0.0002821 | t=30: 18 min |
| Phase separation | 128³ | 0.01233 | 0.00002052 | t=30: 5.0 h |

These are actual Apple M2 Max timings of short runs, with extrapolated full-interval times clearly separated. They include kernel synchronization, copies, and timestep guards. They exclude initial compilation, file export, optics, and rendering; concurrent host work can change throughput. Dense means WaveForcing, Excitable, nutrient Diffusion and Deposition; it does not include recurrent neural state or particles. Phase means PhaseSeparation alone, including its nutrient/activity coupling; the phase profile does not enable deposition by default. The small 12³ phase CPU/Metal comparison after ten steps matched composition, interface, activity, and nutrient exactly at Float32 precision.

Explicit Cahn–Hilliard refinement is the expensive case: its fourth-order operator forces approximately `dt ∝ h⁴`. Doubling resolution makes eight times more cells and can require sixteen times more steps for the same physical time. A semiimplicit conserved solver is the appropriate next computing improvement before committing to full phase 256³ trajectories. Changing `kappa`, shrinking the physical duration, or changing the extent would change the experiment, so those adjustments require explicit new profiles. A very coarse extrapolation from 128³ predicts phase 256³ near 27 days if per-step cost grows eightfold; that is an unmeasured algorithmic estimate, not a GPU benchmark. Dense 256³ could be around 38 minutes under the analogous eightfold per-step and fourfold step-count extrapolation. Measure a capped pilot rather than accepting either prediction as capacity planning.

| Storage / explicit field allocation | 128³ | 256³ |
|---|---:|---:|
| One Float32 scalar | 8 MiB | 64 MiB |
| Dense resident world + scratch (16 scalars) | 128 MiB | 1 GiB |
| Phase resident world + process (20 scalars) | 160 MiB | 1.25 GiB |
| Dense selected capture (3 scalars) | 24 MiB | 192 MiB |
| Phase selected capture (6 scalars) | 48 MiB | 384 MiB |
| Dense full state observation (8 scalars) | 64 MiB | 512 MiB |
| Phase full state observation (12 scalars) | 96 MiB | 768 MiB |

These are counts of explicit arrays, not total process RSS or Metal-driver memory. Full observations omit world scratch arrays but include phase process scratch. Host export can temporarily add copies. Ten intervals plus initial/final observations can consume gigabytes at 256³; the runner prints a conservative raw disk budget before starting. A final optical export adds three or four scalar volumes, and rendered photographs add their own storage. Rendering uses a separate image resolution and ray-sample budget and is deliberately off in the profiles.

Refinement should hold extent, diffusivities, forcing wavelength, interface width coefficients, and initial physical modes fixed. The phase profile uses a common smooth analytic perturbation across grids; independent cell noise would inject different short wavelengths at each resolution. Dense initialization is sampled in physical coordinates. `material_mass` sums cell inventories and therefore grows with cell count for equal concentration; compare `sum(field) × cell_volume`, spatial restrictions, interface area, and spectra when comparing resolutions. If particles are later added, convert their cell-unit inventories and coordinates so that physical cargo and positions remain fixed. Merely reusing the same numerical storage values would change their relation to the refined grid.

For convergence, use a staged 48³/96³/192³ family at common physical parameters and times, with a capped pilot at each stage. Restrict the fine fields to the coarse cells before comparing deposition, composition, interface width and integrated budgets. These stages are a study design; the supplied production profiles are 128³/256³. Enlarging a photograph or interpolation of a coarse volume improves display sampling but introduces no additional resolved waves, interfaces or material dynamics.

| Process family | CPU | Metal |
|---|---|---|
| Dense field processes | Verified | Verified through 128³ |
| Phase separation | Verified | Verified parity and scaling through 128³ |
| Acoustic propagation | Verified | Verified small-grid parity |
| Recurrent neural field | Verified | Verified small-grid parity |
| Thermal writer | Verified | Portable kernels; separate validation report |
| Particles, ecological carriers, bond graph, dielectric growth | CPU implementations | Unsupported |
| GRIN observer | CPU observation | Transfer fields to host explicitly |

A four-step 16³ phase runner smoke test completed with streamed captures and all full-state exports. Its raw files occupied 1,081,344 bytes, below the printed conservative 1,671,168-byte budget. A separate four-step dense16³ Metal runner smoke also passed, including nonzero periodic forcing and exact agreement of scheduled and Float64 world time. The test used a temporary fresh directory; it did not launch either full-resolution profile.
