# Where the note goes

A fresh temporal relative of **A chord held in salt**, preserving the published still and its source. Three authored coherent pressure packets propagate through finite material. Pressure-squared exposure converts finite precursor into slower mineral; a separate constant dissolution returns mineral to precursor. The resulting interfaces change later wave propagation. The camera stays fixed while actual mineral interfaces grow, recede and reconnect; warm light is an explicitly assigned optical surrogate for instantaneous absolute pressure.

## Shared laws and implementation

`studies/wave_memory.jl` calls the production `AcousticField` directly. Its new KernelAbstractions `WaveMemory` process owns a host-support field and cumulative dose, borrowing the acoustic envelope. A single material kernel performs exponential finite-donor conversion, dissolution, dose accumulation and complementary precursor reconstruction. Both transfers use their own pre-step donor, so neither donor can be overspent. Material-plus-precursor is locally fixed. No latent or acoustic energy debit accompanies conversion: this is an exploratory, nondimensional material law, not calibrated crystallization or closed thermodynamics.

The published chord producer's scalar material loop prevented device execution. Moving that loop into a backend-generic kernel is the main structural improvement. A minimal World now holds only nutrient, deposit, activity and pressure, with no unused World scratch fields. Ten persistent Float32 state arrays occupy 40 bytes per voxel. GPU runtime and validation-reduction temporaries and peak memory have not been instrumented; this is not a claim of allocation-free complete steps.

The narrow shared change in `src/acoustics.jl` factors the existing Gaussian initializer through `Interior.seed_acoustic_field!`. The new initializer accepts host or device observations, checks shape/finiteness/timestep, preserves caller input, and uses the original half-acceleration face operator to prepare arbitrary zero-velocity pressure fields. The existing Gaussian convenience function retains its generated pressure and previous-field behavior. No wave timestep algorithm is duplicated in the study.

## Numerical evidence

`scripts/check_acoustic_initialization.jl` checks an independent Fourier half-acceleration solution, Gaussian/general route parity, device-input initialization, input preservation, and invalid shape, timestep and nonfinite inputs. CPU and Metal pass; the Fourier previous-field error is 1.4e-18 and the device/host 20-step fixture is exact.

`scripts/probe_wave_memory.jl` compares complete coupled trajectories against the scalar reference. All six CPU arrays are bit-exact after 100 steps. CPU/Metal maximum error is 3.50e-7. The reversible packet fixture additionally passes with nonnegative donors, local inventory residual 5.96e-8 and zero-wave dissolution error 3.73e-9. These small-grid tests establish portability and finite transfer, not artwork-resolution convergence.

Warm complete 64³ Metal steps take 3.05–3.45 ms across three 40-step repeats. This includes production acoustic validation/reduction, synchronized wave kernel, pressure mirror copy, material kernel and clock update. CPU2 scalar measurement is 17.86 ms; CPU2 kernel measurements vary from 28.66 to 42.73 ms on the shared host. At 48³ the ordering reversed (kernel 7.97–10.22 ms versus scalar 16.48 ms). There is no general CPU speedup claim. These are bounded shared-host observations, not isolated-machine sustained benchmarks. Initialization, export and rendering are excluded from these step timings.

## Actual-state film

The native 128³ run advances 600 steps at Float32 dt=.006 and records 41 observations including initialization. The movie uses the 40 stored observations at steps 15 through 600, eight frames per second; no temporal interpolation or moving camera. These exported observations contain mineral, fixed support, absolute pressure and pressure-squared envelope; they are not full restart checkpoints. Conversion gain is 5, threshold .018 and dissolution rate .25. Pressure has no ongoing source after initialization. The ellipsoidal support, three packet orientations/envelopes/phases and optical cut planes are authored. Subsequent fringes, cavities and changing connections are computed.

The capture takes 11.80 seconds including export; this includes device-to-host copies and duplicate f32/volume files and is not a pure simulation benchmark. All 41 captures have maximum reported local inventory error 5.96e-8. Across consecutive observations, 27,008 voxel crossings recede below mineral isovalue .3 despite net mineral growth. This is a threshold-crossing count, not a count of topology changes.

The renderer uses exactly plane-clipped mineral geometry and an amber back shell, following the reviewed chord composition. Mineral emission strength is .3 times min(abs(pressure),2), sampled trilinearly on the current surface. The fixed cutoff, light transfer, studio lights, camera and isovalue apply to every frame; there is no per-frame normalization. All visual motion comes from computed material and pressure observations.

Private outputs: `outputs/wave-memory-packets-128`, `outputs/wave-memory-film-128`. The earlier moving-point-transducer pilot is retained as a rejected direction because it produced a smooth bean rather than a rich interior. No old artwork was changed or historical run migrated.
