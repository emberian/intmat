# The warm thread remembers

A new optical study for the older Acoustic amber specimen: two thermal writing
histories share a finite mineral precursor; sound propagates through the medium
their intersecting histories have changed. The photograph shows actual Julia
fields inside a deliberately staged amber observation chamber.

One source traverses a continuous 3:7 torus knot. The other follows a rotated 2:5
knot after a delayed start, with a different conversion strength. These paths are
designed inputs, not emergent trajectories. The study supplies positions to the
existing Julia `write_heat_kernel!`; its thermal diffusion, cooling, dose and
conservative conversion laws are unchanged. The writers retain separate
thermal histories but consume the same `World` nutrient and deposit arrays.
`AcousticField` reads that deposit every step and reduces local squared wave
speed using its face coefficient. The scalar wave uses the production periodic
grid; the glass chamber is optical staging. Sound does not change the writers here.

The refined run uses native 384³ cells in extent 2³, 2,500 steps at dt = .0032,
and captures at approximately t = 2.6656, 5.3344 and 8. Source width .026 resolves
Gaussian sigma by 4.992 cells (FWHM 11.75 cells). Maximum source moves per step
are .010404 and .007281, both below half a source sigma. At speed .37 and
frequency 1.7, the unmodified wavelength spans 41.79 cells; even saturated
mineral leaves approximately 10.8 cells per wavelength. The acoustic timestep
limit is .007314 and the thermal limit is .026055. These are model units and a
finite resolution artwork, not calibrated laser physics or a convergence study.

Measured evolution took 521.71 seconds, with a median 197.67 milliseconds per
complete composed step on local Metal. At t = 8 the last source has been off
for .75 model time: maximum temperature is .00524, while 1,621,249 mineral
cell-equivalents remain. Shared inventory suppresses 88,761.96 cell-equivalents
of independent double conversion. The direct composition identity agrees to
1.61e−6; material relative error is −1.28e−10. These are dense composition
timings, not arena timings or an isolated throughput benchmark.

The exact new source-position wrapper and production thermal/acoustic
composition passed a 150-step, 18³ CPU/Metal fixture: maximum pressure error
5.54e−6, material error 1.19e−7 and temperature error 5.96e−8. Broader, faster
sources exercise both paths in the small grid. A separate 24³ causal fixture
checks shared mineral against `1-(1-d₁)(1-d₂)` to 8.94e−7; independent inventories
would claim 373.84 extra cell-equivalents. Removing only mineral-to-speed
feedback changes an identical seeded pulse by a maximum pressure of 1.066.
These fixtures test composition and portability, not artwork convergence.

Optics use native deposit levels .99, .65 and .11 for a copper-gold dense heart,
transparent gel and a thin-film sheath. The open chamber is staged glass, not a
simulated material boundary. Absolute saved pressure, divided by its frame
maximum, thresholded at .08 and gated by mineral times 6, drives cyan volume
emission. Pressure sampled at mineral vertices also tints the gel toward cyan,
normalized by maximum pressure within deposit > .11; amplitude above .12 drives
emission at gain 3. These are chosen optical transfers, not physical
luminescence. There is no invented interior geometry or simulation in Python:
Blender extracts the Julia isosurfaces and photographs them.

Replay, with one local Metal owner at a time:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 WARM_N=384 \
  julia --project=gpu --threads=1 --gcthreads=1 studies/warm_thread.jl FRESH_OUTPUT
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --project=. --threads=1 scripts/check_warm_thread.jl FIXTURE.json
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --project=gpu --threads=1 scripts/check_warm_thread_metal.jl DEVICE_FIXTURE.json
/Applications/Blender.app/Contents/MacOS/Blender -b -t 2 \
  --python scripts/render_warm_thread.py -- FRESH_OUTPUT IMAGE.png \
  --size 1600 --samples 64
```

The current helper pins source, allocates only four used World fields, checks an
8 GiB resident-plus-capture memory estimate, and permits importing its wrapper
for the device fixture. The last two guards were added after the retained native
run; its exact source remains pinned beside its fields under
`outputs/art/runtime-warm-knots-v2-20261016`.

The earlier 256³ pulsed helix/Lissajous run remains intact under
`outputs/art/runtime-warm-thread-20261016`. Its 2,400 steps at dt = .0025 reached
t = 6 in 178.29 seconds of measured evolution; median complete step 66.62 ms,
final mass residual 1.44e−11. Reviewed previews were too lobe-like, motivating
finer continuous sources and more intricate actual trajectories. The failed
source-motion guard at dt = .004 is preserved separately. Native field captures,
VDB files and large original images stay outside Git.
