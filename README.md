# Interior

A Julia studio for matter with an interior life. Local activity moves through fields and particles, spends feedstock to leave deposits, and encounters the changed material on its next update. The intended toolbox extends to waves, electrochemistry, phase change, mechanics, changing interaction graphs and optical feedback.

This checkout contains a working first composition, not all those solvers. The supplied `interior-living-materials.md` and `interior-living-lab.html` are preserved historical design/playback artifacts. Their earlier Python implementation was not supplied.

## Run a specimen

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=4 examples/run.jl examples/choir.toml outputs/choir
python3 scripts/build_viewer.py outputs/choir
```

Open `outputs/choir/viewer.html` for offline recorded 3D playback. Drag to orbit; inspect time and interior sections. It displays captured Julia states, not a separate JavaScript simulation. TOML parameters control the experiment; viewer controls only change observation.

The run saves Float32 fields for every captured frame, full configuration and mass diagnostics in `run.json`, a native `specimen.ppm`, and volume data in `volume/`. On this Mac, FFmpeg can convert the PPM to PNG. Keep interesting intermediate states: this default scene eventually reacts through much of the domain.

```sh
/Applications/Blender.app/Contents/MacOS/Blender -b \
  --python scripts/render_blender.py -- \
  outputs/choir/volume outputs/choir/cycles.png --samples 32 --size 640
```

See [rendering](docs/rendering.md) for the optical model, exports, and Blender/OpenVDB requirements. Native and browser previews use straight-ray transport. Cycles supplies the separate volumetric rendering route; exporting density does not itself implement refractive-index gradients or polarization.

## Compose processes

```julia
using Interior
world = World(Grid(48,48,48); seed=17)  # use Grid(48,48) for 2D
seed_activity!(world)
processes = (
    WaveForcing(amplitude=0.07f0),
    Excitable(diffusion=0.003f0),
    Diffusion(:nutrient, 0.002; material_feedback=3),
    Deposition(rate=0.55f0),
)
for _ in 1:400
    step!(world, processes, 0.01f0)
end
```

A process is a concrete Julia type with `advance!(world, process, dt)`, parameter validation, and `stability_limit(world, process)`. Tuples specify causal update order; this is first-order operator splitting. Fields stay in backend arrays. `ParticleCarrier` adds moving state, inventory and a fixed recurrent rule, and can join the same tuple on CPU. Its latent coordinates are separate from grid activity.

For a neural material, use `examples/neural.toml`: `NeuralField` carries two persistent latent channels, configurable recurrent weights, oscillator/random initialization, and material/nutrient/pressure sensing. Both this process and the analytic grid processes have run on Metal. With the optional GPU environment installed, use `julia --project=gpu --threads=4 examples/run.jl examples/neural.toml outputs/neural-metal metal`. The scene must have `particles.count=0` for this backend.

The current analytic excitation is FitzHugh–Nagumo-like. The standing wave is an imposed kinematic stimulus, not a solved acoustic field. Deposition converts nutrient to solid and dissolution reverses it. Particle uptake and secretion move the same inventory through a third reservoir. Changes in deposit lower nutrient mobility and affect activity, creating a feedback loop.

Grid fields use declared nondimensional per-cell inventory units, and particle storage shares those units. Raw totals conserve within a fixed discretization. A physical concentration model across changing resolutions requires explicit cell-volume conversion; equal raw parameters are not a resolution-independent particle chemistry. The computational box is periodic. The ellipsoidal artwork mask is an observation crop, not a mechanical boundary.

## Verify and extend

```sh
julia --project=. --threads=4 test/runtests.jl
julia --project=. --threads=4 examples/verify.jl
```

Tests cover conservation, positivity, identity under particle reordering, physical movement under grid refinement, transport integration and composition. The verification script records matched ablations and timestep refinement separately from unit tests. See [computing](docs/computing.md) for the optional Metal environment and measured grid throughput; particle exchange is currently CPU-only.

[Process atlas](docs/process-atlas.md) develops seven potential artworks with specific representations, exchange contracts, numerical limits and primary sources. [Project tracker](TRACKER.md) records actual outstanding capabilities. Current recurrent rules are supplied dynamical systems: no learning, evolution or self-reproduction is claimed.
