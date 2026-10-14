# Interior

A Julia studio for matter with an interior life. Local activity moves through fields and particles, spends feedstock to leave deposits, and encounters the changed material on its next update. The intended toolbox extends to waves, electrochemistry, phase change, mechanics, changing interaction graphs and optical feedback.

This checkout contains an executable family of exploratory models. They share state, recording and observation infrastructure, while retaining distinct numerical laws. The supplied `interior-living-materials.md` and `interior-living-lab.html` are preserved historical design/playback artifacts. Their earlier Python implementation was not supplied.

## Continuous active-compartment arena

The live arena couples inherited recurrent inhabitants, separate lifetime plasticity, conserved phase boundaries, signal-controlled permeability, and paid deposition/erosion/remodeling. Food splits into structural cargo and spent inert fuel; energy comes from catabolism. A chemostat supplies and removes material with explicit budgets. The earlier `Ecology` example remains available as a separate simpler model.

Start the simulation, then run the viewer server in a second terminal:

```sh
julia --project=. --threads=4 examples/arena.jl --profile examples/arena/active-2d.toml --output outputs/arena/active-2d --continuous
python3 scripts/serve_arena.py outputs/arena/active-2d --port 8780
```

Open <http://127.0.0.1:8780/arena.html>. The live inspector exposes actual fields, identities, inherited genome signatures, plastic memory, actions and budgets. Pause, resume and checkpoint controls act on the local runtime. Checkpoints occupy a bounded ring; a run directory has an owner marker and a single writer. Restart instructions and limitations are in [the arena guide](docs/arena.md). A published replay will be explicitly recorded playback, since GitHub Pages does not execute Julia.

[The neural ecology](docs/neural-ecology.md), [active medium](docs/compartments.md), [semi-implicit solver](docs/implicit-phase.md), and [matched experiments](docs/arena-experiments.md) describe the laws and observations. Sustained turnover, extinction and loss of lineage diversity are outcomes to study; this system does not establish open-ended evolution.

## Explore the studio

Visit the [published material gallery](https://emberian.github.io/intmat/) for recorded volumes and rendered studies. See [publishing](docs/publishing.md) to refresh its curated assets. Open `outputs/index.html` for the generated local gallery. The local preview server, when running from this checkout, serves it at `http://localhost:8769/outputs/index.html`.

| Direction | Executable experiment | What is implemented |
|---|---|---|
| Active material | `examples/run.jl examples/choir.toml outputs/choir` | Excitable fields, mobile recurrent particles, bounded deposition |
| Neural constitution | `examples/run.jl examples/neural.toml outputs/neural` | Persistent recurrent field, configurable rules, two-way material sensing |
| Electrical growth | `examples/electrical.jl` | Harmonic potential solve, screened stochastic branching, material conversion |
| Beam writing | `examples/writing.jl` | Moving Gaussian heating, diffusion, irreversible dose and written scaffold; living/passive pair |
| Acoustics | `examples/acoustic.jl` | Heterogeneous damped wave propagation, intensity envelope, mineral feedback |
| Phase separation | `examples/phase.jl` | Conservative Cahn–Hilliard demixing, interface activity, phase-dependent nutrient flux |
| Deformation/fracture | `examples/mechanics.jl` | Persistent spring graph, driven deformation, irreversible bond breakage and energy ledgers |
| Evolving ecology | `examples/ecology.jl` | Individual inherited neural rules, birth-only mutations, funded actions, material/energy ledgers |
| Curved-path optics | `examples/grin.jl` | Eikonal rays through baked spatial index; constant-index comparison and refinement |
| Coupled experiments | `examples/compositions.jl` | Thermal/neural/material and acoustic/neural/material ablations |

Run Julia examples with `julia --project=. --threads=4 ...`. Individual examples record finite experiments, not calibrated predictions of physical specimens. Check their documentation for output arguments; full-state records refuse silent replacement, so use a fresh output directory for reruns.

The [composition guide](docs/composition.md) introduces named `Experiment` stages, explicit field publishers, fixed-step observation-independent scheduling, `OpticalMap`, and material diagnostics including carrier-held inventory. [Recording](docs/recording.md) preserves owned neural, particle and geometric state as well as world fields. Records are view-independent samples with source hashes, not a promised restart format.

The [Fraser research](docs/fraser-experiments.md) identifies The MOld, links original posts, and distinguishes author descriptions from observations of the downloaded video. Its compartmental structures suggest an especially useful next synthesis: active boundaries that alter permeability and optics.

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

## Larger experiments

[Measured bigrun profiles](docs/bigruns.md) cover 128³ and 256³ with fixed physical extents, stable timesteps, streamed captures and explicit memory/disk budgets. Preview a budget before allocating the volume:

```sh
julia --project=gpu examples/bigrun.jl examples/bigruns/dense-128.toml outputs/bigruns/my-study --dry-run
```

Remove `--dry-run` for the full configured interval, or add `--steps 100` for a capped pilot. These dense profiles use analytic excitation; neural and particle workloads have separate costs.

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

[Process atlas](docs/process-atlas.md) develops seven potential artworks with specific representations, exchange contracts, numerical limits and primary sources. [Project tracker](TRACKER.md) records actual outstanding capabilities. The fixed-rule choir and neural scenes have no online learning or evolution. `Ecology` is a separate evolutionary model with explicit birth, inheritance and mutation; its outcomes do not establish open-ended intelligence or biological organization.
