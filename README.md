# Interior

**Evolving matter. Julia computes it; renderers observe it.**

The material has state: phase boundaries, permeability, deposits, signals, heat,
pressure and damage. Its inhabitants have state: recurrent brains, plastic memory,
inherited metabolism, ancestry and finite inventories. Each changes what the other
encounters next.

[Gallery](https://emberian.github.io/intmat/) · [Metabolic replay](https://emberian.github.io/intmat/metabolic/arena.html) · [Working state](TRACKER.md)

## Run

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=4 examples/arena.jl \
  --profile examples/arena/metabolic-2d.toml --output outputs/my-arena --continuous
python3 scripts/serve_arena.py outputs/my-arena --port 8781
```

Open <http://127.0.0.1:8781/arena.html>. Inspect fields, individuals, inherited
allocations and budgets. Pause or checkpoint without restarting. Use a fresh output
directory; [resume explicitly](docs/arena.md). One writer per run. Checkpoints retain
full state and check Julia/source versions. The runner is not a reboot service.

`active-2d.toml` selects the earlier single-resource ecology. `metabolic-2d.toml`
selects the multispecies model. `metabolic-96.toml` enlarges the domain at the same cell spacing; these are distinct experiments.

## What evolves

| Mechanism | Contract |
|---|---|
| Recurrent controller | Heritable weights; separate bounded lifetime plasticity; newborn memory resets. |
| Metabolism | Enzyme allocations sum to one. Higher efficiency costs throughput. Assimilation competes with fuel use. |
| Cross-feeding | Partial catabolism exports lower-potential substrate. Exported chemical energy cannot also fund work. |
| Shared boundaries | Agents pay cargo and energy to maintain film. Film alters transport and phase forcing; decay returns structural matter. |
| Selection | Birth, mutation, death and spatial competition. Neither cooperation nor increasing complexity is guaranteed. |

[Chemistry](docs/resources.md) · [Metabolism](docs/metabolic-ecology.md) · [Boundary maintenance](docs/cooperative-boundaries.md) · [Experiments](docs/metabolic-experiments.md)

## Compose

A `World` owns grid fields. Processes own additional state. An ordered tuple defines
causal order; `advance!` advances one process, `step!` advances the tuple.
[Experiments](docs/composition.md) add declared ports, timestep checks and observations.
Borrowed references do not confer ownership of a second inventory.

| Process | Entry point |
|---|---|
| Excitable / recurrent material | `examples/run.jl`, `examples/neural.toml` |
| Thermal writing / acoustics | `examples/writing.jl`, `examples/acoustic.jl` |
| Dielectric branching | `examples/electrical.jl` |
| Conserved phase separation | `examples/phase.jl`; [semi-implicit solver](docs/implicit-phase.md) |
| Deformation / fracture | `examples/mechanics.jl` |
| Evolving inhabitants | `examples/ecology.jl`, `examples/arena.jl` |
| Graded-index light transport | `examples/grin.jl` |

[Process atlas](docs/process-atlas.md) describes the larger toolbox. The supplied
`interior-living-*` files preserve the original design context.

## Measure before scaling

The semi-implicit constant-mobility solver completed a **64³, t=10** phase run in
**19 s**; measured timings imply **84×** explicit throughput over equal simulated
time. This is neither an accuracy guarantee nor a variable-mobility solver.

At **2,400 inhabitants**, reusable neighborhood/controller buffers reduced matched
complete neural steps **19.13→11.20 ms**, allocations **28.42→4.32 MB**. Complete
metabolic steps measured **12.74 ms / 5.10 MB**. Births resize once per step; dead
individual histories are pruned. [Profiling](docs/arena-performance.md)

Metal supports several dense grid processes. The evolving population remains a CPU
implementation. Device-resident controllers beat CPU; per-tick transfers reverse
that result. A real Kaggle T4 ran the dense-field probe; TPU and hbox work have
separate status. [Accelerators](docs/accelerators.md) · [Larger-run budgets](docs/bigruns.md)

Small runs establish mechanisms and expose failures. They do not bound long-run
behavior. Sustained turnover is not open-ended evolution; founder loss is not loss
of all genotype variation. Compare seeds, timescales, capacities and resolutions.

## Invariants and limits

- Material and chemical work have separate ledgers. Structural recycling does not regenerate fuel.
- Inventories use per-cell units. Cross-resolution chemistry needs explicit cell-volume normalization.
- Fields are periodic. Optical crops and polished surfaces are rendering choices, not solved gel walls.
- Observation records preserve sampled state; executable checkpoints are a separate format.
- GitHub Pages replays recorded states. It does not run Julia.

```sh
julia --project=. --threads=4 test/runtests.jl
python3 -m unittest discover -s scripts -p 'test_*.py'
```

[Bake and photograph](docs/rendering.md) · [Record](docs/recording.md) · [Publish](docs/publishing.md) · [Fraser source research](docs/fraser-experiments.md)
