# Interior

**Evolving matter. Julia computes it; renderers observe it.**

The material has state: phase boundaries, permeability, deposits, signals, heat,
pressure and damage. Its inhabitants have state: recurrent brains, plastic memory,
inherited metabolism, ancestry and finite inventories. Each changes what the other
encounters next.

[Gallery](https://emberian.github.io/intmat/) · [Films](https://emberian.github.io/intmat/animations/) · [Laboratory](https://emberian.github.io/intmat/laboratory/) · [Living cells](https://emberian.github.io/intmat/cellular/) · [Architecture](docs/cellular.md) · [Working state](TRACKER.md)

## Material collection

Six families: phase surfaces, written inclusions, dielectric branches, fractured
networks, acoustic/neural contours and 3D metabolic inhabitants. Each photograph
observes a saved simulation through chosen optics. The first series adds 36
specimens at 1600². Refinement now evolves native 384³ phase, writing and waves;
768³ writing is complete. Refined 2400² photographs, recurrent matter, coupled optics and a repairing weave bring the gallery to 79
studies, with corrected crop geometry and resin, wax, frost, bronze and copper
optics. Native 1152³ Julia writing completed in 56.9 minutes with 18.35 GB of state.
Constructive ecology now has actual 64³ trajectories and a 128³ GPU pilot;
mechanics retains its actual deformed bond graph.
[Recipes](docs/collection-20261015.md) · [Spatial refinement](docs/spatial-refinement.md)

**Recurrent matter:** eight coupled internal channels per voxel write and erase
solid; solid changes their communication and thresholds. Native 256³ trajectories
now run entirely on the GPU in Julia. These are designed, untrained dynamics with
material memory. [Photograph and motion](https://emberian.github.io/intmat/recurrent/) · [Model and throughput](docs/recurrent-material.md)

**The nacre fault:** native 256³ recurrent mineral and pressure waves change one
another. Perforated folds are photographed through nacre glass, copper and
refractive pressure fronts: a 2000² still and 65-frame film. The folds are seeded;
their subsequent motion comes from the coupled Julia simulation.
[Study](docs/nacre-fault.md) · [Films](https://emberian.github.io/intmat/animations/)

**The second seam:** damage opens a woven sheet; finite repair stock makes new
connections across the wound. The photographed surface follows the deformed
network. [Mechanism and photograph](docs/second-seam.md)

**Interior instruments:** rotate and cut through actual resonant shells in
[What the bell kept](https://emberian.github.io/intmat/collection-20261015/resonant-reliquary/viewer.html).
The [laboratory](https://emberian.github.io/intmat/laboratory/) maps controller
memory and action-credit learning. Three founder banks improve frozen resource-use
tests after episodic training; autonomous learning remains under development.
[Mechanism and limits](docs/research/action-credit-materials.md)

**Paid construction:** a brief secretion pulse unlocks food and repays its maker
in three matched founder banks; without raw substrate it stays costly. Explore
[who gains and what remains](https://emberian.github.io/intmat/laboratory/secretion-payback/).
This is a measured opportunity for learning, not learned cooperation.

[Construction credit](https://emberian.github.io/intmat/laboratory/construction-credit/) exposes the next constraint: reserve demand can consume precursor before construction. Joint reserve/material control restores an action the learner can change.

[Learning, held in matter](https://emberian.github.io/intmat/laboratory/living-construction/)
follows 24→30 moving cells through six births. Later offspring receive paid copies
of learned readouts; 17 actual captures expose bodies, ancestry and material.
The mechanism works; an advantage from inheritance remains unestablished.

[After the feed](https://emberian.github.io/intmat/laboratory/living-patch/)
follows a more productive material recipe through 33 captures. Compare actual
bonds, retained catalyst and stored food before and after the finite feed ends.
The local Julia workbench can now fork the complete living state—including
pending learning credit—and change what a branch senses.
[Run it](docs/native-workbench-control.md#living-construction-adapter-bounded-local-session).

**Memory that changes construction:** continuously visible resource ramps can
reach the same current balance while predicting opposite future arrivals. A paid
controller learns to allocate from that history; swapping its history reverses
the benefit. Correcting a saturated predictor enables the learning.
[Explore the experiment](https://emberian.github.io/intmat/laboratory/continuous-ramp/) · [Mechanism and transfer](docs/research/continuous-ramp-construction.md).

**The light takes the long way:** a changing 256³ mineral bends light through
three fixed apertures. [Rotate, cut and compare actual paths](https://emberian.github.io/intmat/collection-20261016/refractive-memory/viewer.html)
across five saved states; the 2000² photograph reveals its amber and copper interior.

**What the light persuades** closes that loop: absorbed light changes recurrent
activity, growth and subsequent illumination. [Explore mineral and absorbed power](https://emberian.github.io/intmat/collection-20261016/illuminated-memory/viewer.html).
**Weather inside a lens** photographs finite mineral growth and thermal erosion
that redirect current through opaline folds.

**What the brine left open:** finite reagents precipitate a porous 128³ mineral.
Its deposits redirect diffusion and reshape the interfaces where reaction occurs.
A separate native 384³ epoch now completes with Float64 chemical inventories in
17.3 minutes; its full fields are preserved for photography.
[Execution record](docs/native-brine-producer.md) · [Explore four recorded states](https://emberian.github.io/intmat/collection-20261016/brine-left-open/viewer.html) · [Process](docs/brine-memory.md).

**The chamber the brine kept:** a section through the new native 384³ mineral.
Saved phase and remaining reagents select copper, petrol-dark metal and thin-film
reflections. The cut exposes computed passages; its open ends are photographic
boundaries. [Photograph and optical mapping](docs/brine-interior-384.md).

**A chord held in salt:** three pressure packets write an interference structure
that changes subsequent wave propagation. An exact cutaway reveals its interior.
[Study](docs/chord-salt.md). Its moving relative,
[Where the note goes](https://emberian.github.io/intmat/animations/#wave-memory-motion-128),
records 40 actual states of wave-written mineral growing and receding.

**Inner shore** enters the light-remodeled mineral: a macro view of its actual
mineral walls and absorbed-power sheets. [Photographic study](docs/inner-shore.md).

**The field inside the fracture:** surviving mechanical bonds conduct voltage;
electrical dose weakens later bonds. Gold potential contours reveal the resulting
field inside fractured violet support. [Model and photograph](docs/voltage-pieces.md).

**After the last feeding:** actual cell bonds become a gold scaffold inside
optically stained catalyst and stored food. The photograph observes a native
32³ living state after its finite feed ends. [Study](docs/living-patch-photograph.md).

## Run

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=4 bin/interior.jl arena \
  --profile experiments/arena/cellular-2d.toml --output outputs/my-arena --continuous
python3 scripts/serve_arena.py outputs/my-arena --port 8782
```

Open <http://127.0.0.1:8782/arena.html>. Inspect fields, individuals, inherited
allocations and budgets. Pause or checkpoint without restarting. Use a fresh output
directory; [resume explicitly](docs/arena.md). One writer per run. Checkpoints retain
full state and check Julia/source versions. The runner is not a reboot service.

`cellular-2d.toml` gives cells physical bodies, persistent contacts, internal
chemistry and division-preserved developmental state. Bonds carry force,
metabolites and regulatory context. There is no colony fitness bonus.
[Design](docs/cellular.md) · [Long soaks](docs/cellular-soaks.md)

`catalytic-mosaic.toml` adds paid, inherited construction of persistent catalysts.
They unlock two inaccessible feedstocks into shared food. Cells can change what
neighbors can eat; no group reward or exclusive access is imposed.
[Mechanism and first trajectories](docs/catalytic-habitat.md)

`active-2d.toml` selects the earlier single-resource ecology. `metabolic-2d.toml`
selects the multispecies model. `metabolic-96.toml` enlarges the domain at the same cell spacing; these are distinct experiments.

## What evolves

The cellular model inherits recurrent development and slow memory through division;
physical contacts carry forces and exchange precursors. Cells allocate work among
synthesis, motion and maintenance, with starvation rather than a fixed age limit.

The earlier metabolic model supplies a complementary set of mechanisms:

| Mechanism | Contract |
|---|---|
| Recurrent controller | Heritable weights; separate bounded lifetime plasticity; newborn memory resets. |
| Metabolism | Enzyme allocations sum to one. Higher efficiency costs throughput. Assimilation competes with fuel use. |
| Cross-feeding | Partial catabolism exports lower-potential substrate. Exported chemical energy cannot also fund work. |
| Shared boundaries | Agents pay cargo and energy to maintain film. Film alters transport and phase forcing; decay returns structural matter. |
| Selection | Birth, mutation, death and spatial competition. Neither cooperation nor increasing complexity is guaranteed. |

In the earlier metabolic substrate, a funded three-seed study found inherited
competitive improvement: evolved
descendants won 48/48 reset-state assays, including 24 without lifetime plasticity.
Noisy controls and small evolutionary N limit the claim; cooperation remains unproven.
[Results and component swaps](docs/alife-adaptation-results.md)

[Chemistry](docs/resources.md) · [Metabolism](docs/metabolic-ecology.md) · [Boundary maintenance](docs/cooperative-boundaries.md) · [Experiments](docs/metabolic-experiments.md)

## Investigate capable worlds

Map the clocks, then intervene: erase neural state, plastic marks or retained
food in matched copies; measure their response to the same nutrient pulse.
Distinguish information retained from information used, and lifetime development
from inherited adaptation. Preserve diverse parameterizations and actual fields
for later simulation and rendering.

[Operating map](docs/garden-operating-map.md) · [Causal windows](docs/native-memory-interventions.md) · [Inventory repair](docs/native-inventory-precision.md) · [Research and PDFs](docs/research/learning-materials.md)

[Explore memory and feedback](https://emberian.github.io/intmat/laboratory/) across
recorded controller regimes. [Material receptors](docs/native-material-sensing.md)
connect retained food and construction to recurrent control; their native GPU
pilot costs 1.1% per complete step. The local [Julia workbench](docs/native-workbench-control.md)
now forks complete states, changes receptor scales and compares independent
futures. Browser controls submit finite commands; the worker pauses between them.

## Compose

A `World` owns grid fields. Processes own additional state. An ordered tuple defines
causal order; `advance!` advances one process, `step!` advances the tuple.
[Experiments](docs/composition.md) add declared ports, timestep checks and observations.
Borrowed references do not confer ownership of a second inventory.

Commands below follow `julia --project=. bin/interior.jl`; `--help` lists them.

| Process | Command |
|---|---|
| Excitable / recurrent material | `run experiments/materials/neural.toml outputs/neural` |
| Thermal writing / acoustics | `study writing`, `study acoustic` |
| Dielectric branching | `study electrical` |
| Conserved phase separation | `study phase`; [semi-implicit solver](docs/implicit-phase.md) |
| Deformation / fracture | `study mechanics` |
| Evolving inhabitants | `arena`, `study ecology` |
| Graded-index light transport | `study grin` |
| Budgeted larger runs | `bigrun experiments/bigruns/dense-128.toml outputs/dense-128` |

`src/` is the engine; `bin/` is its command interface; `experiments/` holds recipes;
`studies/` holds reproducible investigations; `scripts/` handles measurement,
rendering and publishing. `outputs/` holds generated state and art. `examples/`
retains compatibility forwards, not the primary interface.

[Process atlas](docs/process-atlas.md) describes the larger toolbox. The supplied
`interior-living-*` files preserve the original design context.

## Measure before scaling

Complete cellular steps now run **9.38–12.80× faster** on four matched hbox workloads,
including 1,024 cells. Numerical specialization and shared contact searches preserve
bit-identical 100-step continuations. No timestep or physics was relaxed.
[Measurements and reproduction](docs/cellular-performance.md)

A second pass adds **1.08–1.27×** on 1,024–4,096-cell workloads and cuts allocated
bytes **23–43%**. The dense recurrent GPU kernel separately improves **2.90×** on
AMD128³: **2.02 billion voxel-steps/s** warm, **1.60 billion** over an actual
10,000-step run including captures/export. These are distinct throughput measures.
[Cellular follow-up](docs/cellular-performance-pass2.md) · [Dense runner](docs/recurrent-material.md)

The semi-implicit constant-mobility solver completed a **64³, t=10** phase run in
**19 s**; measured timings imply **84×** explicit throughput over equal simulated
time. This is neither an accuracy guarantee nor a variable-mobility solver.

At **2,400 inhabitants**, reusable neighborhood/controller buffers reduced matched
complete neural steps **19.13→11.20 ms**, allocations **28.42→4.32 MB**. Complete
metabolic steps measured **12.74 ms / 5.10 MB**. Births resize once per step; dead
individual histories are pruned. [Profiling](docs/arena-performance.md)

**Device-native cells now run on AMD**, including inherited controllers, division,
contacts, transport and paid habitat construction. Their Float32 synchronous law
differs from the serial CPU ecology. Checkpoints retain complete future state;
generic recording copies only explicitly requested observations.
[Device model and measurements](docs/device-cellular.md)

Concrete slab work reduces matched 192³ brine initialization **19.0→0.737 s**
with bit-identical fields. This improves preparation, not sustained stepping.
An exact empty-chemistry branch improves complete living CPU steps **7.4–7.6%**;
GPU whole-step gain is below measurement resolution.
[Initialization](docs/brine-initialization.md) · [Living steps](docs/living-step-performance.md)

Recent decomposition experiments distinguish storage, submission and simulation:

Paid-learning batches now advance 32 independent 8³ worlds **26.7–28.1× faster**
than queued standalone worlds on AMD. Including setup and observation, the
16-step assays gain **6.90–7.77×**, with exact compared states and payments.
This accelerates experiments in parallel; it does not speed one world's clock.
[Measurements and restart contract](docs/learning-world-batches.md)

| Change | Measured result | Boundary |
|---|---|---|
| Axis support reconstruction | 221× less support storage at 256³; 0.230× step throughput | Carry traversal improves the first version; materialized remains faster. |
| HIP graph replay | 1.020× / 1.002× whole-step gain at 32³ / 128³ | Submission shrinks; device execution dominates. |
| Independent-world batching | 2.41× aggregate throughput for four 64-cell worlds | Exact tested states; accelerates ensembles, not one trajectory. |
| Local harmonic correction | 1.22× total time gain over nine 257³ solves | Recorded growth geometry; coupled selection gate pending. |
| Continuous living tissue | 2.78–2.84 billion site-steps/s at 128³ | Different seven-field law; no individual ancestry. |

[Support](docs/device-support-performance.md) · [Graphs](docs/device-graph-replay.md) · [World batches](docs/cellular-world-batches.md) ·
[Local solves](docs/electrical-temporal-reuse.md) · [Tissue](docs/living-tissue.md)

Metal and AMD support dense Julia processes; native 768³ writer steps measured
29.15 ms on hbox. A real Kaggle T4 ran the dense-field probe; TPU remains a separate
unverified route.
[Reproduce](docs/hbox-runtime.md) · [Accelerators](docs/accelerators.md) · [Larger-run budgets](docs/bigruns.md)

Small runs establish mechanisms and expose failures. They do not bound long-run
behavior. Sustained turnover is not open-ended evolution; founder loss is not loss
of all genotype variation. Compare seeds, timescales, capacities and resolutions.
[Overnight observations](docs/alife-overnight.md) · [Controlled adaptation results](docs/alife-adaptation-results.md)

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
