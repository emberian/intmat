# Native dielectric branching refinement

These are experimental study programs, not a replacement for the core `DielectricGrowth` solver. They grow new conducting cells on the requested native lattice. The geometric multigrid solver uses red/black smoothing, restricted residuals and interpolated corrections; the final harmonic field and growth decisions remain on the fine grid. Conductors stay at potential one and the exterior box at zero.

The hierarchy requires cubic grids with `n = 2^k + 1`, at least 5: for example 33, 65, 129 or 257. Both production studies reject unsupported sizes before allocating the simulation or creating an output directory.

## Run

Use a fresh recipe slug: completed and partial output directories are protected. The Metal environment already supplies the required packages; no `JULIA_LOAD_PATH` override is needed.

```sh
julia --project=gpu --threads=1 bin/interior.jl study refinement_branches_metal grow experiments/refinement/branch-thunder-root-257-metal.toml
```

For a bounded batch, set `BRANCH_BATCH_EVENTS`. Each batch stops between complete events and saves a portable checkpoint; ordinary checkpoints are also saved every 25 events.

```sh
BRANCH_BATCH_EVENTS=200 julia --project=gpu --threads=1 bin/interior.jl study refinement_branches_metal grow MY_RECIPE.toml
BRANCH_BATCH_EVENTS=200 julia --project=gpu --threads=1 bin/interior.jl study refinement_branches_metal grow MY_RECIPE.toml outputs/refinement/MY_SLUG-metal/checkpoint.jls
```

An existing partial output requires its own checkpoint explicitly. Resume validates grid, seed, eta, tolerance and extent before writing. To continue a completed specimen to a longer horizon, use a fresh recipe slug and its existing checkpoint. The nutrient/deposit reservoir is reconstructed exactly from occupied sites and the fixed .8 transfer; unused `World.activity` remains zero, while the exported optical activity is a declared deposit transfer.

The CPU reference study has no partial-resume interface:

```sh
julia --project --threads=1 bin/interior.jl study refinement_branches MY_RECIPE.toml
```

The files can also be run directly with `julia --project[=gpu] --threads=1 studies/NAME.jl ...`. GPU benchmarks require the explicit `benchmark` subcommand:

```sh
julia --project=gpu --threads=1 bin/interior.jl study refinement_branches_metal benchmark 65 129 257
```

## Measured specimen and checks

`outputs/refinement/branch-thunder-root-257-metal` contains a fresh 257³ seed831/eta7 run over the same physical extent `(2,2,2)` as the earlier 64³ collection. It reached its declared **800-event horizon**, not ground: 801 conducting sites, 827 exact face-adjacency edges, 192 tips and 203 branching vertices. Its occupied span is 91×48×99 cells. Equal event counts or event-clock times across resolutions do **not** imply equal physical reach or matching morphology: every event adds one native cell.

Growth took **1001.88 seconds** on the tested Metal device, excluding initialization and export. Late solves took about .9 seconds. The final independent core harmonic residual was **4.371e-7**, below the 1e-6 target; potential stayed in `[0,1]`, nutrient plus deposit drift was exactly zero, and the graph was independently checked connected.

At 65³, CPU and Metal multigrid potentials were bit-identical. A 100-event 33³ growth comparison also matched bit-for-bit. Against the existing CPU SOR solver, tight small-domain solves selected the same 100 events, with maximum potential discrepancy 3.22e-6. These checks establish solver parity for the tested cases, not cross-resolution morphological convergence.

Artifacts retain native scalar fields, exact `channels.obj`, full recorded observations, the latest resumable checkpoint, 200-event geometry diagnostics, a render job and validation results. `source/` preserves the exact producer sources and environments used for this specimen, with verified hashes, before subsequent CLI cleanup. Optical filament radius is styling; it does not change simulated channel thickness or graph connectivity. For resumed batches, the current report's elapsed time covers the final invocation only.
