# Growing cellular material

`CellularEcology` makes the inhabitants part of the material: cells occupy space,
grow from imported precursors, push on one another, attach, exchange intracellular
chemistry and divide. The older `NeuralEcology` and `MetabolicEcology` remain
available as different substrates and as the models used by their published studies.

```sh
julia --project=. --threads=2 --gcthreads=1 bin/interior.jl arena \
  --profile experiments/arena/cellular-2d.toml \
  --output outputs/my-cellular-arena --continuous
python3 scripts/serve_arena.py outputs/my-cellular-arena --port 8782
```

Use a fresh output directory. The same runner, snapshots, controls and strict
source-aware checkpoints support this model. `studies/cellular_probe.jl` runs a
bounded matched comparison with junction transport or adhesion removed.

## Why another substrate

Reading the coupled source suggests structural constraints, not a proof that the
previous system cannot do something interesting. Its point inhabitants can learn,
evolve, compete, signal, remodel fields and improve inherited performance. The
adaptation studies establish some of that. They do not exhaust its potential.

Several pathways toward a different kind of organization were simply absent:

- A point agent had no body volume, persistent contact topology or internal
  chemical transport. Building a film changed a common environment but could not
  construct a mechanically connected organism.
- Film lowered both environmental diffusion and the payer's own uptake/motion
  through one permeability factor. Retention could help in particular geometries;
  a general selective advantage for cooperation was not built into that equation.
- Fixed lifespan death, unbounded storage and a global population cap supplied
  much of the old life history. Division reset recurrent and plastic state, which
  is useful for inherited-performance assays but discards a possible developmental
  channel during ordinary growth.
- Cell-local inventories, gradients and daughter offsets changed physical meaning
  when the environmental grid was refined.
- Immediate energy gain modulated all Hebbian updates, without a delayed
  eligibility mechanism or an established credit-assignment interpretation.

The new design takes an analytical step beyond those pathways. A persistent
contact graph can carry chemistry, force and regulatory context at once. Growth
changes that graph, and state can survive its division. Spatial differentiation,
transport networks, fragmentation and coordinated behavior become expressible
without assigning a bonus for colony size or calling a connected cluster an
organism. Whether useful versions evolve remains open; short pilots are clues,
not a boundary on that possibility.

## What the cell does

Each cell has physical position, structural biomass, two intracellular precursor
stores, a finite chemical-work reservoir, recurrent state, slow adaptation,
bounded plastic marks, an inherited weight array and stable identity. A body
radius follows total material volume at a declared density. Cells are soft disks
in 2D and soft spheres in 3D.

Four regulated enzyme allocations sum to one: import A, import B, synthesize
biomass and oxidize precursors. Import draws from a normalized kernel of fixed
physical radius. Both precursors are necessary for synthesis; either can supply
work through oxidation. Oxidation is limited by catalytic capacity and a
regulated reserve target plus expected work. It does not obligatorily burn all
incoming material when the reservoir is full. The exponent on transporter
allocation is an explicit synthetic tradeoff, not a biological universal.

Maintenance and neural operation consume work. Unfunded maintenance erodes
biomass; a sufficiently small cell dies. Shrinking capacity degrades excess
precursor to inert waste. There is no maximum-age death. Biomass over the division
threshold can divide when its controller permits; body material, precursors and
work split equally. Daughter placement breaks symmetry locally, after which
contact forces separate the bodies. Recurrence, adaptation and plastic marks
are inherited as developmental state; sparse mutations alter the daughter's
genotype. These are separate hereditary and developmental channels.

The eleven action channels are movement xyz, four enzyme-allocation logits,
adhesion, junction openness, division and reserve target. Sensors include
physical resource concentrations/gradients, internal inventories, size,
bond degree/tension, contact compression, optional phase/signal context and age.
Neighbor recurrent states enter through open junctions. Local bounded
Oja-like plasticity affects recurrent rows; it is developmental plasticity,
not a claim of successful reinforcement learning.

## Contacts and exchange

Spatial bins find contacts in physical coordinates. Repulsion and tensile bonds
act equally and oppositely. Mutually adhesive neighbors can form bonds without an
ancestry gate. Stable-ID endpoints survive array reordering; withdrawal, excessive
strain and death remove bonds. Division can retain an adhesive parent–daughter
contact within the configured degree bound.

Junctions equalize intracellular concentrations by conservative positive
two-compartment exchanges. The graph update is an ordered splitting, not an exact
simultaneous network solve. Shared chemistry has an actual local consequence:
complementary precursors can reach a cell that otherwise lacks one. Exchange also
lets cells lose useful material. No cooperation reward is added.

Active movement pays drag and positive changes in the coarse contact potential.
Passive relaxation backtracks until that potential does not increase. Adaptive
substeps replace the old fixed cell-displacement clamp. Contacts remain soft;
crowded bodies can overlap. Biomass-derived porosity changes extracellular
transport, while the importer samples local concentrations without multiplying
its capacity by that same permeability a second time.

## Units, accounting and composition

Cell positions and radii use physical lengths. Environmental arrays store mass
per voxel. The cellular TOML builder accepts concentrations and source/sink rate
densities, multiplies them by physical voxel volume, and leaves total source
budgets unchanged. A `CellStencil` uses the same physical support after grid
refinement. Unresolved sub-grid support falls back to a nearest voxel; this is
still a discretization error, not exact cross-resolution equivalence.

Mass is environmental resources + biomass + intracellular precursors + explicitly
retained inert rounding residue. The chemical ledger is environmental potential +
intracellular potential + stored work + spent work + dissipated chemical energy,
adjusted for external supply/removal. Actual Float32 field changes fund Float64
cell transfers. `work_spent` includes synthesis, maintenance and motors; it is not
all called heat. Mechanical potential/relaxation have their own ledger. Changing
body radii and coarse motor work do **not** establish a calibrated, closed
thermomechanical energy model.

The process composition is `(resources, cells)` by default. Optional
`[phase] enabled=true` adds the existing medium and phase stages. Cells borrow
the pool and medium; only the separately scheduled owners advance those fields.
The pool's own `medium` reference remains empty so its permeability can combine
the background medium with body porosity. This is a concrete composition of open
dynamical systems: each connection identifies a flux with equal loss and gain at
its endpoints. An additional categorical abstraction is unnecessary to obtain
that contract.

The population capacity remains a computational guard, with recorded blocking
and deterministic, time-varying keyed priorities among simultaneous divisions.
It is not the intended carrying-capacity mechanism. Long experiments should
leave headroom and inspect blocking. Other deliberate constraints remain: a fixed
neuron count per cell, fixed chemical species and stoichiometry, isotropic bodies,
first-order operator splitting, prescribed resource supplies and a CPU population
implementation. These are invitations for further design, not inferred limits
on artificial life.

## Initial checks and explorations

Focused checks exercise proportional physical-support allocation, body contacts,
funded motion, junction conservation, finite stores, complementary synthesis,
starvation, developmental inheritance, sparse mutation, stable IDs and exact
checkpoint continuation. They test whether the chosen machinery runs as stated.

The first short matched probe grew 24 founders into connected populations with
births and actual junction transfer. Its no-junction and no-adhesion arms followed
different trajectories. Those observations neither identify evolved cooperation
nor show that one architecture is fitter; the conditions and founder priors were
constructed, and the trajectories are short. The study retains complete sampled
states, source hashes and individual metrics for subsequent examination.

On a packed 1,000-cell contact benchmark, broadphase reuse and a less conservative
adaptive bound reduced one mechanical step from about 245 ms / 218 MB allocated
to 52 ms / 18.5 MB. This measures mechanics only, not complete arena throughput.
The physical-support buffers and columnar cell state offer a path to larger
populations, while an accelerator-resident lifecycle is still future work.
