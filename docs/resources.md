# Finite chemical resources

`ResourcePool` holds CPU Float32 concentration fields independently of the
structural `World.nutrient` and `World.deposit` fields. Default species are two
substrates with potentials 4 and 3 and a byproduct with potential 1. These are
model units and an exploratory energetic constraint, not measured chemistry.
Structural nutrient, solid, and agent cargo have zero chemical potential.

Construct with `ResourcePool(world; initial=(1,1,0), medium=medium)`.
`fields[species]` exposes each concentration array. A source or sink has a
per-cell rate, a finite global mass budget, and an optional nonnegative spatial
mask. Supply `source_masks`/`sink_masks` as vectors of grid-shaped arrays;
missing masks are uniformly one. Zero masks disable exchange locally. Initial
fields can be shaped directly before recording the initial inventory.

`allocate_resources!(pool, cells, requests)` accepts a species-by-consumer
request matrix and cell indices. Consumers at a shared cell receive proportional
shares of each species. Returned Float64 amounts reflect actual representable
Float32 withdrawals, and the pool records those same amounts in `withdrawn`.
The consumer must track all withdrawn mass in its stores or byproducts and all
withdrawn potential in usable energy, heat, or retained chemical potential.
`resource_add!` reports actual additions and records `returned` (or `supplied`
with `external=true`); it cannot establish that the caller owned that material.
Direct external additions do not debit the configured automatic source budget.

`resource_convert!` transfers equal mass between two local fields only downhill
in potential. Its `(mass, energy, heat)` result uses the actual field deltas.
Efficiency splits lost potential between work and heat; it cannot change the
amount released. Ordinary Float32 representation gives small mass roundoff;
ledger accounting uses actual field changes. Reverse conversion is rejected.
Crossfeeding through the low-potential byproduct can extract only its remaining
potential; recycling structural material cannot recharge a substrate.

`advance!` applies finite source/sink exchanges and conservative periodic
harmonic-face diffusion using `medium_diffuse!`. Stable explicit substeps ensure
nonnegative transport, including large requested timesteps. Optional medium
permeability is read at this stage, so tuple ordering determines which medium
state is used. The pool never changes structural world mass or clock itself.
`step!` advances the clock once for the complete composition.

`resource_mass`, `resource_potential`, and `resource_statistics` expose inventory
and cumulative exchange ledgers. Potential exchange is species mass times its
potential. Closed conversion satisfies initial potential = current potential +
released work + heat up to field representation error. External source/sink
ledgers and consumer withdrawals/returns must be included for open compositions.
`record_process_fields` includes fields, masks, budgets, rates, and statistics.
