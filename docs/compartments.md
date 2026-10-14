# Active compartments

`ActiveMedium(world, phase)` makes a separately advanced
`ImplicitPhaseSeparation` field into an environment for inhabitants. Its
`membrane` aliases the world's deposit. Inhabitants transfer their own cargo
into that material when building and recover it when eroding: the environment
never creates construction material. Nutrient transport uses periodic symmetric
face fluxes with harmonic face permeability and internal explicit substeps;
transport conserves resource mass and keeps nonnegative concentrations positive.

Interfaces and deposited membrane reduce permeability. Diffusible signals open
local permeability, decay exponentially, and bias the phase chemical potential.
Agents write nonnegative increments into `signal_source` and signed remodeling
increments into `remodeling`; the medium consumes and clears both each stage.
Remodeling accumulates in `forcing`, which relaxes exponentially. Deposit couples
through a quadratic composition energy: its derivative is proportional to
`deposit * phi`, pushing dense deposit sites toward mixed composition. These
forces enter the phase solver under a Laplacian, deforming compartment interfaces
without adding phase mass. Signals are informational fields, not part of the
nutrient/cargo/deposit material budget.

Use a phase constructor with `mineral_pinning=0, nutrient_diffusion=0,
interface_activity=0, activity_relaxation=0`. This chooses constant phase mobility
and assigns nutrient transport to the medium exactly once. Advance processes in
`(ecology, medium, phase)` order to apply inhabitant cues immediately, or
`(medium, phase, ecology)` for a one-frame feedback delay. `advance!(world, medium,
dt)` does not advance the phase or world clock. The current implementation is CPU
only. Arrays available for controllers are `phase.phi`, `phase.interface`,
`permeability`, `signal`, `membrane`, and pending `remodeling`.

`record_state(world, (medium, phase), destination; step=...)` records every
medium array and parameter, labels the membrane/world-deposit alias, and records
the referenced phase separately. Recording is observational and preserves grid
shape and native array values.
