# Coupled arena audit

The material budget counts nutrient and deposit inventories per cell, plus each
inhabitant's cargo and waste. These are inventories rather than continuum
concentrations: grid refinement changes total material unless initial inventory
is rescaled. Phase composition is a separate signed conserved quantity. Signal
and remodeling fields are information and chemical-potential cues; they do not
add material. The neural ecology converts `catabolic_fraction` of each uptake into inert
waste and earns `energy_yield` work per catabolized inventory unit; only the
remainder becomes cargo. Spent work becomes heat. Waste participates in material
mass but cannot be eaten again. Death returns cargo, while waste remains inert.
Consequently, in a closed arena cumulative chemical work is bounded by initial
work plus energy yield times initial usable material. Recycling cargo cannot
create unlimited work. Chemostat supply and removal are external material
flows and must be subtracted from the observed inventory change.

`test/arena_composition.jl` exercises coupled transport and forced phase on a 3D
grid, a thin periodic grid, and the all-singleton grid. It compares transport
through a funded deposited barrier with signal-opened transport, includes a zero
population, and checkpoints the actual coupled FFT-backed phase/environment
rather than only a trivial carrier. Resume must preserve the environment's
membrane/world-deposit alias and its reference to the very same phase process;
subsequent fields, IDs, inherited weights, plastic marks, recurrent state,
positions, and phase must agree. Physical-time cue refinement and zero work
funding are checked independently. A periodic boundary crossing immediately
before death returns cargo to the final occupied cell, preserving the budget. The focused suite passes 57 assertions. An independent writer/reader run in
separate Julia processes also reproduced exact fields, inherited weights,
plastic marks, hidden state, and phase after eight resumed steps. FFTW native
plans are rebuilt by execution rather than serialized as persistent state.

Counterexamples identified during implementation:

- Summing only active axes on the all-singleton grid requires an explicit zero
  reduction identity, both for Fourier Laplacian eigenvalues and transport CFL.
- A local population neighborhood must include z in 3D and deduplicate periodic
  aliases on dimensions shorter than three cells. Repeated aliases inflate
  density sensed by controllers.
- Mutable environment signals, forcing, resource fields, and phase fields need
  validation before coupled mutation. A harmonic face connecting two zero
  permeabilities must yield zero flux rather than `0/0`.
- Signal and remodeling outputs are rates: per-stage increments require a
  timestep factor. Zero requested work must not fund active actions.
- Snapshot limits apply to the whole output, including multiple populations;
  limits must be validated before array/index construction.

The stabilized phase update removes the stiff linear fourth-order CFL bound.
Its nonlinear term is explicit and signed composition is not bounded to [-1,1],
so arbitrary timesteps do not carry an unconditional energy or accuracy promise.
The phase energy diagnostic omits work from externally changing chemical
potential; it must not be interpreted as a monotonically decreasing total
energy during active forcing. Temporal refinement and forced mass conservation
are distinct checks.

`record_state` remains a view-independent observation format, not restart state.
The arena checkpoint stores complete process state and aliases, rejects a Julia
version or implementation-hash mismatch by default, and is a trusted local
artifact. Checkpointing must preserve controller inheritance, within-lifetime
plasticity, hidden state, event counters, and random state independently;
recording effective controller weights alone cannot recover that distinction.

The closed-fuel counterexample runs 500 steps with reproduction, death and cargo
recycling, checks the largest material/energy discrepancies over the trajectory,
and requires extinction followed by a plateau of cumulative work and heat. It
also tests zero catabolic fraction: assimilation retains material but grants no
energy. This refutes the previous retained-matter uptake rule, which recycled
fuel while earning fresh work on every pass. The corrected ledger ties earned
energy exactly to cumulative inert waste rather than repeated uptake totals.
