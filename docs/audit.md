# Numerical and composition audit

Reviewed the Julia implementation on 14 October 2026. Historical Python claims were not used as evidence. Scope: dense fields, particle exchange, scene construction, rendering, and the example runner.

## Reproduced and fixed

`step!(w,(Deposition(),Diffusion(:missing,0.01)),0.01)` previously changed deposit before throwing, while time stayed zero. Process preflight now dispatches on the world as well as the process. Invalid diffusion field names and incompatible particle storage/backend, malformed mutable particle state, malformed rules, and invalid particle-exchange inventories are rejected before the composition runs. Every process checks structural, capability, parameter and timestep preconditions. Dense processes do not repeatedly scan every field for arbitrary numerical corruption; the runner checks health at recording cadence, and extrema needed by numerical bounds are evaluated during those bounds. This does not promise rollback after an arbitrary kernel or custom process fails.

Positive finite real timesteps could underflow to zero or overflow to infinity when converted to Float32. The public step and process check now reject those conversions. `test/audit.jl` has fifteen passing focused assertions covering composition preflight and reaction limits.

Preflight also compares the actual Float32 timestep used by the kernels. A Float64 limit can round upward during conversion; comparing the original value accepted a step that a later process rejected only after earlier processes had changed fields. A phase-process regression now checks rejection before preceding deposition mutates the world.

## Remaining findings sent to implementation owners

- Excitable reaction stepping originally had a fixed 0.05 cap independent of gain and recovery parameters. On a seeded 4×4 world, `Excitable(gain=1f8)` accepted dt=0.01 and produced activity magnitude about 25,034 in one step. A new guard bounds the reaction Jacobian across current activity extrema and the full nonnegative-deposit threshold range, combines its maximum row sum with the diffusion rate, and requires dt ≤ 0.8 divided by that combined rate. Existing diffusion and 0.05 caps also apply. This guards local explicit updates; it is not a global nonlinear stability theorem, nor does it guarantee a preceding custom process will keep activity in the preflight range.
- Scene particle options were all converted to Float32; integer `hidden=4` became 4f0 and failed matrix allocation. Integer constructor options are now preserved.
- Cutaway originally removed primary-ray material while shadow integration still included it. The renderer owner fixed shadow cropping.
- Volume export originally recorded one shared shape without checking that all exported grids match. The renderer owner added shape rejection and a regression test.

## Contracts that are internally coherent

Nutrient, deposit and particle storage use per-cell inventory units, so their raw sums form the implemented conserved budget. No additional cell-volume factor belongs in that budget without changing the representation. Particle speed is physical distance per time and converts through grid spacing; the existing refinement test exercises that conversion. Particle exchange itself is not resolution invariant: increasing resolution changes the total initial inventory and the sensed inventory per support. A future concentration interface must convert explicitly at the exchange boundary.

World seed controls explicit initialization helpers, not an implicit random process at each step. Dense stepping is deterministic given initialized fields. Optical and export coordinates deliberately normalize each domain axis to [-1,1]; an anisotropic physical extent is rescaled visually. Native 2D fields are rendered as a slab spanning that normalized depth. These are declared optical interpretations rather than conservation bugs.

## Next architecture extension

Build one encounter between two particle constitutions sharing nutrient and deposit: one secretes scaffold, the other dissolves deposit into its storage and can redeposit it. A single cell-local exchange allocation should gather both populations' uptake/erosion requests, bound them against shared inventories, and commit matching credits. Keep private hidden state with each carrier and persistent IDs. This extends the existing conservative exchange into a reusable named-port boundary and creates feedback with observable consequences; it does not require a generic solver framework or training infrastructure.

Validate the combined grid-plus-carrier inventory, invariance under carrier/particle storage reordering, and a paired feedback-disabled run. A fair shared allocation must avoid making population tuple order decide who gets the last nutrient.

## New process composition checks

`test/process_composition.jl` exercises actual process seams:

- Thermal writing converts nutrient into scaffold, then material sensing changes a neural field's trajectory relative to the same writing with material feedback disabled. Subsequent dissolution changes deposit while irreversible dose persists. Dose is a history variable, not a second material inventory.
- A propagating acoustic pulse writes pressure and excites activity; deposition and a CPU particle carrier spend the same bounded material budget. Attempting to change the acoustic leapfrog timestep is rejected before earlier stages or particle state mutate. The acoustic implementation does not calculate velocity, so this test establishes excitation/exchange composition rather than wave-driven particle advection.
- Electrical growth changes the quasi-static potential as conductor geometry advances. Four events produce identical geometry and potential whether scheduled with dt=0.1 or dt=0.025. The field residual meets its stated tolerance and material remains conserved after an additional deposition stage. This is a lattice breakdown surrogate with a grounded box and an externally powered conductor, not a charge or energy conservation claim.
- Phase interfaces, a heritable particle ecology and deposition conserve both signed phase composition and their joint nutrient/deposit/storage budget while maintaining nonnegative inventories. Invalid mutated ecology controls reject before preceding deposition changes the world.
- Experiment observation at different capture cadences produces identical trajectories, detached snapshots and one initial/final capture. Shared diagnostics include particle and ecology-held inventory; world grid mass alone would incorrectly report uptake as lost material. Health diagnostics include direct process-owned floating arrays.
- A mechanical graph shares the experiment clock while retaining its own reference geometry and loading history. Optical observation is read-only, and malformed graph state rejects before material conversion. This establishes scheduling/state composition; grid deposition does not yet modify bond constitutive properties or repair topology.
- Written dose can be mapped explicitly into a smooth index field. Tracing a frozen ray through it bends toward the written region and increases optical path length; observation leaves world fields, dose and time unchanged. This is phenomenological dose-to-index interpretation, not simulated pulse-scale refractive modification.

The acoustic source, thermal beam, and electrode boundary conditions are explicit external drives. They belong to different numerical methods: fixed-step wave propagation, explicit diffusion/cooling with a prescribed source, and event-scheduled elliptic solves. Their useful shared boundary is named physical state, conservative material transfers, diagnostics, and capability/timestep preflight. The compositor should retain those differences.

A follow-up review exposed a controller-clock seam in Ecology: the original hidden recurrence advanced once per solver step while physical rates used dt. In a frozen one-agent environment with all physical rates zero and fixed recurrent weights, t=0.1 produced hidden=0.19737533 with one step and 0.66135645 with 100 steps. Coupling to a fourth-order phase solver would therefore accelerate the controller merely by requiring smaller numerical timesteps. This is fixed: Ecology now uses a configurable physical response time with relaxation weight `-expm1(-dt/response_time)`. A timestep-refinement regression exercises the isolated controller, and the finite ecology experiment was regenerated.

Metal verification of the earlier reaction guard completed successfully: extrema reductions work on the device, and the dense CPU/Metal benchmark differs by at most 1.49e-8. This does not establish GPU support for the new electrical process or CPU particle exchange.
