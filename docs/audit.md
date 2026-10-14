# Numerical and composition audit

Reviewed the Julia implementation on 14 October 2026. Historical Python claims were not used as evidence. Scope: dense fields, particle exchange, scene construction, rendering, and the example runner.

## Reproduced and fixed

`step!(w,(Deposition(),Diffusion(:missing,0.01)),0.01)` previously changed deposit before throwing, while time stayed zero. Process preflight now dispatches on the world as well as the process. Invalid diffusion field names and incompatible particle storage/backend, malformed mutable particle state, malformed rules, and invalid inventories are rejected before the composition runs. This checks known preconditions; it does not promise rollback after an arbitrary kernel or custom process fails.

Positive finite real timesteps could underflow to zero or overflow to infinity when converted to Float32. The public step and process check now reject those conversions. `test/audit.jl` has fifteen passing focused assertions covering composition preflight and reaction limits.

## Remaining findings sent to implementation owners

- Excitable reaction stepping originally had a fixed 0.05 cap independent of gain and recovery parameters. On a seeded 4×4 world, `Excitable(gain=1f8)` accepted dt=0.01 and produced activity magnitude about 25,034 in one step. A new guard bounds the reaction Jacobian across current activity extrema and the full nonnegative-deposit threshold range, combines its maximum row sum with the diffusion rate, and requires dt ≤ 0.8 divided by that combined rate. Existing diffusion and 0.05 caps also apply. This guards local explicit updates; it is not a global nonlinear stability theorem, nor does it guarantee a preceding custom process will keep activity in the preflight range.
- Scene particle options were all converted to Float32; integer `hidden=4` became 4f0 and failed matrix allocation. The scene owner is preserving integer constructor options.
- Cutaway originally removed primary-ray material while shadow integration still included it. The renderer owner fixed shadow cropping.
- Volume export originally recorded one shared shape without checking that all exported grids match. The renderer owner added shape rejection and a regression test.

## Contracts that are internally coherent

Nutrient, deposit and particle storage use per-cell inventory units, so their raw sums form the implemented conserved budget. No additional cell-volume factor belongs in that budget without changing the representation. Particle speed is physical distance per time and converts through grid spacing; the existing refinement test exercises that conversion. Particle exchange itself is not resolution invariant: increasing resolution changes the total initial inventory and the sensed inventory per support. A future concentration interface must convert explicitly at the exchange boundary.

World seed controls explicit initialization helpers, not an implicit random process at each step. Dense stepping is deterministic given initialized fields. Optical and export coordinates deliberately normalize each domain axis to [-1,1]; an anisotropic physical extent is rescaled visually. Native 2D fields are rendered as a slab spanning that normalized depth. These are declared optical interpretations rather than conservation bugs.

## Next architecture extension

Build one encounter between two particle constitutions sharing nutrient and deposit: one secretes scaffold, the other dissolves deposit into its storage and can redeposit it. A single cell-local exchange allocation should gather both populations' uptake/erosion requests, bound them against shared inventories, and commit matching credits. Keep private hidden state with each carrier and persistent IDs. This extends the existing conservative exchange into a reusable named-port boundary and creates feedback with observable consequences; it does not require a generic solver framework or training infrastructure.

Validate the combined grid-plus-carrier inventory, invariance under carrier/particle storage reordering, and a paired feedback-disabled run. A fair shared allocation must avoid making population tuple order decide who gets the last nutrient.
