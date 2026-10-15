# A weather archive

An electrical tree heats the material it crosses. Heat spreads and cools; exposure
converts a local nutrient reservoir into persistent mineral. That dose changes
where the next branch can grow. Three deliberately placed, unequal electrode
stubs share one screened potential field, so all three changing conductors affect
later frontier choices.

This is an invented coupled material, not a calibrated discharge, dielectric
breakdown or amber chemistry model. The roots and chamfered susceptibility envelope
are designed geometry. The growing conductor topology and thermal/mineral fields
are simulated in Julia. Optical amber, metal, emission and the observation chamber
are chosen interpretations.

## Implementation

`studies/electrothermal_archive.jl` combines the existing dielectric-growth
multigrid study with `Interior.thermal_exposure`, the same law used by
`ThermalWriter`. Each event:

1. Gathers the current potential, temperature, dose and nutrient at the conductor
   frontier.
2. Selects a face-connected cell using the field drop raised to `eta`, multiplied
   by local nutrient, a bounded dose-dependent susceptibility, and the designed
   specimen envelope. The random variate is the existing deterministic electrical
   event hash.
3. Converts 0.25 nutrient to conductor mineral at that cell and holds it at unit
   potential. The outside electrical box remains grounded.
4. Deposits a Gaussian heat source at the selected site, diffuses/cools
   temperature, accumulates thresholded exposure, and conservatively converts
   nutrient to additional mineral. Heat uses the existing periodic thermal
   stencil; the active specimen is well inside that boundary.
5. Re-solves the harmonic potential for the changed conductor. Saved captures
   also check the independent core residual.

The dose multiplier is `exp(clamp(memory_gain*dose, 0, 3))`; mineral conversion
uses the positive exponential transfer `food * -expm1(-conversion*exposure)`.
Food below the conductor deposit threshold cannot support another electrical
occupation. The explicit thermal positivity bound is checked before evolution.

The study owns five scalar material/thermal arrays rather than allocating a full
World. It preflights a conservative host/device working-set estimate. Gzip
captures preserve native Float32 values exactly; there is no spatial resampling.
Each capture includes exact occupied-cell positions, face edges and birth events,
plus temperature, dose and mineral. These are visualization captures rather than
complete restart checkpoints; final nutrient is also saved.

## Validation and production

`benchmarks/2026-10-15-amber-archive-verification.json` records a 33³, 180-event
coupled fixture. CPU and Metal chose the same 180 cells; maximum field differences
were 2.39e-7, relative material error 1.09e-10, and the independent final harmonic
residual was 1.59e-7. Disabling dose feedback first changed the trajectory at event
7. The fixture scales heat width and diffusion to its coarse grid; it checks the
coupled implementation, not resolution-invariant morphology.

Native production was launched from a standalone source snapshot at
`outputs/art/runtime-amber-archive-20261015`, with all Julia source, project/manifest
files and the recipe hashed before launch. The intended 513³ run has 4096 events,
`dt=0.005`, and 65 captures across t=0–20.48. Its approximately 8.6 GiB preflight
fits a 12 GiB study budget on the 96 GiB local machine. Early solves took 3.2–4.1 s
per event at 4–5 V-cycles; these are observed production timings, not a controlled
benchmark. Production status and actual final timing are in the output run report
when complete. No previous run was migrated, reseeded or replaced.

## Optical interpretation

`scripts/render_electrothermal_archive.py` photographs saved states at a fixed
camera. Conductors use the exact saved graph. Two corner-cutting passes only
soften degree-two centerlines within a voxel; endpoints and junctions are retained,
and no branches or connectivity are invented. Chosen radii make the network
legible. Copper marks older graph edges, while cyan emission marks recent edges;
this is not a solved luminescence law.

Two supported views are being developed from the same captures: a faceted amber
observation chamber, and a carrier-free mineral envelope that allows the actual
thermal writing to constitute the object's surface. The latter is a thresholded
saved-mineral surface. The volume alternative maps saved mineral to scattering
and saved temperature to blue emission. Each rendered image records its exact
thresholds, optics, camera, source-run hash and renderer hash. Display choices do
not feed back into the simulated dynamics.

### Native growth reach

Event count is not a resolution-independent measure of development: each event adds one occupied cell. At native 513³ event 704, new conductors reached at most 0.0573 physical units from the designed roots; the 129³ pilot at 1800 events reached 0.4971. The 4096-event native run therefore does not promise the pilot's physical reach or contact between roots. Its shared harmonic field and local thermal/mineral feedback remain active. A close interior photograph may be the appropriate view of this run. A separately initialized longer 257³ study with closer roots is now running for inter-root thermal encounters; see `docs/amber-confluence.md`. It is a new composition, not a migration or continuation of the pinned 513³ trajectory.
