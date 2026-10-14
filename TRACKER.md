# Interior — working state

2026-10-14 · ChatGPT integration ledger. The destination is an extensible Julia art simulator for phenomena within—and constituting—materials, across grids, particles, graphs and interfaces, with light transport revealing the resulting structures.

## What runs

- Shared Julia `World`, ordered `Experiment` stages, declared ports, explicit publisher conflicts, preflight validation and per-process timestep bounds. Observations are detached from integration; arbitrary runtime errors do not imply transactional rollback.
- Dense excitable and recurrent neural fields; conservative material conversion and deposit-dependent transport. Metal validation exists for dense, neural, writing, acoustic and phase kernels. Moving particles and evolutionary populations currently run on CPU.
- Persistent recurrent particles, resource-funded inherited neural ecology, and a deformable spring/bond graph with irreversible damage. These carriers retain their own state and identities; ecology is distinct from fixed-rule recurrence.
- Thermal writing and irreversible dose, heterogeneous damped acoustic waves, harmonic-field dielectric branching, conservative Cahn–Hilliard phase separation, deformation/fracture, and graded-index ray observation. Each has executable examples and focused numerical checks.
- Actual coupled ablations: writing changes neural deposition; acoustic excitation reaches particle exchange; phase interfaces coexist with finite-material ecology; mechanical history remains distinct within shared experiments.
- Schema-v2 immutable state observations preserve native floating precision, owned process arrays, exact IDs, source hashes and parameters. This is an observation format, not yet an executable restart format.
- Offline WebGL2 recorded-volume viewers, straight-ray previews, and Blender/Cycles photographs of exported fields, actual interfaces and surviving graph geometry. Finished studies include an amber branching conductor, a phase-interface labyrinth, a written inclusion and a fractured network.
- Curated GitHub Pages collection and automated deployment from `dev`. Generated raw output and third-party reference footage stay outside Git. Snapshot commits are pushed as work converges.
- Bounded 32³/64³/128³ CPU/Metal scaling measurements, explicit larger-run profiles and a streaming runner. Full 256³ runs remain a future experiment; small-step benchmarks are not completed long runs.

## Evidence and boundaries

- The integrated suite passed 638 assertions before the long-run clock precision change; subsequent checks are recorded in the work session and test output.
- A completed 32³ phase run evolved to t=54: energy 8210.76 → 5549.54, phase standard deviation 0.06894 → 0.66403, relative conserved-mass drift 4.92e-5. No mean correction or clipping was applied.
- The ecology example produced births, deaths, mutation and surviving lineages with explicit material and energy ledgers. These observations do not establish multicellularity or open-ended evolution.
- The 96³ lace study improved spatial artifacts from the coarser recording. Larger optical image dimensions alone do not improve simulation resolution.
- Periodic numerical fields and the optical object's polished surface/crop are separate choices. No deformable enclosing gel boundary is yet solved. Ecology's current optical export is a shallow extrusion of its 2D population.
- Array inventory units conserve within a fixed discretization; particle chemistry still needs explicit cell-volume normalization for cross-resolution physical equivalence.
- Original supplied HTML/design documents remain historical context. Their earlier NumPy source was not supplied. Current measurements refer to this Julia implementation.

## Next substantial work

1. **Active compartments:** combine phase boundaries, neural activity, permeability, resource-funded populations and optical interfaces in one 3D specimen. Establish explicit conservative surface/volume exchanges and compare isolated versus coupled runs.
2. **Long phase runs:** add a semi-implicit/spectral or equivalent tested integrator. With fixed physical extent/interface parameters, explicit phase timestep scales roughly as h⁴; increasing cubic resolution therefore grows work roughly as n⁷. Preserve the material law when optimizing it.
3. **More neural constitutions:** add local plasticity separately from inheritance; behavioral descriptors and matched perturbation sweeps; investigate organism-scale structure without programming it into rendering.
4. **GPU carriers:** compact-support neighborhoods and conservative particle-grid exchange; then persistent evolving bonds and 3D ecology. Keep exact identity and material ledgers through birth/death/reordering.
5. **Material deformation:** couple mechanics to transport, damage, growth and healing through documented work/material exchanges. The present spring graph is a carrier, not a solved continuum gel.
6. **Optics as a participant:** absorbed light as a thermal/chemical input; interface refraction alongside graded-index transport; specialized wave/polarization optics where a specimen warrants it.
7. **Sustained experiments:** tested restart loader, interrupted-run recovery, multi-seed ensembles, compressed/chunked state output and temporal morphology metrics. Make expensive runs recoverable before routinely running overnight.
8. **Photographic series:** track the same evolving specimen across time, material mappings and camera views. Export a small curated public subset; retain raw state and scene provenance locally.

See `docs/process-atlas.md` for the broader artistic/research map, `docs/fraser-experiments.md` for primary-source research, and `docs/bigruns.md` for measured scaling and concrete run profiles. These next directions remain part of the project, rather than being reduced to the initial specimens.
