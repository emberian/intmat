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
- Bounded 32³/64³/128³ CPU/Metal scaling measurements, explicit larger-run profiles and a streaming runner. The full dense128 study completed 3,414 steps to t=10.001953125 in 100.92 seconds of loop time, with 11 captures and maximum relative material drift 3.51e-07. Full 256³ runs remain a future experiment.

## Evidence and boundaries

- The integrated suite passes 1,049 assertions (1,038 in the integrated run, plus 11 additional initializer/migration checks) after the active-arena convergence, including recording fidelity and composition. Small Metal checks cover wave/acoustic source phase and pulsed thermal writing with the new host clock.
- A completed 32³ phase run evolved to t=54: energy 8210.76 → 5549.54, phase standard deviation 0.06894 → 0.66403, relative conserved-mass drift 4.92e-5. No mean correction or clipping was applied.
- The ecology example produced births, deaths, mutation and surviving lineages with explicit material and energy ledgers. These observations do not establish multicellularity or open-ended evolution.
- The 96³ lace study improved spatial artifacts from the coarser recording. Larger optical image dimensions alone do not improve simulation resolution.
- Periodic numerical fields and the optical object's polished surface/crop are separate choices. No deformable enclosing gel boundary is yet solved. Ecology's current optical export is a shallow extrusion of its 2D population.
- Array inventory units conserve within a fixed discretization; particle chemistry still needs explicit cell-volume normalization for cross-resolution physical equivalence.
- Original supplied HTML/design documents remain historical context. Their earlier NumPy source was not supplied. Current measurements refer to this Julia implementation.

## Next substantial work

1. **Persistent active compartments:** stronger funded feedback changes phase topology in a matched five-run study, but coarsening continues. Test more seeds, longer lifetimes, local interfacial agents, explicit conserved surface/volume exchange, and heritable boundary maintenance. Count connected regions in both phase signs.
2. **Ecological differentiation:** multiple resources, heritable metabolic tradeoffs and spatial niches; compare diversity and extinction against simpler controls. Current long runs sustain turnover but lose founder lineages. Mutation within surviving lineages remains substantial.
3. **Evolving organization:** conservative persistent bonds, division and organism-scale behavioral descriptors; investigate cooperative transport and collective structures without imposing their appearance in rendering.
4. **GPU carriers and 3D:** port conservative fair allocation, local neighborhoods and mutable populations while retaining exact identities and ledgers. The continuous neural arena currently runs on CPU; small 3D tests run, but are not evidence of resolution-independent chemistry.
5. **Material deformation:** couple mechanics to transport, damage, growth and healing through documented work/material exchanges. The present spring graph is a carrier, not a solved continuum gel.
6. **Optics as a participant:** absorbed light as a thermal/chemical input; interface refraction alongside graded-index transport; specialized wave/polarization optics where a specimen warrants it.
7. **Long studies:** use the tested bounded checkpoint ring, source-aware restore and single-writer ownership for multi-seed ensembles. Add chunked trajectory output, lineage sampling, and physical cell-volume normalization before comparing resolutions.
8. **Photographic series:** use the semi-implicit solver for longer 3D boundary dynamics and photograph actual evolving interfaces/inclusions across time. Keep simulation resolution distinct from optical image dimensions.

See `docs/process-atlas.md` for the broader artistic/research map, `docs/fraser-experiments.md` for primary-source research, and `docs/bigruns.md` for measured scaling and concrete run profiles. These next directions remain part of the project, rather than being reduced to the initial specimens.

## Active convergence: compartments and continuous arena

Implemented and integrated: `ImplicitPhaseSeparation`, `ActiveMedium`, `NeuralEcology`, and restartable `Arena`. Inhabitants have inherited recurrent controllers, separate bounded lifetime plasticity, fair proportional resource access, catabolic fuel use, paid signaling/deposition/erosion/remodeling, persistent ancestry and mutation. A closed-fuel test demonstrates extinction and a work plateau rather than free-energy recycling.

The constant-mobility semi-implicit phase solver completed a 64³ t=10 run in 19 seconds; measured per-step timings imply about 84× equal-physical-time throughput versus the explicit reference at its stability limit. Mean composition drift stayed below 2e-10. This is not a guarantee of nonlinear stability or timestep-independent trajectories.

A continuous 64² arena is running locally at `http://127.0.0.1:8780/arena.html`, with measured fields, individual selection, genome/plastic-state inspection, live controls and bounded checkpoints. It is an ordinary local process, not an installed reboot service. The curated Pages replay contains 25 actual snapshots with display fields sampled to 32² and displayed agent floats rounded; it is explicitly recorded playback.

Two coarse fair-share runs reached t=1000 with populations 383/398, no blocked births, and continued birth/death turnover, but only 1/3 surviving founder lineages. Feedback ablation also persisted; this does not establish a diversity benefit from feedback. In the boundary sweep, strong funded feedback retained two substantial positive domains at t=160, while both controls had one component in each phase; interface length still decreased. See `docs/arena-experiments.md`, `docs/arena-experiment-results.json`, and `docs/active-compartment-studies.md` for settings and counterexamples.
