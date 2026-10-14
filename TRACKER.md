# Interior — working state

2026-10-14 · ChatGPT integration ledger. The destination is an extensible Julia art simulator for phenomena within—and constituting—materials, across grids, particles, graphs and interfaces, with light transport revealing the resulting structures. The first specimen is a foundation, not a reduction of that ambition.

## Current work

- Julia package with KernelAbstractions grid kernels, excitable dynamics, imposed standing-wave forcing, conservative state-dependent nutrient diffusion and deposition.
- Mobile recurrent particles carrying identity, hidden state and material inventory; deterministic conservative exchange with the grid. CPU path presently.
- Ordered process composition and explicit timestep restrictions. Refusing incompatible inputs before stepping is being audited; no general rollback or adaptive solver is promised.
- Reproducible TOML scenes, raw 3D snapshots, volume export, native straight-ray volume rendering and actual Blender/Cycles pipeline.
- Active work: neural grid process, Metal validation/benchmark, interactive WebGL playback, numerical ablations, visual scouting, Fraser Price research using Scry.
- Original supplied design/viewer preserved. Their historical NumPy simulator source was not in this directory. New measurements refer to Julia only.

## Decisions

- Physical ports are shared; arbitrary latent channels are private to their rule.
- Activity, stored material memory, plasticity and evolution are distinct. Current rules have no training, inheritance or population evolution.
- This prototype uses a periodic computational box. An ellipsoid in the artwork is an optical crop, not a simulated gel boundary.
- No universal generic PDE solver: specialized wave, elliptic, mechanical and event solvers should compose through explicit exchanges.
- Package environments and manifests record tested versions. GPU-capable grid kernels do not imply GPU particle support.
- Artistically, preserve internal interfaces and traces: a fully reacted uniform sphere is an uninteresting final frame, even if numerically valid.

## Next substantial capabilities

1. Compact-support particle neighborhoods with conservative exchange on GPU, followed by persistent bond graphs.
2. Multiple neural constitutions sharing named chemical/mechanical ports, a behavioral atlas with perturbations and temporal descriptors.
3. Real wave propagation or quasi-static potential solve with deposition feedback, replacing prescribed forcing where the artwork needs it.
4. Material-attached transport and deformation, phase separation and moving optical boundaries.
5. Explicit resource-funded birth, mutation and inheritance for evolutionary experiments; local plasticity independently switchable.
6. Optical feedback (absorbed light as a stimulus), refracting interfaces, then specialized birefringence/wave optics where justified.

Each should have a concrete specimen and a numerical counterexample test, not only an interface definition. Read docs/process-atlas.md for candidate compositions and primary sources.
