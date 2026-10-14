# Interior: a process atlas

Design and research notes, 14 October 2026. The compositions below are proposals. They do not claim that their solvers, training systems, or optical effects are implemented. The Julia core under construction tests grid processes, moving particles, conservative exchange, and volume observation; it is a foundation for this atlas rather than the atlas itself.

`interior-living-materials.md` describes a historical NumPy version 0.0.2, including measurements and tests. That source tree and its verification files were absent from this checkout when this document was written. Treat those paragraphs as historical reports, not verification of the present Julia implementation. The supplied HTML viewer is an artifact, not a substitute for source or reproducible tests.

## Artistic direction

The strongest image here is a material making its own interior: a rhythm becoming mineral bands, a travelling excitation becoming a scar, or a colony building an optical obstacle that changes its next decision. Transparent amber is useful because it makes depth legible, but it should not become a mandatory aquarium around familiar creatures. Opacity, cavities, skins, phase boundaries, and damage should emerge from state.

I would begin with a choir that deposits matter, then a moving population that erodes and repairs those deposits. This creates a consequential encounter between two carriers with a small solver budget. Electrical trees and acoustic writing are excellent subsequent pieces. Full cavitation and wave-optical nanostructure are expensive destinations; decorative glow should not quietly stand in for their physics.

## The exchange contract

Every composition uses named physical ports rather than sharing arbitrary neural hidden coordinates. A process owns its state and proposes a flux, rate, stress, or event. The integrator commits the admissible change. A minimal contract records:

- Quantity, units or declared nondimensionalization, carrier, support, coordinate frame, and ownership.
- Sampling and deposition weights, boundary behavior, sign convention, and source/sink counterpart.
- Update interval, lag, stability restriction, and whether feedback needs iteration.
- Positivity or resource limits, conserved totals, dissipation or external power, and diagnostic residuals.
- Event identity, rule provenance, lineage, and complete rollback state when steps can be rejected.

Particle uptake must subtract exactly the amount credited to particle inventory, including cell volume in the grid conversion. When several particles compete, allocate from the shared available inventory before applying uptake. Independent per-particle clipping can overdraw a cell. Scatter atomics may introduce floating-point order differences; deterministic gather or sorted segmented reduction costs more but is useful for reference runs. Material and charge budgets are different contracts. Preserving one does not establish the other.

## Seven compositions

### 1. Choir in amber

Local excitable or phase oscillators modulate precipitation from a diffusing feedstock. Deposits lower mobility and shift oscillator coupling, so the choir writes a score that alters its own acoustics of communication. The first artwork should show long periods of almost invisible activity interrupted by a mineral front, rather than continuous rainbow boiling.

**State and ports:** grid concentration `c`, immobile deposit `m`, activity `z`; optionally particle oscillators with individual phase and inventory. Activity requests a nonnegative reaction extent, bounded by available `c`; reaction transfers `c → m`. Deposit returns mobility and sensed alteration. Oscillator dynamics may be invented, but chemical bookkeeping is explicit.

**Numerics and project:** start with 2D no-flux finite-volume diffusion and a small fixed recurrent rule, then the same law in a shallow 3D slab. Explicit diffusion has a step restriction proportional to cell-width squared. Resolve fronts before interpreting their thickness. Compare coupled, output-disabled, and feedback-disabled runs with identical forcing. Then add moving oscillators without erasing their identity during reordering.

### 2. A colony that carries its constitution

A particle population moves through gel, senses soluble resources and nearby deposits, and secretes a trail that can become a porous scaffold. Different populations may share the same nutrient while responding differently to the scaffold. Their seam is more interesting than guaranteed combat or cooperation.

**State and ports:** positions, persistent identifiers, inventory and local neural state on particles; concentration and deposit on the grid. A rule requests velocity or force, uptake, and secretion. Pairwise forces need their own momentum convention; an imposed velocity is a kinematic active model. Nutrient enters particle inventory and secretion spends it. Named observations permit distinct controllers without averaging hidden states or weights.

**Numerics and project:** implement a bounded moving-particle study first, then spatial bins and compact-support neighborhoods. Naive all-pairs perception is quadratic; a bounded-density cell list approaches linear work but dense clusters remain expensive. [Neural Particle Automata](https://selforg-npa.github.io/) updates both continuous positions and internal states with a shared learned rule and SPH perception; its custom CUDA training kernels are a precedent, not a drop-in portable Julia implementation. [Flow-Lenia](https://arxiv.org/abs/2506.08569) is an alternative conserved grid formulation with localized parameters. Neither establishes conservation for our independently composed processes.

### 3. Sound writes a slowly changing instrument

An acoustic field drives local deposition or material rearrangement; evolving density and stiffness then retune the field. The specimen becomes an instrument whose resonant geometry is also its memory. First use a slowly varying intensity envelope and acknowledge the approximation. Later show how a node moves when its surrounding material changes.

**State and ports:** pressure/velocity on a staggered grid, density and acoustic response from material state; cycle-averaged intensity or strain drives a slow reaction. Reaction supplies updated constitutive coefficients at a declared refresh interval. Sound power enters through a boundary or source. An intensity-gated reaction is an artistic coupling until calibrated to an actual mechanism.

**Numerics and project:** homogeneous wave propagation and reflection test, then heterogeneous standing waves, then slow feedback. For conventional explicit wave stepping, the step scales with cell width over sound speed; the shortest relevant wavelength controls resolution. Subcycle waves or solve a harmonic field between material updates. [k-Wave's formulation](https://k-wave.org/documentation/kspaceFirstOrder2D.php) provides a verified reference for staggered pseudospectral acoustics; its [stability documentation](https://www.k-wave.org/documentation/checkStability.php) shows why its k-space correction has different restrictions from a simple finite-difference solver. Those guarantees cannot be transferred to our kernel.

### 4. Lightning cultivates a garden

Electrical channels branch inside a cloudy dielectric, then become sites where a second chemistry deposits reflective matter. A subsequent pulse sees the conductive history of earlier pulses. This can produce delicate trees embedded in a changing material rather than a branching animation drawn on top.

**State and ports:** potential and conductivity fields, graph or voxel channel occupancy, soluble precursor, deposited material. A stochastic dielectric-breakdown proposal grows channels from local field strength. Channels return conductivity and reaction support. Electrodeposition is a separate alternative: ion transport, potential, and interfacial kinetics transfer dissolved ions to solid while boundary current supplies charge.

**Numerics and project:** first solve the elliptic potential field around a fixed channel, then add growth events, then deposition feedback. Repeated potential solves dominate; a graph-only random tree misses field screening. The original [dielectric breakdown model](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.52.1033) is a stochastic macroscopic fractal model, not a resolved plasma discharge. [Cogswell's electrodeposition model](https://arxiv.org/abs/1411.6615) couples electrochemical interfacial kinetics to phase-field growth. A glowing dielectric tree and a metal dendrite may look related while requiring different state and balances.

### 5. A wound learns its own repair

A viscoelastic specimen fractures; active particles colonize crack faces and deposit a softer repair phase. The repaired region redirects the next fracture, preserving a visible history of damage. This is a compelling eventual extension because memory has mechanical consequences.

**State and ports:** displacement or velocity, stress and constitutive history; scalar fracture damage; particle inventory; repair phase. Elastic energy drives damage, damage changes stiffness, damaged support recruits particles, and bounded deposition changes the constitutive mixture. Damage irreversibility must be separated from new repair material: lowering the damage field alone can erase history and create artificial healing energy.

**Numerics and project:** fixed damaged geometry and diffusion first, then an elastic equilibrium solver, then a phase-field fracture model. The [phase-field fracture review](https://arxiv.org/abs/2309.03996) explains diffuse cracks without explicit crack-front tracking. The regularization length must be resolved and alters the problem; coupled mechanics and damage usually require iterative solves. This is substantially beyond a local diffusion kernel. Begin in 2D and make repair physically consequential before rendering a full block.

### 6. The beam's writing becomes the next beam's lens

A moving optical focus changes absorption, index, or porosity. Later illumination travels through the written structure and stimulates a mobile active phase. A beam should have a finite focus and a dose history; a cursor painting arbitrary bright voxels loses the material premise.

**State and ports:** absorbed dose, temperature or explicit modification variable, index and absorption, plus active state. A prescribed focus deposits energy; modification changes optical response. A first approximation uses straight-ray attenuation and a declared phenomenological dose response. Continuous index gradients, scattering, polarization, and coherent interference require progressively different transport models.

**Numerics and project:** dose accumulation and persistent modification, then light attenuation feedback, then resolved interfaces or graded-index paths. [Ultrafast laser writing in silica](https://eprints.soton.ac.uk/477361/1/Laser_Photonics_Reviews_2023_Lei_Ultrafast_Laser_Writing_in_Different_Types_of_Silica_Glass.pdf) documents material-dependent modifications driven by nonlinear absorption. Index change, birefringent nanostructure, and void formation are not one universal engraving operation. A macroscopic voxel rule may be artistically useful without resolving the laser pulse or nanoscale mechanism.

### 7. A constellation of collapsing cavities

Bubbles in a liquid pocket or soft specimen rhythmically expand and contract. Some flashes become chemically consequential seeds, leaving pale rings or nucleation sites after the light has gone. Sparse flashes against a dark translucent volume would be stronger than assigning every bubble continuous emission.

**State and ports:** discrete bubble center, radius and radial velocity, gas content; ambient pressure, dissolved gas and slow reaction fields. Bubble dynamics reads pressure and exchanges gas with the carrier. Emission is a separate calibrated or explicitly invented closure; chemical activation consumes a declared photochemical dose.

**Numerics and project:** isolated radial bubble under imposed pressure, then sparse bubbles sampling a wave field, then feedback and moving cavities. Collapse is stiff; substepping and shape-stability limits matter. [Single-bubble sonoluminescence](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.74.425) concerns light from strongly collapsing acoustically driven gas bubbles. It is transient light emission, not subsurface laser engraving or permanent solid modification. A spherical radial model cannot resolve jets, interacting interfaces, or violent collapse in a gel; a threshold flash is an illustrative model, not simulated sonoluminescence.

## Neural rules without premature training infrastructure

[Growing NCA](https://distill.pub/2020/growing-ca/) demonstrated local learned rules trained for growth, persistence and regeneration. Its hidden channels are local computational state. A fixed trained rule need not learn online, and attractive regeneration does not confer a material budget. Interior can initially use a small hand-designed recurrent rule behind the same semantic observation/action interface, then substitute trained rules once the composition is numerically meaningful.

Keep physical variables separate from hidden channels. Training should reward observable behavior across perturbations and environments, not merely a final RGB image. Runtime inference, offline training, local plasticity, and reproductive inheritance are four different projects. Begin with inference; do not require reverse differentiation through adaptive events and particle sorting to make the first artwork.

## Julia and GPU implementation path

[KernelAbstractions](https://juliagpu.github.io/KernelAbstractions.jl/) supports writing kernels for different backends; portability of code does not promise identical performance or floating-point trajectories. Its documentation spans an important version boundary: 0.10 changes the CPU backend to PoCL and reorganizes backend interfaces. Pin the actually tested package versions and read documentation for that version. Do not upgrade a working 0.9 kernel by copying unversioned examples.

[Metal.jl](https://metal.juliagpu.org/stable/) is the local Apple GPU path, while [CUDA.jl](https://juliagpu.org/backends/cuda/index.html) is the NVIDIA path. Float32 is a sensible shared artwork representation. Metal's [`mtl` adaptor](https://metal.juliagpu.org/stable/api/array/) can change floating-point element types; make precision an explicit constructor choice and retain a CPU high-precision reference where needed. Runtime/device capability checks and an actual kernel launch establish availability; package installation alone does not.

Use structure-of-arrays particle storage, preallocated double buffers, explicit Float32 literals in device kernels, and keep intermediate fields on the device. Separate synchronization for measurement from rendering cadence. Compile warm-up should be reported separately from steady-state kernel time. A portable dense stencil is the first shared abstraction; FFTs, sparse linear solves, device sorting, automatic differentiation, and interoperability with graphics all need individual backend verification. KA alone supplies none of those algorithms.

Budget before scaling. One Float32 field at 128³ costs 8 MiB; twenty fields cost 160 MiB before scratch arrays, neural activations, particle bins, checkpoints, and rendering buffers. At 256³ those twenty fields cost 1.25 GiB. Doubling resolution in 3D multiplies cell count by eight and can multiply explicit diffusion steps by four for the same simulated duration: roughly 32 times the stencil work. Wave stepping instead adds approximately a factor of two in steps. These are scaling estimates, not measured performance.

Rendering should consume derived extinction, scattering and emission fields, not overwrite simulation channels. A fast emission/absorption ray marcher is an excellent early observation tool, but scattering, refraction and birefringence remain separate milestones. A 2D section is also a valid artwork and a much faster laboratory; a 3D image must reveal occlusion and interior structure to justify its extra simulation cost.

## Staging toward the larger toolbox

1. **Exchange specimen:** moving stateful particles and grid feedstock with one bounded uptake/deposition law, conserved inventory, deterministic reference cases, and an honest volume view.
2. **Choir and contact:** oscillatory deposition plus two controllers sharing named physical ports; perturbation experiments and a curated parameter atlas. No claim of evolution.
3. **Fields with consequences:** one acoustic-envelope or quasi-static electrical process, coupled both ways to deposited microstructure; compare lagged coupling with refined update intervals.
4. **Moving material:** advection or deformation with material-attached memory. Add a representation because an artwork requires it, not because every process must use every carrier.
5. **Learned constitutions:** train a rule against physical observables and robustness once carrier transfers are reliable; preserve inference provenance and keep online adaptation optional.
6. **Research-grade pieces:** fracture/repair, bubble collapse, and optical feedback acquire specialized solvers and validation. Each should earn its own resolution and approximation study.

The toolbox becomes simpler when exchange, state ownership, observation, and reproducible experiments are shared. It becomes less capable when all mechanisms are flattened into one generic update callback with no units, no time-scale contract, and no notion of an elliptic solve or discrete event. The architecture should share bookkeeping while allowing those numerical differences to remain real.
