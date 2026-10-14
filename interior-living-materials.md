# Interior / living materials

**Design revision, 14 October 2026.** This is a proposed architecture and artistic direction, followed by an exact account of the executable study. It does not claim that the broad architecture is already implemented.

## The material participates

The project is not restricted to neural creatures inside transparent matter. It should also admit neural processes that remodel the matter, neural dynamics that constitute an active material, and activity that propagates over interfaces and networks created by other processes. The outer specimen need not be a glass capsule: it can deform, dissolve, phase-separate, become opaque, grow a skin, or cease to have a useful boundary.

I want the internal computation to become material structure: a rhythm that becomes a sequence of mineral bands, cooperation that becomes a load-bearing membrane, a chemical exchange that becomes an optical interface. These are artistic model proposals, not claims that an arbitrary neural network has biological properties.

There is no need to prescribe eyes, hunger, fighting, or reproduction for every material. A nonreproducing excitable medium, a stateful particle assembly, an evolving population, and a passive crystal should all be admissible. The neural part is another family of processes, not a replacement for the process toolbox.

## Four distinct kinds of change

**Dynamical activity:** positions, activations, concentrations, and phases change under fixed rules. Persistent oscillation does not establish learning.

**Material memory:** plastic strain, deposited material, topology, stored orientation, and other historical variables change the subsequent response.

**Plasticity:** a rule's parameters change during a run in response to experience or a local update law.

**Evolution:** inheritable differences affect reproduction and persistence. Mutation at reproduction differs from independent random perturbation of a controller every tick. A simulation that provides a split action supplies that reproduction primitive, even when subsequent population behavior is emergent.

These should be independent switches and independent data structures. Training a rule offline should not force plasticity or evolution during the artwork's execution. A fixed rule can produce dynamic art; an evolutionary run can also become visually monotonous.

## A material definition is more than a shader

| Component | Examples of controls | Required interpretation |
|---|---|---|
| Carrier and microstructure | Fluid, elastic network, viscoelastic continuum, porous matrix, grains, multiple phases | State representation, coordinate frame, and approximation regime |
| Passive response | Elastic moduli, viscosity, surface tension, transport coefficients, solubility, heat capacity | Values may be spatial fields and state-dependent functions, not only constants |
| Active response | Excitability, recurrence, sensing range, signal delay, latent state size, contractility | Distinguish an invented local rule from a calibrated physical surrogate |
| Resources | Feedstock, chemical fuel, external heating, light, dissolved gas | Define what is actually conserved and where sources and sinks reside |
| Adaptation | Shared rule, per-species rules, local plasticity, birth/mutation/inheritance | Define provenance and transfer of identity; do not silently average arbitrary neural weights |
| Environment and score | Boundaries, nutrient supply, pressure, sound, illumination, periodic or irregular forcing | A changing environment is part of the specimen, not merely a render setting |
| Optical response | Spectral absorption and scattering, phase function, index, emission, interfaces, orientation | Optical properties derive from material state; diagnostic RGB is not a physical material |

Useful exploration coordinates include the ratio of material-remodeling time to neural-response time, resource-replenishment time to consumption time, and interaction length to specimen size. They help compare regimes without treating every numeric slider as independent. These ratios are design coordinates, not a universal similarity theory for all included processes.

A neural policy should usually request admissible fluxes, reaction rates, active stresses, or bond events. It should not overwrite density or create material for free. A positive mobility can be parameterized positively; a reaction extent can be limited by available species; a passive constitutive branch may use a free energy with a separately modeled active power input. For an intentionally impossible material, we may relax such constraints, but the exception should be explicit. Persistent activity in an invented rule is not evidence that its mechanical work has been energetically paid for.

## State, rule, and geometry are separate

A conceptual specimen state is `(q, z, theta, G)`: physical/material state, internal activity state, rule parameters, and interaction geometry/topology. A neural rule reads local observations and proposes changes; another solver integrates physical constraints. This notation is not a demand that all processes be additive ODEs: algebraic constraints, PDEs, discrete events, and implicit solves remain different numerical objects.

Representations must coexist:

- Eulerian grids for continuous fields and wave/transport solvers.
- Particles or material points for moving internal state, inclusions, or active matter.
- Graphs for bonds, membranes, branches, and changing interaction networks.
- Surfaces and interfaces for cavities, boundaries, thin films, and optical discontinuities.

Material-attached state should move with the material. A fluid's neural activity may require both advection and local updating. A coating may evolve on an interface that another solver is deforming. A particle's lineage must survive a spatial reorder. Transfer between grid, particle, and surface domains is itself an explicit modeled operation.

**Compatibility of physical ports does not imply compatibility of latent coordinates.** Two separately trained networks can use their hidden channels in unrelated ways, even with identical architecture. Exchanging a concentration or force has defined semantics; mixing `hidden[3]` or averaging arbitrary weights does not. Shared-rule models, rule-conditioned models, and explicit adapters are possible strategies. The material interface should not require all neural processes to share a latent basis.

## Composition is causal and numerical

A useful loop might be:

`local activity -> secretion -> phase separation -> transport resistance -> local observations -> activity`.

Another could be:

`activity -> bond remodeling -> deformation -> acoustic field -> activity`.

The optical loop is particularly attractive:

`activity -> microstructure -> light field -> absorbed stimulus -> activity`.

These are candidate model networks. Each edge needs a law and an exchange contract; a diagram is not a solver. Strongly coupled components may require iteration or a joint solve. Acoustic cycles, controller updates, mechanical relaxation, and slow mineral growth cannot automatically share one stable time step.

The contract should record inputs/outputs, units, support, frame, validity, transfers, time-step restrictions, state ownership, and validation checks. Rollback must also restore internal activity, evolving rules, events, and RNG state. For stochastic processes, merely seeding a generator by step number does not guarantee equivalence when adaptive time steps change.

## Research seeds

The following are complementary precedents, not interchangeable engines and not an identification of the software in Fraser Price's clips.

1. **Growing Neural Cellular Automata** and **Self-Organising Textures** demonstrate local learned rules for growth/regeneration and dynamically evolving textures. The behavior of a fixed trained rule is distinct from online learning. [Growing NCA](https://distill.pub/2020/growing-ca/) · [Textures](https://distill.pub/selforg/2021/textures/)
2. **Neural Particle Automata, SIGGRAPH 2026** evolves continuous particle positions and internal states through learned local rules, using SPH-based perception. This is particularly relevant to making activity travel with matter. [Paper](https://arxiv.org/abs/2601.16096) · [Authors' demo](https://selforg-npa.github.io/)
3. **Pandora: Stateful Particle Lenia**, updated July 2026, attaches a changing state vector to particles, which modulates their interactions. Its author supplies WebGPU and CUDA/PyTorch implementations. It is an experimental project, not a validation of arbitrary general claims about the inevitability of life. [Technical article](https://pajouheshgar.github.io/blog/posts/pandora-artificial-life/)
4. **Flow-Lenia** adds mass conservation and localized rule parameters, enabling multiple kinds of patterns in one world. This is a useful precedent for representing different constitutions without giving every species an isolated universe. [2025 paper](https://arxiv.org/abs/2506.08569)
5. **MeshNCA** evolves dynamic textures on mesh vertices and can produce multiple appearance attributes. It is relevant to skins, membranes, and internal boundaries, without making the entire volume a mesh. [Project](https://meshnca.github.io/)
6. **From Cells to Pixels, SIGGRAPH 2026** separates a coarse evolving NCA from a local implicit appearance decoder. This suggests a useful rendering-resolution boundary. Decoded fine appearance is not automatically physically simulated fine structure. [Project](https://cells2pixels.github.io/)
7. **ALIEN** provides a closely related particle-based artificial-life environment with neural control and specialized cell functions. Its existence is useful implementation context; this project has not established it as the source of the supplied clips. [Author's description](https://chrxh.itch.io/alien)

Reuse of any external code or pretrained weights requires examining its own license. References here are not bundled implementations.

## Observation should reveal structure

Diagnostic modes should expose activity, concentrations, pressure, stresses, rule identity, and causal contributions. The artwork view should instead transport light through geometry and constitutive optical properties. A colored neural state plot can be informative without being the desired final image.

Blender can consume OpenVDB volume data and render absorption/scattering/emission, together with meshes for sharp interfaces. Its official documentation describes these capabilities: [volume objects](https://docs.blender.org/manual/en/latest/modeling/volumes/introduction.html), [volume materials](https://docs.blender.org/manual/en/latest/render/materials/components/volume.html). Mitsuba provides heterogeneous media and a research-oriented light-transport environment: [participating media](https://mitsuba.readthedocs.io/en/stable/src/generated/plugins_media.html).

Neither exporting a density grid nor selecting a renderer automatically implements stress birefringence, thin-film interference, continuously varying-index ray paths, or wave-scale structural coloration. Those require the relevant optical model and solver. Fine appearance synthesized by a decoder may remain render-only; feeding it back into transport needs an explicit resolved or homogenized coupling.

## A family of artworks, not one effect

**A choir that writes its score into matter.** Coupled local oscillators gate deposition in a yielding matrix. Deposits change communication and transport; rhythms can outlast their original excitation as a structural record. Initially this can be an invented coupled law, not a simulation of a known organism.

**A material that remodels its illumination.** Activity changes an optical microstructure. Light redistributed by that microstructure stimulates subsequent change. This is an ambitious light/material feedback experiment, not just a luminous shader.

**Two constitutions meeting.** Two active materials share a transport field and mechanical boundary but have different local rules. Their contact could produce a persistent seam, synchronization, segregation, or dissolution. The goal is to discover the interaction, not guarantee cooperation or anthropomorphize the result.

I would curate a behavioral atlas rather than optimize a single "interestingness" score. Candidate descriptors include temporal spectrum, coherence length, interface area, material turnover, persistence of memory, spatial localization, and recovery after a perturbation. Noise or population explosion should not win just because it produces more activity. Selection should use multiple seeds and perturbations, with rule provenance and simulation settings retained.

## What version 0.0.2 actually implements

The executable extension is still a small 3D fixed-grid NumPy coupling test. It adds six bounded neural state fields, an immutable recurrent rule, local neural-state diffusion, a neural modulation of an abstract reaction rate, and a material-alteration input to the neural rule. The carrier parameters are separate JSON-configurable transport, thermal, reaction, and alteration coefficients. The initial mask and chemical arrangement are still the existing ellipsoidal toy scene.

The default network contains three deliberately supplied leaky-tanh oscillator pairs plus seeded random perturbations. A separate random-rule initialization is available. Neither is trained or evolved. There is no reproduction, selection, motion, mechanical constitutive solver, learned surrogate, or energetic cost for latent computation. Hidden state is part of the explicit world state, so failed steps do not secretly advance an RNN or RNG.

The experiment runs three versions with identical initial data and weights: all connections; removal of alteration-to-neural sensing only; and removal of neural modulation of reaction. It exports full final 3D arrays and weights, plus an interactive offline HTML viewer of central slices. The viewer is a diagnostic, not a photorealistic render. Only diagnostic pixel encoding is clipped/quantized; the simulation does not silently clip invalid physical state.

At 20 cubed cells and final nondimensional time 8, the removed-sensor run produced about 1.70% less integrated reaction product, and the removed-neural-output run about 1.86% less, than the coupled run. These are measurements of the supplied toy laws, not physical predictions or measures of emergent intelligence. Time-step refinement and residuals are recorded under `examples/living-verification.json`; pointwise neural activity changes more than the integrated chemical observable. Twenty-four tests passed in the supplied environment. The viewer controls were checked in Chromium.

The next architecture milestone should introduce a second state carrier and meaningful two-way exchange, rather than indefinitely enlarging this dense-grid example. Neural state moving with particles or a deforming interface, while exchanging a conserved species with a grid, would test something the current kernel genuinely does not establish.
