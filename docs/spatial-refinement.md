# Spatial refinement

2026-10-15. Native simulation cells and photograph pixels are separate budgets.

Four times the cells **along each axis** means 64 times the volume samples;
eight times means 512 times. These studies initialize and evolve fresh fields at
fixed physical extent. They do not enlarge a saved low-resolution volume.

## Completed simulations

| Process | Native grid | Physical duration | Evolution time | Evidence |
|---|---:|---:|---:|---|
| Thermal writing | 768³ | 5; 1600 steps | 651.77 s, Metal | 453 million cells; all three native fields finite. |
| Thermal writing | 384³ | 5; 500 steps | 26.41 s, Metal | Finite native deposit, temperature and dose. |
| Acoustic deposition | 384³ | 6; 1500 steps | 207.75 s, Metal | Finite pressure, accumulated intensity and deposit. |
| Conserved phase | 384³ | 60; 240 steps | 825.7 s, CPU | Mean drift 5.48e−13; energy 64263→45650. |

Timings exclude compilation, allocation, export and photography. This shared
host slowed sustained GPU steps relative to short pilots; use complete-run
measurements when budgeting. Recipes live in `experiments/refinement/`; native
fields and run manifests in `outputs/refinement-20261015/`.

Writing at 96³ and 384³ uses the same path, width, extent, duration and dt=.01.
Its volume-weighted deposit integral changes .0129232→.0129752 (+.402%).
The .018 source width spans less than one cell per standard deviation at 96³,
about 3.5 at 384³. Refinement resolves the writing process itself.

Phase at 96³ and 384³ uses identical physical modes and dt=.25. At t=60,
matched-lattice composition differs by RMS .00794, maximum .0330; discrete
energy differs .386%. These comparisons support this specimen's refinement;
they do not establish convergence for every parameterization.

## Memory is architecture

The default `World` allocates eight fields and eight scratch volumes.
`ThermalWriter` needs two of those fields and owns three more arrays. The study
uses the existing parametric constructor to allocate exactly those five arrays:
**9.06 GB at 768³, versus 34.43 GB** for the default carrier plus writer.

The physics kernel is unchanged. A 120-step 96³ comparison is bitwise equal in
nutrient, deposit, temperature, dose and clock. This specialized world supports
the writer; processes or generic diagnostics requiring absent fields need a
different allocation contract. Raw field exports are not executable checkpoints.

At 768³, explicit thermal diffusion requires a smaller timestep. The recipe
uses dt=.003125, below .0033405, for 1600 steps to the same physical time 5.
The complete run took 651.77 s of evolution (409.73 ms median step). Its three
native field exports occupy 5.44 GB. The deposit integral is .0128622 versus
.0129752 at 384³; this comparison changes both spatial and temporal resolution.

## Correct the photograph too

The old renderer multiplied a binary spherical crop into a scalar field before
extracting its surface. That introduced artificial staircase rings. New
photographs extract the actual field first, then cut the mesh with analytic
geometry. A controlled 96³ comparison changes only that operation and removes
the rings: this improvement is independent of simulation refinement.

Native 384³ phase meshes use a declared .01 surface-meshing adaptivity to remove
redundant near-planar polygons. The scalar field is not downsampled or smoothed.
Surface-only rendering loads the required native channel through a memory map,
avoiding unrelated optical-array copies. The new target is 2400² photographs
with 128 maximum adaptive Cycles samples.

Material presets change actual optical behavior: resin and quartz absorption;
wax subsurface scattering; frosted transmission; thin-film interference;
anisotropic bronze; copper with a chosen patina; porcelain and basalt.
Patina is a shader choice, not simulated chemical growth. Transparent materials
need sufficient optical thickness and suitable lighting to remain legible.

## Ecology needs a different refinement

The old ecology dots were recorded agents, displayed as colored markers. New
material views hide those markers and expose the actual paid film and deposits.
The six diagnostic exports remain native **32³**. Apparent connectivity changes
strongly with the contour threshold; they do not show robust evolved anatomy.

Simply allocating 128³ would change the biology. A 32³/64³ resource probe kept
total inventory and funding equal using cell-volume scaling, but an unchanged
agent uptake request still changed .12025→.06050. Comparable refinement needs:

- Concentration-based sensing and saturation.
- Uptake and deposition over fixed physical supports.
- Physical interaction radii and offspring offsets.

Those changes remain outstanding. Mechanics already stores continuous node
positions and actual bonds; increasing a display voxel grid would not refine
its mechanical degrees of freedom. Dielectric growth is being moved to a
residual-controlled multigrid solver for a fresh 257³ specimen; equal event
counts at different lattice spacings do not represent equal physical growth.

[Collection recipes](collection-20261015.md) · [Alife evidence](alife-adaptation-results.md)
