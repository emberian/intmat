# Seeing inside the specimen

`render_volume((density=carrier, mask=boundary, activity=activation, structure=deposit), "frame.ppm"; settings=RenderSettings())` integrates through a three-dimensional scalar field. It returns a binary PPM with exposure mapping and sRGB encoding. Calling without a path returns linear RGB in `(3,width,height)` order. `orbit` and `elevation` control the camera; `cutaway` removes points whose specimen z coordinate exceeds its value, exposing deeper layers. All fields occupy `[-1,1]^3` with cell-centered samples. A mask provides the geometric specimen boundary; absent a mask the specimen is the full box.

Density and structural deposit increase spectrally colored absorption. Positive activity emits cyan light; the carrier scatters a directional warm light. The ray marcher uses an exact constant-segment Beer–Lambert update, trilinear field interpolation, and an optional twelve-sample light transmittance approximation. Rendering is threaded by image row. `scattering=0` disables shadow work. This is a fast straight-ray volume preview: no refraction, interface Fresnel terms, caustics, multiple scattering, or physical wavelength calibration. Render emission is an optical design choice; it does not pay the simulation's energy balance.

## Blender / Cycles

`export_volume(fields, "examples/volume")` exports each three-dimensional field as Mitsuba VOL v3 and raw float32, with a `volume.json` manifest. The raw layout is x fastest (Fortran order), little endian; bounds and shape are explicit. On a little-endian Julia host, no extra Julia packages are needed.

Run:

```sh
/Applications/Blender.app/Contents/MacOS/Blender --background --python scripts/render_blender.py -- examples/volume examples/cycles.png --samples 64 --size 640
```

The script requires `numpy` and `openvdb` or `pyopenvdb` inside Blender's Python. Blender 5.0.1 on this workstation supplies both numpy and openvdb. Other installations may need a compatible OpenVDB Python build. The script reads actual exported grids, masks them, writes `specimen.vdb`, saves `scene.blend`, and renders absorption, scattering and activity emission with Cycles. No external generated imagery is used. Cycles' shader coefficients differ from the native preview; these are two artistic optical interpretations of the same simulation state, rather than a renderer equivalence test.

Density is required. Supply density as carrier plus deposit and use structure for the deposited phase. Activity is optional. OpenVDB conversion preserves the full grid; emission does not manufacture filaments absent from activity. The Blender scene has no invented outer glass shell: sharper physical interfaces need an explicit surface geometry model.

The `.vol` grids are also directly usable by Mitsuba's `gridvolume` plugin. They carry scalar values and coordinates, not a full optical material or ready-made scene.

For the choir sample, `cycles-fronts.png` selects saved frame 0300 (t=6), with density `0.06 + 0.8*deposit` and the original viewing mask. The Blender controls `--activity-threshold .65 --activity-power 2 --emission 2` emphasize excitation crests using `max((activity-.65)/.35,0)^2`. This monotone optical mapping reveals high-activity regions; it is not a new simulated field or evidence of additional dynamics. The final t=20 specimen is much more uniformly active and naturally reads as a glowing cloud. Use lower thresholds or power one to see weaker activity.

`--interfaces` provides a second artistic transfer: activity emission follows the spatial gradient magnitude, while amber density is `0.025 + 0.7*edge_strength(deposit)`. Each edge strength is normalized by its 98th percentile inside the viewing crop, capped at 1.5, then raised to power 1.4. Derivatives use the normalized specimen coordinates. This exposes moving excitation fronts and accumulated phase boundaries; it does not claim microscopic interfaces or calibrated optical constants. The raw `.f32` files remain simulation state. The converted VDB contains the chosen optical interpretation; `optical-mapping.json` records controls and normalization scales for reproducibility.

`--smooth N` optionally applies N separable three-tap `[0.25,0.5,0.25]` smoothing passes to each scalar field before calculating interface derivatives. It is an optical homogenization of the coarse grid, recorded in the mapping file; it does not increase the simulation resolution. High powers of gradient magnitude can exaggerate grid artifacts, so compare smooth and unsmoothed views before interpreting fine contour patterns.

The final finer specimen is `outputs/lace-fine`: 96³ actual simulation cells at t=7. Render `cycles.png` with `--interfaces --emission .35 --activity-threshold .55 --activity-power 3 --orbit .3`, without optical smoothing. A second camera (`cycles-above.png`) uses `--orbit .9 --elevation 1.15` on the same baked volume. Both are 640px, 48-sample Cycles images. Comparing 40³ and 96³ versions greatly reduces the apparent fine contour banding; the remaining broad hollow fronts are the stable visual feature. The baked float fields, settings, and VDB are retained for subsequent photography.

## Photographing different baked phenomena

The Blender entry point supports `--palette amber-cyan`, `amber-gold`, `glass-ice`, and `acoustic-pearl`. These set absorption/scattering/emission colors; they do not change the simulation. `--emission-field FIELD` selects an actual exported scalar (for example `deposit`, `temperature`, or `activity`), while `--emission`, `--absorption`, `--scattering`, and `--density-scale` control optical coefficients. Deposited lightning branches can use direct deposit emission, avoiding gradient-based contour artifacts. Thermal tracks can distinguish residual deposited material from current temperature by choosing the source explicitly.

Every new render writes a sibling `<image>.optical.json` and `<image-stem>.blend`, preserving camera and material controls independently for multiple photographs of the same baked arrays. The volume directory retains the most recently converted optical VDB and mapping; raw simulation arrays are unchanged.

### Electrical tree

`outputs/electrical/tree/cycles-filaments.png` photographs the actual occupied conducting sites' face adjacency. `--conductive-graph --graph-radius .007` draws each positive-axis adjacency as an eight-sided tube between recorded voxel centers. The radius is an explicit optical styling choice, smaller than the grid cell; it is not a predicted physical channel diameter. These edges are conductive connectivity, not a uniquely retained branch lineage. Their lattice angles reflect the six-neighbor growth model. Volumetric emission is disabled in this mode so branch width is not dominated by voxel interpolation. The amber carrier still uses the exported volume. The `amber-gold` palette and material's faint emission aid visibility; no additional branches are created.

### Thermal writing

`outputs/writing/{scaffold,living}/cycles.png` uses the same glass-ice material and camera for a direct morphology comparison. Emission source is the actual permanent `deposit` field (gain2); final temperature has mostly cooled, so the track's persistent pale glow is an artistic material interpretation, not current thermal radiation. `cycles-above.png` photographs the living specimen from a second angle at gain.65. The prescribed toolpath is a seeded scaffold; recurrent feedback changes material deposition, rather than originating the path.

## Cut and polished optical specimens

`--polished-carrier` adds a closed refractive viewing surface fitted to the support bounds of the supplied mask (level .01, with half-cell allowance). The default is an ellipsoid; `--carrier-shape box` provides a beveled box for a block specimen. This is an optical extraction from the periodic simulation, not a new mechanical wall or change to its dynamics. IOR1.48 and roughness.035 are chosen artistic constants, overridable by `--ior` and `--surface-roughness`. Cycles transports refraction at this surface, unlike the native straight-ray preview. A dark studio ground and area-light reflections make the cut legible. The fitted axes are stored in each image's optical metadata. These are optical interpretations, not calibrated predictions of real amber or glass.

Tree `polished.png` uses sharply drawn actual conductor connections and bulk scattering .004. Writer `polished-inclusion.png` additionally extracts `structure=.25` via OpenVDB `convertToPolygons` (adaptivity .05), then maps the resulting cell-centered coordinates to specimen coordinates. Its pale frosted material is selected for the deposited phase, with much weaker volumetric emission (.08). This surface derives entirely from the deposited scalar; the isovalue and generated mesh counts are recorded. The interpolation and threshold do not establish a simulated microscopic interface.

`--bond-graph OBJ` renders the vertices and surviving line segments of a supplied actual graph. For the mechanics photograph, `frame-0010.obj` supplies the deformed unbroken bonds; the carrier in `volume-stage` is a constant optical staging field, not a claimed mechanical density. `--graph-radius .0018` controls visible bond thickness. Broken bonds are absent because the simulator omitted them from that snapshot. Graph photographs therefore show actual surviving connectivity and positions, with chosen optical material and cut.

Acoustic `cycles-pressure.png` uses actual instantaneous `pressure` as its source, taking positive compression above .06 with gain8; the carrier is diluted for visibility. This reveals the suspended wave pattern rather than the broad, partly saturated deposited field. `cycles-energy.png` instead observes the actual recorded running squared-pressure intensity with square-root optical transfer. Neither emission model claims that acoustic pressure naturally emits visible light.

New photographs retain an independent sibling `.vdb`, `.blend`, and `.png.optical.json`. Preferred preexisting hero photographs were also archived this way after detecting that the original shared `specimen.vdb` could be replaced by another optical mapping. Copy the scene and its sibling VDB together when moving a photograph. The volume filepath is relative for new scenes.

`outputs/acoustic/polished-wave.png` photographs the same positive-compression crest transfer through a polished refractive carrier. The recorded acoustic field is transient; treating it as visible optical emission is a chosen observation law, not deposited solid geometry. The companion pressure image exposes its morphology without the cut's refraction. Tree `polished-above.png` provides a second camera on the same baked conductor and material as `polished.png`.

### Ecology: a shallow population observation

`outputs/ecology/population-glass.png` shows the current final population: 132 recorded agents, colored by their six surviving founder lineages, at their actual periodic x/y positions. Spheres of radius 0.015 are observation glyphs, not predicted cell shapes or sizes. No edges or organism geometry are added. The recorded deposit is a 64×64×8 extrusion of the two-dimensional simulation. `--z-scale .125` gives the optical block a 2×2×0.25 aspect ratio, with isotropic displayed grid-cell spacing. Every population point remains at z=0. Ecological dynamics remain two-dimensional; depth belongs to the observation. Population IDs and lineage records remain in `population.json`, while optical metadata records the source, count, marker radius and extrusion scale.

## Immutable captures from larger runs

`scripts/export_capture.jl RUN OUTPUT [TIME=5] [CARRIER=.08] [DEPOSIT_GAIN=.7]` exports the nearest complete recorded `frames.jsonl` snapshot into the optical volume contract, reading only saved scalar files. It never advances or recreates state. It validates file sizes, exports the original deposit/activity/nutrient plus optical aliases and viewing mask, and saves `capture-source.json` with actual capture time, source prefix, grid shape, simulation extent, and chosen density/crop interpretation. The current helper targets cubic dense profiles; rendering coordinates are normalized as documented above.

`outputs/bigruns/dense-128-study/photograph-1200.png` is a1200px,64-sample photograph of actual128³ frame00001707 at t=5.0009765625. Its radial interior remains smooth; this dense run demonstrates numerical resolution and the bake/photograph pipeline, rather than claiming intricate morphology. The separate phase, tree, and mechanics specimens provide stronger structural artworks.
