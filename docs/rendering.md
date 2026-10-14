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
