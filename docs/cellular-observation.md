# Cells in motion and material

[Interactive observations](https://emberian.github.io/intmat/cellular/) contain
18 actual archived states and 360 short-continuation frames from the first two
cellular soaks. A pair of 2000² Cycles photographs shows early connected and later
dispersed mosaic populations through the same optical treatment.

At the captured 86-minute mark, broad had 212 cells, 3,828 births and 3,680 deaths;
mosaic had 219 cells, 4,817 births and 4,662 deaths. Neither had blocked births.
Contacts were sparse: zero and nine bonds in those particular latest snapshots.
Earlier reconstructed frames contain 345 and 431 bonds, with large connected
assemblies. Initial adhesive founders and rich initial stocks helped create that
regime; these images do not establish evolved multicellular organization.

Later small groups remain dynamic. In a separate broad continuation from model
time 5596.5 through 5636, 9 births and 10 deaths accompany 5.27 material units of
junction transfer. Twenty-three of its 71 starting ID edges appear in every saved
frame. The mosaic continuation has 10 births, 12 deaths and .28 transferred units;
two of its 20 starting edges occur throughout. Sampled persistence does not exclude
breaking and reforming between frames. These are observations, not a comparative
fitness assay or a limit on what the model can express.

The browser can color measured bodies by metabolic allocation, energy, adhesion,
or recurrent state. It retains physical positions/radii, stable identities and
recorded bonds. Precursor fields use six decimal places for display; cell state
retains exported precision. The chemistry contrast is fixed across each selected
trajectory. It does not interpolate motion between samples.

## Reproduce

The capture lives under `outputs/alife-observation-20261015/`. Its extraction
manifest records original checkpoint hashes and loaded Julia source. Archived
states are historical observations. Early and mature animations are explicitly
labeled continuations of the original checkpoints with unchanged code and
parameters, computed separately from the live jobs. The early continuations
match independently saved live counters and exact bond sets at sampled steps.

```sh
python3 scripts/curate_cellular_observations.py outputs/alife-observation-20261015/extracted
/Applications/Blender.app/Contents/MacOS/Blender -b -t 4 \
  --python scripts/render_cellular_art.py -- SNAPSHOT.json IMAGE.png \
  --size 2000 --samples 96 --device METAL
python3 scripts/build_pages.py
```

The material renderer maps measured planar radii to shallow lenses. Lens height,
pearlescence, gel depth and lighting are optical staging, not inferred anatomy.
Teal-to-amber body color maps relative import-A allocation. The bath uses actual
precursor concentrations, and filaments use actual periodic ID bonds. No nuclei,
branches or procedural microtexture are invented. Each photograph records its
source hash, mapping and geometry policy beside the retained Blender scene.

The larger artistic destination is still activity that changes what the material
does: deposits that redirect communication, moving interfaces that retain memory,
and internal dynamics that become structure. A readable picture of cell contacts
is a useful observation, not completion of that destination.
