# What the bell kept

A new independent descendant of **Memory rings**, whose smooth sliced shell hid
most of its interior. This specimen starts from three eccentric, corrugated
shells with unequal open mouths and branching fissures. The starting shape is
composed deliberately. The artwork asks what a resonant material retains inside
that shape, rather than claiming the shape spontaneously evolved.

`studies/resonant_reliquary.jl` supplies that initial condition to the existing
interfacial-choir runner. Eight recurrent channels write solid from finite
nutrient and return eroded solid to nutrient. Scalar pressure changes the fast
recurrent pair; solid changes acoustic propagation and recurrent communication.
The actual calculation is Julia on Metal, with a native 128³ grid over a 64³
coordinate extent, 960 steps of .025, and five retained full states through t24.
This is a bounded artistic run with fixed designed dynamics, not an ecological
experiment or calibrated acoustic model.

The source and recipe are retained beside the saved fields in
`outputs/collection-20261015/resonant-reliquary/`. The ordinary capture contract
retains native solid, activity and pressure; no rendering operation advances or
resimulates them. The nacre and warm-metal surfaces are optical isovalues of
solid, with activity/pressure controlling the selected optical material.

The selected photograph records **t12**, solid contours .12/.40, 1600² pixels and
96 Cycles samples. Three pilots compared reflective glass, quieter porcelain,
and a later t24 composition. The quieter t12 view retains the clearest nested
openings. The final PNG, scene, optical metadata and review stay local; its
reviewed WebP is 163,836 bytes. The unchanged simulation advanced 960 steps in
55.339 seconds including five captures; final material error was +0.00015645
on an initial 118817.35795 inventory (1.32e-9 relative). Every capture was finite,
with nonnegative nutrient. This inventory check does not establish a biological
or calibrated physical interpretation.

`python3 scripts/build_reliquary_viewer.py outputs/collection-20261015/resonant-reliquary`
reuses `build_viewer.py` to package t0/t12/t24. Nonoverlapping 2³ arithmetic means
reduce native128³ fields to64³; each field then uses one fixed linear8-bit scale
across all three times. The public self-contained viewer is about3.16MB.
Source/capture hashes and the reduction are embedded. No frames are interpolated
or resimulated. The surface mode raymarches the actual reduced solid at .12;
its chosen tint follows recorded activity. Volume, signed-pressure/activity/solid
slices, rotation, zoom and a z cut remain available. Browser optics are a scalar
observation, not Cycles refraction. Fine fissures soften at display resolution.

Desktop rendering, cut and rotation pixel changes, signed pressure scale, actual
t24 selection and390px mobile layout were checked. The existing gallery card
finds the co-located `viewer.html` without a new publication path. Full staging
with the laboratory addition measured159.42MB, below the unchanged160MB cap.
