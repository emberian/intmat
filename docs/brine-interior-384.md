# The chamber the brine kept

A close observational section of the completed native 384³ brine material exposes
its arches and copper-lined chambers. The earlier 128³ photograph and all full
states remain preserved. This is a separate finer-grid run and a new photograph,
not an enlargement of the earlier image or its 64³ point-sampled replay.

The shared Julia producer completed 5760 steps at dt=.0025, reaching
14.399999678 model time units. Phase/forcing/permeability are Float32; chemical
A, B and mineral are Float64; FFT arithmetic is Float64/ComplexF64. Extent remains
24³ and spacing is .0625. The original bounded reaction and both mineral feedback
edges are retained. The frozen epoch records core commit `aae69fe`, exact executed
study hashes and recipe SHA256
`cb23920f1f50695a43665ed7074324439ed8cc17a6cea6166e92fe2a2784fd70`.
The actual brine study is hash
`5ab3ad1ba93f07383af119a59b76d8ccf5b390c807d0dc141bb9d517b5287639`.

The final mineral mass is 716.8644887, chemical relative residual −9.90e−16 and
signed phase mean error −9.36e−11. Producer evolution took 839.42 seconds, with
median warm complete step .14445 seconds; total wall time was 1036.94 seconds.
These are one completed dense trajectory, without a claim of spatial convergence
or demographic/ecological behavior. Finer resolution and changed precision do
not by themselves establish converged pore geometry. The phenomenological law
still does not close thermodynamic energy.

All four full fields are read directly at their recorded native types and Julia
column-major order. The renderer uses read-only Float64 chemical memmaps and an
OpenVDB DoubleGrid for mineral extraction at concentration .30, the same physical
isovalue as the 128³ photograph. It does not downsample or first cast chemical
state to Float32. Mesh adaptivity .005 is an optical surface simplification,
recorded independently of the native scalar grid. The full surface has 503,255
vertices and 557,107 faces; the chosen section has 297,467 vertices and 328,078
faces. The four field digests match the independently verified epoch inventory.

The camera observes a deliberate uncapped planar cut: retain normalized x≤.10,
equivalent to model x≤13.2 with the original 0–24 grid coordinate. No source scalar
is changed and no invented cap or interior mesh is supplied. Dark oval openings
at that plane are observation boundaries, not a claim of naturally hollow tubes.
Photographic coordinates center the domain and divide by half-extent 12.

The phase-metal transfer uses actual signed phase to mix navy and umber, select
metallic response, and assign a chosen 180–400nm thin-film interference thickness.
Residual A/(A+B) adds copper. Actual remaining reagent A+B modulates brightness
and roughness. These are explicit artistic optics, not an oxidation solve or
simulated film thickness. There is no texture noise or invented microscopic
geometry. In the observed section, signed phase ranges approximately
−.5454…+.4885, reagent ratio .0088….9813, and remaining reagents .6891…1.0068.
The photograph recipe records every transfer and camera/cut setting.

Replay from the immutable local epoch, coordinating the local Metal slot:

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 TBB_NUM_THREADS=2 \
/Applications/Blender.app/Contents/MacOS/Blender -b -t 2 \
  --python scripts/render_brine_interior.py -- NATIVE_EPOCH FINAL.png \
  --device METAL --size 2000 --samples 128 --mineral-iso .30 --no-wet \
  --mode interior --distance 2.45 --elevation .48 --orbit 1.1 \
  --target -.12 0 0 --cut-x .10 --optics phase-metal
```

The renderer is `scripts/render_brine_interior.py`; optical recipe is
`experiments/art/brine-interior-384.toml`. Full fields, scenes and high-resolution
PNG remain private artifacts outside Git. Publication should retain the native
384³ label, source precision, explicit cut and chosen optical interpretation.

Root approved the exact 2000²/128-sample final. The Q95 publication candidate is
190,574 bytes and was visually checked after compression. Photograph/render
metadata records 39.34 seconds including source-file hashing after rendering; this
is not a native simulation timing. The original run report is byte-identical to
`benchmarks/2026-10-16-native-brine-384-f64-complete.json`, allowing the
photograph source hash to be checked against compact tracked evidence.
