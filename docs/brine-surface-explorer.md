# Brine: an interior you can open

The native 384³ explorer makes the completed mineral surface navigable: orbit,
zoom, move an uncapped section through X/Y/Z, switch which side remains, and
compare saved phase, reagent share and remaining reagent on the same surface.
“Chambers” opens the photographic section; “Whole surface” restores the complete
contour. A chosen copper material provides the initial optical reading.

This is one actual saved state at t14.399999678, step 5760, from the completed
384³ Julia brine epoch. It is not a live simulation or a reconstruction of the
older 128³ trajectory. The earlier replay and photographs remain separate.

## Source and display

Builder: `scripts/build_brine_surface_explorer.py`. Template:
`web/brine_surface_explorer.html`. Candidate:
`outputs/brine-surface-explorer-384-20261016/`.

All four native field files are hashed against the verified epoch inventory and
read in Julia x-fast order. Mineral/A/B remain Float64; phase remains Float32.
OpenVDB DoubleGrid extracts mineral=.30 directly from 384³. There is no scalar
downsampling. Adaptive surface polygons at 0.02 yield 275,297 vertices,
84,112 triangles and 233,223 quads; triangulation yields 550,558 display triangles.
Quads retain their original vertices and deterministic diagonals. Display
adaptivity is separate from native grid resolution; its geometric approximation
is not a measured continuum error bound.

Saved phi/A/B are sampled trilinearly at native extracted vertices, clamped at
voxel-center domain boundaries. Shader interpolation continues those surface
samples across display triangles. Normals are area-weighted vertex normals.
Positions, normals and scalar samples use 16-bit display encodings; triangle
indices use Uint32. Maximum coordinate quantization errors in model units are
[.000149715,.000134846,.000101118], versus native voxel spacing .0625. This bounds
coordinate encoding alone, not extraction/adaptivity or numerical physics error.
Maximum scalar encoding errors are[7.889e-6,7.439e-6,5.883e-6] for phi/A/B.

The complete surface is present in the payload. Plane clipping discards fragments
only; it adds no cap and does not alter source fields or geometry. Model domain
coordinates span 0…24 on each axis. The initial section keeps x≤13.2.

Optical mode maps actual phi, A/(A+B), and A+B to chosen colors, roughness and
metallic response under studio lights. This is an artistic shader, not the
photograph's thin-film renderer or an oxidation simulation. The other modes show
saved signed phi, reagent A share, A+B and neutral material. Color ranges remain
fixed to the full display surface while sections change.

## Payload and checks

The gzip geometry is 7,229,415 bytes; metadata and viewer bring the exhibit to
about 7.25 MB, inside the explicit 8 MB build budget. Decoded buffers occupy
11,562,048 bytes (plus alignment recorded in metadata). `mesh.json` records source
and binary SHA256 digests, source/runtime pins, contour settings, counts,
quantization bounds and array layouts. The browser checks both compressed and
decoded digests, lengths, attribute layouts and triangle index bounds. WebGL2
and browser gzip decompression are required; no remote library is loaded.

Native extraction took 128.99 s on the local CPU with thread limits 2. This is mesh
preparation, not simulation performance. No local Metal rendering was used.
Desktop 1440 px and mobile 390 px reviews exercised all five colors, whole/section
presets, axis/side/plane controls, zoom, malformed-layout/index rejection and
unchanged source metadata. WebGL reported no errors; mobile has no overflow.
Tiny Blender tests independently check x-fast trilinear sampling, boundary
clamping and display quantization bounds.

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 TBB_NUM_THREADS=2 \
/Applications/Blender.app/Contents/MacOS/Blender -b -t 2 --python-exit-code 1 \
  --python scripts/build_brine_surface_explorer.py -- \
  outputs/art/native-brine-384-f64-20261016 NEW_DESTINATION \
  --adaptivity .02 --max-bytes 8000000

/Applications/Blender.app/Contents/MacOS/Blender -b -t 2 --python-exit-code 1 \
  --python scripts/test_brine_surface_math.py
PLAYWRIGHT_MODULE=PATH_TO_PLAYWRIGHT node scripts/test_brine_surface_explorer.cjs \
  NEW_DESTINATION REVIEW_DESTINATION
```

Root reviewed desktop sections, reagent colors and both portrait presets. The
published route is `/laboratory/brine-interior/`, linked from its photograph.
The existing 128³ replay and the photograph's image bytes remain unchanged.
