# Kept, released, remade

A compact planar motion window into both completed stronger-feedback constructive
continuations, seeds 211 and 307, from model time 3000 to 3100. The exhibit retains
all 201 recorded states and all living cells in each recording, with actual
physical radii, periodic bonds, catalyst and retained-food fields. It adds no
simulation or temporal interpolation. Hidden-state traces, adaptation/controller
readouts, internal stocks and paid interval flows are selectable for any recorded
cell. The sibling-window button translates the whole periodic domain to center
the featured pair; the whole-interval button restores the original view.

Per-cell secretion, construction work and motor work come from deterministic
physical-step replay sidecars. Every replayed physical snapshot matches the
original recording exactly. Removed-cell contributions and floating summation
residuals remain explicitly unassigned rather than being attributed by location,
ancestry or controller output. These are descriptive observations, not a new
adaptation or cooperation assay. [Regime context](constructed-food-fixed-development.md).

Build and verify with standard Python:

```sh
python3 scripts/build_constructed_motion.py
python3 scripts/verify_constructed_motion.py --report outputs/constructed-motion-review-20261015/encoding-audit.json
```

The four deployable files in `outputs/constructed-motion-20261015` are `index.html`,
`recording.json.gz`, `provenance.json`, and the browser-captured `preview.webp`.
The complete window is approximately 2.96 MB. The template is the separate
`web/constructed_motion.html`; the shared arena viewer is unchanged. Browsers use
native gzip decompression, and the loader also accepts transparently decompressed
HTTP responses.

Compact display encoding uses fixed per-field 8-bit scales shared across both
recordings, 16-bit physical geometry and body/battery scalars, 8-bit bounded
controller states, and Float32 internal stocks/paid flows. Field byte differences
are lossless relative to that declared display quantization. The audit checks all
402 frames, 49,909 cell observations, 1,646,592 field samples, 93,270 bond endpoints
and 198,592 assigned paid-flow values against the originals. Maximum observed
position error is 0.00012208 physical units; neural error is below 0.00393701;
paid-flow relative error is below 5.96e−8. Full original observations and exact
Float64 replay sidecars remain separately retained.
