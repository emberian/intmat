# Published material gallery

The gallery is served at <https://emberian.github.io/intmat/>. GitHub Pages deploys on pushes to `dev`, or through the **Publish material gallery** workflow's manual run button.

The published collection contains rendered simulation images, sanitized recording metadata and self-contained WebGL2 viewers. The viewers replay recorded fields; the browser does not run Julia or Blender. Documentation and source links point to this repository.

## Refresh the collection

After generating local specimens in `outputs/`, run:

```sh
python3 scripts/build_pages.py --refresh
python3 -m http.server 8770 --directory _site
```

Review <http://localhost:8770>, then commit the updated `site/collection/` assets with the source change. `site/manifest.json` explicitly selects specimens and files; add new images there before refreshing. Missing optional images are omitted, so new renders can be added progressively. The existing snapshot is regenerated only with `--refresh`.

`python3 scripts/build_pages.py` stages the checked-in collection into `_site/` without reading simulation outputs. This is the command used by CI and requires only Python's standard library. No simulations, cloud GPU runtimes, external image services or frontend dependencies are required.

The curated bundle has a 160 MB budget, expanded for the reviewed material films and constructive motion exhibit. `site/arena/` also contains a recorded neural arena: 25 measured frames, agent inspection, field selection and playback, with no live backend. Its compact `replay.json` and `latest.json`, HTML viewer and preview image are curated from `outputs/arena/replay-fairshare-20261014/`; JSON paths are sanitized through the publication builder’s `public` helper. Refresh this explicit arena snapshot separately when publishing a new recording. Raw `.f32` fields, VDB volumes, Blender scenes, browser captures and third-party reference footage are outside the publication manifest. Metadata and embedded viewer data are sanitized; staging fails if a local filesystem or localhost reference remains.

## Repository setup

Enable GitHub Pages using GitHub Actions as its source. For this repository the equivalent authorized CLI command is:

```sh
gh api --method POST repos/emberian/intmat/pages -f build_type=workflow
```

If Pages already exists, use `--method PUT` to update its build type. The workflow uploads `_site/` and deploys into the `github-pages` environment with Pages and identity-token write permissions.

## Metabolic arena recording

`site/metabolic/` publishes eight evenly spaced actual snapshots from the sixteen-frame `outputs/arena/replay-metabolic-96-20261014/` recording, including both endpoints. All 96² spatial fields and all captured agents in each selected frame are preserved. Arena snapshots can cap the captured/displayed agent list (3,000 agents in this recording), while the measured population may be larger; the viewer reports both counts. `public_curation` metadata records original snapshot indices and the selection policy. `latest.json` is the final retained frame. This is recorded playback with resource, boundary and metabolic measurements; no live backend is connected. The original 25-frame neural arena remains available separately.

## Material photographs

The collection opens with **Amber window**, a selected 0.82 composition contour from a recorded 32³ demixed field, photographed inside a transparent amber display matrix. Its small provenance and optical metadata accompany the approved image; raw fields and Blender scenes stay outside the public bundle.

**Metabolic gel** follows a recorded 96² arena at simulation time 1698, using a shallow optical extrusion of the phase and shared-film fields. The image includes the 3,000 captured agents from a measured population of 5,599. Public assets contain the approved photograph and small provenance/optical metadata; the exact observation and Blender scene remain local. Both arena replays remain available.

## Larger artistic collections

The gallery reads nested `collection-YYYYMMDD/<specimen>/run.json` metadata, with searchable titles/descriptions and family filters. Existing photographs and both arena replays remain available.

After the renderer has finished and visually reviewed the selected images, curate them using Python with Pillow:

```sh
python3 scripts/curate_collection.py outputs/collection-20261015
python3 scripts/build_pages.py
```

The curator reads the batch’s `collection-*-selected.json` lists. Each published specimen requires a completed `final-status.json`, a matching approved `art-reviewed.json`, and the original final image hash. Preview approvals do not qualify. Optional `--ids ID ...` limits an incremental publication to selected specimens. The current phase selection uses the aged specimens, excluding superseded pilots.

Approved 1600px finals become quality-95 WebP photographs; full PNG originals remain in local outputs. Public assets include small sanitized run, recipe, render-job and optical metadata. Titles and descriptions follow the reviewed render job; simulation metrics remain in run metadata. Raw states, scenes and private filesystem references are excluded. `site/collection-batches.json` records the curated assets without removing older publication history. CI stages the tracked files and requires no Pillow or rendering runtime.

## Native-grid refinements

A `refinement-YYYYMMDD` batch publishes only records in its explicit `refinement-publication-selected.json`; diagnostic comparison directories are excluded. Each record supplies `id`, `directory`, `native_shape`, `size`, and `samples`, with optional `comparison_notes`. Final status, review hashes and image dimensions must match the selected render job. Images retain the job's full pixel dimensions (for example 2400px / 128 samples), with a minimum of 1600px / 96 samples.

```sh
python3 scripts/curate_collection.py outputs/refinement-20261015
python3 scripts/build_pages.py
```

The declared native grid must match the recorded source grid. Gallery facts distinguish that simulation grid from photograph pixels; newer refinements lead the gallery while earlier photographs remain browsable. Native solver resolution is never inferred from image dimensions or an enlarged render.

## Recorded material films

Render numbered PNG frames from actual sampled states. Keep a fixed optical setup
when motion is intended to show the dynamics alone. In that output directory,
write `sequence.json` with `id`, `title`, `description`, `fps`, `source` (grid,
simulation/source hashes and optical interpretation), and ordered `frames`:
`[{"file":"frame-0000.png","time":0.0,"sha256":"…"}, …]`.
Times must be finite and strictly increasing; all frames must share dimensions.
An optional zero-based `poster_frame_index` selects an actual reviewed frame for
the film poster; omission uses the midpoint. Set it before approving the manifest.

After inspecting the sequence, `art-reviewed.json` records `status: "approved"`,
`reviewed_by`, and `sequence_sha256` matching the entire sequence manifest.
The curator checks every frame hash, review, encoded frame count and duration:

```sh
python3 scripts/curate_animation.py outputs/my-film/sequence.json
python3 scripts/build_pages.py
```

Use a Python environment with Pillow and system `ffmpeg`/`ffprobe`. The curator
encodes H.264 MP4 at CRF18, up to1280px wide, with two encoder threads and a12MB
per-film publication budget. `--max-width` and `--budget-bytes` are explicit
overrides; the whole gallery budget still applies. Source PNGs stay untouched.
One video frame corresponds to one simulation sample: no frame interpolation or
invented dynamics. The film metadata maps each playback frame to its simulation
time. The poster retains original render dimensions.

Publication creates a fresh `site/animations/ID/` atomically, after all checks;
an existing version is never overwritten. Commit that directory and the film
viewer. The builder discovers its metadata and adds it to the material films
catalog. Inspect playback and seeking locally, then publish through the ordinary
gallery workflow. Optical animation and camera motion, when used, belong in the
source interpretation rather than being described as simulation motion.

For a local gallery with working video scrubbing, use the byte-range server:

```sh
python3 scripts/serve_gallery.py --directory _site --port 8770
```

It binds to loopback and serves partial video responses. Python's basic
`http.server` can play a movie but may reset seeks to the beginning. Film URLs
include their content hash so old cached server responses cannot persist across
new published assets.
