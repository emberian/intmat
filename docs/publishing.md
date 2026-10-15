# Published material gallery

The gallery is served at <https://emberian.github.io/intmat/>. GitHub Pages deploys on pushes to `dev`, or through the **Publish material gallery** workflow's manual run button.

The published collection contains rendered simulation images, sanitized recording metadata and three self-contained WebGL2 viewers. The viewers replay recorded fields; the browser does not run Julia or Blender. Documentation and source links point to this repository.

## Refresh the collection

After generating local specimens in `outputs/`, run:

```sh
python3 scripts/build_pages.py --refresh
python3 -m http.server 8770 --directory _site
```

Review <http://localhost:8770>, then commit the updated `site/collection/` assets with the source change. `site/manifest.json` explicitly selects specimens and files; add new images there before refreshing. Missing optional images are omitted, so new renders can be added progressively. The existing snapshot is regenerated only with `--refresh`.

`python3 scripts/build_pages.py` stages the checked-in collection into `_site/` without reading simulation outputs. This is the command used by CI and requires only Python's standard library. No simulations, cloud GPU runtimes, external image services or frontend dependencies are required.

The curated bundle has a 128 MB budget. `site/arena/` also contains a recorded neural arena: 25 measured frames, agent inspection, field selection and playback, with no live backend. Its compact `replay.json` and `latest.json`, HTML viewer and preview image are curated from `outputs/arena/replay-fairshare-20261014/`; JSON paths are sanitized through the publication builder’s `public` helper. Refresh this explicit arena snapshot separately when publishing a new recording. Raw `.f32` fields, VDB volumes, Blender scenes, browser captures and third-party reference footage are outside the publication manifest. Metadata and embedded viewer data are sanitized; staging fails if a local filesystem or localhost reference remains.

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
