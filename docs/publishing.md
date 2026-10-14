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

The curated bundle has a 80 MB budget. `site/arena/` also contains a recorded neural arena: 25 measured frames, agent inspection, field selection and playback, with no live backend. Its compact `replay.json` and `latest.json`, HTML viewer and preview image are curated from `outputs/arena/replay-fairshare-20261014/`; JSON paths are sanitized through the publication builder’s `public` helper. Refresh this explicit arena snapshot separately when publishing a new recording. Raw `.f32` fields, VDB volumes, Blender scenes, browser captures and third-party reference footage are outside the publication manifest. Metadata and embedded viewer data are sanitized; staging fails if a local filesystem or localhost reference remains.

## Repository setup

Enable GitHub Pages using GitHub Actions as its source. For this repository the equivalent authorized CLI command is:

```sh
gh api --method POST repos/emberian/intmat/pages -f build_type=workflow
```

If Pages already exists, use `--method PUT` to update its build type. The workflow uploads `_site/` and deploys into the `github-pages` environment with Pages and identity-token write permissions.

## Metabolic arena recording

`site/metabolic/` publishes eight evenly spaced actual snapshots from the sixteen-frame `outputs/arena/replay-metabolic-96-20261014/` recording, including both endpoints. All 96² spatial fields and all captured agents in each selected frame are preserved. Arena snapshots can cap the captured/displayed agent list (3,000 agents in this recording), while the measured population may be larger; the viewer reports both counts. `public_curation` metadata records original snapshot indices and the selection policy. `latest.json` is the final retained frame. This is recorded playback with resource, boundary and metabolic measurements; no live backend is connected. The original 25-frame neural arena remains available separately.
