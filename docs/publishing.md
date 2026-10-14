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

The curated bundle has a 30 MB budget. Raw `.f32` fields, VDB volumes, Blender scenes, browser captures and third-party reference footage are outside the publication manifest. Metadata and embedded viewer data are sanitized; staging fails if a local filesystem or localhost reference remains.

## Repository setup

Enable GitHub Pages using GitHub Actions as its source. For this repository the equivalent authorized CLI command is:

```sh
gh api --method POST repos/emberian/intmat/pages -f build_type=workflow
```

If Pages already exists, use `--method PUT` to update its build type. The workflow uploads `_site/` and deploys into the `github-pages` environment with Pages and identity-token write permissions.
