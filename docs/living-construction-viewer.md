# Living construction: inspect actual captures

`build_living_construction_viewer.py` packages `living-construction-v1` Julia
observations into a self-contained WebGL2 page. It also rebuilds from its compact
`living-construction-viewer-v1` artifact. Neither path advances the simulation.

```sh
python3 scripts/build_living_construction_viewer.py RECORDING.json outputs/patch-candidate \
  --max-bytes 6000000 --packed-output outputs/patch-candidate/packed-source.json \
  --title 'After the feed: a living scaffold'
```

The default exhibit budget is 1 MB. A larger budget must be explicit; the builder
does not drop captures, reduce the grid or invent intermediate states to fit it.
The 33-capture patch candidate is about 4.30 MB. Its source and exhibit are
separate from the previously published 17-capture recording; do not overwrite
that historical artifact when staging a new run.

The approved patch exhibit is reproducible from the committed compact source:

```sh
python3 scripts/build_living_construction_viewer.py docs/figures/native-living-patch.json \
  site/laboratory/living-patch --max-bytes 6000000
```

## Read the material

- Scrub exact captures or play forward. Playback holds a capture until the next
  recorded time; slider positions are capture indices and labels show model time.
- Choose construction, bound A/B, raw/free A, or catalyst A/B fields. Both colors
  use maxima over the entire recording, including while comparing captures.
  Opacity, spatial filtering and ray sampling affect display only.
- Pin the current capture, move elsewhere, then **Show reference** / **Return to
  current**. The camera, selected ID and field scales stay fixed. The table keeps
  explicit reference/current times and reports sums of recorded voxel amounts.
- The funding strip reports the displayed capture's cumulative external supply
  against the declared source budgets. It distinguishes external funding from
  energy already stored in cells and material; a spent feed is not an empty system.

Domain dimensions and voxel volume are in model units. Source fields are amounts
per voxel, not concentrations. The browser's packed fields are Float32; scalar
ledgers retain their exported precision. Native positions and radii set body
geometry. The finite box depicts a periodic domain, not a physical container.

## Read the cells and network

Select an actual resident ID. Selection persists when it is absent or not yet
observed. Parent ID, lineage and generation come from the export; no trajectory
matching or cell identity is inferred from proximity.

**Selected ID trail** displays unconnected positions from completed captures up
to the displayed time. These small markers are historical samples, not extra
cells or a reconstructed continuous path. Energy history likewise uses separate
captured points.

Optional `frame.bonds` records contain `a_id`, `b_id`, `rest_length`, `a_position`
and `b_position`. The builder requires distinct active endpoints, one undirected
record per pair, positive finite rest length, and endpoint positions exactly
matching their captured cells. The viewer draws minimum-image periodic edges,
clipped at domain seams. It does not add inferred proximity edges. Network
component counts and largest-component size describe these links at this capture;
they do not establish multicellular individuality or learning benefit.

Birth records retain copying work, copied coefficient norm, funding outcome and
adjusted pending samples. Requested and funded secretion totals appear when the
producer exports them. These observations expose mechanisms without assuming
that inherited coefficients improve performance.

## Integrity and checks

The compact artifact includes the original report's SHA256, source metadata and
a SHA256 for decoded field bytes. Rebuilding verifies field encoding, exact
bytecount, frame offsets, fixed scales, finite fields, cell identity continuity
and bond endpoints. Gzip decode is bounded by the declared size. Field order is
`x + nx*(y + ny*z)`; geometry does not transpose it.

```sh
python3 scripts/test_living_construction_viewer.py
node scripts/test_living_construction_geometry.cjs
PLAYWRIGHT_MODULE=/path/to/playwright node scripts/test_living_construction_viewer.cjs outputs/patch-candidate
```

The Python checks cover raw/compact rebuild equivalence, invalid fields and
identities, bond endpoint corruption, compressed-artifact corruption and explicit
artifact budgets. Geometry fixtures cover periodic seams, physical edge lengths,
shifted origins, missing IDs and native components. Browser checks cover actual
capture timing, retained selections, reference flips, fixed scales, trails,
funding, WebGL errors and narrow layouts. Publication still goes through the
reviewed gallery staging workflow; this builder changes no shared index.
