# Lossless recording delivery

The four curated WebGL volume viewers store their existing base64 field arrays
in one deterministic gzip stream. This changes transport only: every decoded
8-bit voxel, fixed scale, recorded time, camera, shader and control remains the
same. The Julia Float32 recordings and all actual frames are preserved.

`build_viewer.py` uses `viewer_storage.compress_html` after assembling the viewer.
The helper wraps initialization in an async function, decompresses the inline
payload, and restores the original `data.fields` before the existing renderer
starts. `build_pages.py --refresh` also upgrades historical HTML after its
public-metadata sanitization, so an explicit refresh keeps the compact format.
Custom adapters, including the reliquary's contour shader, keep their
source transformations. There are no fetches, CDN dependencies, or required
server configuration: downloaded HTML still opens offline.

The implementation follows the WHATWG [Compression Streams standard](https://compression.spec.whatwg.org/#interface-decompressionstream)
using `Blob.stream().pipeThrough(new DecompressionStream('gzip'))` and
`Response.text()`. Gzip decompression is [available across current browsers
since May 2023](https://developer.mozilla.org/en-US/docs/Web/API/DecompressionStream).
Older browsers receive a visible upgrade message. The format intentionally uses
gzip rather than newer, less widely supported compression formats.

## Controlled migration (2026-10-16)

Only these existing curated HTML files were upgraded; no collection refresh or
source simulation was needed. Original HTML was backed up privately under
`/private/tmp/intmat-viewer-originals-20261016`.

| Recording | Previous bytes | Compressed bytes |
| --- | ---: | ---: |
| Electrical | 3,343,034 | 41,565 |
| What the bell kept | 3,161,071 | 1,230,988 |
| Choir | 5,392,776 | 2,436,867 |
| Neural | 4,880,126 | 1,694,931 |
| Total | 16,777,007 | 5,404,351 |

This frees **11,372,656 bytes** without removing art or reducing quality. A
complete staged gallery has 68 specimens and occupies **148.05 MB**, passing the
existing 160 MB cap. This is an on-disk/static delivery measurement; a web host
may apply its own HTTP compression, so it is not a claim of proportional network
speedup.

All **12,535,296 decoded quantized bytes** match the prior viewers exactly;
metadata and scales also compare equal. Eight adapter/storage tests and the
optional-gallery-viewer regression pass. The storage tests check deterministic
output, idempotence, byte preservation and refusal of unknown bootstraps.

For future explicit upgrades:

```sh
python3 scripts/viewer_storage.py --backup /private/tmp/viewer-backup PATH/TO/viewer.html
python3 scripts/test_viewer_storage.py
python3 scripts/test_build_viewer.py
```

The CLI only changes explicitly named viewers, writes original backups with
exclusive creation, and replaces each viewer atomically. It never deletes old
outputs or caches.

An illustrative second gzip pass (as an HTTP server might apply) produces
4,035,575 bytes for all four old HTML files versus 4,064,282 bytes for the new
ones: a 0.71% increase in network transfer under HTTP gzip. The main benefit
is freeing the gallery's uncompressed artifact budget and avoiding parsing huge
uncompressed inline JSON on arrival, rather than promising faster cold loads.

The optional `scripts/test_viewer_storage_browser.cjs PRIVATE_BACKUP_DIRECTORY`
harness launches an isolated browser, compares old/new framebuffer hashes and
exercises frame scrubbing, slices, rotation, zoom, playback and mobile layout.
Set `PLAYWRIGHT_MODULE` to the installed Playwright module if needed. Its cold
ready times include software rendering and shader compilation; they should not
be treated as stable load-performance measurements.

Observed old/new first-render times in the isolated Chromium software-rendering
run were Electrical 17.8/30.3 s, reliquary 39.1/50.5 s, Choir 62.9/23.7 s and
Neural 59.8/49.9 s. The large variation and shared-host contention preclude a
cold-load speed claim. All four first-render framebuffer SHA256 hashes match
exactly between old and compressed files, with no WebGL errors. Desktop controls
and 390×844 responsive layout were exercised; this is browser emulation rather
than a physical mobile GPU benchmark.

Neural's final full-volume mobile screenshot exceeded the initial 30 s capture
timeout under software rendering. A separate 390×844 slice-view run passed
mobile scrubbing/playback, zero overflow, no script errors and `GL_NO_ERROR`,
and saved a reviewed screenshot. The reusable harness allows a longer capture
timeout; no rendering law was changed to accommodate this host limitation.
