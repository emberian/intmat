# Exact static recording transport

The publication build now packages eight legacy JSON replays as gzip assets.
Repository recordings, images, films, numerical values and captured frames remain
unchanged. Only the publication staging tree uses the compressed transport.

The byte survey found no exact duplicate assets larger than 5 KB. The largest
avoidable cost was plain JSON: arena and metabolic replays plus six cellular
recordings. Images and films are already distinct art assets; recompressing those
lossily is unnecessary for this change.

| Recording | Frames | Original bytes | Published gzip bytes |
| --- | ---: | ---: | ---: |
| Neural arena | 25 | 14,443,191 | 3,998,524 |
| Metabolic arena | 8 | 19,754,343 | 5,083,582 |
| Broad early | 100 | 8,163,080 | 3,157,101 |
| Broad mature | 80 | 6,107,140 | 2,273,057 |
| Broad history | 9 | 694,008 | 237,548 |
| Mosaic early | 100 | 8,470,808 | 3,281,577 |
| Mosaic mature | 80 | 6,459,165 | 2,403,097 |
| Mosaic history | 9 | 713,261 | 246,727 |
| Total | 411 | 64,804,996 | 20,681,213 |

The complete reviewed stage falls from **159,397,735 to 115,281,578 bytes**.
Net saving is **44,116,157 bytes**, including 7,626 bytes of shared loader,
manifest and rewritten-reference overhead. The unchanged 160 MB cap has
44,718,422 bytes of headroom in that stage. These are publication artifact
sizes, not a claim of equal network acceleration: a host may already compress
plain JSON over HTTP. Local first-load timings varied by roughly 0.1 seconds;
the extra checksum/decompression work is not presented as a decoding speedup.

## Shared publication path

`build_pages.py` calls `static_recording_transport.stage_recordings` after copying
the complete public tree. The helper handles an explicit eight-file set, verifies
that gzip roundtrips to identical bytes and parsed states, and refuses unfamiliar
loader/catalog layouts. It checks private/local references in decoded JSON before
compression so binary packaging cannot hide them from the publication gate.

Only freshly copied staging files are replaced with `.json.gz`. The original
`site/` JSON files remain available in the repository. The two arena loaders,
the cellular loader, its default selection and all six catalog entries are
rewritten together. A single `recording-transport.js` serves all three pages.
There are no additional remote dependencies. The existing compact embedded
living-construction recordings keep their own codec and are not rewritten.

The staged `recording-transport.json` records original/compressed paths, sizes,
SHA256 digests and frame counts. The browser verifies decoded byte count and
SHA256 before returning a JSON response to the original UI. An HTTP failure
retains the UI's original unavailable-recording behavior. A missing browser
`DecompressionStream` produces a visible upgrade message when decompression is
needed; missing secure-context hashing also produces a visible message.

The loader inspects the fetched body's gzip magic bytes. If a static host sends
a `.gz` asset directly, the loader uses native `DecompressionStream('gzip')`.
If the host sends `Content-Encoding: gzip`, browser networking has already
decoded the body: the helper skips decompression and checks those original JSON
bytes directly. A real HTTP test with `DecompressionStream` removed passed this
already-decoded path. Thus no double decompression depends on server MIME or
Content-Encoding conventions.

Relative URLs were tested under a non-root `/published/` prefix, matching the
structure of a GitHub Pages project base path. The helper resolves its integrity
manifest against its own script URL and each recording URL against the document
base. All three page routes and all six catalog options loaded successfully.

## Verification

```sh
python3 scripts/test_static_recording_transport.py
python3 scripts/build_pages.py --destination /private/tmp/recording-transport-review
```

Three narrow tests cover exact bytes/parsed states/checksums/catalog rewrites,
private-reference rejection before mutation, and refusal of unknown loaders.
The actual eight-file build independently passed the same byte and state gates.

Browser comparisons of each original/staged loader produced identical SHA256
hashes for its full parsed recording and initial rendered canvas. Scrubbing,
playback, field/color controls, all cellular catalog choices and 390×844 mobile
layouts passed without script errors or overflow. The metabolic comparison loads
the full original 19.75 MB replay, preserving all eight published states and
96² fields. No original selection is further reduced by this transport change.
