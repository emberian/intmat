# Native garden observatory

The local observatory shows completed observations of the native three-dimensional constructive garden. It does not run a simulation, read checkpoints, resume a run, or generate intermediate cell positions.

```sh
python3 scripts/observe_native_garden.py RUN_DIRECTORY outputs/native-garden-observatory --follow --max-wall-seconds 1800 --artifact-budget-bytes 256000000
python3 -m http.server 8769 --bind 127.0.0.1
```

Open `http://127.0.0.1:8769/outputs/native-garden-observatory/`. For a completed run omit `--follow`. Julia must be available with this project's dependencies; the detached frame exporter uses Serialization, JSON3 and SHA, without loading the simulation or its device backend. The observer stops at the earliest of its wall bound, the recorded run deadline, completion, an export failure, or its artifact budget. Existing observations remain preserved; no rolling deletion occurs.

Only safe relative `frames/*.jls` paths explicitly listed in `run.json` are admitted, with exact recorded size and SHA256. The exporter checks SHA again before and after reading. Atomic `observatory.json` updates expose only completed exports. Partial manifest reads retain the last good index; the browser likewise retains its last good frame. A failed export remains diagnostic and is never published in the frame index. Restart with the same source and destination to admit subsequent observations. The destination is bound to the resolved source path and first recorded frame SHA; a different path or replaced run is rejected before exports. Legacy indexes are bound only after matching every preserved frame SHA to the supplied source manifest. Previous session stop reasons are archived on restart; the new session evaluates its bounds afresh. Use a fresh destination for a different run.

Drag to rotate and scroll to zoom. Cells use their actual physical radii and positions; clicking a visible sphere or selecting its ID opens identity, internal resources, hidden state, adaptation, actions and catalytic expression. Bonds connect actual surviving IDs using the shortest periodic displacement and are clipped to the observation box. The box is a periodic view boundary, not a membrane. Spheres crossing it can extend beyond the displayed field volume.

Amber and cyan channels are actual field concentrations (or dimensionless permeability). The volume transfer function is a display choice: opacity and display maxima change visibility, not material quantities. Opaque cell bodies occlude fields behind them; fields in front attenuate cell visibility. A Z slice is also available. Integrated field stocks use actual voxel volume. Field scales can be locked across observations for comparison.

Native 64³ samples and 32³ concentration block means are explicitly labelled; the resolution menu can separate them. Both retain the full actual cell and bond records. Playback advances physical model time and selects recorded observations without temporal interpolation. Last-step paid secretion is the realized amount during the final physical step only, not the amount over the capture interval. This schema does not provide interval-integrated per-ID work.

The small real smoke recording is available at `outputs/native-garden-observatory-smoke-20261015/`: eight actual cells, native 8³ and averaged 4³ fields at physical times 0 and 0.002. It validates display and export plumbing, not long-run ecology.

Checks: `python3 scripts/test_observe_native_garden.py` covers hash/size/path admission, partial manifest retention, budget admission failed export retention, source identity rejection and restart bounds. Isolated Chromium additionally exercised actual volume rendering, rotation, identity selection, native/coarse switching, slices, time navigation and mobile layout.
