# The shared weather

Three nearby electrical roots write a common mineral interior. Every new conductor changes the harmonic field seen by all three; its heat diffuses, leaves persistent exposure, and converts local feedstock into mineral. Later growth responds to both that memory and the remaining feedstock. The initial roots and susceptibility envelope are designed; the subsequent network and mineral field evolve through the coupled Julia process in `studies/electrothermal_archive.jl`.

This is a separate composition from the native 513³ amber archive. A growth event occupies one lattice cell, so a fixed event count does not represent equal physical development at different resolutions. The new recipe uses native 257³, three closer roots, and 24,576 requested events. It preserves the original 513³ run unchanged.

## Numerical and causal checks

The final source uses explicit face-connected rasterization of each designed root, including samples that cross all three axes. At 257³ the three initial components contain 24, 22 and 21 cells. This corrects a possible diagonal gap in the earlier seed rasterizer for new runs only; earlier source pins and outputs remain intact.

The 33³ coupled CPU/AMD fixture matched all 180 growth choices. The largest field difference was `2.3842e-7`, the relative material-balance error `7.99e-11`, and the final independent harmonic residual `1.1921e-7`. Removing dose feedback first changed the chosen trajectory at event 7. The coarse fixture scales the heat width and diffusion to its grid; it establishes portable coupled execution and causal memory feedback, not resolution-independent morphology. Results are in `benchmarks/2026-10-15-amber-confluence-v3-verification.json`.

The solver uses the previously verified AMD workgroup `(32,4,2)` and `DeviceMGNorm` workspace. These change execution layout and the residual reduction workspace, not the electrical law. A 512-event native257 pilot of the preceding seed raster completed in 33.80 seconds including captures; its median warm complete event took 51.6 ms and its 95th percentile 63.7 ms. These are production observations rather than a controlled cross-machine speedup. The final connected-root production has a distinct source pin and trajectory.

If all frontier growth propensities vanish, the producer now records an arrested run with its actual event count, final capture and termination reason. It does not fabricate remaining frames or throw away the completed trajectory. A zero-propensity fixture exercised this complete path.

## Production and optics

The final immutable runtime is `outputs/art/runtime-amber-confluence-v3-20261015`, based on committed core `02fbeff6363c7546995566900e813e62078819b7` with the separately hashed study and recipe. Its source manifest includes every source file and the verified remote AMD environment. The working-set estimate is 1.62 GB, within a 6 GiB study budget and an 8 GiB host process limit. Captures every 384 events preserve all native Float32 temperature, dose and mineral values losslessly, plus exact graph positions, face edges and birth events. Final feedstock is saved as well. These are observation captures, not complete restart checkpoints.

Production completed all 24,576 events in 1,538.36 seconds, including captures. The final state contains 24,643 conductive cells, maximum dose 2.0844, maximum mineral fraction 0.73810 and minimum nutrient fraction 0.26190. Its recorded relative material-balance error is 4.64e-10 and independent harmonic residual 6.78e-6. The final `run.json` records actual event count, termination, elapsed time, per-capture invariants and source hashes; the independent saved-state validator passed all 65 captures and 803,075 conductor observations. It checks persistent birth identities, exact face adjacency and attachment to older sites. The final pointwise food-plus-mineral error is 6.91e-6 and global relative error 4.64e-10; source hashes match both the immutable runtime and exported source copies. The report is `benchmarks/2026-10-15-amber-confluence.json`. The selected photograph and fixed-camera film use those actual states. Amber, copper, cyan emission, contour thresholds and camera placement are optical choices; they are not a solved optical feedback model.

The selected midpoint is event 11,520 (time 28.8). Its more open branching is preferred for the still; the complete sequence retains the denser later growth. The three electrical components remain separate throughout the saved history, and midpoint mineral contours at 0.007, 0.02 and 0.04 also have three components. Their shared harmonic screening is real, but this run does not demonstrate fusion of their thermal envelopes. Clear mineral on a dark photographic background keeps the actual copper-colored conductor paths legible; warmer absorption and brighter cove variants are retained as unselected optical studies.

The reviewed 2000², 128-sample photograph is now curated in the
[gallery](https://emberian.github.io/intmat/#collection-20261015/amber-confluence).
Its original PNG SHA256 is
`bd3f6b0c6e667ee9001107db6301e958480091175975615affd28288988f9bcf`.
The 65-frame film is complete and visually reviewed, using the same actual captures, fixed camera and optical mapping at 960²/48 samples. All frames use Blender 5.2.2 on hbox CPU cores 14–15, following a reviewed host-portability trial. Per-frame images, consumed source fields, timestamps and renderer hashes were checked before the exact sequence manifest was approved. Its SHA256 is `c72975dbd2ce97a09161cb18d98feae8ce14290cda97d761c0a24b73869fed68`; the film follows time 0–61.44 at eight recorded observations per playback second, without interpolated or reversed frames. Compact movie publication follows the shared animation curator.
