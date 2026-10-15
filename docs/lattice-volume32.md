# A larger native lattice volume, 2026-10-15

The frozen lattice law completed a fresh **32³, spacing-one** specimen to t=320. This is eight times the volume and four times the duration of the earlier 16³ pilot. It retains all 42 actual 3D captures, including the explicitly imposed central spherical injury at t=160, radius 6. It demonstrates executable funded turnover and inherited spatial variation; it does not establish adaptation, anatomical organization or GPU performance.

The recipe is `experiments/alife/lattice-volume32.toml`: seed 11, full condition, eight founder codes, 4,096 shared policy codes, initial occupancy 0.1, precursor 0.6, dissolved foods 0.25 and finite complementary reserves 8, dt=0.1, 3,200 steps. Chemistry and physical spacing match the earlier pilot. Increasing domain size increases total inventory; this is a larger specimen, not a fixed-domain resolution-convergence experiment.

| Measurement | Completed result |
|---|---:|
| Final population | 13,899 |
| Births / deaths | 236,432 / 225,778 |
| Mutation events / living inherited codes | 9,376 / 523 |
| Maximum living generation | 40 |
| Injury-region bodies: before / after / final | 392 / 0 / 394 |
| Injury-region matrix: before / after / final | 300.608 / 0 / 176.089 |
| Complete-step compute | 42.424 s |
| Wall time including captures, excluding startup/warm-up | 43.677 s |
| Startup-inclusive resource-record to final READY | 67.660 s |
| Complete-step throughput | 2.472 million voxel steps/s |

The last timing is timestamp-derived from the resource record and final log modification time; it includes Julia startup, the logged 5.8-second package precompile and kernel warm-up, but excludes wrapper setup before the resource record. It is not a separately instrumented whole-command timer. Systemd recorded 109.636 CPU seconds and a 563.2M memory peak. The enforced job used CPU8–9, two Julia threads, GC1, BLAS1, nice15, MemoryMax=8 GiB, swap0 and a 1,800-second timeout. Resource controls and the systemd journal are retained. No GPU was used; this is not a controlled speed comparison against moving agents or the separate packing law.

All sampled states were finite and nonnegative, within capacity, and satisfied lifecycle accounting. Maximum absolute sampled material/energy residuals were 0.289908 / 0.809605 against initial inventories 561,371.20 / 1,893,390.40: approximately 5.2×10⁻⁷ and 4.3×10⁻⁷ relative. Independent reconstruction from each native frame matched recorded material and accounted energy within 2.33×10⁻¹⁰. Heat comes from the recorded ledger because capture frames contain the named environmental, cell and hidden channels. Population counts matched the native integer fields exactly. Two independent source-pinned final checkpoint restorations followed by 17 further complete steps agreed bit-for-bit across every dynamic array, tick and clock. This proves tested restart determinism, not real-number exact conservation.

All 42 captures were visually reviewed using body and matrix fields in three orthogonal central sections, with a common mass scale and no interpolation. The body field is granular; the matrix nearly spans the volume. The imposed sphere is a **central void**, not a clipped face. Bodies recolonize rapidly; the lower-matrix region remains visible through the final frame. That contrast is a useful native material subject, but the spherical geometry was imposed rather than evolved. One inherited bit gives diagnostic code geography, not a functional phenotype. Final neighboring live sites agree on their full code 19.26% of the time, versus 2.78% under independent abundance-matched draws; local inheritance makes such clustering unsurprising and does not establish cooperation.

Review assets:

- `outputs/alife-20261015/lattice-volume32-review/native-sections.png`: initial, pre-injury, post-injury and final body/matrix/code sections.
- `outputs/alife-20261015/lattice-volume32-review/all-captures-sections.png`: all 42 captures, three orthogonal central sections per field.
- `outputs/alife-20261015/lattice-volume32-review/field-contract.json`: all capture paths, hashes, dimensions, times, scalar names and source channel offsets. Native `.f32` arrays are `(32768,15)`, column-major; `.i32` arrays are `(32768,3)`. Each scalar column reshapes to `(32,32,32)` in Fortran order, extent `(32,32,32)`, spacing one. No optical transfer or contour geometry is baked into these fields.

The immutable hbox run is `/othersys/tmp/intmat-lattice-volume32-20261015`; its complete downloaded copy is `outputs/alife-20261015/lattice-volume32`. The original artifact directory was 121,884,773 apparent bytes, including copied source/environment, native frames, checkpoint, recording and resource evidence, below the approved 150 MB. Later tiny replay/review evidence is retained separately or alongside it. No original history was overwritten. Machine-readable summary: [lattice-volume32-results.json](lattice-volume32-results.json).

The runner adds an explicit `--coordinated-cpu-volume` opt-in limited to exactly 32³ at extent32³, one seed/full arm, at most 3,200 steps, 42 captures and two Julia threads. Without that flag the original 4,096-voxel guard remains. Eight focused guard assertions passed, including rejected oversized grids, extra seeds/arms and excessive capture schedules. Model SHA256 remains `e6d6cab2b0624965a5e7ec66b7f9f9ce44a01f9ded5653db660de61368cde9b5`; the extended runner archived with this run is `5ff3bad263ef4642438a315e3fee59e2ece6b83e9b8f56e73c52bf501641c87f`. Earlier pilots retain their earlier exact runner copies. The prepared AMD gate remains unchanged and awaits its allocation.
