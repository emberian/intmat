# The second seam

A curved, layered spring cloth is pulled, sheared and twisted from its two opposite edges. Original connections break irreversibly. After a delay, a finite construction stock at their endpoints may pay for a new, softer connection. The new connection begins at its present separation, so it adds no instantaneous spring energy, but it transmits subsequent loads and can break again. A copper seam is a second history, not a reset of the first.

`studies/repair_weave.jl` composes Interior's existing `BondNetwork` mechanics with this study-local construction rule. The initial curved lattice, weakened half-plane notch and grip motion are designed. All subsequent node positions, strains and failures come from the mechanical integrator. There is no self-contact model; the specimen is a spring network, not a calibrated continuous cloth.

## Construction contract

An eligible broken edge waits the configured delay and must remain within 1.8 times its own rest length. Both endpoints must afford half of `0.18 × current separation`. This is a normalized construction budget, not a calibrated mass law; node inertia stays fixed. Candidate edges are funded in persistent ID order. Each new edge receives a new ID, records its broken parent, takes 70% of the parent's stiffness, and uses the current separation as its rest length. At most two repair generations are allowed. Parent damage, rest lengths and IDs never change.

External grip work, kinetic and elastic energy, damping and fracture loss remain in the existing mechanical ledger. Repair stock has a separate balance. The producer checks nonnegative stock, finite motion and the stock balance at every capture; it checks zero initial extension and broken ancestry at every insertion. `scripts/check_repair_weave.py` independently checks saved identities, ancestry, irreversible damage, finite motion, stock and energy residuals over the captured history.

## Runs and evidence

The 1,440-node twisting pilot (`outputs/art/repair-weave-pilot-v2`) ran two six-unit loading cycles. It made 1,632 new bonds; 360 of these later broke. Its final stock residual was `1.71e-13`, and its mechanical energy residual was `0.00258` against `7.9853` external work. The otherwise matching zero-stock control made no new bonds and ended with `3.3013` external work. These runs establish that repair changes later load transfer and failure; they do not establish that repair improves every outcome.

The scalar-loop mechanics optimization was replayed against this complete pilot: all 61 geometry files and all history metrics matched exactly. Its source hash is recorded in `docs/bond-network-performance-results.json`. The source-pinned finer run uses that implementation, 5,472 nodes, 52,928 original bonds, a deeper initial fold, opposing 0.55-radian grip twists and 65 captures over 24,000 steps. This is a finer discrete specimen, not a demonstrated continuum-convergence result. Its runtime is `outputs/art/runtime-second-seam-20261015`; all original runtimes and pilot states are retained.

## Optical interpretation

`render_repair_weave.py` photographs the saved positions and surviving straight bonds. Original fibres and funded repair fibres have different chosen materials. A triangular optical panel appears only if all three connections supporting it survive; a panel supported by any repair edge receives the copper material. This makes actual topological openings visible. Panels, fibre radii, colors and strain-dependent thin-film thickness are optical choices and add no forces. The camera is fixed across a sequence, with no interpolation or animated deformation added by the renderer.

The finer run completed in 301 seconds and made 1,075 repair bonds, of which 377 later failed. Across 65 captures, the independent validator checked 3,486,175 bond observations. The maximum stock residual was `1.14e-13`; the maximum absolute mechanical energy residual divided by `max(1, |external work|)` was `0.0004215`. The measured report is `benchmarks/2026-10-15-second-seam.json`.

The reviewed 2000², 128-sample photograph at step 15000 is published in the
[material gallery](https://emberian.github.io/intmat/#collection-20261015/second-seam).
Its triangular panels expose surviving connections; copper distinguishes funded
repairs. The completed fixed-camera film uses a finer 193-frame replay, capturing every 125 steps through time 0–12. Its 960², 64-sample images were rendered consistently in Blender 5.0.1 on Metal, with no interpolation or reverse padding. The film plays at 24 fps and uses time 7.5 as its poster. Every original geometry file, image, timestamp, renderer and run hash passed validation; artist and root reviewed all 193 frames. Reviewed sequence SHA256: `27580dd8ef61dffcc6be05d89e14488b572ac95ddbe0352626e1924bd8f8f385`.
Full PNGs, scenes and recorded mechanics remain in local outputs.
