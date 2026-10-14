# A material that bends the view through it

The separate CPU graded-index observer follows curved rays through a frozen scalar refractive-index volume. This is an executable geometric-optics branch alongside the ordinary volume renderer and Cycles observer. A background chart makes positional distortion visible; the same specimen can be observed again without advancing its thermal evolution.

The research connection is direct: [Pediredla et al., Path Tracing Estimators for Refractive Radiative Transfer](https://imaging.cs.cmu.edu/rrte/) shows heterogeneous-index rendering and ultrasonic virtual lenses. Their [author manuscript, section 3](https://imaging.cs.cmu.edu/rrte/assets/Renderer_AuthorCopy.pdf) presents nonlinear ray tracing and refractive transport; it also explains why connecting two prescribed points requires solving a curved-path boundary-value problem. Their renderer uses cubic B-spline index interpolation and specialized path connection machinery. Interior's preview implements an initial-value ray observer with trilinear interpolation, not those path-tracing estimators. No external implementation was copied.

For geometric arc length `s` and unit direction `v`, Fermat's ray equation `d(nv)/ds = ∇n` gives:

```
dr/ds = v
dv/ds = (∇n - v (v⋅∇n))/n
```

The projection removes the direction-parallel gradient, keeping the exact direction derivative perpendicular to `v`. We evaluate the analytic gradient of the same trilinear interpolant used for `n`, integrate with explicit midpoint, and renormalize direction. A plane-crossing fraction limits the final segment. This integrator is suited to a bounded preview, not a guarantee of long-distance accuracy; gradient discontinuities at voxel-cell boundaries mean an ideal second-order convergence rate is not promised. Step refinement measures the resulting endpoint error.

The index array must be a finite positive CPU strided array with nonempty dimensions. GPU fields require an explicit `Array` transfer. Rendering validates the field once before tracing rays; independent ray calls validate their field. A finite arc-length budget makes trapped or backwards rays unresolved. The image paints such pixels magenta and records their count rather than inventing a background hit.

## The first frozen specimen

Run `julia --threads=4 --project=. bin/interior.jl study grin` from the project. It advances a 64³ `ThermalWriter` for 140 steps of 0.01, then freezes temperature and deposit. The optical interpreter maps temperature to an exaggerated positive index perturbation, multiplied by a taper that reaches ambient index 1 near every box boundary. The background is a fixed luminous grid/check/circle chart; an orthographic camera observes it through the slab.

This thermo-optic law is invented and deliberately strong so the bend is visible. It is not a calibrated glass or amber coefficient. The taper is consequential: there is no unimplemented Snell/Fresnel index jump at the bounding box. Supplying arbitrary nonambient values on a boundary does not add a valid physical interface. The `GradedIndex` constructor cannot know whether the caller intends an infinite medium or a physical specimen; the caller must choose a smooth ambient transition for this observer.

`outputs/grin/constant-index.ppm` and `thermal-index.ppm` use the same frozen specimen and chart. Only the index observer changes. `observation.json` records the thermal state, optical interpretation, index-array hash, camera/chart settings, selected curved paths, unresolved counts, and two successive endpoint refinements. Exported temperature, deposit, and refractive-index grids are view independent; the image is not the stored specimen.

The executed specimen reached time 1.399999 with index range **1.0–1.0639205**. The 320² image had **zero unresolved rays**. At 64², the maximum background-hit coordinate difference was **0.00123356** for steps 0.012 versus 0.006 and **0.000616251** for 0.006 versus 0.003. This decreases but does not show ideal second-order convergence on this trilinear field. The largest hit-coordinate difference between constant and varying index was **0.261756** in specimen coordinates. PNG previews were visually inspected: the straight chart remains square, while the frozen thermal tracks bend its lines into localized ripples and distort its circle. These are geometric-optics preview results, not physical predictions for a glass sample.

## What the observer establishes

Focused tests check a straight ray and optical length in a constant medium, exact interpolation of a linear index field, bending toward higher index, decreasing positional error under refinement, and the conserved longitudinal optical momentum `n*v_z` in a stratified `n(x)` field. Invalid fields and finite path budgets are exercised. The example reports endpoint changes between step sizes 0.012, 0.006 and 0.003 separately from the chart image, whose sharp edges can magnify tiny positional differences.

The focused test run passed **22/22** assertions. This establishes the listed invariants and rejection cases; it does not validate a calibrated optical constitutive law or a general refractive radiative-transfer solver.

This observer has no scattering, absorption, emission, polarization, interference, spectral dispersion, interfaces, or point-to-point light connection solver. The camera and background endpoints are in ambient index 1, so pure background radiance can be sampled without an endpoint refractive radiance correction. Adding internal emission or scattering requires the relevant refractive transport weights; it should not reuse ordinary source integration silently.

The next strong artistic study is an acoustic pressure field making a transient lens around a permanent written structure. A chart or caustic target reveals the transient lens, while another view reveals the deposited memory. The pressure-to-index law, attenuation, and observer must remain separate. Later a real light/material feedback loop should feed absorbed stimulus into material evolution, rather than feed camera pixels into it.
