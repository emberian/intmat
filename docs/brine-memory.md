# What the brine left open

A finite precipitation process leaves an open copper and pale verdigris cage
inside an evolving, conserved phase material. Mineral begins at zero. The image
photographs its native concentration surface, rather than prescribing pores.

The Julia study composes Interior's stabilized semi-implicit Cahn–Hilliard phase
solver with conservative harmonic-face diffusion of two reagents. A bounded
reaction consumes A and B equally and deposits twice that amount as immobile
mineral. Signed phase is never clipped or mean-corrected after initialization.
Reaction vanishes below physical gradient magnitude .35 per model length and
reaches full rate at .833333. These thresholds remain fixed across grids.
There is no basal reaction in the selected recipe.

Mineral changes subsequent transport through resistance 14 and changes phase
chemical potential through `4 * mineral * phi`, favoring local mixing. The
interfaces that gate reaction therefore respond to prior deposition. This is a
phenomenological constitutive law, without a calibrated chemistry or closed
thermodynamic energy model.

The selected initialization uses twelve fixed physical Fourier modes at RMS .12,
with A in a finite bent lens of thickness 4 and B outside. It supplies neither
mineral nor a pore mask. Boundaries remain periodic. The 128³ grid covers extent
24³: spacing .1875; phase characteristic length sqrt(.55)=.742 spans 3.96 cells,
and reservoir thickness spans 21.3 cells. Fine local fronts can remain sensitive
to discretization.

The bounded final reaches t=14.4 with 720 steps at dt=.02. Evolution took 119.61 s,
median complete step 168.61 ms, on two CPU threads. It ends with mineral mass
723.67859, phase range [-.92131,.90670], nonnegative reagents, and signed phase
mean error 8.57e-11. Chemical error is −.0065683 against initial inventory
12860.98533 (relative −5.11e-7).

An independent raw-field auditor verifies both A+mineral/2 and B+mineral/2,
not merely their sum. Their final absolute errors are −.00004188 against
1677.54223 and −.00652642 against 11183.44309, respectively (relative −2.50e-8
and −5.84e-7). Mineral is nondecreasing in every cell across saved frames.

Matched 40³ controls remove each feedback edge separately. Removing only
transport resistance changes mineral by maximum .16272; removing only phase
forcing changes phase by maximum .59247. The closed mineral inventory is
758.77582 versus 1360.98501 without either feedback. A below-threshold interface
produces exactly zero reaction. An oversized reaction step preserves inventory
and nonnegativity. Halving dt .04→.02 on the same 40³ grid changes mineral mass
by .2762%, phase RMS by .01334 and mineral RMS by .00344. The coarse40 and final128
mineral inventories differ about 4.5%; these checks establish bounded sensitivity,
not spatial convergence or exact pore geometry. Earlier basal/lamellar pilots
produced weak pads and remain preserved outside the public gallery.

The final extracts mineral concentration .30 directly, without cropping,
deformation, invented meshes, or a low-concentration enclosing film. Residual
A/(A+B) trilinearly sampled at vertices selects copper-to-verdigris color and
metallic response .12–.62. Saved absolute phase sets roughness .40–.62. These
chosen optical transfers do not simulate oxidation. The floor and lights are
photographic staging. Native geometry has 78,330 vertices and 83,728 faces.

The producer ran from an isolated source stage of commit `365d5d5` plus the
study and recipe hashes recorded in run.json. Subsequent repository changes do
not alter that run. Native fields and high-resolution photographs stay outside
Git; compact published evidence links them through SHA256 digests.

Replay into fresh output directories:

```sh
JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 julia --project=. --threads=2 --gcthreads=1 \
  studies/brine_memory.jl experiments/art/brine-porous-128.toml FRESH_OUTPUT
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. --threads=1 \
  scripts/check_brine_memory.jl FRESH_CHECK_OUTPUT experiments/art/brine-porous-pilot.toml
julia --project=. --threads=1 scripts/audit_brine_memory.jl FRESH_OUTPUT/closed AUDIT.json
/Applications/Blender.app/Contents/MacOS/Blender -b -t 2 \
  --python scripts/render_brine_memory.py -- FRESH_OUTPUT/closed FINAL.png \
  --device METAL --size 1800 --samples 96 --mineral-iso .30 --no-wet
```

Coordinate the local Metal slot before rendering. All simulation algorithms are
Julia; this run uses the existing CPU FFT phase backend, without GPU speed claims.
