# Weather inside a lens

A new preserved variant of the gallery's static electrical-tree idea: a folded
material that changes which parts carry current. The previous electrical and
fracture works remain unchanged.

Three folded precursor lamellae, their small initial mineral texture, and the
opposing constant-voltage electrode planes are designed. The potential, local
Joule power, temperature, deposited mineral, erosion and recovering insulating
damage are computed in Julia. The electrodes remain fixed: the bright regions
are not moved with an imposed toolpath or traveling voltage pattern.

Conductivity increases quadratically with mineral and decreases exponentially
with damage. Moderate temperature converts precursor into mineral; sufficiently
high temperature returns mineral to precursor and creates damage. This opens a
negative feedback path that redirects later current. Precursor and mineral
exchange locally without creating material. Electrical source/cooling are open;
temperature and kinetics are nondimensional, and conversion does not debit
latent heat or chemical free energy. This is an exploratory material law, not
calibrated electrochemistry or a closed thermodynamic model.

The thermal update reuses Interior's periodic Laplacian. The electrical problem
has fixed x electrode planes and no-flux y/z faces. The native 65³ and 129³ grids
sample the same continuous geometry/texture over a two-unit span with node
spacing 2/(n−1). Physical sheet widths, diffusivity and all local rates stay fixed. The periodic thermal graph wraps all n nodes, giving ring length n*h: 2.03125 at 65³ and 2.015625 at 129³. This boundary approximation precludes a strict periodic-domain convergence claim; matched controls at each grid share it.
At 129³, dt=.01 gives thermal positivity factor 0.204608 < .9. The earlier 32³ pilot
used voxel-random texture and cell-centered samples; it is a distinct artistic
pilot, not grid-convergence evidence.

The variable-conductivity study solver follows the existing branch prototype's
nested hierarchy, full weighting and trilinear correction. Its harmonic face
conductances and physical-h residual extend the operator to evolving material
and no-flux faces. Cached face weights and strided color sweeps avoid repeatedly
recomputing conductance in hot loops. Separate concrete thermal/material and
Joule functions prevent captured, swapped arrays from becoming boxed. Shared
production source is untouched.

Four discrete manufactured 17³/33³ fixtures (constant and smooth heterogeneous
conductivity) converge below 1e-6 physical residual/maximum conductivity, with
maximum solution error below 5e-8. Independent positive-edge enumeration checks
flux divergence, source power and integrated Joule power at every captured
material state; it does not call the multigrid residual routine.

The completed 65³ coupled run took 91.95 seconds for 500 steps on two CPU threads.
All 501 solves converged below 1e-4; the independently calculated maximum residual
is 9.986e-5. Maximum source/Joule power defect is 4.00e-5. Global material residual
is zero at the reported Float64 sum. It transfers 56604.37 units into mineral and
29403.44 back into precursor; these are nodal inventory sums, not physical grams.

A matched 65³ control freezes conductivity at its initial field while retaining
the same material/thermal laws, geometry and forcing. Between t = 3.75 and t = 5, the
control retains 100% of its former top 2% Joule hotspots and their mean power is
unchanged. The coupled material retains 26.4%, and the former hotspots' mean
power falls to 27.3%. The control reaches Tmax=.07694 with no erosion or damage;
the coupled material reaches Tmax=1.45016. This supports feedback-driven
redistribution in this bounded study; it does not establish robust oscillation
or a calibrated traveling-wave law.

The optical study photographs the actual mineral level set. Actual local Joule
power drives chosen pink emission, and damage drives roughness/thin-film optics.
The blue absorbing polished prism, illumination and optical transfer functions
are photographic staging. The carrier is not claimed as simulated geometry.

Reproduction uses fresh output directories:

```
JULIA_NUM_THREADS=2 julia --project=. studies/weather_lens_solver_check.jl
JULIA_NUM_THREADS=2 julia --project=. studies/weather_lens.jl experiments/art/weather-lens-65.toml outputs/weather-lens-65-new
JULIA_NUM_THREADS=2 julia --project=. studies/weather_lens.jl experiments/art/weather-lens-65-fixed.toml outputs/weather-lens-65-fixed-new
JULIA_NUM_THREADS=2 julia --project=. studies/weather_lens.jl experiments/art/weather-lens-129.toml outputs/weather-lens-129-new
python3 scripts/inspect_weather_lens.py outputs/weather-lens-65-new
```

Each run retains its source snapshots, recipe, scalar fields and independent
checks privately under outputs. The earlier slow profiling attempts are
preserved as partial runs. Moving hot loops behind concrete function barriers
preserved all eight exported fields bit-for-bit through 200 common steps.

The completed 129³ coupled study ran 500 steps in 562.6 seconds. All 501 electrical
solves converged, with independent maximum residual 9.975e-5 and maximum
source/Joule power defect 3.17e-5. Global material residual is 5.82e-11 and maximum
local residual 3.77e-15. The same-resolution frozen-conductivity control ran 19.19
seconds by reusing its converged immutable potential/current field. Its hotspots
remain fixed (100% retention, power ratio1), while the coupled study retains
32.7% and the former hotspots' power falls to 40.4%. The coupled study transfers
267083 nodal units back into precursor; the control has no erosion or damage.
See `weather-lens-controls.json` for the matched measurements at both grids.
These are bounded behavior comparisons, not an oscillation or grid-convergence
claim. Timings describe complete study workloads on a shared host, not controlled
hardware benchmarks.

Optical scalar exports are rounded Float32 samples of the native Float64
simulation. Independent operator and inventory checks use native Float64 state
before publication. The 129³ optical photograph selects the actual t = 4 snapshot;
it does not synthesize an intermediate state.

The root-approved final is1600×1600 at64 Cycles samples. The reviewed curation
origin is `outputs/collection-20261016/weather-lens-129`; its final-status and
art-reviewed hashes both identify the original PNG. The collection selection
explicitly requests minimum_samples=64 with the root's visual-review reason;
the default96-sample curator policy remains intact. The curator publishes the
canonical quality95 WebP and compact control summary. No sample count is relabeled.
The exact renderer used for the finished photograph is retained privately as
`render-source.py` in that reviewed origin; its recorded SHA256 is unchanged.

For future photographs, `scripts/render_weather_lens.py` adds optional
`--device METAL` and `--threads N` arguments, retaining CPU/two-thread defaults.
Read-only enumeration verifies Blender5.0.1 can select the local38-core M2 Max
Metal device. No new GPU benchmark was run. GPU ownership and additional CPU
threads must follow the current shared-host coordination; the hardcoded hbox
CPU quota is not a local Blender capability limit.
