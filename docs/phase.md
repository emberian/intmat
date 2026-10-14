# Conserved phase separation

`PhaseSeparation(w)` adds a signed, process-owned composition φ. A homogeneous
mixture demixes into positive and negative domains, then reduces interface area
by coarsening. The mean composition selects approximately balanced connected
labyrinths (mean 0) or dispersed droplets (moderately negative mean). These are
regimes of this equation, not scripted geometry.

The periodic finite-difference model is

    F = Σ [bulk (φ²−1)²/4 + κ |∇forward φ|²/2] cell_volume
    μ = bulk (φ³−φ) − κ Δφ
    ∂t φ = div(Mface grad μ)
    Mface = mobility / [1 + mineral_pinning (deposit_left+deposit_right)/2].

Every face transfer is antisymmetric. No composition clipping, normalization,
or mass correction occurs during evolution. Float32 rounding limits observed
mass conservation. `phase_energy(w,p)` evaluates the matching discrete energy.
`stability_limit` bounds the explicit fourth-order operator with a margin on the
current cubic slope; it is conservative, state-dependent, and not a proof of
unconditional energy decay. The tested small-step regime decreases energy. The caller supplies finite,
nonnegative world deposit/nutrient fields; structural backend and size validation
preflights process ownership, while the timestep guard scans composition extrema.
Physical-grid refinement requires smaller steps proportional to spacing⁴.
For fast exploration use spacing near one; κ=0.7 then interfaces span a few cells.

Matter participates in the surrounding simulation. Deposit slows demixing by
reducing mobility. The positive phase reduces nutrient permeability according
to `1/(1+phase_barrier*(1+tanh(φ))/2)`; symmetric face means conserve nutrient.
The interface gradient drives activity toward a bounded target with linear
relaxation. Activity can feed `Deposition`, which pins later motion. Parameters
can disable either coupling for isolated thermodynamic experiments. Nutrient
and activity feedback are constitutive artistic extensions, not a calibrated
binary-alloy model. Existing unrelated `Diffusion` processes do not automatically
see composition: use the phase process's nutrient diffusion or explicitly add
another diffusion channel. Tuple ordering remains consequential.

`phi`, `chemical_potential`, `interface`, scratch, and every effective parameter
are exposed to `record_process_fields`; initial state can be recorded before the
first step. `phase_optical_fields` derives density from a smooth positive-phase
fraction and emission from actual interface gradients. Optical tanh changes only
the observation, never the conserved composition.

Run `julia --project examples/phase.jl` (set `INTERIOR_PHASE_OUTPUT` to a fresh
directory when repeating). It records initial, demixed, and coarsened
32³ states, metrics and a rendered view, rather than representing depth by
stacking a 2D pattern. This model lacks hydrodynamics, elasticity, thermal noise,
solidification and a validated physical timescale. Lamellar equilibrium would
require extra interactions or confinement; ordinary Cahn–Hilliard coarsens.

The core equations follow the conserved phase-field formulation documented by
[Idaho National Laboratory](https://tmap8.inl.gov/modules/phase_field/Phase_Field_Equations.html)
and the derivation and discrete conservation discussion in
[the University of Münster numerical methods notes](https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ss2016/num_methods_ii/lecture-cahn-hilliard.pdf).

The bundled 32³ run evolves 18,000 steps of 0.003 (nominal final time 54;
Float32 world-time accumulation reports 53.9897). Composition standard deviation
rises from 0.06894 to 0.66403; discrete energy falls from 8210.76 to 5549.54.
Mean interface gradient falls from 0.37758 at time 27 to 0.37006 at time 54,
showing modest late coarsening as growing deposit slows motion. The final signed
composition mass drifts by 0.19354 from −3932.15973, or 4.92×10⁻⁵ relative;
this measured Float32 limitation is exposed, not repaired in the state export.
`metrics.json` records the observations. These nondimensional outcomes show this
implementation demixes, and do not establish a calibrated material prediction.
