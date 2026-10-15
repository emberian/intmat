# Spatial modes after the turnover bake

The native 128³ specimen is a useful transient but coarsens. These separate
parameter tests use the same `5377470` turnover law; they do not alter that
specimen, its rates or its source pin. They distinguish spatial organization
from the number of activity crossings.

Three 64² tests increased matrix feedback, delayed memory feedback or both.
All retained repeated material turnover, but all eventually concentrated their
excitation spectrum at box scale. Memory delayed coarsening without preventing
it. Their recipes and `living-tissue-echo-results.json` retain these outcomes.

A low-gain test (gain 1.8, capacity .32, inhibitory diffusivity 4) stabilized the
uniform neural mode while admitting a weakly growing finite wavelength around
5.5 once biomass increased. The nonlinear 64² test retained a spectral peak
around 4.8 through t100, with little high-frequency power, but excitation was
weak (final std .014) and no strict material reopening occurred. This checks a
mechanism rather than meeting the intended art/material behavior.

A stronger recipe keeps gain 2.8 and capacity 1.2 and changes inhibition time to
.5 and diffusivity to 4. Its frozen reaction has a stable uniform mode and a
growing finite band around 3.0–3.8 for B above approximately .12. The actual 64²
test produced strong excitation and material reopening, but 30–45% of excitation
fluctuation power lay beyond .7Nyquist. Halving dt from .005 to .0025 preserved that
high-frequency content: at t20 the fractions were .426/.420, and at t40
.305/.309. This is not an acceptable resolved art result.

The separate 128², dt .00125, t25 probe completed with the same physical start:
restricting initial fine N/B/x cell averages to 64² differs by at most
5.96e-8/2.98e-8/7.45e-9. Its dominant radial wavelength stays 2.667 through t5–25,
about 14 fine cells, while the coarse t20 peak was .96. Thus finite-band structure
survives refinement and the coarse short-wave peak was spatially biased. Fine
power beyond its own .7 Nyquist still grows from .00124 at t10 to .0794 at t20 and
.1535 at t25; late sharpening has not been shown converged. At t20, halving dt
on 64² changes spent material 186.55→185.83, whereas 128² gives 178.14 and changes
biomass 201.93→226.33. Time refinement alone would have missed this spatial bias.

The fine job used unchanged immutable source/environment 5377470 and a separately
hashed read-only recipe, 120s/8GiB/no swap/CPU10–11/nice15. Its 20k steps and 25 later
captures took 1.7783s, excluding initial capture, JIT/warmup and final
checkpoint/record/maps. Artifacts total 16.42MB. Exact parameters, hashes, initial
restriction errors and compared spectra/inventories are in
[`living-tissue-mode-results.json`](living-tissue-mode-results.json).

## Frozen reaction bound

This result concerns the continuous three-signal reaction with material and
payment factors held fixed. It omits derivatives of material transport,
metabolism, renewal and phase work. It is not a full material equilibrium,
discrete-time stability theorem or nonlinear pattern guarantee.

Let s=B/(.12+B)<1, v=s·paid·N/(.15+N), and let z be the tanh slope at the chosen
signal operating point. Then 0≤v≤s and 0≤z≤1. For the fast recipe, define

```
a = v(2.8z−1) − (1−s),   b = 2.3vz,   c = .3vz
J = [ a  −b  −c;  2  −2  0;  .025  0  −.025 ]
```

Its characteristic polynomial is λ³+Aλ²+Bλ+C, with

```
A = 2.025−a                         ≥ .225
B = 2(b−a)+.025(2−a+c)              ≥ v+2(1−s)+.005
C = .05(b−a+c)                     > 0
C                                   ≤ .05[v+(1−s)]
AB−C                                ≥ .175v+.4(1−s)+.001125 > 0
```

These bounds follow from a≤1.8, b−a=v(1−.5z)+(1−s), and
b−a+c=v(1−.2z)+(1−s). They imply that every eigenvalue has negative real part:
positive coefficients exclude a nonnegative real root. If a complex pair had
real part α≥0 and the remaining root were −γ, positive A gives γ>2α; direct
factorization gives `AB−C=−2α[(α²+β²)+γ(γ−2α)]≤0`, contradicting the bound.

The simulation's finite Float32 coefficients have ample margin for the actual
concentration range. The support assumption s<1 and fixed-payment derivative
are part of the statement, rather than assertions about every representable
floating-point state. `studies/living_tissue_mode_hypothesis.jl` exports actual
axis/diagonal stencil dispersion for both 64² and 128² and 2D/3D directions.
Spatial diffusion can destabilize a finite band even though the local reaction
is stable; the failed coarse nonlinear test demonstrates why that fact alone
does not establish resolved evolving organization.
