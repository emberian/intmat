# Acoustic amber

A localized transducer launches waves through the specimen. Their running intensity stimulates the excitable activity field, which gates conversion of nutrient to mineral. Deposited mineral slows subsequent waves, changing interference and trapping pressure near the source. This is an invented acoustic/material feedback law, with externally supplied transducer power; it is not calibrated sonochemistry or a full elastic solid.

The pressure-like scalar evolves by

```
p_tt + 2γ p_t = div(C(deposit) grad p) + S(x,t)
C_face = c₀² / (1 + β (deposit_left + deposit_right)/2)
S = amplitude × Gaussian(x-center,radius) × sin(2π frequency t)
```

The face coefficient is symmetric and positive for nonnegative deposits. A centered second-order leapfrog stencil advances pressure and its previous value, with centered damping. The multidimensional Courant limit is `dt ≤ 0.9 / (c₀ sqrt(sum(1/hᵢ²)))`; mineral only reduces wave speed, so this is an upper-speed bound. This method and its Courant restriction follow standard wave finite differences, as described in [MIT's wave-equation finite difference notes](https://ocw.mit.edu/courses/18-325-topics-in-applied-mathematics-waves-and-imaging-fall-2015/resources/mit18_325f15_appendix_b/). [deal.II's wave tutorial](https://dealii.org/9.3.3/doxygen/deal.II/step_23.html) provides another primary explanation of wave propagation and time discretization.

The intensity variable is an exponential running average of squared pressure. Activity receives `gain × intensity/(scale + intensity)`, an explicitly bounded artistic stimulus. Nutrient plus deposit remains conserved because conversion uses the existing limited reaction. Acoustic energy is dissipated by damping, supplied by the source, and altered by changing material coefficients; total energy across acoustics and mineralization is not tracked. Boundary conditions are periodic. This means returning waves are part of the modeled domain, rather than reflections from a glass boundary. The pressure scalar is copied into the world's diagnostic pressure; this process does not solve velocity or strain.

Run the specimen with `julia --project=. bin/interior.jl study acoustic`. It uses a 48³ grid and records its parameters and measurements in `outputs/acoustic/run.json`, full pressure/intensity/deposit/activity Float32 arrays, an optical volume, and a native volume render. The selected time is 6: the longer run at time 14 converted over 90% of the resource, obscuring the spatial wave memory. Selecting an intermediate specimen preserves the experimentally generated structure rather than implying a guaranteed long-term pattern.

Thirteen targeted tests establish propagation from a localized pulse, finite bounded pressure, Courant rejection, fixed-timestep rejection, mineral-to-wave feedback, acoustic-to-deposition feedback, conserved feedstock, and strong damping reducing discrete mechanical energy. A 10×10×8 CPU/Metal comparison after twenty steps produced maximum pressure difference 2.91×10⁻¹¹ and activity difference 7.11×10⁻¹⁵. This validates the same propagation kernel on the actual Apple GPU; it does not establish grid convergence or physical calibration.

`AcousticField` owns pressure, previous pressure, scratch, and intensity envelope arrays. Its fixed timestep is part of its state. Construct it for the chosen world/backend; do not reuse it with a different grid. `seed_acoustic!` initializes a localized displacement with zero initial velocity using a half-acceleration Taylor expansion. Use a fresh world when seeding this initial-value experiment.
