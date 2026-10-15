# The nacre fault

Three perforated folds share a changing mineral inventory. A pressure field passes
through their displaced openings; eight recurrent channels per voxel deposit and
erode the folds, and the resulting mineral slows subsequent pressure propagation.
The artwork observes this exchange through nacre glass, warm copper and refractive
pressure fronts.

The selected specimen is `collection-20261015/nacre-fault-openwork`. Its still
observes **t = 8**, rendered at **2000² / 128 samples**. Its film uses **65 actual
256³ samples**, t = 0–16 at intervals of 0.25, with a fixed camera and lights.
Playback at eight frames per second lasts 8.125 seconds; no intermediate geometry,
camera orbit or reversed playback supplies its motion. The film is a finite
excerpt, not a seamless cyclic simulation.

## What is simulated

`studies/interfacial_choir.jl` composes the existing Julia `AcousticField` and
`RecurrentMaterial` on Metal. The domain is periodic, extent 64 in each direction,
with spacing 0.25 and a fixed Float32 timestep of 0.025. The 640-step run captures
three full native scalar volumes every ten steps: mineral, fast recurrent activity
and pressure. Full final latent channels, nutrient and acoustic history are retained
locally alongside the source snapshot and captured fields.

The three folded, perforated sheets and coherent latent phases are deliberately
seeded geometry. They are not claimed to arise from noise. Nutrient and mineral
initially share this support. The recurrent law transfers nutrient into mineral and
back, and nutrient diffuses. A localized sinusoidal acoustic source supplies external
drive. Acoustic face coefficients depend on current mineral. After acoustic and
recurrent advances, a bounded pressure signal relaxes the first recurrent pair
with strengths 0.75 and 0.15; it affects the next material update. These are designed,
untrained dynamics, not a calibrated biological or acoustic material.

The selected production trajectory remained finite, with nonnegative nutrient.
Its largest observed material balance error was **6.67×10⁻¹⁰ of initial inventory**.
The approximately 56-second bake includes capture transfers and writes and is not
a hardware benchmark. Each production run used a separate immutable source and
GPU-environment snapshot; none of the long-running ecological experiments was
changed.

A matched 24³ / 160-step fixture compared CPU and Metal and independently disabled
each new feedback direction. Maximum CPU/Metal differences were 1.79×10⁻⁷ for
mineral and nutrient, 1.08×10⁻⁶ for latent state and 4.03×10⁻⁶ for pressure. Removing
pressure-to-memory coupling changed mineral by as much as **0.337**; removing
mineral-dependent wave speed changed pressure by as much as **0.824**. Thus both
connections have a measurable effect in this fixture. This is a coupling check,
not resolution convergence or evidence of learned behavior. See
[the fixture record](../benchmarks/2026-10-15-nacre-fault-verification.json).

## What the optics choose

`scripts/render_interfacial_choir.py` meshes the recorded mineral at 0.16 and 0.60.
The lower contour is translucent nacre; the higher contour has a rough copper
material. The activity sampled at each surface vertex chooses its tint. Saved
pressure chooses a thin-film thickness of `380 + 170 tanh(4 pressure)` nanometres
in the shader. That thickness is an artistic transfer, not a simulated film.

A third surface follows recorded positive pressure at 0.12. It is displayed as a
weakly emitting refractive film, cropped to radius 0.86 and |y| < 0.58 in normalized
observation coordinates. A cubic opacity fade across the outer 0.07 radial and
0.05 depth widths hides the discretized face-cut rim. This is a
pressure visualization, not a claim that sound creates a separate glass phase or
visible light. Cycles transports light through the chosen surfaces; neither optical
feedback into the material nor physical structural coloration is solved. There is
no procedural displacement or generated image substitution.

## Reproduce

Use a fresh output directory. The existing Metal project is described in
[the accelerator notes](accelerators.md).

```sh
julia --project=gpu --threads=2 studies/interfacial_choir.jl \
  experiments/art/nacre-fault-openwork.toml outputs/my-nacre
julia --project=gpu --threads=2 studies/interfacial_choir.jl \
  --verify outputs/nacre-verification.json
/Applications/Blender.app/Contents/MacOS/Blender -b -t 4 \
  --python scripts/render_interfacial_choir.py -- \
  outputs/my-nacre outputs/my-nacre/still --steps 320 \
  --size 2000 --samples 128 --orbit .55 --elevation .36 --save-scene
/Applications/Blender.app/Contents/MacOS/Blender -b -t 4 \
  --python scripts/render_interfacial_choir.py -- \
  outputs/my-nacre outputs/my-nacre/animation \
  --size 960 --samples 48 --orbit .55 --elevation .36
```

The normal Blender Metal environment succeeds. Forcing
`CYCLES_METAL_ADAPTIVE_COMPILE=0` and `CYCLES_METAL_SPECIALIZATION_LEVEL=0` failed
shader compilation in this installation; those overrides are not part of the
recipe. The still also runs through the collection queue's `interfacial-choir`
renderer dispatch. The ordinary reviewed-hash curator publishes its photograph;
[the animation curator](../scripts/curate_animation.py) accepts its reviewed frame
sequence and exact simulation timestamps. Raw volumes, Blender scenes and failed
or unselected studies stay outside the public bundle.

[Watch the film](https://emberian.github.io/intmat/animations/#nacre-fault) · [View the still](https://emberian.github.io/intmat/#collection-20261015/nacre-fault-openwork).
