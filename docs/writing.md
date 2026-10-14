# Focused thermal writing

`ThermalWriter(world; ...)` owns temperature, scratch temperature, and irreversible dose arrays. A Gaussian beam follows a deliberate stationary, three-turn helix, or closed 3D Lissajous toolpath. The source raises temperature; diffusion and cooling spread and remove it. Temperature above a threshold accumulates dose and conservatively converts local nutrient to deposited scaffold. Cooling cannot erase dose. Subsequent dissolution processes may erase deposited material while the exposure record remains.

This is phenomenological, nondimensional laser writing. The smooth `pulse_rate` envelope is a controllable heating modulation, not resolved ultrafast laser pulses. Neither optical breakdown nor sonoluminescence is simulated. Rendered glowing loops are an interpretation of written inclusions, not computed optical caustics.

`width`, `heat`, `diffusion`, `cooling`, `threshold`, and `conversion` control source and material response. `center`, `radius`, `height`, `period`, `start`, and `stop` control the toolpath and shutter. The heat solver timestep must satisfy the combined diffusion/cooling bound returned by `stability_limit`. Path time is simulation time; render capture timing does not change the writer.

Compose `(writer, neural, transport, deposition)` to let the persistent neural field sense the written scaffold through its deposit input. Beam writing changes the neural material's environment rather than replacing its recurrent rule. Material mass is nutrient plus deposit and is conserved by writing.

Run `julia --project=. --threads=4 examples/writing.jl` to generate paired 48³ specimens, with and without neural dynamics, under `outputs/writing`. Each run uses a fresh directory so complete state exports cannot overwrite earlier observations. An optional first argument specifies the output directory. Each specimen contains all physical and process-owned fields as Float32 snapshots, a volume export, run metadata, and a native volume rendering.
