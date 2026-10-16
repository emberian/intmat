# The fault remembers — stored motion

This is a new bounded capture study using the same seed, law and final time as
the approved still. It does not resume or modify an historical trajectory.
Forty actual states are captured at rounded uniform step intervals from step 0
to 1,450 (physical time 0 to 1.015). Capture instrumentation performs no state
update. The final arrays match the approved still's full stored state exactly:
position, reference, bonds, damage, peak strain, dose, temperature and voltage.
The source state hash is
`a08b764f3071b452e2e65614dd6bd1d7bcfdea0659517f3aef96c49c3d521936`.

The recorded graph retains full Float64 coordinates and fields privately in
`outputs/fault-remembers-motion/states`. `capture-verification.json` retains
individual source hashes, steps, times and damage/dose progression. Its checks
require finite fields, voltage in [0,1], irreversible damage, nondecreasing
per-node dose, and exact final-state equality. The simulation remains Julia.

Motion optics use one fixed camera and a fixed global dose maximum of
0.04419005268498425. The intact teal laminae gradually acquire the bronze record
and then expose their changing interior. Each photograph uses the actual
stored deformed face complex and surviving ties. No coordinates, cracks or
fields are interpolated between stored states. Optical material thickness and
lighting remain the choices described in `fault-remembers.md`.

The three-keyframe CPU prototype (640², 16 Cycles samples, two threads) took
51.1, 43.4 and 27.7 seconds. Those exact reviewed frames are reused in the full
sequence. The full sequence has forty actual photographs, played at eight
frames per second; physical time and playback time are distinct. Encoding uses
the existing `curate_animation.py` validation and publication workflow.

```
INTMAT_FAULT_CAPTURE=1 INTMAT_FAULT_OUTPUT=outputs/fault-remembers-motion JULIA_NUM_THREADS=2 julia --project=. studies/fault_remembers.jl
python3 scripts/verify_fault_motion.py outputs/fault-remembers-motion outputs/fault-remembers/state.json
/Applications/Blender.app/Contents/MacOS/Blender -b -t 2 --python scripts/render_fault_motion.py -- outputs/fault-remembers-motion outputs/fault-remembers-motion/prototype --prototype --size 640 --samples 16
/Applications/Blender.app/Contents/MacOS/Blender -b -t 2 --python scripts/render_fault_motion.py -- outputs/fault-remembers-motion outputs/fault-remembers-motion/frames --reuse outputs/fault-remembers-motion/prototype --size 640 --samples 16
```

Full graph dumps and rendered PNG sequences remain private outputs. Publication
contains only the reviewed compact film, poster and source/time mapping.

The completed film contains forty reviewed 640×640 photographs at 8 fps,
duration 5.000 seconds. Publication is 646,083 bytes including the film,
poster and metadata. The full render measured 819 seconds excluding the
three-keyframe prototype (122 seconds); host contention makes these workflow
costs, not a hardware benchmark. FFprobe confirms forty decoded frames at
640×640 and the expected five-second duration. Decoded initial, middle and
final keyframes were reviewed after H.264 encoding.

`prepare_fault_film.py` verifies the complete rendered/source-state mapping,
prepares four contact sheets and writes the sequence manifest. Visual approval
covers the manifest hash in `art-reviewed.json`; the shared curator then checks
all forty PNG hashes, encodes the film, verifies its decoded frame count and
publishes atomically. No raw graph arrays or PNG sequence are published.

```
python3 scripts/prepare_fault_film.py outputs/fault-remembers-motion outputs/fault-remembers-motion/frames
# Review the contact sheets and actual keyframes, then record their sequence hash.
python3 scripts/curate_animation.py outputs/fault-remembers-motion/frames/sequence.json --destination site/animations --max-width 640 --budget-bytes 1000000
```

The last command preserves existing publications and requires a fresh version
ID if the publication directory already exists.
