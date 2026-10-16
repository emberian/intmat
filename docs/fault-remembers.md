# The fault remembers

A fresh, bounded material artwork; historical fracture states remain untouched.
`studies/fault_remembers.jl` constructs 5,670 nodes in seven laminae, with
63,209 central-force bonds. The shared Julia `BondNetwork` implements actual
loaded dynamics, damping, measured tensile rupture and an energy ledger.
Lamination changes both stiffness and conductance. Weak, offset flaw seeds
compete during tensile and shear grip loading; a small spatially varying
out-of-plane load curls the plates.

Every ten mechanical advances, 24 damped Jacobi sweeps relax a voltage field
on the surviving graph. The grips impose zero and unit potential. Current
squared times resistance supplies a nonnegative local Joule record; its
accumulated dose reduces failure thresholds. Failed bonds cease conducting.
This exploratory coupling gives material memory a causal role in subsequent
rupture. It is an artistic law with an approximate relaxed voltage field,
not calibrated dielectric breakdown, a fully converged circuit solve or a
thermodynamically closed electrothermal model.

The photographed geometry is a face complex: a lamina quad survives only if
all four actual perimeter bonds survive, and uses their deformed coordinates.
Actual surviving axial ties reveal connections across layers. Thin plate
thickness, beveled edges, studio illumination and dose/strain-to-metal optical
transfer are photographic choices. No crack curves are drawn in the renderer.

Run with the existing Julia project:

```
JULIA_NUM_THREADS=2 julia --project=. studies/fault_remembers.jl
/Applications/Blender.app/Contents/MacOS/Blender -b -t 2 --python scripts/render_fault_remembers.py -- outputs/fault-remembers/state.json outputs/fault-remembers/final.png --size 1400 --samples 64
```

The generator validates the final mechanical state and retains measured energy,
source hash and bounded duration in `outputs/fault-remembers/run.json`. The
private full graph state stays outside the published gallery; the finished
photograph and compact run/optical records suffice for viewing.

For the completed seed-7126 study: 1,450 steps at dt=0.0007, 6,383
broken bonds, maximum dose 0.0441900527. The mechanical energy/work residual is
0.00348623 against work 3.85618 (0.0904%); this ledger excludes the explicitly
open electrical supply and phenomenological threshold law. Final fields are
finite, dose and temperature are nonnegative, the relaxed voltage lies in
[0,1], and exported damage counts agree with the run record.
