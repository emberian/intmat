# A yielding web remembers its tear

`BondNetwork` is a third carrier: the material is a deforming graph, rather than
fixed grid cells or migrating agents. A 10×10×10 disordered network has local
edge, face-diagonal and body-diagonal springs. It is a modest mass–spring
specimen, not a calibrated finite-element solid. The lattice gives construction
coordinates; geometry and velocities subsequently live in independent arrays.
No periodic links connect opposite faces.

For bond e connecting nodes i,j, let r=xⱼ−xᵢ, l=|r| and reference length L.
Its elastic energy is ½(1−d)k(l−L)², and force on i is
(1−d)k(l−L)r/l; j receives its negative. Each node evolves by Newton's
second law. This is the standard energy-based spring model described in
[Baraff, Witkin and Kass's particle dynamics and energy notes](https://www.cs.cmu.edu/~baraff/pbm/).
Our failure rule is a deliberately simple bespoke constitutive law: the first
observed positive extension strain exceeding a stored threshold sets d=1
forever. The current bond's elastic energy moves into the damage reservoir.
Compression does not initiate fracture, broken bonds exert no force, and the
maximum historical strain remains recorded. Damage is neither erased on unloading
nor inferred from an image.

A weakened half-plane seeds a notch. It changes local failure thresholds, not
initial shape, opacity or topology. In the example, opposite x faces move apart
with a slight transverse shear. Interior nodes accelerate, redistribute tension,
and lose connectivity. The surviving spring web deforms around the opening.
Damage endpoint emission is an optical scar, not mechanical repair. There is no
nutrient uptake, healing, volume preservation, contact or self-collision yet;
those require additional constitutive mechanisms rather than relabeling springs.

Coordinates are physical world units centered at zero. k has force/length units,
m has mass units, and damping is inverse time. `stretch_rate` and `shear_rate`
are the total relative speed of the two grips, length/time; each moves at half
that speed. `load[:,i]` is an external force. Velocity Verlet integrates free
nodes with exact exponential damping half steps. Grip nodes are prescribed
positions with stored zero velocity (they are constraints, not inertial nodes).
The exact elastic energy change from imposing grip displacement is external work;
constant nodal load contributes F·Δx. Elastic energy removed on fracture and
kinetic energy removed by damping have separate reservoirs. `total -
external_work` should stay near its initial value; its residual is integration
error, not a claim of exact mechanical conservation.

`stability_limit(w,p)` takes a conservative tangent stiffness Gershgorin estimate
including transverse stiffness of compressed springs, then uses 0.35/√bound.
It is a local explicit integration guard, not a nonlinear global stability theorem.
Very collapsed springs shrink the admissible step. `step!` validates against it;
the example uses 0.0005 time units. For force-symmetry tests use `bond_forces`,
which excludes loads and constraints.

`record_process_fields` includes exact integer node/bond IDs, reference geometry,
current positions/velocities, edge endpoints, rest lengths, thresholds, peak strain,
damage, grips, loads, all parameters and energy ledgers. Broken edges stay in the
graph with damage=1; IDs never depend on which bonds remain visible. This lets
recordings inspect graph memory without interpreting a rendered picture.

`mechanical_optical_fields` samples the real current live segments into volume
fields and samples broken endpoints as a scar. It does not render the original
rest lattice. Fine bonds can alias on coarse grids; `radius` is an optical sampling
width. `export_bond_graph` writes current OBJ vertices and live `l` segments for
a cylinder-based external renderer. Example graph JSON retains broken edges too.

Run `julia --project bin/interior.jl study mechanics`. The example saves eleven actual
simulation snapshots, graph geometry, volume fields and renders under
`outputs/mechanics`, with energy and fracture counts in `run.json`.
