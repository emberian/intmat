# The field inside the fracture

A new crystalline interpretation of **Fault script** and **Ice cleavage**. A genuinely three-dimensional spring material is driven in tension and shear. Its surviving bonds conduct a quasi-static voltage field; Joule dissipation accumulates a cooling dose that lowers subsequent tensile failure thresholds. Failure stops conduction on that bond, redistributing the next electrical field. The photograph exposes this computed electrical structure inside fractured support.

The gold curtains are **potential iso-surfaces**, not additional physical sheets or simulated deposits. Six tetrahedra subdivide each reference-grid cell; a tetrahedron is displayed only while all six of its actual mechanical bonds survive. Voltage is interpolated linearly inside those surviving tetrahedra and observed at nine fixed levels. Actual deformed node coordinates determine their shape. Normals consistently point toward increasing voltage; this avoids false alternating optical highlights without smoothing geometry or fields. A .004 optical foil thickness and a 16-shade accumulated-dose palette make the contours visible. These are declared optical assignments.

The etched-violet back surface is the actual exposed boundary of surviving tetrahedra, exactly clipped at y=.18 as an optical section. There is no smooth enclosing box, procedural noise, painted crack mask or particle fill. Quenched node disorder, a sinusoidal conductivity grain field, the shared lower-half-plane weak seam and two additional small notches are authored inputs. Subsequent failures, deformations and voltage redistribution are computed. The conservative six-live-edge display criterion can omit regions still connected by individual surviving bars; it is an observation complex, not an independently simulated continuum volume.

## Processes and approximations

The study calls production `BondNetwork` for central-force motion, damping, irreversible strain failure, grip loading and its mechanical energy ledger. No shared simulation source changes are needed. The electrical solve is study-only, matrix-free diagonally preconditioned CG on the actual surviving graph. Conductances use reference bond lengths, fixed designed grain coefficients and failure state; they do not change continuously with elastic stretch. Floating components get an arbitrary constant zero gauge and produce no spurious current or Joule heat.

Every solve independently enumerates edge currents to check the normalized node residual, with a 2e-6 rejection gate. The same edge conductance defines Joule power. Electrode power and summed bond Joule power are compared. The dose is a dimensionless low-pass power/exposure surrogate, not Kelvin temperature. Prescribed electrodes and grips supply power; no finite electrode reservoir, latent-heat debit or closed coupled thermodynamic model is claimed. Electrical and threshold updates occur at finite cadence, with mechanics advanced between updates.

## Bounded evidence

The small matched 24×18×10 pilot has 4,320 nodes and 48,800 bonds. After 750 steps, electrical weakening gives 3,570 broken bonds versus 220 without weakening. Final electrical power is 166.67 versus 245.94 under identical unit electrode potentials. The coupled mechanical work-ledger residual is .0421%, versus .000575% for the control. Both arms use the same retained pilot source. The overly late 1,300-step image is retained privately as a rejected state: most voltage drop had crowded into tiny remaining necks and its broad contour interior vanished.

An independent dense linear-system fixture checks the graph solver to maximum solution error 3.41e-12, with relative source/Joule power defect 2.05e-16. It explicitly includes an isolated floating node and checks the maximum principle.

The selected fresh 48×36×20 run has 34,560 nodes, 419,228 bonds and 28,606 broken bonds. It advances 1,554 steps over nominal physical duration .525 in 106.02 seconds on CPU2. Final independent normalized residual is 1.47e-7, relative source/Joule power defect 6.64e-8, and mechanical energy-ledger residual .0672% of external work. One detached component has a constant gauge. The observed complex contains 89,230 surviving and 98,300 omitted tetrahedra, and 10,370 potential-contour faces.

The spacing ratio relative to the pilot is .482873. Edge stiffness scales with this ratio, node mass and distributed load with its cube, and conductance cross-section with its square before division by reference bond length. Per-node heating is divided by the volume ratio. The timestep scales with spacing and is rounded downward to reach the same duration; physical loading rates and conductivity grain lengths remain fixed. The finer realization resamples quenched disorder and is **not a demonstrated continuum-converged solution**. The scaled variant also uses actual elapsed electrical-update intervals; the retained pilot/control pair uses its original fixed interval approximation.

Sources: `studies/voltage_between_pieces.jl`, `scripts/check_voltage_pieces.jl`, `scripts/render_voltage_pieces.py`. Compact evidence: `docs/voltage-pieces-controls.json`. Private native state, exact source copy and reviewed photography are under `outputs/voltage-pieces-48-scaled`. Original gallery pieces remain unchanged.

## Subsequent boundary handling

After the approved photograph, the study solver gained explicit zero identities for an empty free-node residual and an empty electrode-power sum. A fully isolated eight-node graph now returns residual and power zero; an all-electrode graph returns residual zero and matching Joule/electrode power 16. Both require zero CG iterations. The independent dense-system fixture retains its 3.41e-12 error. Evidence is `docs/voltage-pieces-boundary-check.json`.

This is subsequent boundary handling only. The artwork was not rerun and its provenance was not changed. Its exact recorded study remains privately retained as `outputs/voltage-pieces-48-scaled/source-study.jl`, SHA256 `e952444f2dcf35b5f25d8af75b035c19ad6c4902299469d62076385f7227c4c3`.
