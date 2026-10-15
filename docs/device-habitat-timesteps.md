# Fixed-duration habitat timestep study

The paid 128³ pilot used dt=.001. This is an experimental discretization choice,
not the largest permitted step. With physical extent4, spacing1/32 and diffusion
.001, the explicit diffusion bound is .146484375. The cellular bound also caps dt
at .25 and response_time/2; the composed process applies the raw-pool bound too.
Those bounds alone do not establish accuracy. The unchanged mechanics stage
additionally requires `2*force_row_sum*(dt/16) <= .45` and
`speed*(dt/16) <= .2*radius` for every moving cell at runtime. The recipe
response_time=.8 gives a neural bound of .4; the global .25 cap and diffusion
therefore dominate the static128³ limit.

`studies/device_habitat_timesteps.jl` creates the same recipe independently for
each timestep. Its default reference is Float32(.0005), then exact powers-of-two
multiples .001/.002/.004/.008, with integer step counts and identical actual
rounded physical duration. Geometry, physical radius, rates, seed and16 mechanics
substeps remain fixed. It compares full free/raw/catalyst/bound fields, neural
states, metabolites, periodic positions, identities and accounting residuals.
Reported simulated-time per wall second is distinct from steps per second.
Construction and detached final snapshots are outside the timer; status checks
are inside it. First-use compilation is included unless separately warmed, so
initial CPU gates are correctness checks rather than throughput evidence.

The first short interval should precede births. A longer interval requires
population/lineage/bond and field-distribution comparisons: divergent lifecycle
ordering makes individual-coordinate errors an inadequate ecological criterion.
No larger timestep is selected for production from the analytic bound alone.
The old pilot and immutable source archives are unchanged.

The initial CPU-KA gate completed all five stepsizes at native8³ with two founders,
capacity32 and duration .00800000038. No births occurred and identities matched.
At dt=.008, relative L2 error against dt=.0005 was9.73e-5 for free fields and
8.44e-4 for catalyst; maximum periodic position difference was1.14e-5. Material
residuals stayed below4.44e-8 and potential/work/heat residuals below1.73e-5.
This tiny coarse-grid check establishes execution and comparison machinery,
not accuracy at128³ or a recommended production timestep. Values are retained in
[device-habitat-timestep-cpu.json](device-habitat-timestep-cpu.json).

## Native AMD 128³ gate

The frozen d4 numerical source with only the carry helper from7d3caeb and study
frome93da97 completed all five candidates at fixed duration .01600000076,512
founders and capacity8192. All stayed finite with unchanged identities/population.
The study body took49.19s including construction, kernel compilation and full
reads; Julia startup, backend imports and package precompilation precede that
timer. This is not process-launch-to-exit timing.
The first .0005 reference includes cold compilation, so its rate is not comparable
to the subsequent candidates as a timestep speedup.

| dt | Free relative L2 | Catalyst relative L2 | Bound relative L2 | Model time / wall second |
| --- | ---: | ---: | ---: | ---: |
| .001 | .000553 | .000511 | .0238 | .00767 |
| .002 | .00154 | .00137 | .0721 | .01546 |
| .004 | .00359 | .00317 | .176 | .03130 |
| .008 | .00804 | .00677 | .480 | .06339 |

These short timers include batch status synchronization but exclude construction
and full observations. Larger timesteps buy more model time per wall second, yet
the newly accumulating bound pool is sensitive to splitting: at .008 its relative
L2 discrepancy is48%, despite maximum absolute voxel error only1.04e-8. That is a
reason to retain a bound-pool accuracy criterion rather than approve .008 from
free-field errors or stability alone. Maximum periodic position error was6.10e-5
for all four coarser candidates. Material residuals were below4.83e-6 and
potential/work/heat residuals below3.87e-4. No production timestep is changed.
A longer fixed-duration study, a finer reference and repeated warm timings are
needed before selecting an ecological production step. Exact metrics, runtime
and source hashes are in[device-habitat-timestep-amd.json](device-habitat-timestep-amd.json).

The bound-food discrepancy needs physical context. At time .016 the reference
contains only .0005101 total bound mass; dt=.008 contains .0006077, a difference
of .00009759 mass units (19.1% of this very early pool, versus48% field L2).
The largest bound concentration discrepancy is .0003398, after dividing voxel
mass by native voxel volume. Catalyst inventory changes from .315548 to .315321;
work from1.690139 to1.689185 and heat from18.08504 to18.08405. These observations
are consistent with an early splitting/startup bias and do not by themselves
establish consequential ecological error. They motivate the longer comparison,
not an automatic rejection of a larger step.

The next prepared study uses fixed duration .256, reference dt=.00025 and
candidates .0005/.001/.002/.004/.008 with the same materialized support strategy
throughout. A tiny coupled warmup compiles the path; each repetition constructs
fresh same-seed states and times complete trajectories. Factor order reverses
on alternating repetitions. Full concentration and total-mass errors, position
over exchange radius, and work/heat differences accompany relative L2 metrics.
There are2016complete steps per repetition; the repeat count and wall budget
will be chosen from the matched dense coupled profile before launching.
