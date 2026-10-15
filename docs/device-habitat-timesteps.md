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

## Matched materialized support at time .256

The longer native128³ study completed11 trajectories: one .00025 reference and
two fresh trajectories at each coarser timestep, all with materialized support.
All remained finite with the same512 identities. The study body took184.20s;
external launch-to-exit time was202.60s, including Julia/package startup. The
first failed driver launch is retained separately and contributed no simulation.

The reference now contains .434539 bound mass and4.95666 catalyst mass. Bound
food is no longer merely a near-zero startup pool:

| dt | Bound mass difference | Bound field relative L2 | Free field relative L2 | Median model time / wall second |
| --- | ---: | ---: | ---: | ---: |
| .00025 | reference | reference | reference | .00636 |
| .0005 | +.00187 (0.43%) | .00727 | .00456 | .01225 |
| .001 | +.00551 (1.27%) | .0170 | .00407 | .02589 |
| .002 | +.01290 (2.97%) | .0409 | .00424 | .04818 |
| .004 | +.02881 (6.63%) | .102 | .00521 | .10374 |
| .008 | +.08197 (18.9%) | .330 | .00796 | .20710 |

Relative bound-inventory percentages use the reference pool. Timers include
complete physical steps and batch status checks; construction and final detached
reads are outside individual timers. The two .002 runs took4.93s and5.76s,
so its rate has visible variation. The other paired durations were closer.

Compared with .001, .002 provided1.86× model-time throughput while changing the
bound inventory by2.97% relative to the fine reference (4.09% field L2). That makes
.002 a plausible exploratory tradeoff for a declared few-percent chemical
accuracy target, rather than a demonstrated equivalent ecological trajectory.
The .004/.008 bound-pool deviations are substantial despite modest free-food
errors. A production timestep has not been changed.

The .00025 reference is not proved converged: its potential/work/heat residual
is .02583 over initial potential4085.76 (6.32e-6 relative), and material residual
.000226 over1456.64. Finer Float32 accumulation increases rounding exposure.
Position differences are not monotone with timestep; coordinate increments can
approach Float32 resolution. At .001 the maximum periodic difference is .00174,
about0.97% of exchange radius; at .008 it is .00235 (1.31%). These are comparison
errors, not evidence that the fine trajectory is exact. The test interval has
no births, so longer ecological/morphological equivalence remains unestablished.

[Exact compact evidence](device-habitat-timestep-materialized-amd.json) retains
per-species mass differences, mass L1 errors, maximum concentration errors,
work/heat differences, runtime/source pins and both timing samples. Numerical
source remains d4d152d with only the7d3caeb support helper; study/test/fixture
were overlaid fromd958d79 with a separately hashed driver parse/CLI correction.
