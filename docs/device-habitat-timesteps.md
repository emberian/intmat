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
