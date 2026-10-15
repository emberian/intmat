# Opening developmental feedback in constructed food

Our next search emphasizes two responsive loops: inherited development that can
move among internal states, and paid construction that appreciably changes the
habitat its producers inhabit. These are promising routes, not necessary
conditions for all persistent organization; stable differentiated states could
also matter. The present law contains both loops, but their operating scales can
make either effectively weak.

There is a concrete reason to investigate gain: the retained seed-211 full
constructive-refuges snapshot at time 3000 has mean absolute hidden activity
about .983, with about 94% of hidden components above .95 in absolute value.
This is the endpoint in
`outputs/constructive-refuges-20261015/211-full/frame-0030000.json`, not a
measurement of how often each individual switched state through its lifetime.

This is an analysis of the implemented equations and a reduced-mode experiment,
not a finding that a whole ecological trajectory oscillates. It changes no law,
checkpoint or running configuration. The companion
[helper](../studies/cellular_feedback_analysis.jl) and
[24 reduced trajectories](cellular-feedback-analysis.json) are reproducible with
`JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=.
studies/cellular_feedback_analysis.jl FRESH_OUTPUT.json`.

## Developmental gain and delayed inhibition

The cellular controller applies `tanh` to inherited weights plus plastic marks.
Hidden state approaches that output, with inhibitory adaptation subtracted inside
the activation; adaptation then approaches the updated hidden state. Actions use
the original controller outputs. Thus stronger adaptation changes later actions
through hidden state; it is not a direct penalty on secretion, growth or motion.
Plasticity updates hidden-output rows only. Action rows and the separate catalytic
readout are inherited rather than trained by this rule.

For one frozen recurrent mode, ignoring the extreme activation cap,

```
dx/dt = [tanh(w*x + I - g*a) - x] / tau_h
da/dt = (x-a) / tau_a
```

At an equilibrium `x=a`, let `s=1-x^2`. The local Jacobian has
`trace=(s*w-1)/tau_h-1/tau_a` and
`det=[1+s*(g-w)]/(tau_h*tau_a)`. A possible oscillatory instability crosses
`s*w=1+tau_h/tau_a`, provided `s*g>tau_h/tau_a`. This is a local condition for
this reduction, not a criterion for the entire network or a proof of a limit
cycle. Saturation matters: at `abs(x)>.95`, activation susceptibility `s<.0975`.
The actual implementation also caps the inverse activation of extreme Float32
outputs; that cap is retained in the reduced numerical helper.

Bounded plastic marks need not imply weak collective feedback. For constant
activity `x` and input `u_j`, an uncapped mark approaches
`M_j=eta*x*u_j/(lambda+eta*x^2)`. Its contribution summed over fixed sensory
inputs is `x*eta*sum(u_j^2)/(lambda+eta*x^2)`: repeated sensory activation can
reinforce the same hidden mode across many columns. At the default
`eta=.002, lambda=.02`, mark relaxation takes roughly 45–50 time units when
`abs(x)<=1`. Individual normalized-input equilibrium marks remain below .1,
yet their combined drive can matter well before the per-mark .3 cap.

The helper projects this law onto one recurrent mark and the aggregate marks
along constant external sensors, retaining discrete update ordering. It omits
chemistry, spatial contacts, neighbor modes, division and mutation; its aggregate
external-mark equation assumes individual caps stay inactive. With inherited
gain `w=1.2`, weak bias `.02`, response time `.8`, adaptation time `8` and
external sensor squared norm `8`, the last 400 of 1,000 time units give:

| Adaptation strength | Plasticity rate | Hidden range | Hidden standard deviation | Zero crossings |
| --- | --- | --- | --- | --- |
| .3 | 0 | −.413 to .561 | .346 | 18 |
| .3 | .002 | .922 to .922 | below .001 | 0 |
| 1.5 | 0 | −.490 to .512 | .363 | 54 |
| 1.5 | .002 | −.484 to .533 | .368 | 52 |

This makes a specific hypothesis worth testing: slow sensory reinforcement can
lock a recurrent mode that otherwise varies, and stronger delayed inhibition can
release it. The reduced baseline can already oscillate; no claim that strength
1.5 universally creates oscillation follows. Full networks have off-diagonal
couplings, moving sensory inputs, inherited variation and changing contact graphs.
All 24 reduced trajectories are finite and bounded, and zero-plasticity and
zero-sensor invariants pass. Halving the step to .05 preserves locking versus
cycling in the two plastic `R=8, tau_a=8` cases (the cycling extrema change by
about .01); these are exploratory dynamics rather than precision cycle estimates.

## Construction must reach an ecological operating range

In the constructed-food recipe, each raw feedstock is supplied at concentration
rate `.006`; soluble food has no continuing source. Catalyst release is
`kcat*e*F/(K+F)` in concentration units, with `kcat=6`, `K=.2` and initial raw
concentration near `F=2` before spatial masks. Matching local supply then requires
only `e≈.0011`. This estimate holds raw concentration fixed and ignores uptake,
transport and niche variation; accumulated or depleted raw food changes it.

Matrix permeability is multiplied by `1/(1+20*e_total)`, so total catalyst
concentration `.05` halves transport. If both species match supply, their total
is about `.0022`: efficient food unlocking can therefore operate with roughly
23 times less catalyst than a strong transport effect requires. Sorptive
capacity is six times the matching enzyme concentration: capacity near the
supply-matching `.0011` is only `.0066`, versus initial soluble concentration
`.3`. Low catalyst abundance can be a viable efficient strategy while producing
little transport geometry. That is a parameter-scale hypothesis, not an
impossibility result or evidence that evolution found an optimum.
The inspected time-3000 snapshot actually contains 7.06 catalyst mass units over
area 256, a mean total concentration about .0276. Its transport feedback is
therefore not absent; spatially concentrated catalysts can make it stronger still.
Local raw depletion and the observed distribution must be measured before using
the supply-matching estimate to explain a trajectory.

The clocks are also separated. Hidden response is `.8`, adaptation 6–8, bare
diffusion across exchange radius 1 is about `1/.2=5`, and free desorption takes
`1/.12≈8.3`. Bound-food equilibration depends additionally on local free food and
binding capacity. Mark relaxation is 45–50; catalyst lifetime is `1/.008=125`.
An enzyme-rich region can therefore retain a chemical history much longer than
one developmental cycle. Matrix resistance lengthens transport further. Actual
paid motion is bounded by speed `.3`, and can be much slower than its nominal
3.3-unit exchange-radius crossing time. These clocks permit delayed local
feedback, but do not guarantee that neighboring cells remain there to exploit it.

## Small discriminating interventions

First retain the existing constructed-food and strong-adaptation comparison.
For a clean follow-up, fix adaptation time at 8 and cross strength `.3/1.5` with
plasticity `0/.002`, keeping funding, founder genotypes and all chemistry matched.
This separates delayed inhibition from the slow reinforcing memory. The current
two recipes change adaptation time as well as strength, so their comparison
tests the combined regime rather than strength alone.

For the ecological loop, a separate four-arm contrast crosses catalytic rate
`6/.3` with matrix transport barrier `20/100`, holding the developmental regime
fixed. At raw concentration 2, rate `.3` requires each catalyst near `.022` to
match supply. With both species at that estimate, total catalyst `.044` gives
barrier-100 permeability near `.19`, versus `.82` for fast catalysis at total
`.0022`. This moves useful construction into a
range where it can reshape transport. It also makes failure to establish food
release plausible: retain extinction and bootstrap exhaustion as outcomes.
These are proposed new runs, not changes to existing pins.

Track living-cell hidden saturation, temporal variance, and same-ID sign changes
alongside actual secretion, catalytic concentration, free/bound food, local
permeability, births, deaths and material/work ledgers. Population-averaged
oscillation alone cannot distinguish individual cycling from replacement of
different cell types; variation among surviving cells alone cannot establish
persistent differentiation. Follow identities and descendants through divisions,
and compare local environmental phase lags before making either claim. A useful
result could be responsive development, sustained construction, or a clear
failure mechanism that identifies the next intervention.
