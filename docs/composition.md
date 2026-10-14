# Ordered experiments

`Experiment` names processes and uses the existing `step!` tuple order. Every
process sees the state left by its predecessor. The world clock advances once
per complete step. Construction and each run preflight all processes with the
existing parameter, backend and timestep checks before observation or mutation.
This does not prove global stability: earlier processes can alter the state on
which a later process's dynamic timestep guard depends.

```julia
w = World(Grid(24,24))
e = Experiment(w, (writer=ThermalWriter(w), neural=NeuralField(w),
    reaction=Deposition(), transport=Diffusion(:nutrient,.001)); dt=.01)
r = run!(e; steps=100, capture_every=20,
    observe=(w, ps, k)->copy(ps[1].dose))
```

The integration schedule is fixed independently of observation cadence. Captures
include step zero and the final step, with no duplicate final sample. Returned
observations and diagnostics are copied, so references to evolving arrays cannot
rewrite earlier samples. Callbacks must be read-only; Julia cannot prevent a
callback from mutating a world. `diagnostics=nothing` disables default host
reductions when they are expensive. A repeated run continues the world clock;
report step numbers start at zero for that invocation. End time reflects the
existing Float32 clock's repeated additions, rather than an artificial reset to
an exact mathematical duration.

`process_ports(p)` describes reads, writes, published fields, owned state, drives,
memory and numerical law. It is documentation, not a conservation or backend
certificate. Custom processes default to unspecified semantics. Two exclusive
publishers require `allow_overwrite=true`. For example, an Excitable followed by
a NeuralField both replaces activity; this explicit opt-in makes the neural
publication the final activity each step. Additive acoustic activity forcing and
material exchange may share fields without that opt-in. Read the process order
in `r.metadata.process_order` when interpreting a result.

Rendering is a separate observation. An `OpticalMap((w,ps)->...)` explicitly
selects physical fields or process-owned arrays for the renderer. It does not
blend arbitrary latent arrays or advance a simulation. Recording can use the
same boundary: `observe=(w,ps,k)->record_state(w,ps,joinpath(out,"step-$k");step=k)`.
Each destination must be new; repeated runs need separate output roots.

Run `julia --project bin/interior.jl study compositions` for writing → neural → deposition
→ nutrient transport and acoustic → neural → deposition → transport. Each
example includes a coupling ablation. Writing uses deposit feedback to change
persistent neural state; acoustics uses pressure feedback. Those differences are
observable causal checks of these configurations, not evidence for a universal
composition guarantee or biological realism.
