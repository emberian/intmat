# Shared recurrent training on the verified AMD runtime

Use the existing KernelAbstractions stack for the first general learner. The
implemented runtime trains the entire input/recurrent/head matrices of a vanilla
tanh RNN by truncated backpropagation. Action, value and prediction heads share
one learned hidden representation. The first experiment explicitly shares
parameters across residents/worlds while keeping hidden state private per ID.
This is an ecological/training choice, not an implementation trick that silently
replaces independent policies.

`studies/shared_recurrent_runtime.jl` contains the reusable core. Physics,
observations, action transforms, rewards and head objectives remain outside it.
No simulation derivative, handcrafted temporal basis or fixed three-feature
readout is substituted for recurrent learning. Initialization uses an orthogonal
recurrent matrix with declared gain .95; global gradient clipping and Adam/SGD
are explicit training choices. Width and horizon are generic.

## Why this route now

The verified hbox package pin is AMDGPU 2.8.0 commit
`c77460e860b35f21277d8f1253e37304bd8b8993`, GPUCompiler 2.14.1, LLVM 10.2.1,
KA 0.9.45, Julia 1.12.4, standalone ROCm 6.4.3, gfx1030 compatibility override on
RX 6750 XT. It runs native KA kernels, but its startup logs explicitly report
rocBLAS and MIOpen unavailable. This is also recorded in `docs/hbox-runtime.md`.
Read-only inspection found batched and strided-batched GEMM bindings in the pinned
AMDGPU source; bindings do not establish that the vendor library exists or works.
No package or runtime installation was performed for this investigation.

A separate pinned training project remains useful for automatic-differentiation
cross-checks or larger models. The current [Flux GPU guide](https://fluxml.ai/Flux.jl/dev/guide/gpu/)
documents AMD support through AMDGPU/ROCm and MIOpen. [NNlib](https://fluxml.ai/NNlib.jl/dev/)
loads AMD rules as an extension; its [batched multiplication](https://fluxml.ai/Flux.jl/dev/reference/models/nnlib/)
represents independent matrix slices and has derivative rules. Generic dense
matrix multiplication on ROCArray dispatches to [rocBLAS](https://amdgpu.juliagpu.org/stable/libraries/linalg).
Installing Flux or Lux alone does not cure the missing vendor runtime.

Lux's current [GPU management guide](https://lux.csail.mit.edu/stable/manual/gpu_management)
labels native AMD support experimental and recommends Reactant for AMD. Reactant
would introduce another compiler/runtime path that is not the verified Interior
stack. It should be a separate experiment, rather than an implicit dependency
upgrade for the first living learner. [Lux parameters/state separation](https://lux.csail.mit.edu/stable/manual/migrate_from_flux)
is convenient for an oracle model, while [Flux recurrent cells](https://fluxml.ai/Flux.jl/stable/guide/models/recurrence/)
also implement the same simple recurrence.

Read-only compatibility inspection, not a successful dependency solve:

| Source inspected | Relevant declared ranges |
| --- | --- |
| [Lux1.31.4 project](https://raw.githubusercontent.com/LuxDL/Lux.jl/main/Project.toml) | Julia1.10+, LuxLib1.15.1+, NNlib0.9.33+, Zygote0.7.5+ |
| [LuxLib1.15.9 project](https://raw.githubusercontent.com/LuxDL/Lux.jl/main/lib/LuxLib/Project.toml) | AMDGPU2.2+, KA0.9.39+, NNlib0.9.33+ |
| [MLDataDevices1.17.10 project](https://raw.githubusercontent.com/LuxDL/Lux.jl/main/lib/MLDataDevices/Project.toml) | AMDGPU2.2+, GPUArrays11.3.2+ |
| [Flux0.16.12 project](https://raw.githubusercontent.com/FluxML/Flux.jl/master/Project.toml) | AMDGPU1/2, GPUArrays11.2/12.1, NNlib0.9.48+ |
| [NNlib0.9.45 project](https://raw.githubusercontent.com/FluxML/NNlib.jl/master/Project.toml) | AMDGPU1/2, KA0.9.2+ |

These are fetched project versions, not a proposed resolved manifest; upstream
main pages can be cached at different revisions. Their broad bounds admit the
verified AMD/KA versions, but that does not prove all transitive dependencies or
GPU primitives. If needed, create a new private project from a copied verified
manifest, preserve every existing package/compiler pin, and resolve only the
additional library there. Start with a CPU gradient oracle. A future GPU library
route needs a separate matched rocBLAS/Tensile runtime overlay and a small real
GEMM/gradient gate. Do not mutate the working producer environment. Avoid relying
on MIOpen fused recurrence: an explicit tanh matrix recurrence needs only dense
linear algebra and elementwise operations once BLAS is available.

## Runtime contract and numerical checks

Parameters are `Wx[H,D]`, `Wh[H,H]`, `b[H]`, `Wy[O,H]`, `by[O]`.
`step!(hnext,y,theta,x,hprev,continuation)` advances a batch of private states.
`hnext` must not alias `hprev`. Block inputs are `x[D,B,T]`, `h0[H,B]` and
`continuation[B,T]`; `forward!` produces hidden states and `y[O,B,T]`.
Continuation is zero at an ID boundary and nonzero within one identity. It resets
the incoming hidden state and cuts reverse credit across that boundary. Bool,
integer and native Float32 zero/one masks are accepted.

The objective supplies `dY[O,B,T]`. `backward!` computes the exact derivative of
that unrolled tanh network, including input/recurrent/head weights, and exposes
`dh0`. It returns summed gradients with no implicit averaging: loss normalization
belongs to the objective. TBPTT detachment ignores `dh0` while carrying the final
hidden state into the next block. It does not reset memory at every update.
Observations, sampled actions, returns and prediction targets are data in this
VJP; the runtime does not differentiate the physical environment.

The cache preallocates h, d-preactivation, y, initial-state gradient, all parameter
gradients, and deterministic partial reductions. Each partial covers 128 samples
in batch-fastest order; the final reduction visits partials in increasing order.
There are no floating gradient atomics. SGD/Adam use reusable moments, a scaled
sum-of-squares norm to clip large finite gradients safely, and a compact failure
flag for nonfinite gradients. Shape/backend/dtype contracts and representable
optimizer options are checked. Call `check_status!` at an update/rollout boundary.

Callers must keep inactive observations finite/zero and set invalid dY to zero;
a reset mask is not a loss mask. Forward must precede backward with the same
weights/data/masks. Keep weights fixed during the recorded rollout. One fresh
on-policy block is supported first; a replay store does not by itself justify
off-policy recurrent targets or stale hidden-state reconstruction.

`test/shared_recurrent_runtime.jl` passes 127 checks: independent matrix-expression
forward oracle, finite differences for every parameter and h0, reset boundaries,
causal one-step/block agreement, deterministic multi-partial gradients,
Float32/Float64 agreement, SGD/Adam/clipping, large finite gradients, nonfinite
rejection, workspace mismatches and mask compatibility.

The bounded native gate compares Float32 AMD forward, BPTT, h0 gradients and Adam
against a Float64 CPU reference. Reset initial-state gradients are exactly zero.
No rocBLAS/MIOpen or new packages are required:

| H | Batch | Horizon | Maximum scaled error | Component block median | Core device arrays |
| --- | ---: | ---: | ---: | ---: | ---: |
|16|7|8|3.52e−7|2.17 ms|.039 MiB|
|64|512|8|8.39e−7|5.63 ms|4.25 MiB|
|128|512|32|1.34e−6|46.99 ms|31.37 MiB|

The component boundary includes forward+BPTT+Adam+compact status synchronization,
with 400 host-allocated bytes per warmed call, zero GC and no timed compilation.
It excludes simulation, rollout collection, objective, observations and transfers;
these are not arena throughput results. The successful job costs 27.50s wall,
30.56s user+1.02s system CPU, maxRSS 1,096,324KiB and zero swaps. Pin: exclusive
RX 6750 XT, CPU10–11 SMT siblings of one physical P-core, Julia2/OpenBLAS1,
nice15, 8GiB/no-swap control. Frozen source/environment:
`/othersys/tmp/intmat-probe-20261014/{source,environment}-shared-recurrent-20261016-v1`.
Report: `benchmarks/2026-10-16-shared-recurrent-native.json`.
The later host validation change accepting Float32 masks leaves device kernel
bodies unchanged; its two compatibility tests are included above.

## Width, memory and batching envelope

The current first actor has three controls, one value and twelve forecasts:
D25/O16. For a tanh RNN, parameter count is
`P = H*(D+H+1) + O*(H+1) = H²+42H+16`.
With shared weights, parameters/gradient/Adam moments cost16P bytes regardless
of resident count. Private parameters would multiply that term by512; they are
not silently substituted for this shared training experiment.

| H | Parameters | Shared theta | Private theta+gradient+Adam, capacity 512 |
| --- | ---: | ---: | ---: |
|16|944|3.69 KiB|7.38 MiB|
|64|6,800|26.56 KiB|53.13 MiB|
|128|21,776|85.06 KiB|170.13 MiB|

A conservative complete training-buffer estimate includes shared theta/gradient/
Adam, h and d-preactivation histories, y and dY, h0/dh0 plus double-buffered live
state, fixed gradient partials and one trajectory ring. The ring stores inputs,
three raw actions, twelve forecast targets, reward/discount/behavior logp and
Int32 ID/validity. It is 180 bytes per capacity-time slot; another retained block
adds only that ring, not a second optimizer. Physical matter and observers are
separate. Generic AD tapes or future vendor workspaces are not included.

| H | Horizon8 envelope | Horizon32 envelope | Forward FLOPs/decision at B512 | Forward/backward block estimate, T32 |
| --- | ---: | ---: | ---: | ---: |
|16|1.96 MiB|7.41 MiB|.934 million|.085 GFLOP|
|64|4.64 MiB|16.74 MiB|6.881 million|.642 GFLOP|
|128|9.19 MiB|32.78 MiB|22.151 million|2.089 GFLOP|

These are explicit arithmetic/storage estimates, not measured time predictions.
The block estimate allows head recomputation and omits scalar activations,
optimizer work and framework overhead. `studies/recurrent_runtime_costs.jl`
reproduces the six rows in
`benchmarks/2026-10-16-recurrent-runtime-envelope.json`. The earlier full 13-action
proposal has O26 and `P=H²+52H+26`; at H128 its private four-buffer term would be
180.20 MiB, while shared theta alone is 90.10 KiB. The distinction concerns parameter
sharing, not whether action/prediction share a representation.

Shared weights turn each time slice into ordinary GEMM over agents/worlds and
let the parameter gradient sum across all valid samples. Independent per-cell
weights instead require strided-batched/grouped matrix-vector products (or small
per-cell GEMMs over replay sequences); a single shared-policy GEMM changes that
model. Never build a giant dense block-diagonal matrix or launch one tiny library
call per resident. Pack world and resident into the batch axis, while identity
masks stop hidden/return continuity at each world's births/deaths. Do not pool
physical resources or transfer hidden state across worlds.

The next throughput gate must include native physical advances, live inference,
actions, rewards, rollout collection, objective, BPTT, optimizer and hidden carry.
Report simulated time/s, native physical steps/s, valid decision samples/s,
updates/s and peak device memory, with observer/serialization overhead separately.
Only then choose more worlds, larger batches, tiling or a BLAS overlay. The
RX 6750 XT's [vendor peak 432 GB/s and 13.31 TFLOP/s](https://www.amd.com/en/products/graphics/desktops/radeon/6000-series/amd-radeon-rx-6750-xt.html)
are not achieved rates or timing forecasts. Training and optical simulation must
share one coordinated GPU reservation; no concurrent independent render should
contaminate measurements or silently steal memory from a sustained run.

## Wind-down checkpoint: whole-pipeline profile remains pending

The user requested wind-down before the replacement full-pipeline job. No v2
job was launched. The first v1 attempt crossed a new adapter identity-isolation
fix and a trainer freeze request; only its own Julia PID 1337290 was terminated
before evidence collection. Its immutable source, environment and log remain:

- `/othersys/tmp/intmat-probe-20261014/source-recurrent-pipeline-20261016-v1`
- `/othersys/tmp/intmat-probe-20261014/environment-recurrent-pipeline-20261016-v1`
- `/othersys/tmp/intmat-probe-20261014/recurrent-pipeline-v1.log`

Do not resume that old adapter pin as the corrected learner. There is no measured
complete-pipeline bottleneck or learning-benefit result yet. The finished native
runtime evidence above is valid and separate. CPU/GPU reservation is released.

Current ownership: performance lane owns `shared_recurrent_runtime.jl`, its
focused tests, native component gate, arithmetic envelope and pending
`recurrent_pipeline_profile.jl`; design owns `native_recurrent_agent.jl` and the
native controller override seam; frontier owns `recurrent_training_objective.jl`,
`recurrent_training.jl` and `train_native_recurrent.jl`. Their current adapter
identity regression and objective/driver validation results should be read before
freezing a new source stage. Root owns snapshot/commit/push and next reservation.

Next bounded check: freeze current validated source into a fresh v2 stage and
copy the verified environment with only its Interior source path changed; do not
install packages. Reuse the v1 controlled launcher template with v2 source/env,
output and log names. Its executable payload is:

```
SWARM_MEM_MAX=8G /home/hbox/bin/swarm-build /usr/bin/time -v \
  timeout --kill-after=15s 300s taskset -c 10,11 \
  "$PROBE_ROOT/julia-1.12.4/bin/julia" --startup-file=no \
  --project="$FRESH_PINNED_ENV" \
  "$FRESH_SOURCE/studies/recurrent_pipeline_profile.jl" "$FRESH_REPORT"
```

Use the existing launcher exports (Julia2/OpenBLAS1, isolated depot with trailing
colon, standalone ROCm 6.4.3, gfx1030 override and opt-level1). The harness already
selects 64³ parallel support, H64/T8 and H128/T32, two warm blocks and three measured
blocks, phase timing and a CPU/native objective VJP check on real rollout data;
final observation is timed separately. It has not yet completed this composed
execution path. Coordinate GPU ownership first. Only after that evidence choose
replica batching, reusing existing constructive/sidecar batch kernels rather
than duplicating physical laws. No sustained or unbounded run is authorized by
this restart note.

Before stopping, design reported the final adapter identity regression26/26
passing. Corrected adapter SHA256:
`e53c3d056c03747b1207531c9e4356db6cca022c1ade5b6da74ab7aa0fc8deeb`;
its core model source:
`be9511b68465fd3e4b1c5aef1d012d73acc19bbd9a16600c5b397564d490d08e`.
This is a CPU integration result; no new native run followed it.
