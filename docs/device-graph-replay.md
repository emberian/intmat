# Experimental AMD graph replay

The cellular macrostep launches roughly 110 fixed-capacity kernels at 16 mechanical substeps. HIP graph replay could reduce host dispatch and driver gaps while retaining the existing kernels and their order. This is a research prototype, not a measured acceleration result.

The pinned hbox Manifest selects `/othersys/tmp/intmat-probe-20261014/AMDGPU-c77460e860b35f21277d8f1253e37304bd8b8993` (version2.8.0). Its `src/hip/graph.jl`, SHA256 `082d96e3a1fe2fe39de63e3f552062504c972d9b0a052310d8b12880bb75dd89`, provides `AMDGPU.HIP.capture`, `instantiate`, and `launch`. The [official graph API](https://amdgpu.juliagpu.org/stable/api/graphs) describes the same sequence. Capture uses the current stream and temporarily disables garbage collection. We explicitly avoid `@captured`: that macro can execute a warm-up body and retries after capture invalidation, which is inappropriate for a stateful biological step unless its effects are separately restored.

## Clock and ownership

Capturing `device_cellular_step!` directly is incorrect. It increments a host `attempted_steps` field once, and embeds that number into lifecycle selection, mutation, daughter placement and inherited catalyst regulation. Replaying it would repeat a random stream.

The prototype `DeviceGraphStepIndex` borrows device `metadata`. Inside existing lifecycle/random kernels its multiplication reads `UInt32(metadata[6])+1`, before the existing last `dc_complete_step!` node increments the completed-step count. On successful runs this equals the ordinary attempted-step index. A latched failure suppresses all later kernels, including completion; its partial state remains unresumable. Graph capture must not change host counters or `World.time`; only actual launches count as attempts. A replay batch must reject completed-counter overflow, preexisting failure, replaced buffers, changed parameters/dt, or mismatched ownership. Restore allocates new buffers and therefore requires a new graph.

The graph owner must retain the process, hook, captured graph, executable graph, and all borrowed arrays for their full lifetime. Capture occurs only after warming the exact kernel specializations outside capture on a disposable/restored fixture. No compilation, synchronization, status download, array allocation, or observation belongs inside capture.

## Staged checks

`studies/device_graph_clock.jl` is an isolated device-index adapter. `studies/device_graph_enqueue.jl` is a private copy of the existing host launch schedule, calling only shared simulation kernels. After exact private parity, production `device_cellular_enqueue!` now contains that unchanged schedule and accepts a scalar or device index. The ordinary API retains its current host counters, checks, launch order and arithmetic. The graph clock and HIP owner remain experimental files outside production includes.

`test/device_graph_clock.jl` checks random values against ordinary UInt32 indices, then compares every full snapshot over a birth/mutation trajectory and strict restore. `scripts/probe_amd_graph_capture.jl` is an unlaunched compatibility probe for a granted exclusive GPU window: it checks capture has no execution side effects, repeated graph launches advance the random stream, and a failure latch suppresses replay. Passing that two-node gate would establish API compatibility only.

The required performance experiment remains a complete-loop alternating comparison from the same checkpoint: ordinary queued steps versus graph replay, identical timestep/substeps/capacity/physics, compilation excluded, synchronized wall time and host enqueue time recorded separately, followed by primary-array, identity, RNG and checkpoint continuation gates. A small-kernel graph result alone is insufficient evidence of faster cellular ecology.

CPU gates passed 170 exact clock/RNG checks plus 42 full-snapshot and strict-restore assertions comparing the frozen pre-extraction schedule, ordinary API, and new device-clock enqueue over 20 three-dimensional birth/mutation steps. The callable 3D constructive habitat gate also passes after extraction. `scripts/benchmark_device_graph.jl` stages alternating whole-loop comparisons, host enqueue time and HIP event stream intervals; no GPU graph result has been measured yet. The pinned graph file was also byte-identical to the older depot copy; the Manifest-selected source above is the authority.

The final CPU gate passed 227 assertions, including coherent cell-buffer/ownership replacement, habitat binding replacement, failed/incomplete state preflight, replay attempt bookkeeping, counter overflow, and empty-population audit. The graph owner retains the ownership objects captured at instantiation and rejects a coherently replaced process before launching the stale graph.
