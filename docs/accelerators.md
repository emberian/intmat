# Accelerator execution routes

The continuous arena currently runs its ecological sensors and recurrent controllers on CPU arrays. Its implicit phase solver explicitly requires CPU and uses FFTW with ComplexF64 workspaces. Moving an unchanged arena to a machine with a GPU therefore does not accelerate it. Dense field kernels already use KernelAbstractions and can be tested on an available Julia GPU backend.

On 14 October 2026, SSH to `hbox` worked. Vulkan enumerated **AMD Radeon RX 6750 XT**, with 12,868,124,672 bytes of VRAM reported by sysfs, and GPU activity was 0% at inspection. The host has an i9-12900, 24 logical CPUs and 123 GiB RAM. It was busy (load about 26, swap nearly full), so its CPU is not spare capacity to assume. Julia was absent from PATH and the usual install locations; no ROCm/HIP libraries were registered with ldconfig. The historical H1 Python environment still has a library directory but its Python executable is gone. Historical successful ROCm runs do not establish a working runtime today. Vulkan enumeration establishes hardware presence, not Julia compute availability.

The subsequent isolated bootstrap restored Julia/ROCm compute without system changes. A native gfx1031 HIP arithmetic kernel passed exact readback; AMDGPU.jl 2.8.0 then reported `functional() == true` with HIP 6.4. AMDGPU automatically applies its gfx1030 compatibility setting for this Navi 22 device. See [the hbox runtime record](hbox-runtime.md) and `scripts/bootstrap_hbox_probe.sh` for the exact fresh directory, resource controls and Interior parity outcome. The minimal runtime omits optional BLAS/FFT libraries, so this does not provide a GPU FFT solver.

The first Interior probe failed before parity/timing: AMDGPU could not compile the tuple-valued `extrema` reduction used by `Excitable` validation (`InvalidIRError` in `partial_mapreduce_device`). The runtime and native arithmetic kernel passed; the dense Interior Julia path needs that reduction repaired. `docs/hbox-accelerator-results.json` records the actual failure and leaves simulation timings null.

Keep new hbox artifacts under `/tank` or a fresh project-specific directory in `/othersys`; respect the existing `swarm-build` resource controls for builds. Do not replace its H1 environments, remove caches, or alter system packages to establish a benchmark. At inspection, `/tank` had 423 GiB free and `/othersys` had 97 GiB free.

## Bounded Julia probe

`scripts/benchmark_accelerator.jl BACKEND N STEPS` accepts `cpu`, `metal`, `cuda`, or `amd`, with bounds of 4–128 cells per side and 1–100 timed steps. The selected optional package must be installed in the active environment. It checks an 8³ CPU/device comparison, warms the exact requested shape, then reports synchronized dense-field time as JSON. It does not advance the live arena, write its checkpoint, or benchmark the ecological workload.

```sh
julia --project scripts/benchmark_accelerator.jl cpu 32 10
julia --project=gpu scripts/benchmark_accelerator.jl metal 32 10
```

CUDA environments need CUDA.jl; hbox would need a working AMDGPU.jl/ROCm environment. [AMDGPU's runtime requirements](https://amdgpu.juliagpu.org/stable/faq) require ROCm and a supported GPU. Availability must be established with a kernel launch and parity check.

## Kaggle access and launch

The existing account is `emberian64`. The private bearer-token file `~/.kaggle/access_token` is present, and an authenticated CLI status query for the previous T4 CPT notebook returned COMPLETE on 14 October. This confirms API access; it does not measure today's remaining quota or imply that old TPU training permissions cover a new long simulation.

The previous workflow used private Python script kernels, uploaded with `uvx --from kaggle kaggle kernels push`, and retrieved results via the CLI/API. Kaggle notebooks are jobs rather than an existing SSH host. Historical GPU and TPU access used T4 and TPU v5e-8; current accelerator names must be selected explicitly. [Kaggle's CLI documentation](https://github.com/Kaggle/kaggle-cli/blob/main/docs/kernels.md) documents the `--accelerator` and `--timeout` flags.

```sh
python3 scripts/prepare_kaggle_probe.py /tmp/intmat-kaggle-probe
uvx --from kaggle kaggle kernels push -p /tmp/intmat-kaggle-probe \
  --accelerator NvidiaTeslaT4 --timeout 1500
uvx --from kaggle kaggle kernels status emberian64/intmat-julia-cuda-bounded-probe
uvx --from kaggle kaggle kernels output emberian64/intmat-julia-cuda-bounded-probe \
  -p /tmp/intmat-kaggle-results --force --file-pattern '^measurements.jsonl$'
```

The prepared private notebook contains only Project.toml, Julia source and the bounded benchmark. It downloads Julia 1.12.4, creates an isolated CUDA environment, limits Julia and precompilation to two CPU threads, records hardware/runtime information, and compares CPU/CUDA on 32³ and 64³ with ten timed steps each. Installation time is separate from warmed step time. The total job is capped at 25 minutes, including installation. Never push a new version over a running notebook: the historical API workflow lost access to running output that way. Use a distinct notebook ID for a distinct simultaneous job.

A TPU is a different path: this Julia KernelAbstractions code has no TPU backend. A JAX port of a pure batched controller can test TPU throughput, but ecological sensor gathering, births, deaths, positions, state transfers and FFT phase evolution still need implementation and timing. Controller-only matrix throughput is not an end-to-end arena speedup.

The private CUDA probe completed successfully on 14 October. Kaggle supplied two Tesla T4 devices with 15 GiB each; the Julia probe used the default single device. Julia 1.12.4 with CUDA runtime 12.8 and CUDACore 6.4.2 completed installation/precompilation in 425.7 seconds, then passed parity with maximum absolute error 3.73e-9. Ten warmed dense-field steps took:

| Grid | Kaggle CPU | Single T4 | CPU/T4 ratio |
| --- | ---: | ---: | ---: |
| 32³ | 54.48 ms | 8.21 ms | 6.64× |
| 64³ | 355.78 ms | 14.44 ms | 24.63× |

Machine-readable results are in `gpu/cuda-verification.json`. These establish a working remote Julia GPU route and dense-kernel execution; they do not measure the continuous arena. The first notebook version accidentally included its Julia installation in output artifacts; download only the named measurement file. The preparation script now stores runtime/source/environment under `/tmp` so future notebook outputs stay small.

The same bounded 64³ ten-step command on the local M2 Max Metal backend took 312.88 ms, compared with 14.44 ms on the T4 (21.66× lower elapsed time). Both passed parity. This is a single comparison on different shared hosts, with host load uncontrolled; it supports using Kaggle for this dense workload, not a universal GPU ranking.

## Measured local controller component

`julia --project=gpu scripts/benchmark_controller.jl metal` benchmarks the actual fused controller with 21 output rows, 44 input columns and independent Float32 weights/marks per organism. On the local M2 Max, twenty warmed repeats on 14 October produced:

| Population | CPU controller | Metal, state resident | Upload all weights/marks/inputs and download output |
| --- | ---: | ---: | ---: |
| 2,400 | 3.320 ms | 0.890 ms | 11.765 ms |
| 10,000 | 9.512 ms | 2.167 ms | 23.094 ms |

Maximum absolute CPU/Metal output difference was 1.19e-7. These are a single bounded component comparison, excluding sensors, phase FFT and ecological changes. Resident execution improved controller time by 3.7–4.4×; staging all arrays every tick was slower than CPU. Keeping the controller state on the device and updating only the required inputs is necessary before expecting an integrated benefit. This result is not evidence that the entire arena runs on Metal.
