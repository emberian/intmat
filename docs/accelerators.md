# Accelerator execution routes

The continuous arena currently runs its ecological sensors and recurrent controllers on CPU arrays. Its implicit phase solver explicitly requires CPU and uses FFTW with ComplexF64 workspaces. Moving an unchanged arena to a machine with a GPU therefore does not accelerate it. Dense field kernels already use KernelAbstractions and can be tested on an available Julia GPU backend.

On 14 October 2026, SSH to `hbox` worked. Vulkan enumerated **AMD Radeon RX 6750 XT**, with 12,868,124,672 bytes of VRAM reported by sysfs, and GPU activity was 0% at inspection. The host has an i9-12900, 24 logical CPUs and 123 GiB RAM. It was busy (load about 26, swap nearly full), so its CPU is not spare capacity to assume. Julia was absent from PATH and the usual install locations; no ROCm/HIP libraries were registered with ldconfig. The historical H1 Python environment still has a library directory but its Python executable is gone. Historical successful ROCm runs do not establish a working runtime today. Vulkan enumeration establishes hardware presence, not Julia compute availability.

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
  -p /tmp/intmat-kaggle-results --force
```

The prepared private notebook contains only Project.toml, Julia source and the bounded benchmark. It downloads Julia 1.12.4, creates an isolated CUDA environment, limits Julia and precompilation to two CPU threads, records hardware/runtime information, and compares CPU/CUDA on 32³ and 64³ with ten timed steps each. Installation time is separate from warmed step time. The total job is capped at 25 minutes, including installation. Never push a new version over a running notebook: the historical API workflow lost access to running output that way. Use a distinct notebook ID for a distinct simultaneous job.

A TPU is a different path: this Julia KernelAbstractions code has no TPU backend. A JAX port of a pure batched controller can test TPU throughput, but ecological sensor gathering, births, deaths, positions, state transfers and FFT phase evolution still need implementation and timing. Controller-only matrix throughput is not an end-to-end arena speedup.
