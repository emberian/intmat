# hbox isolated AMD probe (2026-10-14)

The Radeon RX 6750 XT is reachable through its existing kernel driver. An isolated
ROCm 6.4.3 runtime enumerated `gfx1031` with `KERNEL_DISPATCH`; HIP allocation,
device memset, synchronization, and a 4096-byte device-to-host copy all returned
status zero, and every copied byte matched. HIPRTC then compiled a native
`gfx1031` kernel assigning `x[i] = 2*i + 1` to 256 Float32 elements. Module load,
kernel launch, synchronization, and readback all succeeded; all 256 results
matched exactly. No architecture override was used. This establishes real GPU
compute, but does not yet establish Interior kernel parity or a speedup.

## Runtime and isolation

Host: Ubuntu 24.10, Intel i9-12900 (24 logical CPUs), about 123 GiB RAM, RX 6750 XT
12 GB. The hbox account already belongs to `render` and `video`, and can access
`/dev/kfd` and the DRM render devices. No Julia was found in the original host
PATH, and the system library cache did not contain HIP/HSA/rocFFT.

The new directory `/othersys/tmp/intmat-probe-20261014` contains Julia 1.11.7, a
copied Interior source tree, an isolated Julia depot/environment, and extracted
official ROCm packages. It is outside the root filesystem. Nothing was installed
with apt, no kernel/driver settings were changed, and previous environments and
caches were retained.

Julia archive:
`https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.7-linux-x86_64.tar.gz`.
ROCm package index:
`https://repo.radeon.com/rocm/apt/6.4.3/dists/noble/main/binary-amd64/Packages.gz`.
Packages were downloaded from that index and extracted with `dpkg-deb -x`:
`comgr`, `hip-runtime-amd`, `hsa-rocr`, `rocm-core`, `rocm-device-libs`, `rocminfo`,
and `rocprofiler-register`. The extracted root is
`rocm-debs/opt/rocm-6.4.3`. These Ubuntu 24.04 binaries load on this 24.10 host;
that observation is not a claim that the OS/GPU pairing is officially supported.

Current [AMDGPU installation documentation](https://amdgpu.juliagpu.org/stable/)
requires Julia 1.10+ and ROCm 6+. Old v0.5 instructions about bundled runtime
artifacts should not be applied to current AMDGPU. Installed AMDGPU 2.8.0 supplies
compiler and device-library artifacts, while HIP/HSA come from the external ROCm
runtime. The isolated extraction also provides ROCm device libraries.
[Installation tips](https://amdgpu.juliagpu.org/stable/install_tips) describe the
architecture override. The previous hbox training environment used
`HSA_OVERRIDE_GFX_VERSION=10.3.0`; the HIP memory probe above succeeded without
that override. AMDGPU 2.8.0 itself detects `navy_flounder` in Linux KFD topology
and sets `HSA_OVERRIDE_GFX_VERSION=10.3.0` during initialization, so its Julia
device report uses `gfx1030`. The separate HIPRTC test used native `gfx1031`.

## Resource controls and replay

`scripts/bootstrap_hbox_probe.sh` runs install, diagnose, and the bounded
`amd 32 10` parity benchmark. Its `verify` mode prints AMDGPU runtime details and
runs parity in one Julia process, avoiding two uncached module loads. Set and
export `PROBE_ROOT` to the isolated
directory before invoking it. It requires the already provisioned Julia/ROCm
files and copied source tree.

Each Julia command uses hbox's existing `/home/hbox/bin/swarm-build`, 8 GiB RAM
limit, no swap, nice 15, two pinned CPUs (0,1), two Julia threads, one precompile
task, and a 900-second timeout. The existing shared slice remains the aggregate
memory boundary. The host is busy: sampled cores 0–15 were 98–100% occupied.
Both the probe scope and shared slice reported unlimited CPU quota with zero
throttled periods, so slow startup compilation reflects contention, not a new
quota failure. Existing jobs were not modified.

## Interior result

Julia 1.11.7 with AMDGPU 2.8.0 installed successfully. `AMDGPU.versioninfo()`
found HIP 6.4.43484 and the Radeon device with 11.984 GiB; `AMDGPU.functional()`
returned true. Optional rocBLAS, rocSOLVER, rocSPARSE, rocRAND, rocFFT, hipTENSOR,
and MIOpen libraries are absent from this minimal runtime. In particular, this
setup does not enable an implicit GPU FFT solver.

The bounded `amd 32 10` Interior probe failed before parity or timing in the
validation call `extrema(w.fields.activity)` for `Excitable`. AMDGPU's tuple
extrema reduction generated an `InvalidIRError` in
`AMDGPU.partial_mapreduce_device`, with unsupported dynamic `pointerset` and
`pointerref` invocations through GPU hostcall allocation/boxing. This is a Julia
GPU reduction/compiler failure, not evidence that the GPU or HIP runtime is
unavailable. No Interior throughput number was produced. Full remote logs are
`install.log`, `diagnose.log`, and `probe.log` in the isolated directory; structured
results are in `docs/hbox-accelerator-results.json`.
