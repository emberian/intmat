# hbox isolated AMD probe (2026-10-14)

The Radeon RX 6750 XT is reachable through its existing kernel driver. An isolated
ROCm 6.4.3 runtime enumerated `gfx1031` with `KERNEL_DISPATCH`; HIP allocation,
device memset, synchronization, and a 4096-byte device-to-host copy all returned
status zero, and every copied byte matched. HIPRTC then compiled a native
`gfx1031` kernel assigning `x[i] = 2*i + 1` to 256 Float32 elements. Module load,
kernel launch, synchronization, and readback all succeeded; all 256 results
matched exactly. No architecture override was used. This establishes real GPU
compute. The tested Julia AMDGPU environments did not establish Interior parity
or throughput: the released stack failed reductions, and a newer upstream stack
failed a minimal KernelAbstractions arithmetic kernel before Interior. The
Interior AMD path is not verified usable on hbox.

## Runtime and isolation

Host: Ubuntu 24.10, Intel i9-12900 (24 logical CPUs), about 123 GiB RAM, RX 6750 XT
12 GB. The hbox account already belongs to `render` and `video`, and can access
`/dev/kfd` and the DRM render devices. No Julia was found in the original host
PATH, and the system library cache did not contain HIP/HSA/rocFFT.

The new directory `/othersys/tmp/intmat-probe-20261014` contains Julia 1.11.7 and
1.12.4, a
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

The default reproduces the released Julia 1.11.7 diagnostic baseline. For the
retained upstream environment, set `PROBE_JULIA_VERSION=1.12.4` and
`PROBE_PROJECT=$PROBE_ROOT/environment112`. Its exact AMDGPU source commit is
`c77460e860b35f21277d8f1253e37304bd8b8993`, chosen because it retains KA 0.9 while
using LLVM 10. `install` explicitly pins the released baseline and should not be
used to replace the upstream environment. `prefetch` installs the declared
LLVMDowngrader lazy artifact without importing the affected JLL. During the probe,
that 2.47 MB artifact was downloaded and SHA256-verified before extraction; no JLL
source was changed.

Replay now uses `JULIA_DEPOT_PATH=$PROBE_ROOT/depot:`: the trailing colon includes
bundled system depots while excluding the normal user depot. The earlier probes
used only the fresh depot and consequently rebuilt/JIT-loaded standard libraries.
Normal compiled modules are the replay default; `PROBE_COMPILED_MODULES=existing`
reproduces the earlier uncached-loading policy. One controlled normal Pkg cache
precompile completed in 128.875 seconds. Writes stay in the isolated first depot.

Each Julia command uses hbox's existing `/home/hbox/bin/swarm-build`, 8 GiB RAM
limit, no swap, nice 15, two pinned CPUs (0,1), two Julia threads, one precompile
task, and a 900-second timeout. The existing shared slice remains the aggregate
memory boundary. The host is busy: sampled cores 0–15 were 98–100% occupied.
Both the probe scope and shared slice reported unlimited CPU quota with zero
throttled periods, so slow startup compilation reflects contention, not a new
quota failure. Existing jobs were not modified.
Setup-plus-compile pipelines used a 1200-second cap; replay can select that with
`PROBE_TIMEOUT_SECONDS=1200`. No Julia probe jobs remain, and the retained directory
occupies approximately 2.0 GiB.

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
`pointerref` invocations through GPU hostcall allocation/boxing. Separate scalar
`minimum`/`maximum` reductions failed too. Their temporary source change passed
26 CPU core checks, finite/nonnegative state guards, and Metal parity, but did not
repair AMD compilation, so the original `Excitable` source was restored.

Julia 1.12.4 with the same released AMDGPU/GPUCompiler 2.6/LLVM 9 stack also failed
`minimum(vec(ROCArray))`. The latest upstream commit required KA 0.10 and could not
resolve with Interior's KA 0.9 constraint. The compatible commit above resolved
with GPUCompiler 2.14.1, LLVM 10.2.1, and backend JLL 23.1.3. After fixing its lazy
artifact bootstrap, its runtime initialized, but a tiny KA arithmetic kernel
failed with `InvalidIRError` in `_pointerset`, boxing and hostcall allocation,
before reductions or Interior. This broader failure means an extrema-only patch
or a custom KA reduction is not demonstrated to solve the problem.

The final normal-cache control requesting released AMDGPU 2.8.0 plus GPUCompiler
2.10.0 and LLVM 9.13.2 was rejected by the resolver: GPUCompiler requires backend
JLL 23.1.1+, while released AMDGPU requires 22.1.8. GPUCompiler 2.10 never executed
a kernel in that control. No further version permutations were attempted.

No Interior throughput or parity number was produced. Remote logs include
`probe.log`, `probe-after.log`, `probe112.log`, `probe-upstream.log`,
`probe-compatible.log`, `probe-final.log`, and `probe-stable.log`. Structured
results, exact package versions, source/tree hashes, and both retained manifest
SHA256 hashes are in `docs/hbox-accelerator-results.json`. No physics source change
remains from this investigation.
