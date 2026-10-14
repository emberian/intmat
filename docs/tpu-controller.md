# Bounded JAX TPU controller component

`scripts/tpu_controller_probe.py` measures one component of the neural ecology:
Float32 controllers with 12 hidden values, 32 sensors, 21 outputs and 44 inputs,
a periodic six-channel local grid gather, exponential hidden-state relaxation,
and bounded motion with supplied funding, permeability and flow. The remaining
26 sensor values are fixture inputs. Resource competition, neighbor search,
plasticity, births, deaths, field evolution and the arena are outside this probe.
The 100-step scan evolves hidden state and positions with fixed weights, fields,
funding and remaining sensor cues. Timings cannot establish whole-arena speedup.

The CPU reference sums columns sequentially, following `neural_controller!`.
A local 64-agent fixture was also evaluated by Julia's actual helper; the
NumPy reference agreed within Float32 tolerance. JAX checks output, hidden,
position and sensor errors against that reference before recording timings.
`Precision.HIGHEST` requests Float32-quality products on TPU; the probe fails
above a maximum absolute error of `3e-5`. This is tolerance parity, not bitwise
identity. The current probe uses one default device; it reports all visible
devices but does not shard across the eight-device TPU slice.

Run locally with an environment containing NumPy and JAX:

```sh
python scripts/tpu_controller_probe.py --sizes 64,1024,8192 --output results.json
```

Create a private code-only Kaggle bundle and request a bounded TPU run:

```sh
python3 scripts/prepare_kaggle_tpu_probe.py /tmp/intmat-kaggle-tpu-probe
uvx --from kaggle kaggle kernels push -p /tmp/intmat-kaggle-tpu-probe --accelerator TpuV5E8 --timeout 600
uvx --from kaggle kaggle kernels status emberian64/intmat-jax-tpu-controller-component-probe
uvx --from kaggle kaggle kernels output emberian64/intmat-jax-tpu-controller-component-probe -p /tmp/intmat-tpu-output
```

Authentication stays with the local Kaggle CLI. The bundle has no datasets,
credentials, transcripts, package upgrades or network downloads. TPU execution
is required; CPU fallback raises an error. The six-hundred-second platform cap
bounds compilation and measurement together. Report separates initial transfer,
first compilation/execution, warm component dispatch, and warm fused scan time.

On 2026-10-14, version 1 was accepted with the explicit `TpuV5E8` request and
600-second execution cap. Latest observed status at 21:46 UTC was
`KernelWorkerStatus.QUEUED`. The job remains submitted; there are no real TPU
timings yet. Acceptance does not establish TPU availability or performance;
retrieve the actual results first. Local JAX CPU checks at 64 and 1024 agents
had maximum output errors of `1.19e-7` and `1.49e-7`; hidden errors were below
`3e-8`, and gathered sensors and positions agreed exactly. Julia versus the
same NumPy sequential fixture at 64 agents differed by at most `5.96e-8`.
A separate 100-step recurrent CPU check differed by at most `1.19e-7` in hidden
state and `9.54e-7` in position, with exact gathered sensors.
