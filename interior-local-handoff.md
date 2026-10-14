# Interior — local continuation

**14 October 2026 · handoff h1 on the unchanged 0.0.2 reference**

## The project to carry forward

Interior is an artistic research instrument for **interacting phenomena in matter**. A specimen combines initial material, local laws, a history of interventions, and an observation setup. The interesting subject is how processes alter the material and thereby change what can happen next.

ember explicitly wants ambition and composability, not a narrow effect generator. The initial reference image suggested transparent or translucent matter with complicated internal structures, but amber, a glass sphere, lightning, a neural automaton, and even transparency are not requirements. Neural activity may inhabit material, remodel it, constitute it, or live on interfaces made by another process. It is one family alongside reaction, transport, mechanics, phase change, fields, growth, waves, and other processes.

The user invited ChatGPT to take artistic and technical initiative, with ember assisting on implementation and compute. Take that invitation seriously: make concrete experiments, form and revise artistic judgments, and preserve their evidence. Do not turn each reversible parameter choice into a question for ember. Do not pretend a new local session has memories beyond the supplied record.

The latest design question was whether this exceeds Cycles, or whether each evolving state can be baked into a renderable instant. The answer adopted here is **renderer-independent material history, followed by explicit optical interpretation**. Cycles can be one observer, without dictating the simulator's laws or the full optical contract.

## Start here on the machine

Unpack the archive and open a terminal in this directory. Use Python 3.10 or newer, with dependency versions compatible with that Python. The supplied exact test environment is recorded in `HANDOFF-STATUS.json`; the broad package dependency ranges are not a reproducibility guarantee.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest -q

# Existing three-run neural/material study and diagnostic slice viewer:
python -m interior.living_demo --n 20 --until 8 --out local-living

# New full-volume, view-independent sample recorder (one selected variant):
python -m tools.bake_sequence \
  --n 20 --until 8 --max-dt 0.008 --every-steps 25 \
  --out local-samples
```

Windows PowerShell activation instead: `.venv\Scripts\Activate.ps1`. Run the `tools` module from this checkout root: it is a handoff utility, not part of the installed `interior` package discovery.

Open `local-living/viewer.html` for the original diagnostic interface. `local-samples/manifest.json` describes the new full 3D samples. A recorder output directory must not already exist; choose a new run name rather than overwriting another experiment.

The full samples are much larger than slice previews. Start small, check bytes per sample and runtime, and only then increase grid resolution, duration, or capture frequency. Do not launch a large parameter sweep as the first action.

## What actually exists

The original implementation is a nondimensional, fixed-grid NumPy reference. It has snapshot-based additive explicit coupling with atomic state commits and rejected-step retries, conservative masked-grid transport, an abstract exothermic reaction, an imposed energy pulse, an irreversible alteration field, and a six-channel recurrent neural field coupled to the reaction and alteration.

The default neural weights contain deliberately supplied oscillator structure plus seeded perturbations. There is no learned or evolved controller, reproduction, selection, mechanical solver, moving carrier, acoustic wave solver, or completed optical renderer. The alteration field is not a simulated fracture. Unit tags are string checks, not dimensional inference. No energetic cost for latent computation is established.

The original 24 tests were rerun successfully. This handoff adds 11 recorder checks, for **35 passing tests** in the recorded environment. It includes an executed eight-sample `8^3` recording at model time 0.25. This is a smoke test, not artwork, physical validation, or a performance benchmark.

All 51 files from the supplied 0.0.2 archive are retained unchanged. Their hashes are in `handoff-evidence/upstream-sha256.json`. The Python package version remains 0.0.2: h1 is an additive handoff, not a new universal solver release.

### New recorder

`tools/bake_sequence.py` observes successfully committed steps of the original engine. It records the initial state, every Nth accepted step, and the final state. There is only one `Engine.run` over the interval. The recorder does not choose a camera frame rate or insert extra integration stops.

Each sample contains full precision 3D material and latent fields, the mask, spacing, actual time and step count, and the existing artistic optical-preview fields. The manifest records coordinates, nondimensional unit tags, field specifications, rule provenance and weights, source hashes, actual sample times, per-file hashes, and numerical audits. No camera or viewing light is baked into the specimen.

It explicitly says that the optical preview's varying refractive-index field requires an appropriate optical interpretation; a backend must not silently treat its presence as support for curved rays. It also states that spectra, a phase function, extracted interfaces, and a render backend are absent.

These are **state samples, not a tested restart format**. There is no loader that reconstructs and resumes a complete world, no temporal interpolation, and no VDB or renderer exporter. Atomic rename gives file visibility, not power-loss durability of a transaction. On I/O failure, the simulation can already have committed a step; the recorder aborts and attempts to mark the partial recording failed.

### A real trap in the old viewer

`interior.living_demo.run_variant` calls `engine.run` separately for display timestamps. This can insert shorter steps and change the numerical trajectory when the requested frame count changes. A small measured example is in `handoff-evidence/old-viewer-cadence-probe.json`.

The new recorder avoids that coupling. Its tests compare capture strides 1, 3, and 100 against one uninterrupted original `Engine.run`: the tested final arrays and full numerical audit are exactly equal within this environment. It does **not** promise exact equality to the old viewer run with its extra frame stops, nor bitwise equality across platforms.

## Read in this order

`AGENTS.md` gives the continuation mandate and acceptance criteria. `OPTICS-BOUNDARY.md` gives the optical architecture, limits of snapshot baking, backend escalation, and freshly checked primary references. `README.md` describes the historical runnable studies. `DESIGN.md` and `LIVING-MATERIALS.md` preserve the wider project direction. `HANDOFF-STATUS.json` distinguishes executed evidence from proposed work.

The `references/` directory carries the five user-supplied images, with provenance and interpretation limits. These are inspiration, not simulator outputs.

The inherited research list is historical design context, not a claim that every external repository or paper was reverified during this handoff. Inspect licenses and current compatibility before copying code or weights. No external research implementation is bundled here.

## The next working session

Build two narrow, demonstrable slices of the broad architecture rather than a giant framework.

**Observation slice:** take a saved specimen, define one optical interpretation explicitly, and render it with a real backend. First validate a simple absorbing slab, then render the same specimen with two cameras or lighting arrangements without rerunning material evolution. Record every optical approximation. Make one compelling still early; do not mistake a shader-tuning exercise for simulated structure.

**Composition slice:** add a genuinely different state carrier, initially particles or a deforming interface with persistent identity and internal state. Exchange one conserved species with a grid in both directions. Establish conservation and reorder invariance before pursuing evolution, complex controllers, or a large ecology. This remains the next important simulation milestone, not displaced by the optical work.

The first slice should not make Cycles a dependency of the simulation. The second should not require a completed custom renderer. Together they test the claimed architecture much more strongly than adding more latent channels to the existing voxel demo.

## Initial artistic direction

A promising first study is **a rhythm that writes a material record**: activity gates deposition, the record changes subsequent transport or communication, and lighting reveals the resulting interfaces or inclusions. It can use deliberately invented laws, but the dependencies must actually exist in the code.

Keep a second study available that stresses the optical boundary: **a changing material that redirects its illumination**. Start with a prescribed index field or ordinary interfaces, then introduce feedback only after the optical response is understood. It need not begin as a complete acousto-optical system.

Do not promise that any arbitrary controller produces life, persistent novelty, cooperation, or beauty. Retain simple passive specimens and fixed-rule active specimens alongside adaptive or evolutionary ones. Let unexpected regimes alter the artistic plan.

## Prompt to give a local coding/research agent

> Continue Interior from this directory. Read START-HERE.md, AGENTS.md, OPTICS-BOUNDARY.md, README.md, DESIGN.md, LIVING-MATERIALS.md, and HANDOFF-STATUS.json, then inspect the implementation. Treat the broad architecture as proposed and the tests as limited empirical evidence. Take artistic and engineering initiative with ember; preserve composability across phenomena, material representations, and optical models. First reproduce the tests and a small sample recording, inspect the actual machine and installed tools, and write a short concrete execution plan. Then implement the smallest real end-to-end observation slice and begin the grid/second-carrier exchange slice. Do not merely produce another grand design document. Do not silently drop unsupported optics, bind simulation timesteps to display frames, hard-code a glass sphere, claim evolution from fixed dynamics, or port everything to a new language before profiling. Keep a running record of decisions, actual commands, results, approximations, and the next blocked edge. Ask ember only for consequential artistic choices, missing capabilities, or compute beyond agreed limits.
