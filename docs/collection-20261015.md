# First collection · 2026-10-15

The first series is complete: thirty-six reviewed photographs across six material
families, published alongside the sixteen earlier studies. Folded inclusions,
branching conductors, fractured filaments, pressure surfaces and living niches
provide six photographs per family. This historical replay guide records their
original simulations and 1600² photography; later native-grid refinements are
separate specimens and preserve this series. The final selection is the union of
`outputs/collection-20261015/collection-*-selected.json`, excluding pilots and
rejected phase attempts.

| Family · six specimens | Actual simulation discretization | What makes the form |
| --- | --- | --- |
| Phase · coral, strata, spinodal, nested interfaces, pearls, open veil | 96³ periodic cells, extent 64³, spacing 2/3; 48³ pilots at the same extent | Conserved semi-implicit Cahn–Hilliard, constant mobility, t40–120. Spinodal develops from weak perturbations; other initial geometries are prescribed. Coral has explicit traveling chemical-potential forcing. |
| Writing · counterspiral, engulfed script, nested script, pearl nuclei, pulse rosary, trefoil | 96³ cells, extent 2³, spacing 1/48; 40³ pilots | Prescribed Gaussian heat paths, diffusion/cooling and threshold nutrient-to-solid conversion; optional fixed recurrent neural activity. |
| Branches · copper thicket, fault script, field comb, ice antlers, three veins, thunder root | 64³ sites; extent 2³ except field comb (3.8, 1.6, 2.4) | Harmonic potential and field-weighted face-connected dielectric growth. Nutrient gates sheet/fork paths; physical anisotropy changes field selection. |
| Mechanics · golden sail, ice cleavage, shear veil, copper lace, parted choir, held breath | Actual spring-node lattices: 20×20×3, 16³, 20×18×4, 20×16×5, 18×18×8, 18×18×10 | Disordered central-force bonds under grip loading, with irreversible strain failure. Thin specimens have finite thickness. The 24³ constant optical carrier is not the mechanical discretization. |
| Waves/neural · crossed choir, pearl cavity, memory rings, oblique fronts, neural scroll, neural lattice | 96³ cells, extent 2³, spacing 1/48 | Four damped scalar-wave/excitable-deposition compositions; two fixed recurrent tanh fields. Pressure and deposited memory are distinct observations. |
| Ecology · gyroid gardens, tidal lamellae, pearl niches, root canals, finite orchard, crossfeed foam | 32³ medium cells, extent 24³, spacing .75; continuous 3D agent positions | Metabolic agents, finite substrates, paid boundary maintenance and phase remodeling in prescribed niches; 1,000 steps. Markers show surviving agents, not voxel samples. |

All six studies simulate on CPU. From the repository root, these are the actual
CLI selectors; environment limits keep these bakes bounded:

```sh
export JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
julia --project=. bin/interior.jl study collection_writing
julia --project=. bin/interior.jl study collection_branches
julia --project=. bin/interior.jl study collection_mechanics 1,2,3,4,5,6
WAVES_N=96 julia --project=. bin/interior.jl study collection_waves 1,2,3,4,5,6
ECOLOGY_STEPS=1000 julia --project=. bin/interior.jl study collection_ecology 1,2,3,4,5,6
julia --project=. bin/interior.jl study collection_phase \
  experiments/collection/phase/phase-{driven-coral-aged,folded-strata-aged,gyroid-cathedral-aged,nested-interfaces-aged,pearl-orchard-aged,perforated-veil-open}.toml
```

Writing also accepts recipe-name substrings, e.g. `collection_writing trefoil`;
`--pilot` selects its 40³ bake. Mechanics accepts `pilot`. Branches and phase
accept explicit recipe paths. Waves, ecology and mechanics use their study's
numbered specification order, not alphabetical filename order. Their generated
TOMLs describe runs; they are not accepted as CLI inputs. Output paths are fixed
under this batch: phase/branches skip completed `run.json`; other studies rewrite
the selected outputs. Use a separate checkout for a fresh comparison.

Each specimen directory retains `recipe.toml`, `run.json`, recorded process state,
`volume/volume.json` and actual Float32 fields, plus `render-job.json`. Mechanics
adds deformed surviving `channels.obj` and `graph.json`; ecology adds
`population.json`. Phase adds quarter-time fields/slices, preserved study source
and `selection-evidence.json`. Run metadata records parameters, metrics and source
hashes; the phase batch evidence quantifies initial-to-final change and mass drift.

The photograph is a separate optical interpretation. Palette, filament radius,
inclusion threshold, surface smoothing, section cuts, carrier, camera and light
are styling. Branch display coordinates normalize its anisotropic lattice;
mechanics displays actual deformed bond segments. Emissive pressure/activity is
an artistic transfer, not simulated light emission. Neither neural family claims
training, and ecology makes no adaptation or multicellularity claim.

The resumable queue uses Blender Cycles; final optics can use Metal independently
of CPU simulation. On the configured macOS Blender installation:

```sh
python3 scripts/render_collection.py --selected --stage preview --workers 1 --threads 2
python3 scripts/render_collection.py --selected --approved --stage final \
  --device METAL --workers 1 --threads 2
```

Final jobs request 1600² pixels and at least 96 maximum samples with adaptive
sampling. `final-status.json` distinguishes completion from a pending job;
`final.png.optical.json` records transfer settings and device. The retained
`final.blend`, VDB and log support optical replay. Publication through
`scripts/curate_collection.py outputs/collection-20261015` requires a matching
visual review and writes photographs under `site/collection/collection-20261015/`.

All thirty-six selected finals completed visual review with matching source-image
hashes before publication. The gallery retains their quality-95 WebP photographs,
small run/recipe/optical metadata and six family filters. Raw simulation fields and
full PNG originals remain under the local output batch.
