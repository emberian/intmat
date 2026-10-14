# Fraser Price's neural material experiments

Research checked 14 October 2026 through public web search, authenticated read-only Scry queries of the public X archive, and timestamped browser inspection of the October recording. Scry returned author timeline text, not simulation code or video analysis. Direct X web fetching was blocked; the original video media could be played directly in a browser. The archive extent ended on 12 October; absence from it does not establish absence from X.

## What the author actually describes

The July experiment is called **The MOld** in [Fraser's 8 July post](https://x.com/fraserpricee/status/2074902548763258948). It is his configurable world of agents with neural controllers, resources, reproduction and mutation. It should not be attributed to ALIEN, Lenia, Growing NCA, or the 2026 Neural Particle Automata paper without further evidence.

The architecture also changed across posts. [1 July](https://x.com/fraserpricee/status/2072295607222325587) describes approximately 64 genes, including aggression, affinity and reproduction parameters, plus randomly initialized MLP controllers, local senses/actions, random spawns, and a Gaussian energy field. [7 July](https://x.com/fraserpricee/status/2074454162146787791) describes a per-particle RNN whose local senses determine energy allocation among configurable action primitives. These are evidence of variants, not a single frozen specification.

The [8 July action description](https://x.com/fraserpricee/status/2074902542039859507) lists movement, energy harvesting, splitting, signaling, sharing, fighting, attraction and repulsion. Actions have energy costs and inefficiencies described as heat losses; zero energy causes death. [The following post](https://x.com/fraserpricee/status/2074902543553929266) says agents choose when to split and can spend more energy on reproductive fidelity: lower spending leads to more mutation, including both weights and architecture. This is an unusually interesting ingredient: the cost of preserving a constitution competes with the cost of surviving now.

[9 July](https://x.com/fraserpricee/status/2075243081222275080) describes a million neural controllers with mutation and no backpropagation or training objective; color approximately represents genome. A [related post](https://x.com/fraserpricee/status/2075243093046010217) describes agents spending more on movement in rapidly shifting environments while drag and environmental forcing are changed. These are author reports; they do not independently establish convergence, physical thermodynamics, or open-ended evolution.

September includes a **fly-connectome variant**. [15 September](https://x.com/fraserpricee/status/2099865719823339912) lists sensing, signaling, movement, splitting, sharing and death under connectome-based control; [another post](https://x.com/fraserpricee/status/2099884910580281347) describes 100,000 fly connectomes, while [an earlier post that day](https://x.com/fraserpricee/status/2099653615359037529) describes one million fly brains. No accessible specification found in this research establishes the neuron reduction, dynamics, biological calibration, or how those scale claims relate. Preserve the distinction between a connectome-derived controller and a validated simulation of an actual fly brain.

October returns to randomly mutated neural/RNN particles in a heat bath. [6 October](https://x.com/fraserpricee/status/2107577589216784671) describes structures resembling membranes and nuclei; [11 October](https://x.com/fraserpricee/status/2109077352940536254) describes hierarchical multicellular behavior. The resemblance and the author's interpretation should not be promoted into a biological finding without lineage, functional organization and perturbation evidence.

## Footage and reproducibility

- [8 July full run from random noise](https://x.com/fraserpricee/status/2074902550969467292/video/1). The author says larger/high-energy agents shrink and die as basal costs increase, lattices form, and energy injection changes during the recording.
- [9 July compilation and thread](https://x.com/fraserpricee/status/2075243081222275080). The [last post](https://x.com/fraserpricee/status/2075243095550017754) says the clips are from one simulation.
- [6 October full simulation](https://x.com/fraserpricee/status/2107577593289748553/video/1), paired with the membrane/nucleus description above.
- [15 September connectome-controlled agents](https://x.com/fraserpricee/status/2099865719823339912), a separate useful reference.

The July and October video URLs were confirmed by the original `t.co` redirects. Source code, a downloadable reproducible configuration, a standalone playable demo, exact integration equations, and a benchmark protocol were not established. The July performance report of 500,000 agents with up to 10,000 parameters each in real time on one GPU is an author claim, not a benchmark we ran. Do not size our first implementation from it: weight storage alone would reach 20 GB if every agent actually held 10,000 distinct Float32 weights.

Public FxTwitter metadata independently returned those two original post IDs and video media: the July recording is 1,110.29 seconds (about 18.5 minutes), the October recording 2,836.9 seconds (about 47.3 minutes). [Fraser's public GitHub account](https://github.com/fraserprice), cross-linked via his [Hugging Face model page](https://huggingface.co/fraserprice/Ternary-Bonsai-2-27B-vllm), had 34 repositories in the [public repository inventory](https://api.github.com/users/fraserprice/repos?per_page=100); none was identified as this implementation. [Another video](https://x.com/fraserpricee/status/2107119547660566771) and its [continuation](https://x.com/fraserpricee/status/2107119550181294163) refer viewers back to the descriptive thread rather than a standalone interactive app. This is a bounded search result, not a claim that code does not exist.

## Actual visual observations

The [6 October full recording](https://x.com/fraserpricee/status/2107577593289748553/video/1) was inspected directly in the browser at separated playback times. After ember confirmed the download was complete, its already-local 720×720 H.264 file was checked with ffprobe (2,836.9 seconds, 625,490,609 bytes) and frames at 15:00, 30:00, 40:00, 40:05 and 46:00 were inspected at native resolution. Its local path is `/Users/ember/dev/intmat/.playwright-mcp/QA834CEv-VlI6UpF.mp4`. Chrome's download history records the corresponding temporary UUID `52d5f47a-d0f8-4395-a64c-6d47a4d6e2c5`; that temporary path no longer exists. No new media download was needed for the local inspection.

These samples cannot establish individual interactions or lineage. Playback timestamps are recording time, not calibrated simulated time. The overlay shows zoom **4.2× at 30:00** versus **1.0× at 40:00 and 46:00**, so apparent scale changes across these samples are partly observational. At 40:00 and 46:00 the displayed population is `999999`; that display does not identify whether a simulation capacity was reached.

| Recording time | Visible structure |
| --- | --- |
| 00:30 | Numerous fine dots over a mottled green/gray background, with broad soft green patches and two pale rounded concentrations on the left. |
| 05:08 | A much sparser field of separated small colored dots over a comparatively even gray/olive background. |
| 15:00 | A dense multicolored congregation occupies the lower-left part of the view, while the upper region remains much sparser. |
| 30:00 | A roughly circular dark rim encloses a dense mixture of tiny differently colored points near the lower center. Larger colored dots occur inside and outside. The apparent diffuse violet interior in the low-resolution browser view resolves into a fine particle congregation. |
| 40:00 and 40:05 | Several open or near-closed dark rings and curved chains coexist with brighter, larger colored disks and smaller scattered dots. Dark boundaries have granular thickness, and some ring interiors contain compact yellow aggregates distinct from the large flat-colored disks. Fine curved tracks are visible between structures. The arrangements differ substantially between the samples; these images do not by themselves distinguish translation, reorganization, or camera movement. |
| 46:00 | A broad yellow-filled band with a dark boundary spans the top, accompanied by smaller dark-rimmed yellow compartments near the bottom. Granular texture is visible inside the band. |

The most striking observed transition is from an almost unstructured swarm to boundaries enclosing visually distinct interiors. A large enclosed band and small compartments coexist at the late timestamp. That difference of scale is valuable for Interior: local activity can form a boundary, then make a new domain where different processes occur. The images do not prove that these are membranes, nuclei, or reproducing organisms; those functional interpretations remain the author's claims.

Native-resolution inspection strengthens a more specific observation: the mature structures have **different spatial roles**, with densely granular interiors, thick granular rims, compact inclusions, and sparse surrounding particles. The desired extension is therefore more than cohesive moving clusters. Interior should admit interface-attached activity and material properties that differ between rim, core and exterior, with transport between those domains. Do not bake those roles into a creature template; make them possible consequences of secretion, binding, resource depletion and permeability.

An especially useful artistic extension is a boundary that changes permeability and optical response. It could gather particles, trap a reagent, sustain an interior oscillation, then rupture and release a cloud. That makes the observed compartmental morphology consequential without requiring a predetermined creature body or a claim of biological identity.

## What to carry into Interior

The valuable direction is **a world whose economics shapes its organization**. Different energy costs, spatial forcing and exchange opportunities make a controller's choices consequential. We can borrow that direction without installing every anthropomorphic action. A material should be able to spend resources on transport, secretion, contraction, signaling, or maintaining a boundary. In a first piece, uptake, movement and deposition already create a substantial ecological loop.

Five particularly useful ingredients:

1. **Internal state travels with matter.** A particle retains its recurrent state and identity while it changes neighbors. Moving it should move its memory, not resample an unrelated grid controller.
2. **An externally driven world has a visible geography.** A shifting supply field can cause migration; a slow supply front can become bands in amber. Record forcing in the simulation state so it is part of the piece's score.
3. **Every action competes for a finite inventory.** Movement and secretion cannot independently spend the same resource. A shared allocation stage connects neural choices to conserved exchange.
4. **Sharing and affinity allow collective structure.** First make pair exchange symmetric and bounded, then let a rule request it. Stable clusters can create material interfaces; an explicit split primitive still remains supplied by the model.
5. **Reproduction fidelity has a price.** Later, distinct inherited rules can trade accurate copying against mutation. This requires lineage and explicit inheritance semantics; mutation of the controller every frame is a different phenomenon.

My artistic preference is to translate genome-colored lattices into structures with different refractive or absorptive properties. Keep genome colors as a diagnostic mode. Let the specimen's enduring appearance come from what the agents consumed, deposited, bound, and damaged. A recurring contraction that writes a thin opaque seam into otherwise clear matter would connect these experiments to Interior more strongly than reproducing their flat dot colors.

The immediate Julia project should therefore prove moving recurrent state plus bounded grid exchange and one irreversible material trace. The next experiment can compare mobile and fixed controllers under changing resource geography. Evolution becomes a subsequent selectable process, once persistence, copying costs and rule provenance have actual state structures. This preserves the inspiration's ambition while giving it a composable material implementation.
