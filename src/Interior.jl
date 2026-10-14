module Interior
using KernelAbstractions
using Random
using LinearAlgebra
using JSON3
using TOML
using SHA
export Grid, World, CPU, Excitable, Deposition, WaveForcing, Diffusion, advance!, step!, seed_activity!, material_mass, stability_limit
include("core.jl")
include("fields.jl")
include("neural.jl")
export NeuralField
include("particles.jl")
export ParticleRule, ParticleCarrier, reorder_particles!, total_particle_nutrient, particle_step!
include("recording.jl")
export record_state, record_process_fields
include("acoustics.jl")
export AcousticField, seed_acoustic!
include("writing.jl")
export ThermalWriter, beam_position
include("electrical.jl")
export DielectricGrowth, solve_potential!, potential_residual, electrical_optical_fields
include("phase.jl")
export PhaseSeparation, phase_energy, phase_optical_fields
include("mechanics.jl")
export BondNetwork, bond_forces, mechanical_energy, break_bonds!, mechanical_optical_fields, export_bond_graph
include("ecology.jl")
export Ecology, ecology_mass, reorder_ecology!, ecology_step!, ecology_statistics
include("experiments.jl")
export Experiment, ExperimentReport, ProcessPorts, process_ports, OpticalMap, run!, experiment_diagnostics, owned_material_mass
include("render.jl")
include("grin.jl")
include("scenes.jl")
export RenderSettings, render_volume, export_volume, volume_sample, transport_segment, build_scene, optical_fields
end
