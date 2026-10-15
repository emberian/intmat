module Interior
using KernelAbstractions
using Random
using LinearAlgebra
using JSON3
using TOML
using SHA
using FFTW
using Serialization
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
include("implicit_phase.jl")
export ImplicitPhaseSeparation
include("compartments.jl")
export ActiveMedium, medium_observe!
include("mechanics.jl")
export BondNetwork, bond_forces, mechanical_energy, break_bonds!, mechanical_optical_fields, export_bond_graph
include("ecology.jl")
export Ecology, ecology_mass, reorder_ecology!, ecology_step!, ecology_statistics
include("neural_batch.jl")
export neural_controller!, neural_controller_backend!
include("neural_ecology.jl")
export NeuralEcology
include("neural_workspace.jl")
export NeuralWorkspace
include("resources.jl")
export ResourcePool, allocate_resources!, resource_add!, resource_convert!, resource_mass, resource_potential, resource_statistics
include("cooperative_boundaries.jl")
export BoundaryMaintenance, maintain_boundary!, boundary_material_mass, local_boundary_benefit
include("metabolic_ecology.jl")
export MetabolicEcology
include("cellular_transport.jl")
include("cellular_contacts.jl")
include("cellular.jl")
export CellularEcology, CellBonds, CellStencil, cellular_statistics, cellular_radii
export add_cell_bond!, cellular_pairs, cellular_mechanics!, exchange_cell_metabolites!
export rebuild_cell_stencil!, sample_concentrations, sample_concentration_gradients, allocate_cell_resources!, deposit_cell_resource!
include("experiments.jl")
include("catalytic_habitat.jl")
export CatalyticHabitat, catalytic_statistics
include("sorptive_matrix.jl")
export SorptiveMatrix, matrix_statistics
include("constructive_cells.jl")
export ConstructiveEcology, anabolic_statistics
export Experiment, ExperimentReport, ProcessPorts, process_ports, OpticalMap, run!, experiment_diagnostics, owned_material_mass
include("recurrent_material.jl")
export RecurrentMaterial
include("dose_writer.jl")
export DoseThermalWriter, writer_material!
include("render.jl")
include("grin.jl")
include("scenes.jl")
include("arena.jl")
include("cellular_arena.jl")
include("catalytic_arena.jl")
include("device_cellular.jl")
include("device_cellular_protocol.jl")
export DeviceCellular, device_cellular_step!, device_cellular_status!, device_cellular_statistics
export device_cellular_initialize, device_support_descriptor, device_support_edges
export device_cellular_snapshot, device_cellular_audit, device_cellular_checkpoint, device_cellular_restore
include("device_habitat.jl")
export DeviceCatalyticHabitat, DeviceConstructiveEcology
export device_habitat_snapshot, device_habitat_audit, device_habitat_checkpoint, device_habitat_restore
export Arena, Chemostat, arena_step!, arena_snapshot, checkpoint!, restore_arena, run_arena!, build_arena
export RenderSettings, render_volume, export_volume, volume_sample, transport_segment, build_scene, optical_fields
end
