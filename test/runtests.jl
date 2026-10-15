using Interior, Test, Random
include("core.jl")
include("particles.jl")
include("render.jl")
isfile(joinpath(@__DIR__,"scenes.jl")) && include("scenes.jl")

include("neural.jl")
include("audit.jl")
include("clock.jl")

include("recording.jl")
include("acoustics.jl")
include("writing.jl")
include("dose_writer.jl")
include("electrical.jl")
include("phase.jl")
include("mechanics.jl")
include("mechanics_allocations.jl")
include("ecology.jl")
include("experiments.jl")
include("process_composition.jl")

include("grin.jl")

include("implicit_phase.jl")
include("compartments.jl")
include("neural_ecology.jl")
include("arena.jl")
include("arena_composition.jl")

include("neural_batch.jl")
include("neural_workspace.jl")
include("resources.jl")
include("cooperative_boundaries.jl")
include("metabolic_ecology.jl")
include("metabolic_composition.jl")
include("cellular_transport.jl")
include("cellular_contacts.jl")
include("cellular.jl")
include("cellular_storage.jl")

include("cli.jl")

include("recurrent_material.jl")
include("recurrent_runner.jl")

include("catalytic_habitat.jl")
include("sorptive_matrix.jl")
include("constructive_cells.jl")
include("matrix_mechanics.jl")
include("cellular_rounding_audit.jl")

module ElectricalMultigridTests
include("electrical_multigrid.jl")
end

module DeviceCellularTests
include("device_cellular.jl")
include("device_support_carry.jl")
include("cellular_initialization.jl")
include("device_habitat.jl")
include("device_habitat_sorption.jl")
include("device_habitat_support.jl")
include("device_habitat_audit.jl")
include("device_cellular_protocol.jl")
end

module DeviceResourceTests
include("device_resource_exchange.jl")
include("device_resource_aliases.jl")
end

module DeviceEnergyTests
include("device_energy_audit.jl")
end
