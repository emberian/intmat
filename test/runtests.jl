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
include("device_inventory_precision.jl")
include("device_material_sensing.jl")
include("device_cellular_protocol.jl")
end

module DeviceResourceTests
include("device_resource_exchange.jl")
include("device_resource_aliases.jl")
end

module DeviceEnergyTests
include("device_energy_audit.jl")
end

module DeviceLifecycleTests
include("device_lifecycle_prefix.jl")
include("device_lifecycle_workspace.jl")
include("device_lifecycle_habitat_aliases.jl")
end

module GardenOperatingMapTests
include("garden_operating_map.jl")
end

module ControllerCapabilityTests
include("controller_capability_atlas.jl")
end

module ConstructiveBatchTests
include("constructive_world_batches.jl")
include("device_material_sensing_batches.jl")
end

module ProportionalDebitTests
include("cellular_proportional_debit.jl")
end

module NativeActionCreditTests
include("native_action_credit.jl")
end

module NativeWorkbenchTests
include("native_workbench_worker.jl")
end

module NativeLifetimeCreditTests
include("native_lifetime_credit.jl")
end

module FeedbackFactorizationTests
include("feedback_factorization.jl")
end

module HabitatActuatorRequestTests
include("habitat_actuator_requests.jl")
end

module LearningWorldBatchTests
include("learning_world_batches.jl")
end

module NativeConstructionCreditTests
include("native_construction_credit.jl")
end

module NativeJointConstructionCreditTests
include("native_joint_construction_credit.jl")
end

module ImplicitPhaseRuntimeTests
include("implicit_phase_runtime.jl")
end

module NativeConstructionLagTests
include("native_construction_lag.jl")
end

module ConstructionLagWorldBatchTests
include("construction_lag_world_batches.jl")
end

module LocalTemporalFeatureTests
include("local_temporal_features.jl")
end

module NativeCohortConstructionCreditTests
include("native_cohort_construction_credit.jl")
end

module NativeLivingConstructionTests
include("native_living_construction.jl")
end

module CohortConstructionWorldBatchTests
include("cohort_construction_world_batches.jl")
end
