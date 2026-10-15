# Callable native reduction gate; caller supplies a complete small habitat.
# Mutated scratch is restored before return. Not a simulation backend.
function device_habitat_observation_gate(w,p)
    expected=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p)).material
    observed=Interior.owned_material_mass(p)
    isapprox(observed,expected;atol=1e-8,rtol=1e-12) || error("scalar material reduction differs from detached audit")
    Interior.process_state_finite(p) || error("valid habitat failed backend finiteness reductions")
    expression=Array(p.habitat.arrays.expression)
    p.habitat.arrays.expression.=NaN32
    !Interior.process_state_finite(p) || error("habitat NaN expression escaped backend finite reduction")
    copyto!(p.habitat.arrays.expression,expression)
    catalyst=Array(p.habitat.arrays.enzyme_bits)
    p.habitat.arrays.enzyme_bits.=reinterpret(UInt32,NaN32)
    !Interior.process_state_finite(p) || error("NaN catalyst bits escaped backend finite reduction")
    copyto!(p.habitat.arrays.enzyme_bits,catalyst)
    Interior.process_state_finite(p) || error("restored observation fixture is not finite")
    observation=Interior.record_process_fields(p)
    observation.habitat_weights==Array(p.habitat.arrays.weights) || error("detached named observation differs")
    original=Array(p.habitat.arrays.weights)[1]
    observation.habitat_weights[1]=123f0
    Array(p.habitat.arrays.weights)[1]==original || error("named observation aliases live state")
    (;owned_material=observed,expected_material=expected,material_scalar_agrees=true,
        finite=true,nan_expression_detected=true,nan_catalyst_detected=true,
        detached_named_observation=true,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN)
end
