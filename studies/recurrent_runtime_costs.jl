#!/usr/bin/env julia
# Explicit arithmetic/memory envelope, not a hardware throughput benchmark.
using JSON3
function recurrent_runtime_costs(out;capacity=512,D=25,actor=3,predictions=12)
    B=capacity;O=actor+1+predictions;rows=Any[]
    for H in (16,64,128),Tsteps in (8,32)
        P=H*(D+H+1)+O*(H+1);blocks=cld(B*Tsteps,128)
        parameter_bytes=4P
        # Ring: x, raw actions, predictive targets, reward, discount, behavior logp;
        # Int32 identity + validity. History mask is included in validity storage.
        ring_bytes=B*Tsteps*(4*(D+actor+predictions+3)+8)
        # Two full hidden trajectories (h and da), output and dY, h0/dh0 +
        # double-buffered live state, fixed partial gradients, theta/g/Adam(m,v).
        training_bytes=16P+8H*B*Tsteps+8O*B*Tsteps+16H*B+4P*blocks+ring_bytes
        forward_flops=2B*(H*D+H*H+O*H)
        block_flops=B*Tsteps*(6H*H+4H*D+8O*H) # includes head recomputation allowance
        push!(rows,(;H,Tsteps,D,O,capacity=B,parameter_count=P,shared_parameter_bytes=parameter_bytes,
            private_parameter_bytes=parameter_bytes*B,shared_adam_gradient_bytes=4parameter_bytes,
            private_adam_gradient_bytes=4parameter_bytes*B,partial_gradient_bytes=4P*blocks,
            training_envelope_bytes=training_bytes,ring_bytes,forward_flops,block_flops))
    end
    report=(;schema="recurrent-runtime-envelope-v1",assumptions="Float32; shared theta; full capacity; no generic AD tape; deterministic128-sample gradient partials; one retained block; extra replay blocks not included; physical simulation/observers/vendor workspace separate",rows)
    open(io->JSON3.write(io,report),out,"w");report
end
abspath(PROGRAM_FILE)==(@__FILE__) && recurrent_runtime_costs(ARGS[1])
