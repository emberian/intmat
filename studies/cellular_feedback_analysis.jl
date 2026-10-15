#!/usr/bin/env julia
# Reduced-mode analysis, not an Interior ecology trajectory or a fitted model.
using JSON3, Statistics

function reduced_feedback(;strength=.3,adaptation_time=6.,plasticity=.002,
        sensor_norm2=8.,recurrent=1.2,bias=.02,response_time=.8,decay=.02,
        dt=.1,horizon=1000.,late_start=600.)
    x=.02; a=0.; recurrent_mark=0.; sensory_drive=0.
    values=Float64[]; capped=0; crossings=0; previous=x
    alpha=-expm1(-dt/response_time);beta=-expm1(-dt/adaptation_time)
    for step in 1:round(Int,horizon/dt)
        # Match activation, adaptation ordering and bounded mark law. Aggregate
        # external marks along fixed sensors with squared norm sensor_norm2.
        output=tanh(Float32((recurrent+recurrent_mark)*x+sensory_drive+bias))
        target=tanh(atanh(clamp(output,-.999999f0,.999999f0))-strength*a)
        old=x;x+=alpha*(target-x);a+=beta*(x-a)
        recurrent_mark=clamp(recurrent_mark+dt*(plasticity*(x*old-x*x*recurrent_mark)-decay*recurrent_mark),-.3,.3)
        # No individual external mark is clipped here: this reduction is exact
        # for collinear frozen sensors only while their individual caps are inactive.
        sensory_drive+=dt*(plasticity*x*sensor_norm2-(plasticity*x*x+decay)*sensory_drive)
        if step*dt>=late_start
            push!(values,x);capped+=abs(output)>.999999f0
            crossings+=signbit(x)!=signbit(previous)
        end
        previous=x
    end
    (;strength,adaptation_time,plasticity,sensor_norm2,recurrent,bias,response_time,decay,
      dt,horizon,late_start,late_min=minimum(values),late_max=maximum(values),
      late_mean=mean(values),late_std=std(values),late_sign_crossings=crossings,
      late_activation_cap_fraction=capped/length(values),
      final_recurrent_mark=recurrent_mark,final_sensory_mark_drive=sensory_drive)
end

function feedback_analysis(output)
    ispath(output)&&error("use a fresh output path")
    runs=[reduced_feedback(;strength,adaptation_time,plasticity,sensor_norm2)
        for strength in (.3,1.5),adaptation_time in (6.,8.),
            plasticity in (0.,.002),sensor_norm2 in (0.,4.,8.)]
    result=(kind="analytical reduced-mode exploration, not ecological evidence",
        assumptions="one recurrent mode; constant external sensors; fixed inherited gain; no growth, chemistry, bonds, mutation or division; external mark caps inactive",
        equations=(hidden="dx/dt=(tanh(w*x+b+m*x+I-g*a)-x)/tau_h",
            adaptation="da/dt=(x-a)/tau_a",
            sensory_mark="db/dt=eta*x*R-(lambda+eta*x^2)*b",
            recurrent_mark="dm/dt=eta*x^2-(lambda+eta*x^2)*m"),
        note="Discrete integration follows updated hidden state and old input for recurrent plasticity, as in CellularEcology; the written equations give its continuous-time limit.",
        runs=vec(runs))
    open(io->JSON3.pretty(io,result),output,"w")
    result
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS)==1||error("usage: cellular_feedback_analysis.jl FRESH_OUTPUT.json")
    feedback_analysis(ARGS[1])
end
