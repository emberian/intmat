# Frozen neural coefficients only: a finite-wavelength hypothesis, not a
# stationary solution or nonlinear guarantee for the evolving material law.
using Interior, LinearAlgebra, SHA
using Interior: JSON3
include("living_tissue_turnover.jl")
using .LivingTissueTurnover

function mode_hypothesis(output;regime="slow")
    ispath(output) && error("fresh output required")
    regime in ("slow","fast","oscillatory") || error("regime must be slow, fast or oscillatory")
    if regime=="oscillatory"
        q=turnover_base_parameters(neural_gain=2.8,capacity=.9,
            excitation_diffusion=.08,inhibition_gain=2.3,inhibition_time=.5,
            memory_feedback=3.5,memory_time=8.,growth=0.,biomass_diffusion=0.)
        biomasses=(.1,.15,.26,.5,.9)
    else
        q=turnover_base_parameters(neural_gain=regime=="fast" ? 2.8 : 1.8,
            capacity=regime=="fast" ? 1.2 : .32,excitation_diffusion=.08,
            inhibition_gain=2.3,inhibition_time=regime=="fast" ? .5 : 4.5)
        biomasses=regime=="fast" ? (.1,.15,.26,.5,1.2,5.) : (.26,.28,.30,.32)
    end
    phase=turnover_parameters(inhibition_diffusion=4.)
    cases=Any[]
    for biomass in biomasses
        support=biomass/(q.tissue_saturation+biomass)
        vitality=support*.6/(q.half_saturation+.6)
        a=vitality*(q.neural_gain-1)-(1-support)
        reaction=[a -vitality*q.inhibition_gain -vitality*q.memory_feedback;
            inv(q.inhibition_time) -inv(q.inhibition_time) 0;
            inv(q.memory_time) 0 -inv(q.memory_time)]
        zero=eigvals(reaction)
        curves=Any[]
        for side in (64,128),dimensions in (2,3),direction in ("axis","diagonal")
            spacing=24/side
            samples=Any[]
            for k in range(0,pi/spacing;length=501)
                symbol=direction=="axis" ? 4sin(k*spacing/2)^2/spacing^2 :
                    4dimensions*sin(k*spacing/(2sqrt(dimensions)))^2/spacing^2
                jacobian=copy(reaction)
                jacobian[1,1]-=q.excitation_diffusion*support^2*symbol
                jacobian[2,2]-=phase.inhibition_diffusion*symbol
                values=eigvals(jacobian);dominant=values[argmax(real.(values))]
                push!(samples,(k=k,growth=real(dominant),frequency=abs(imag(dominant))))
            end
            peak=samples[argmax([s.growth for s in samples])]
            push!(curves,(side,dimensions,direction,peak,
                wavelength=peak.k>0 ? 2pi/peak.k : nothing,samples))
        end
        push!(cases,(biomass,support,vitality,reaction_jacobian=reaction,
            zero_mode=[(real=real(z),imag=imag(z)) for z in zero],curves))
    end
    result=(model_sha256=LivingTissueTurnover.TURNOVER_SHA,
        source_sha256=bytes2hex(sha256(read(@__FILE__))),base_parameters=q,
        regime,inhibition_diffusion=phase.inhibition_diffusion,cases,
        convention="Freeze N=.6, B at listed value, M=x=y=q=0 and paid fraction1. Actual B/N/M/work/reserve dynamics are omitted; the frozen state need not be an equilibrium. Inhibition and excitation use their actual physical diffusivities and discrete axis/diagonal stencil symbols. A stable zero mode and growing finite mode motivate a separate nonlinear test, not a material organization or convergence claim.")
    open(io->JSON3.pretty(io,result),output,"w")
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (1,2) || error("fresh-output.json [slow|fast|oscillatory]")
    mode_hypothesis(ARGS[1];regime=length(ARGS)==2 ? ARGS[2] : "slow")
end
