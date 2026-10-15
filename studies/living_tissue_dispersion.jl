# Frozen-coefficient neural dispersion of the continuum turnover law. Tissue,
# fuel and matrix evolve in actual runs; this is not an equilibrium proof.
using Interior,LinearAlgebra,SHA
using Interior:JSON3
include("living_tissue_turnover.jl")
using .LivingTissueTurnover

function tissue_dispersion(output)
    ispath(output)&&error("fresh output required")
    curves=Any[]
    for gain in (2.5,2.8,3.1),biomass in (.12,.26,.5),excitation_diffusion in (.03,.08)
        q=turnover_base_parameters(neural_gain=gain,excitation_diffusion=excitation_diffusion,
            inhibition_gain=2.3,inhibition_time=4.5)
        n=.6;s=biomass/(q.tissue_saturation+biomass);v=s*n/(q.half_saturation+n)
        # x=y=memory=matrix=0, paid=1. Tissue and fuel are held fixed here.
        a=v*(q.neural_gain-1)-(1-s)
        reaction=[a -v*q.inhibition_gain -v*q.memory_feedback;
            inv(q.inhibition_time) -inv(q.inhibition_time) 0;
            inv(q.memory_time) 0 -inv(q.memory_time)]
        dx=q.excitation_diffusion*s*s;dy=.15
        measurements=Any[]
        for side in (64,128),direction in ("axis","diagonal")
            h=24/side;maximum_k=direction=="axis" ? pi/h : sqrt(2)*pi/h
            samples=Any[]
            for k in range(0,maximum_k;length=401)
                symbol=direction=="axis" ? 4sin(k*h/2)^2/h^2 : 8sin(k*h/(2sqrt(2)))^2/h^2
                jacobian=copy(reaction);jacobian[1,1]-=dx*symbol;jacobian[2,2]-=dy*symbol
                values=eigvals(jacobian);dominant=values[argmax(real.(values))]
                push!(samples,(k=k,laplacian_symbol=symbol,growth=real(dominant),frequency=abs(imag(dominant))))
            end
            peak=samples[argmax([x.growth for x in samples])]
            push!(measurements,(side,direction,peak,wavelength=peak.k>0 ? 2pi/peak.k : nothing,samples))
        end
        push!(curves,(gain,biomass,excitation_diffusion,inhibition_diffusion=dy,parameters=q,
            tissue_support=s,vitality=v,effective_excitation_diffusion=dx,reaction_jacobian=reaction,
            zero_mode_eigenvalues=[(real=real(z),imag=imag(z)) for z in eigvals(reaction)],measurements))
    end
    open(io->JSON3.pretty(io,(model_sha256=LivingTissueTurnover.TURNOVER_SHA,
        base_sha256=LivingTissueTurnover.LivingTissueStudy.MODEL_SHA,curves,
        interpretation="Frozen-coefficient3-channel neural linearization at zero hidden state, no matrix, precursor.6, paid fraction1. Tissue/fuel/reserve are held fixed despite evolving in the actual law. Axis and diagonal curves use the discrete stencil symbol. This supplies a wavelength/underresolution hypothesis; it is not a nonlinear pattern, stability, or biological-organization guarantee.")),output,"w")
end
if abspath(PROGRAM_FILE)==@__FILE__;length(ARGS)==1||error("fresh-output.json");tissue_dispersion(only(ARGS));end
