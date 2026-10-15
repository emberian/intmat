# Matched continuum-tissue feedback contrasts, not an adaptation assay.
using Interior, TOML, SHA, Statistics
using Interior: JSON3
include("living_tissue.jl")
using .LivingTissueStudy

function tissue_components(mask)
    nx,ny=size(mask);seen=falses(nx,ny);queue=CartesianIndex{2}[];sizes=Int[]
    for start in CartesianIndices(mask)
        (!mask[start] || seen[start]) && continue
        empty!(queue);push!(queue,start);seen[start]=true;head=1
        while head<=length(queue)
            i,j=Tuple(queue[head]);head+=1
            for neighbor in (CartesianIndex(mod1(i-1,nx),j),CartesianIndex(mod1(i+1,nx),j),
                CartesianIndex(i,mod1(j-1,ny)),CartesianIndex(i,mod1(j+1,ny)))
                if mask[neighbor] && !seen[neighbor];seen[neighbor]=true;push!(queue,neighbor);end
            end
        end
        push!(sizes,length(queue))
    end
    sort!(sizes;rev=true)
end
function tissue_topology(state)
    size(state,3)==1 || throw(ArgumentError("2D topology observation required"))
    mask=state[:,:,1,2].+state[:,:,1,3].>.15f0
    components=tissue_components(mask);voids=tissue_components(.!mask)
    # These are periodic grid components at an explicit threshold. A separate
    # background component is not interpreted as a physical surrounding wall.
    (threshold=.15,occupied_components=length(components),occupied_sizes=components,
        void_components=length(voids),void_sizes=voids)
end
function tissue_ablation(recipe,output)
    c=TOML.parsefile(recipe);ispath(output)&&error("fresh output required")
    mkpath(output);cp(recipe,joinpath(output,"recipe.toml"));cp(@__FILE__,joinpath(output,"ablation-source.jl"))
    write(joinpath(output,"study-source.jl"),LivingTissueStudy.MODEL_SOURCE)
    cp(joinpath(@__DIR__,"living_tissue_runner.jl"),joinpath(output,"living_tissue_runner.jl"))
    for (name,source) in LivingTissueStudy.SHARED_SOURCES;write(joinpath(output,name),source);end
    q=tissue_parameters(;NamedTuple(Symbol(k)=>v for (k,v) in c["model"])...)
    dims=(128,128,1);grid=Grid(dims...;extent=Tuple(c["extent"]));dt=.025f0;steps=3200
    arms=(("full",q),("barrier-off",merge(q,(barrier=0f0,))),
        ("activity-control-off",merge(q,(activity_control=0f0,))),
        ("both-off",merge(q,(barrier=0f0,activity_control=0f0))))
    results=Any[];finals=Dict{String,Array{Float32,4}}();matched_gate=.5f0
    for (name,parameters) in arms
        parameters.activity_control==0 && (parameters=merge(parameters,(constant_activity=matched_gate,)))
        w=tissue_world(grid;seed=c["seed"]);p=LivingTissue(w;parameters,pattern=c["seed_pattern"])
        initial=tissue_statistics(w,p);history=Any[];probe=Tuple(argmax(view(p.state,:,:,:,2)))
        previous=p.state[probe...,5];crossings=0
        elapsed=@elapsed for k in 0:steps
            if k>0
                step!(w,p,dt)
                x=p.state[probe...,5];crossings+=x*previous<0;previous=x
            end
            if k%100==0
                state=Array(p.state);s=tissue_statistics(w,p;state)
                weighted_gate=sum(Float64,max.(0f0,view(state,:,:,:,5)).*view(state,:,:,:,2))
                push!(history,merge((step=k,topology=tissue_topology(state),biomass_weighted_positive_gate=weighted_gate,
                    material_residual=s.material-initial.material),s))
            end
        end
        name=="full" && (matched_gate=Float32(sum(h.biomass_weighted_positive_gate for h in history)/sum(h.biomass for h in history)))
        finals[name]=Array(p.state);tissue_checkpoint(w,p,joinpath(output,"$name.jls"))
        push!(results,(name,parameters,elapsed_seconds=elapsed,fixed_site_crossings=crossings,history))
        println(name," ",elapsed,"s ",last(history));flush(stdout)
    end
    baseline=finals["full"]
    differences=Dict(name=>(biomass_rms=sqrt(mean(abs2,view(state,:,:,:,2).-view(baseline,:,:,:,2))),
        matrix_rms=sqrt(mean(abs2,view(state,:,:,:,3).-view(baseline,:,:,:,3))),
        excitation_rms=sqrt(mean(abs2,view(state,:,:,:,5).-view(baseline,:,:,:,5))),
        maxabs=maximum(abs.(state.-baseline))) for (name,state) in finals)
    report=(dimensions=dims,extent=grid.extent,dt,steps,duration=steps*dt,arms=results,constant_gate=matched_gate,
        final_differences_from_full=differences,
        study_source_sha256=LivingTissueStudy.MODEL_SHA,
        source_sha256=bytes2hex(sha256(read(@__FILE__))),
        interpretation="Same initial state and parameters except the named interventions. Activity-control-off replaces the dynamic positive-phase gate by the full arm's biomass-weighted mean positive gate at 100-step observations. This controls approximate mean gating, not realized total construction or work (trajectories differ). It tests feedback consequences; not adaptive benefit, organisms, or evolution.")
    open(io->JSON3.pretty(io,report),joinpath(output,"ablation.json"),"w")
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2||error("recipe fresh-output");tissue_ablation(ARGS...)
end
