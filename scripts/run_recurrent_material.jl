#!/usr/bin/env julia
using Interior,Statistics
using Interior: JSON3
const loaded_model_source=read(joinpath(@__DIR__,"..","src","recurrent_material.jl"),String)
isdefined(Interior,:RecurrentMaterial) || Base.include(Interior,joinpath(@__DIR__,"..","src","recurrent_material.jl"))
using Interior: RecurrentMaterial
using Metal
Metal.functional() || error("Metal GPU unavailable")
function main()
    side=length(ARGS)>0 ? parse(Int,ARGS[1]) : 64
    steps=length(ARGS)>1 ? parse(Int,ARGS[2]) : 3000
    root=length(ARGS)>2 ? ARGS[3] : "outputs/recurrent-material"
    side>0 && steps>0 || throw(ArgumentError("side and steps must be positive"))
    profile=get(ENV,"RECURRENT_PROFILE","all")
    profile in ("all","braided","slow-memory","threshold-feedback-off") || throw(ArgumentError("unknown profile"))
    ispath(root) && throw(ArgumentError("output directory already exists; use a fresh destination"))
    milestones=Set(max.(1,round.(Int,steps.*(.25,.5,.75,1.))))
    mkpath(dirname(abspath(root)))
    mkdir(root)
    write(joinpath(root,"loaded-model-source.jl"),loaded_model_source)
    cp(@__FILE__,joinpath(root,"runner-source.jl");force=true)
    for (name,kwargs) in (("braided",(coupling=.12,feedback=1.2,writing=.25,erosion=.06,contrast=.7)),
                          ("slow-memory",(coupling=.28,feedback=2.,writing=.5,erosion=.035,contrast=1.1)),
                          ("threshold-feedback-off",(coupling=.12,feedback=0.,writing=.25,erosion=.06,contrast=.7)))
        profile!="all" && name!=profile && continue
        w=World(Grid(side,side,side;extent=Float32.((side,side,side)));backend=MetalBackend(),seed=31)
        p=RecurrentMaterial(w;kwargs...); dir=joinpath(root,name);mkpath(dir)
        before=sum(Float64,Array(w.fields.nutrient))+sum(Float64,Array(w.fields.deposit));samples=Any[]
        open(joinpath(dir,"volume.json"),"w") do io
            JSON3.pretty(io,(shape=w.grid.size,extent=w.grid.extent,dtype="float32",order="Julia-column-major",dt=.05,capture_every=100,latent_channels=8))
        end
        # Compile once, restore initial condition, then time only evolution.
        advance!(w,p,.05); Metal.synchronize()
        w=World(w.grid;backend=MetalBackend(),seed=31);p=RecurrentMaterial(w;kwargs...)
        started=time()
        for k in 1:steps
            advance!(w,p,.05);w.time+=Float64(Float32(.05))
            if k%100==0 || k==steps
                activity=Array(w.fields.activity);solid=Array(w.fields.deposit);latent=Array(p.state)
                push!(samples,(step=k,time=w.time,activity_std=std(activity),deposit_mean=mean(solid),deposit_std=std(solid),latent_std=std(latent),material_error=sum(Float64,Array(w.fields.nutrient))+sum(Float64,solid)-before,finite=all(isfinite,latent)))
                # Compact full 3D fields, including eight latent channels at four milestones.
                for (field,array) in (("activity",activity),("deposit",solid))
                    open(joinpath(dir,"$field-$k.f32"),"w") do io;write(io,array);end
                end
            end
            if k in milestones
                record_state(w,p,joinpath(dir,"state-$k");step=k)
            end
        end
        elapsed=time()-started
        result=(side=side,steps=steps,dt=.05,backend="Metal",seconds=elapsed,voxel_steps_per_second=side^3*steps/elapsed,settings=kwargs,samples=samples)
        open(joinpath(dir,"metrics.json"),"w") do io;JSON3.pretty(io,result);end
        println(name," ",side,"³ ",round(result.voxel_steps_per_second/1e6;digits=2)," M voxelsteps/s final ",last(samples));flush(stdout)
    end
end
main()
