#!/usr/bin/env julia
"""Bounded two-factor living patch probe; all evolution uses native composition."""
module LivingPatchRecruitment
using Interior, KernelAbstractions, TOML, SHA
using Interior: JSON3
include("native_living_construction.jl")
const N=NativeLivingConstruction

function spatial(p,h,frame)
    x=Interior.device_habitat_snapshot(p);s=x.cells.state;g=x.cells.geometry
    q=x.cells.parameters;V=prod(Int.(g.dims));volume=Float64(g.volume)
    active=findall(==(1),s.active)
    ea=Float64.(reinterpret(Float32,x.arrays.enzyme_bits[:,1]))
    eb=Float64.(reinterpret(Float32,x.arrays.enzyme_bits[:,2]));e=ea+eb
    total=sum(e);near=0.;saturated=0.;occupied=0
    for site in 1:V
        pos=Interior.dc_coords(site,g.dims)
        overlap=false
        for i in active
            r2=0.
            for d in 1:3
                dx=abs((pos[d]-.5)*g.spacing[d]-s.position[i,d])
                dx=min(dx,Float64(g.extent[d])-dx);r2+=dx*dx
            end
            overlap|=r2<=Float64(q.exchange_radius)^2
        end
        near+=overlap ? e[site] : 0.;occupied+=overlap
        for k in 1:2
            raw=Float64(x.feedstock.fields[site,k]);enz=k==1 ? ea[site] : eb[site]
            saturated+=enz*raw/(Float64(x.parameters.half_saturation)*volume+raw)
        end
    end
    a=Interior.device_habitat_audit(x)
    hdim=Int(g.hidden)
    (;time=frame.time,population=length(active),births=frame.births,deaths=frame.deaths,
        max_generation=maximum(s.generations[active];init=0),
        total_raw=sum(x.feedstock.fields[:,1:2]),total_free=sum(x.cells.environment.fields[:,1:2]),
        total_bound=sum(x.arrays.bound),total_catalyst=total,
        total_usable=sum(N.L.local_value(s,q,i) for i in active),
        converted=sum(a.solubilized),funded_construction=sum(frame.funded_secretion),
        requested_construction=sum(frame.cumulative_requests),
        catalyst_fraction_near_living=total>0 ? near/total : 0.,
        living_support_fraction=occupied/V,
        catalyst_effective_volume_fraction=total>0 ? total^2/(sum(abs2,e)*V) : 0.,
        catalyst_weighted_raw_saturation=total>0 ? saturated/total : 0.,
        material_receptor_max=[maximum(abs,s.inputs[active,2hdim+15+k];init=0f0) for k in 1:3],
        material_residual=frame.ledger.material_residual,potential_residual=frame.ledger.potential_residual,
        supplied=frame.ledger.supplied,learning_work=sum(c.learning_work for c in frame.cells))
end

function run(profile,out;duration=8.,capture_time=2.,seed=211)
    ispath(out)&&error("fresh output required");mkpath(out)
    sourcefiles=[@__FILE__,joinpath(@__DIR__,"native_living_construction.jl"),joinpath(@__DIR__,"native_cohort_construction_credit.jl"),
        joinpath(@__DIR__,"native_joint_construction_credit.jl"),joinpath(@__DIR__,"native_construction_credit.jl"),joinpath(@__DIR__,"native_lifetime_credit.jl")]
    hashes=Dict(basename(f)=>bytes2hex(sha256(read(f))) for f in sourcefiles)
    for f in sourcefiles;write(joinpath(out,"source-"*basename(f)),read(f));end
    base=TOML.parsefile(profile);runs=Any[]
    for raw_diffusion in (0.,.006),decay in (.2,.05)
        name="raw$(raw_diffusion)-decay$(decay)";cfg=deepcopy(base);cfg["seed"]=seed
        cfg["feedstock"]["diffusion"]=[raw_diffusion,raw_diffusion];cfg["catalysis"]["decay_rate"]=decay
        recipe=joinpath(out,name*".toml");open(io->TOML.print(io,cfg),recipe,"w")
        f=N.fixture(recipe);w,p=f.w,f.p;dt=cfg["native"]["dt"]
        h=N.sidecar(p;dt,horizon=2.,act_time=duration,rate=.003,learning=true,seed=6201,inheritance_fraction=.5)
        origin=N.C.A.audit(p);fr=N.frame(w,p,h,origin);samples=Any[spatial(p,h,fr)];started=time()
        for _ in 1:round(Int,duration/capture_time)
            N.L.advance!(w,p,h,capture_time;dt)
            fr=N.frame(w,p,h,origin);push!(samples,spatial(p,h,fr))
            println((;arm=name,samples[end]...));flush(stdout)
        end
        open(io->JSON3.write(io,fr),joinpath(out,name*"-final.json"),"w")
        push!(runs,(;raw_diffusion,decay,samples,elapsed_seconds=time()-started))
    end
    report=(;schema="living-patch-recruitment-v1",seed,duration,capture_time,source_hashes=hashes,
        scope="new independent worlds, same finite funding, moving/dividing paid policy; no ecological conclusion from prefix",runs)
    open(io->JSON3.write(io,report),joinpath(out,"summary.json"),"w");report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    LivingPatchRecruitment.run(ARGS[1],ARGS[2];duration=length(ARGS)>2 ? parse(Float64,ARGS[3]) : 8.)
end
