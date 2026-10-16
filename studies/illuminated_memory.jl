# Separate optical-feedback epoch; reuses the original producer initialization,
# recurrent material and GRIN tracer without changing their source.
include("refractive_memory.jl")
using Interior.KernelAbstractions

@inline function optical_cell(p,n)
    q=ntuple(d->clamp((p[d]/2+.5)*n+.5,1.,Float64(n)),3)
    lo=ntuple(d->floor(Int,q[d]),3);hi=ntuple(d->min(lo[d]+1,n),3)
    lo,hi,ntuple(d->q[d]-lo[d],3)
end
function optical_launches(grid)
    result=NamedTuple[];z=-.65
    for strand=0:2,j=1:grid,i=1:grid
        u=2*((i-.5)/grid-.5);v=2*((j-.5)/grid-.5);u*u+v*v>1&&continue
        phase=2pi*strand/3+2.6z;cx=.27cos(phase)+.08z;cy=.27sin(phase)
        push!(result,(origin=(cx+.095u,cy+.095v,z),direction=(-.27*2.6sin(phase)+.08,.27*2.6cos(phase),1.)))
    end
    result
end
"""Return absorbed power per material volume, conserving input ray weights.
Optical arc length is dimensionless (box width 2); material box width is 64.
Alpha is per optical length per unit solid, NOT per material length.
"""
function absorbed_power(solid;grid=13,ds=.002,alpha=2.,power=1000.,contrast=.65)
    n=size(solid,1);field=GradedIndex(indexfield(solid,contrast,550.));launches=optical_launches(grid)
    deposition=zeros(Float64,size(solid));emergent=0.;unresolved=0.;absorbed=0.;ends=NTuple{3,Float64}[]
    for launch in launches
        ray=trace_grin(field,launch.origin,launch.direction;ds=ds,max_length=7.,target_z=1.05,_validated=true)
        weight=power/length(launches)
        for segment=2:length(ray.path)
            a=ray.path[segment-1];b=ray.path[segment];p=ntuple(d->(a[d]+b[d])*.5,3)
            any(x->abs(x)>1,p)&&continue
            lo,hi,w=optical_cell(p,n);density=0.
            @inbounds for k=0:1,j=0:1,i=0:1
                f=(i==0 ? 1-w[1] : w[1])*(j==0 ? 1-w[2] : w[2])*(k==0 ? 1-w[3] : w[3])
                density+=solid[i==0 ? lo[1] : hi[1],j==0 ? lo[2] : hi[2],k==0 ? lo[3] : hi[3]]*f
            end
            segment_length=sqrt(sum((a[d]-b[d])^2 for d=1:3));loss=weight*(-expm1(-alpha*density*segment_length));weight-=loss;absorbed+=loss
            @inbounds for k=0:1,j=0:1,i=0:1
                f=(i==0 ? 1-w[1] : w[1])*(j==0 ? 1-w[2] : w[2])*(k==0 ? 1-w[3] : w[3])
                deposition[i==0 ? lo[1] : hi[1],j==0 ? lo[2] : hi[2],k==0 ? lo[3] : hi[3]]+=loss*f
            end
        end
        ray.status==:reached ? (emergent+=weight) : (unresolved+=weight)
        push!(ends,ray.position)
    end
    voxel_volume=(64/n)^3
    deposited=sum(deposition);abs(deposited-absorbed)<1e-8*max(power,1)||error("lost deposition")
    abs(power-absorbed-emergent-unresolved)<1e-8*max(power,1)||error("lost ray energy")
    dose=Float32.(deposition./voxel_volume)
    dose,(incident=power,absorbed=absorbed,emergent=emergent,unresolved=unresolved,residual=power-absorbed-emergent-unresolved,deposited=deposited,voxel_volume=voxel_volume,max_density=maximum(dose),rays=length(launches)),ends
end
@kernel function photodrive!(latent,dose,dt,gain)
    i,j,k=@index(Global,NTuple)
    @inbounds begin
        relax=-expm1(-gain*dose[i,j,k]*dt)
        latent[i,j,k,1]+=relax*(1f0-latent[i,j,k,1])
    end
end
function illuminated(config,out)
    c=TOML.parsefile(config);n=c["resolution"];dt=Float32(c["dt"]);steps=c["steps"];refresh=c["refresh_steps"];mode=c["mode"]
    mode in ("coupled","frozen","dark")||error("mode");ispath(out)&&error("Refusing existing output");mkpath(out);cp(config,joinpath(out,"recipe.toml"))
    mkpath(joinpath(out,"source"));cp(joinpath(REPO,"src"),joinpath(out,"source","src"))
    for f in ("refractive_memory.jl","illuminated_memory.jl");cp(joinpath(@__DIR__,f),joinpath(out,"source",f));end
    w=World(Grid(n,n,n;extent=(64,64,64));backend=MetalBackend(),seed=618)
    p=RecurrentMaterial(w;coupling=.23,feedback=1.4,writing=.8,erosion=.34,contrast=.65)
    dt<=Interior.stability_limit(w,p)||error("dt unstable")
    d,f,l=initial_lenses(n);copyto!(w.fields.deposit,d);copyto!(w.fields.nutrient,f);copyto!(p.state,l)
    mass0=sum(Float64,d)+sum(Float64,f);dose=similar(w.fields.deposit);fill!(dose,0);ledger=Any[];frames=Any[];history=Any[];absorbed_energy=0.;incident_energy=0.;emergent_energy=0.;unresolved_energy=0.
    function capture(step)
        dir=joinpath(out,"frames",lpad(step,6,'0'));mkpath(dir)
        dd=Array(w.fields.deposit);ff=Array(w.fields.nutrient);aa=Array(w.fields.activity)
        all(isfinite,dd)&&all(isfinite,ff)&&minimum(dd)>=0&&minimum(ff)>=0||error("bad state")
        for (name,v) in (("solid",dd),("nutrient",ff),("activity",aa),("dose",Array(dose)));write(joinpath(dir,name*".f32"),v);end
        push!(frames,(step=step,time=w.time,directory=relpath(dir,out)))
        push!(history,(step=step,mass_error=sum(Float64,dd)+sum(Float64,ff)-mass0,solid_mean=mean(dd),solid_range=extrema(dd)))
    end
    capture(0);current=nothing
    elapsed=@elapsed for step=1:steps
        if step==1 || ((step-1)%refresh==0&&mode=="coupled")
            host=Array(w.fields.deposit)
            q,current,ends=absorbed_power(host;grid=c["ray_grid"],ds=c["ray_ds"],alpha=c["alpha"],power=mode=="dark" ? 0. : c["power"],contrast=c["contrast"])
            copyto!(dose,q);push!(ledger,(step=step-1,time=w.time,power=current,endpoints=ends));println("REFRESH ",step-1," ",JSON3.write(current));flush(stdout)
        end
        advance!(w,p,dt)
        photodrive!(w.backend,256)(p.state,dose,dt,Float32(c["gain"]);ndrange=w.grid.size)
        w.time+=Float64(dt)
        absorbed_energy+=current.absorbed*dt;incident_energy+=current.incident*dt;emergent_energy+=current.emergent*dt;unresolved_energy+=current.unresolved*dt
        (step%c["capture_every"]==0||step==steps)&&capture(step)
    end
    Metal.synchronize()
    sourcehash=Dict{String,String}()
    for (dir,_,files) in walkdir(joinpath(out,"source")),f in files
        path=joinpath(dir,f);sourcehash[relpath(path,joinpath(out,"source"))]=bytes2hex(sha256(read(path)))
    end
    report=(id=c["id"],title=c["title"],description=c["description"],family="optics",shape=[n,n,n],extent=[64,64,64],time=w.time,steps=steps,dt=dt,elapsed_seconds=elapsed,initial_material_mass=mass0,frames=frames,history=history,params=c,source_sha256=sourcehash,ledger=ledger,energy=(incident=incident_energy,absorbed=absorbed_energy,emergent=emergent_energy,unresolved=unresolved_energy,residual=incident_energy-absorbed_energy-emergent_energy-unresolved_energy),units="Power is arbitrary energy/time; alpha per optical length per solid; deposition divided by material voxel volume (64/n)^3. Gain has material volume/energy units. Optical width2 corresponds to material width64.",law="Beer–Lambert absorption along actual GRIN paths drives recurrent channel1 toward+1; existing finite-stock recurrent process changes mineral. Frozen optical power density held between explicit refresh clocks. No temperature model or optical scattering.")
    open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w");println("COUPLED_READY ",out)
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    illuminated(ARGS[1],ARGS[2])
end
