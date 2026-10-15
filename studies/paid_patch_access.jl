#!/usr/bin/env julia
#= Paid, stationary patch assay using the native field kernels on one CPU thread.

This is an intervention on a finite funded builder, not an evolving population.
After its paid construction pulse, physiology and position are frozen. Read-only
virtual sensors measure availability; they do not remove food or receive rewards.
=#
using Interior, KernelAbstractions, Serialization, SHA, TOML
using Interior: JSON3
isdefined(Interior,:DeviceCellular) || Base.include(Interior,joinpath(dirname(pathof(Interior)),"device_cellular.jl"))
isdefined(Interior,:DeviceConstructiveEcology) || Base.include(Interior,joinpath(dirname(pathof(Interior)),"device_habitat.jl"))

@kernel function pay_patch!(s,e,arg,q)
    i=@index(Global,Linear)
    a,b,mass,energy,work,heat,build=Interior.device_anabolic_build!(arg,i,s.metabolites[i,1],
        s.metabolites[i,2],s.mass[i],s.energy[i],0f0,q.growth_work,q.potentials[1],q.potentials[2])
    s.metabolites[i,1]=a;s.metabolites[i,2]=b;s.mass[i]=mass;s.energy[i]=energy
    Interior.dc_ledger_add!(s.ledger,s.ledger_low,i,1,work)
    Interior.dc_ledger_add!(s.ledger,s.ledger_low,i,2,heat)
end

function fixture(cfg,raw_diffusion,matrix_on;side=32)
    spacing=.0625;extent=side*spacing;volume=spacing^3
    grid=Grid(side,side,side;extent=(extent,extent,extent))
    w=World(grid,CPU(),NamedTuple(),NamedTuple(),0.,211)
    food=ResourcePool(w;names=[:a,:b,:waste],initial=[0.,0.,0.],potentials=[4.,4.,0.],diffusion=[.006,.006,.001])
    cell=CellularEcology(w,fill(extent/2,3,1),food;capacity=1,hidden=8,seed=211,biomass=1.6,
        energy=1.,metabolites=.12,density=128.,exchange_radius=.18,speed=.06,basal_rate=.008,
        adaptation_strength=1.5,adaptation_time=8.)
    # Match the garden's source concentration rate, with a finite 80-time budget.
    supply=1/64;rate=supply*volume;budget=supply*extent^3*80
    raw=ResourcePool(w;names=[:raw_a,:raw_b],initial=[2volume,2volume],potentials=[4.,4.],
        diffusion=[raw_diffusion,raw_diffusion],source_rates=[rate,rate],source_budgets=[budget,budget])
    habitat=CatalyticHabitat(w,cell,raw;Dict(Symbol(k)=>v for (k,v) in cfg["catalysis"])...)
    mq=Dict(Symbol(k)=>v for (k,v) in cfg["matrix"])
    if !matrix_on;mq[:capacity_per_enzyme]=0.;mq[:transport_barrier]=0.;end
    matrix=SorptiveMatrix(w,habitat;mq...)
    c=Interior.DeviceCellular(w,cell;capacity=1,mechanical_substeps=16,bin_width=.25,energy_type=Float64)
    p=Interior.DeviceConstructiveEcology(w,c,habitat;matrix)
    w,p
end

function weights_at(g,position,R)
    weights=zeros(Float64,prod(Int.(g.dims)))
    for site in eachindex(weights)
        index=Interior.dc_coords(site,g.dims);r2=0.
        for d in 1:3
            dx=abs((index[d]-.5)*g.spacing[d]-position[d]);dx=min(dx,g.extent[d]-dx);r2+=dx^2
        end
        r2<R^2 && (weights[site]=(1-r2/R^2)^2)
    end
    weights./=sum(weights);weights
end

function row(w,p,near,far,exposure)
    c=p.cells;h=p.habitat;g=c.geometry;e=c.environment
    audit=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    concentrations=Float64.(e.fields[:,1:2])./g.volume
    center=g.extent./2
    second_moment=0.;free_mass=0.
    for site in axes(e.fields,1)
        index=Interior.dc_coords(site,g.dims);r2=0.
        for d in 1:3
            dx=abs((index[d]-.5)*g.spacing[d]-center[d]);dx=min(dx,g.extent[d]-dx);r2+=dx^2
        end
        mass=Float64(e.fields[site,1])+Float64(e.fields[site,2]);second_moment+=mass*r2;free_mass+=mass
    end
    (time=w.time,audit,free_mass,near_free=vec(near'*concentrations),far_free=vec(far'*concentrations),
        near_raw=[sum(near.*h.feedstock.arrays.fields[:,k])/g.volume for k in 1:2],
        near_bound=[sum(near.*h.arrays.bound[:,k])/g.volume for k in 1:2],
        near_matrix_permeability=sum(near.*h.arrays.permeability),near_combined_permeability=sum(near.*e.permeability),
        free_rms_distance=free_mass>0 ? sqrt(second_moment/free_mass) : 0.,exposure=copy(exposure))
end

function frame(path,w,p)
    c=p.cells;h=p.habitat;dims=Int.(c.geometry.dims)
    serialize(path,(time=w.time,shape=dims,extent=c.geometry.extent,
        free=reshape(copy(c.environment.fields[:,1:2]),dims...,2),raw=reshape(copy(h.feedstock.arrays.fields[:,1:2]),dims...,2),
        enzyme=reshape(copy(reinterpret(Float32,h.arrays.enzyme_bits)),dims...,2),bound=reshape(copy(h.arrays.bound),dims...,2),
        permeability=reshape(copy(c.environment.permeability),dims...),matrix_permeability=reshape(copy(h.arrays.permeability),dims...),
        volume=c.geometry.volume,position=copy(c.state.position),radius=copy(c.state.radii)))
end

function assay(profile,output;horizon=60.,dt=.02,side=32)
    ispath(output)&&error("fresh output required");horizon>0&&dt>0&&side>=16||error("invalid bounds")
    mkpath(output);cfg=TOML.parsefile(profile);results=Dict{String,Any}()
    report=Dict("kind"=>"paid stationary patch accessibility intervention; not an evolution/fitness assay",
        "source"=>Interior.DEVICE_CELLULAR_SOURCE_PIN,"habitat_source"=>Interior.DEVICE_HABITAT_SOURCE_PIN,
        "driver_sha256"=>bytes2hex(sha256(read(@__FILE__))),"profile_sha256"=>bytes2hex(sha256(read(profile))),
        "horizon"=>horizon,"dt"=>dt,"shape"=>[side,side,side],"spacing"=>.0625,
        "dose_per_catalyst"=>.006,"sensor_offset"=>.4,"frozen"=>"builder position, biomass and post-pulse intracellular state; no upkeep, motion, uptake, division or further construction",
        "sensors"=>"read-only normalized native-radius kernels; exposure is concentration*time, not captured food or fitness",
        "arms"=>results)
    for (name,D,matrix) in (("immobile_full",0.,true),("mobile_full",.006,true),("immobile_catalytic",0.,false),("mobile_catalytic",.006,false))
        dir=joinpath(output,name);mkdir(dir);w,p=fixture(cfg,D,matrix;side);c=p.cells;h=p.habitat;g=c.geometry
        Float32(dt)<=Interior.stability_limit(w,p)||error("step exceeds diffusion limit")
        initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
        h.arrays.requests[1,:].=.006f0
        pay_patch!(CPU(),1)(c.state,c.environment,h.argument,c.parameters;ndrange=1)
        # Rebuild physical porosity after material leaves the frozen builder.
        Interior.device_cellular_launch!(c,Interior.dc_clear_occupancy!,c.environment;ndrange=prod(g.dims))
        Interior.device_cellular_launch!(c,Interior.dc_clear_heads!,c.environment;ndrange=length(c.environment.heads))
        Interior.device_cellular_launch!(c,Interior.dc_geometry!,c.state,c.environment,c.parameters,g;ndrange=1)
        Interior.device_cellular_launch!(c,Interior.dc_porosity_scatter!,c.state,c.environment,g;ndrange=1)
        Interior.device_cellular_launch!(c,Interior.dc_porosity_commit!,c.environment,c.parameters;ndrange=prod(g.dims))
        body=copy(c.environment.permeability)
        Interior.device_cellular_launch!(c,Interior.dh_permeability!,h.arrays,c.environment,g,h.parameters;ndrange=prod(g.dims))
        near=weights_at(g,g.extent./2,.18);far=weights_at(g,g.extent./2 .+(.4,0.,0.),.18)
        near_sites=findall(>(0),near);far_sites=findall(>(0),far)
        exposure=zeros(2,2);rows=Any[];push!(rows,row(w,p,near,far,exposure));frame(joinpath(dir,"frame-00000.jls"),w,p)
        started=time();delta=Float32(dt);steps=round(Int,horizon/Float64(delta));capture=max(1,round(Int,1/Float64(delta)))
        for step in 1:steps
            Interior.device_anabolic_begin!(h,w,c,delta)
            Interior.device_cellular_launch!(c,Interior.dc_diffuse!,c.environment,c.parameters,g,delta;ndrange=(prod(g.dims),3))
            Interior.device_cellular_launch!(c,Interior.dc_commit_resources!,c.environment;ndrange=(prod(g.dims),3))
            copyto!(c.environment.permeability,body)
            Interior.device_anabolic_finish!(h,w,c,delta)
            w.time=step*Float64(delta)
            for k in 1:2
                exposure[1,k]+=Float64(delta)*sum(i->near[i]*Float64(c.environment.fields[i,k]),near_sites)/g.volume
                exposure[2,k]+=Float64(delta)*sum(i->far[i]*Float64(c.environment.fields[i,k]),far_sites)/g.volume
            end
            if step%capture==0||step==steps
                Interior.device_cellular_status!(c);push!(rows,row(w,p,near,far,exposure))
            end
            if step%round(Int,10/Float64(delta))==0||step==steps
                frame(joinpath(dir,"frame-"*lpad(step,5,'0')*".jls"),w,p)
            end
        end
        last=rows[end].audit
        mass_error=last.material-initial.material-last.supplied+initial.supplied
        potential_error=last.potential-initial.potential+last.work-initial.work+last.heat-initial.heat-last.supplied_potential+initial.supplied_potential
        @assert last.finite && all(>=(0),c.environment.fields) && all(>=(0),h.arrays.bound)
        @assert abs(mass_error)<.002 && abs(potential_error)<.008
        result=(raw_diffusion=D,matrix,initial,rows,material_error=mass_error,potential_error,seconds=time()-started)
        results[name]=result;open(io->JSON3.pretty(io,result),joinpath(dir,"result.json"),"w")
        open(io->JSON3.pretty(io,report),joinpath(output,"report.json"),"w")
        println((name,time=w.time,converted=last.solubilized,bound=last.bound,near_exposure=vec(exposure[1,:]),far_exposure=vec(exposure[2,:]),mass_error,potential_error));flush(stdout)
    end
    report
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)>=2||error("PROFILE FRESH_OUTPUT [HORIZON=60] [SIDE=32] required")
    assay(ARGS[1],ARGS[2];horizon=length(ARGS)>2 ? parse(Float64,ARGS[3]) : 60.,side=length(ARGS)>3 ? parse(Int,ARGS[4]) : 32)
end
