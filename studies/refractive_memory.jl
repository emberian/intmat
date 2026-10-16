# A recurrent, finite-stock material observed by Interior's existing GRIN rays.
# The material never receives an optical energy source; tracing is read-only.
using Interior, Metal, SHA, TOML, Statistics
using Interior: JSON3, RecurrentMaterial
const REPO=dirname(@__DIR__)
function initial_lenses(n)
    solid=zeros(Float32,n,n,n); food=similar(solid); latent=zeros(Float32,n,n,n,8)
    for k=1:n,j=1:n,i=1:n
        x,y,z=ntuple(d->2f0*((i,j,k)[d]-.5f0)/n-1,3)
        stock=0f0
        # Three soft helical lenses are an authored initial material geometry.
        for strand=0:2
            phase=2f0*pi*strand/3+2.6f0*z
            cx=.27f0*cos(phase)+.08f0*z;cy=.27f0*sin(phase)
            envelope=exp(-((z/.74f0)^8))
            stock+=.80f0*exp(-((x-cx)^2+(y-cy)^2)/.032f0)*envelope
        end
        stock=min(stock,1.15f0)
        phase=8x-6y+5z+1.1f0*sin(6y+4z)+.8f0*cos(7x-3z)
        fraction=.48f0+.15f0*sin(phase)
        solid[i,j,k]=stock*fraction;food[i,j,k]=stock*(1-fraction)
        for pair=1:4
            a=phase+pair*.8f0+(.5f0+.2f0*pair)*sin(4z-5y+pair)
            latent[i,j,k,2pair-1]=.65f0*sin(a);latent[i,j,k,2pair]=.65f0*cos(a)
        end
    end
    solid,food,latent
end
function produce(config,out)
    c=TOML.parsefile(config);n=c["resolution"];steps=c["steps"];dt=Float32(c["dt"])
    ispath(out)&&error("Refusing existing run");mkpath(out);cp(config,joinpath(out,"recipe.toml"))
    mkpath(joinpath(out,"source"));cp(joinpath(REPO,"src"),joinpath(out,"source","src"));cp(@__FILE__,joinpath(out,"source","refractive_memory.jl"))
    sourcehash=Dict{String,String}()
    for (dir,_,files) in walkdir(joinpath(out,"source")),f in files
        path=joinpath(dir,f);sourcehash[relpath(path,joinpath(out,"source"))]=bytes2hex(sha256(read(path)))
    end
    w=World(Grid(n,n,n;extent=(64,64,64));backend=MetalBackend(),seed=618)
    p=RecurrentMaterial(w;coupling=.23,feedback=1.4,writing=.8,erosion=.34,contrast=.65)
    dt<=Interior.stability_limit(w,p)||error("Unstable dt")
    d,f,l=initial_lenses(n);copyto!(w.fields.deposit,d);copyto!(w.fields.nutrient,f);copyto!(p.state,l)
    initialmass=sum(Float64,d)+sum(Float64,f);d=nothing;f=nothing;l=nothing;GC.gc()
    history=Any[];frames=Any[]
    function capture(step)
        dir=joinpath(out,"frames",lpad(step,6,'0'));mkpath(dir)
        d=Array(w.fields.deposit);f=Array(w.fields.nutrient);act=Array(w.fields.activity)
        all(isfinite,d)&&all(isfinite,f)&&all(isfinite,act)&&minimum(f)>=0||error("invalid fields")
        for (name,v) in (("solid",d),("nutrient",f),("activity",act));write(joinpath(dir,name*".f32"),v);end
        metric=(step=step,time=w.time,mass_error=sum(Float64,d)+sum(Float64,f)-initialmass,solid_range=extrema(d),solid_mean=mean(d),nutrient_min=minimum(f))
        push!(history,metric);push!(frames,(step=step,time=w.time,directory=relpath(dir,out)));println(JSON3.write(metric));flush(stdout)
    end
    capture(0)
    elapsed=@elapsed for step=1:steps
        advance!(w,p,dt);w.time+=Float64(dt)
        (step%c["capture_every"]==0||step==steps)&&capture(step)
    end
    Metal.synchronize()
    report=(id=c["id"],title=c["title"],description=c["description"],family="optics",shape=w.grid.size,extent=w.grid.extent,time=w.time,steps=steps,dt=dt,elapsed_seconds=elapsed,initial_material_mass=initialmass,frames=frames,history=history,params=c,source_sha256=sourcehash,initialization="Three soft helical mineral/nutrient lenses; phase-textured eight-channel recurrent initialization. Initial geometry is authored, not emergent from noise.",simulation="Existing RecurrentMaterial conservative writing/erosion and diffusion, fixed designed recurrence. Optical observation does not feed back on material.")
    open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w")
end
function loadfield(root,frame,name,n)
    path=joinpath(root,frame["directory"],name*".f32");filesize(path)==n^3*4||error("field size")
    reshape(reinterpret(Float32,read(path)),n,n,n)|>Array
end
function indexfield(solid,contrast,wavelength)
    n=size(solid,1);index=similar(solid)
    # Weak explicitly invented Cauchy-like dispersion; the outside is ambient1.
    gain=contrast*(.94+.06*(550/wavelength)^2)
    for k=1:n,j=1:n,i=1:n
        x,y,z=ntuple(d->2f0*((i,j,k)[d]-.5f0)/n-1,3)
        smoothfade(v)=begin t=clamp((.94f0-abs(v))/.12f0,0f0,1f0);t*t*(3f0-2f0*t) end
        taper=prod(smoothfade(v) for v in (x,y,z))
        index[i,j,k]=1+gain*solid[i,j,k]*taper
    end
    index
end
function observe(root,step;contrast=.65,ds=.003,visual_grid=17,screen_grid=128)
    run=JSON3.read(read(joinpath(root,"run.json"),String),Dict{String,Any});n=run["shape"][1]
    frame=only(filter(f->f["step"]==step,run["frames"]));tag=get(ENV,"OPTICS_TAG","");occursin(r"^[a-z0-9-]*$",tag)||error("invalid tag");out=joinpath(root,"optics",lpad(step,6,'0')*(isempty(tag) ? "" : "-"*tag));ispath(out)&&error("Refusing existing optical observation");mkpath(out)
    solid=loadfield(root,frame,"solid",n);rays=Any[];screens=Any[];checks=Any[]
    wavelengths=(460.,550.,630.);launch_span=1.15;beam_mode=get(ENV,"BEAM_MODE","parallel");beam_mode in ("parallel","helical")||error("beam mode")
    target_z=beam_mode=="parallel" ? 1.45 : 1.05;screen_span=2.4;screen_pixels=256
    function launches(grid)
        result=NamedTuple[]
        if beam_mode=="parallel"
            for j=1:grid,i=1:grid
                x=((i-.5)/grid-.5)*launch_span;y=((j-.5)/grid-.5)*launch_span
                push!(result,(origin=(x,y,-1.08),direction=(0.,0.,1.)))
            end
        else
            # Three small apertures are aligned to the authored initial lens
            # tangents. Later fields do not move or retarget these launchers.
            z=-.65
            for strand=0:2,j=1:grid,i=1:grid
                u=2*((i-.5)/grid-.5);v=2*((j-.5)/grid-.5);u*u+v*v>1&&continue
                phase=2pi*strand/3+2.6z;cx=.27cos(phase)+.08z;cy=.27sin(phase)
                push!(result,(origin=(cx+.095u,cy+.095v,z),direction=(-.27*2.6sin(phase)+.08,.27*2.6cos(phase),1.)))
            end
        end
        result
    end
    display_launches=launches(visual_grid);dense_launches=launches(screen_grid)
    straight(launch)=ntuple(d->launch.origin[d]+(target_z-launch.origin[3])*launch.direction[d]/launch.direction[3],3)
    total_elapsed=@elapsed for wavelength in wavelengths
        index=indexfield(solid,contrast,wavelength);field=GradedIndex(index);write(joinpath(out,"index-$(Int(wavelength)).f32"),index)
        trace(launch,h;save=true)=trace_grin(field,launch.origin,launch.direction;ds=h,target_z=target_z,max_length=7.,save_path=save,_validated=true)
        # All launched rays are kept, including bounded unresolved paths.
        for launch in display_launches
            ray=trace(launch,ds)
            displaypath=ray.path[1:2:end];displaypath[end]!=ray.path[end]&&push!(displaypath,ray.path[end])
            push!(rays,(wavelength=wavelength,launch=launch.origin,launch_direction=launch.direction,constant_index_position=straight(launch),status=ray.status,position=ray.position,direction=ray.direction,optical_length=ray.optical_length,arc_length=ray.arc_length,path=displaypath))
        end
        # Dense equal-weight probe rays for a geometric flux-density screen.
        hist=zeros(Int,screen_pixels,screen_pixels);unresolved=0;outside=0;deflections=Float64[]
        hits=Vector{Any}(undef,length(dense_launches))
        Threads.@threads for q in eachindex(dense_launches)
            launch=dense_launches[q];rr=trace(launch,ds;save=false);hits[q]=(reference=straight(launch),ray=rr)
        end
        for hit in hits
            rr=hit.ray
            if rr.status!=:reached;unresolved+=1;continue;end
            push!(deflections,hypot(rr.position[1]-hit.reference[1],rr.position[2]-hit.reference[2]))
            ix=floor(Int,(rr.position[1]/screen_span+.5)*screen_pixels)+1;iy=floor(Int,(rr.position[2]/screen_span+.5)*screen_pixels)+1
            if 1<=ix<=screen_pixels&&1<=iy<=screen_pixels;hist[ix,iy]+=1;else;outside+=1;end
        end
        coarse=Float64[];fine=Float64[];fixture=Any[]
        for q in unique(round.(Int,range(1,length(display_launches);length=16)))
            launch=display_launches[q]
            a=trace(launch,ds*2;save=false);b=trace(launch,ds;save=false);c=trace(launch,ds/2;save=false)
            err(u,v)=sqrt(sum((u.position[d]-v.position[d])^2 for d=1:3))
            push!(coarse,err(a,b));push!(fine,err(b,c));push!(fixture,(launch=launch.origin,statuses=(a.status,b.status,c.status),coarse_difference=coarse[end],fine_difference=fine[end]))
        end
        push!(checks,(wavelength=wavelength,index_range=extrema(index),total=length(dense_launches),reached=length(deflections),unresolved=unresolved,outside=outside,binned=sum(hist),mean_transverse_deflection=mean(deflections),maximum_transverse_deflection=maximum(deflections),max_coarse_difference=maximum(coarse),max_fine_difference=maximum(fine),fixture=fixture))
        sum(hist)+outside+unresolved==length(dense_launches)||error("lost ray weight")
        push!(screens,(wavelength=wavelength,counts=vec(hist),total=length(dense_launches)));println(JSON3.write(checks[end]));flush(stdout)
    end
    # Constant-index control uses exactly the same established integrator.
    uniform=GradedIndex(ones(Float32,4,4,4));constant_errors=Float64[]
    for launch in display_launches[1:max(1,length(display_launches)÷16):end]
        rr=trace_grin(uniform,launch.origin,launch.direction;ds=ds,target_z=target_z)
        push!(constant_errors,maximum(abs.(rr.position.-straight(launch))))
    end
    maximum(constant_errors)<1e-10||error("straight-ray control")
    observation=(shape=run["shape"],step=step,time=frame["time"],source_run_sha256=bytes2hex(sha256(read(joinpath(root,"run.json")))),source_solid_sha256=bytes2hex(sha256(reinterpret(UInt8,vec(solid)))),tracer_sha256=bytes2hex(sha256(read(joinpath(REPO,"src","grin.jl")))),observer_source_sha256=bytes2hex(sha256(read(@__FILE__))),contrast=contrast,index_law="n(lambda)=1+contrast*(.94+.06*(550nm/lambda)^2)*solid*C1 smoothstep boundary taper; invented Cauchy-like constitutive law, no calibration",coordinate_map="Material extent64^3 mapped affinely to optical box[-1,1]^3",ds=ds,beam_mode=beam_mode,constant_index_max_error=maximum(constant_errors),target_z=target_z,launch_span=launch_span,screen_span=screen_span,screen_pixels=screen_pixels,visual_grid=visual_grid,screen_grid=screen_grid,elapsed_seconds=total_elapsed,rays=rays,screens=screens,checks=checks,limitations="Frozen-field geometrical ray initial-value observer, no interference/polarization/scattering/Fresnel or optical feedback. Luminous path display is an observation glyph; histogram is equal-weight geometric flux, not a wave caustic solver.")
    open(io->JSON3.pretty(io,observation),joinpath(out,"observation.json"),"w")
    cp(joinpath(REPO,"src","grin.jl"),joinpath(out,"grin-source.jl"));cp(@__FILE__,joinpath(out,"observer-source.jl"));println("OPTICS_READY ",out)
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    ARGS[1]=="produce" ? produce(ARGS[2],ARGS[3]) : foreach(step->observe(ARGS[2],step;contrast=parse(Float64,get(ENV,"INDEX_CONTRAST",".65")),ds=parse(Float64,get(ENV,"RAY_DS",".003")),visual_grid=parse(Int,get(ENV,"RAY_GRID","17")),screen_grid=parse(Int,get(ENV,"SCREEN_GRID","128"))),parse.(Int,split(ARGS[3],',')))
end
