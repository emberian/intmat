using Interior, JSON3, TOML, SHA, Statistics, LinearAlgebra, FFTW
BLAS.set_num_threads(1); FFTW.set_num_threads(2)
const ROOT=dirname(@__DIR__)
function initial(kind,x,y,z)
    # All arguments are periodic physical coordinates on [0,64)^3.
    a,b,c=2pi*x/64,2pi*y/64,2pi*z/64
    g=sin(3a)*cos(3b)+sin(3b)*cos(3c)+sin(3c)*cos(3a)
    if kind=="spinodal"
        return .14*(sin(5a)*cos(4b)+sin(4b)*cos(5c)+sin(5c)*cos(4a))+.035sin(7a+3b-5c)
    elseif kind=="pearls"
        return tanh(2.3*(cos(3a)+cos(3b)+cos(3c)-1.8))
    elseif kind=="gyroid"
        return tanh(1.5*g)
    elseif kind=="membrane_open"
        return tanh(2.5min(.58-abs(sin(2c+.35sin(2a)*sin(2b))),cos(3a)+cos(3b)+.6))
    elseif kind=="membrane"
        return tanh(2.0*(.42-abs(sin(2c+.35sin(2a)*sin(2b)))))*(.5+.5tanh(2.5*(cos(3a)+cos(3b)+.9)))-.55
    elseif kind=="strata"
        return tanh(2.0*(sin(4c+.85sin(2a)+.6cos(3b))+.24g))
    elseif kind=="shells"
        return tanh(2.1*(.28-abs(cos(2a)+cos(2b)+cos(2c))))
    else
        return tanh(1.6*(g+.48cos(2a-3b)+.35sin(4c)-.35))
    end
end
function bake(recipe)
    cfg=TOML.parsefile(recipe); slug=cfg["slug"]; out=joinpath(ROOT,"outputs","collection-20261015",slug)
    isfile(joinpath(out,"run.json")) && return
    mkpath(out); cp(recipe,joinpath(out,"recipe.toml");force=true); cp(@__FILE__,joinpath(out,"study-source.jl");force=true)
    traces=Any[]
    tend=get(cfg,"t_end",6.0); finaldt=get(cfg,"dt",.1)
    for (label,n,dt,steps) in (("pilot",48,finaldt,round(Int,tend/finaldt)),("final",96,finaldt,round(Int,tend/finaldt)))
        w=World(Grid(n,n,n;extent=(64,64,64));seed=cfg["seed"])
        p=ImplicitPhaseSeparation(w;noise=0,kappa=cfg["kappa"],mobility=get(cfg,"mobility",1.0),mineral_pinning=0,nutrient_diffusion=0,interface_activity=.65,activity_relaxation=.3,stabilization=3)
        for I in CartesianIndices(p.phi)
            x,y,z=Tuple(I).-1
            p.phi[I]=initial(cfg["initial"],64x/n,64y/n,64z/n)
        end
        Interior.phase_observe!(w,p)
        phi0=copy(p.phi); m0=sum(Float64,p.phi); e0=phase_energy(w,p)
        function temporal(s)
            label=="final" || return
            mkpath(joinpath(out,"temporal"))
            open(joinpath(out,"temporal","phi-step-$s.f32"),"w") do io; write(io,p.phi); end
            open(joinpath(out,"temporal","slice-step-$s.ppm"),"w") do io
                write(io,"P6\n$n $n\n255\n")
                for j in 1:n,i in 1:n
                    u=(1+tanh(2p.phi[i,j,n÷2]))/2
                    write(io,UInt8(round(255u)),UInt8(round(255sqrt(u))),UInt8(round(255(1-u))))
                end
            end
        end
        temporal(0)
        label=="final" && record_state(w,(p,),joinpath(out,"state-initial");step=0)
        history=Any[(step=0,time=0.,energy=e0,mass=m0,range=extrema(p.phi))]
        elapsed=@elapsed for s in 1:steps
            if cfg["active"]
                for I in CartesianIndices(p.phi)
                    x,y,z=(Tuple(I).-1).*(2pi/n)
                    p.chemical_forcing[I]=.42sin(3x+.55sin(w.time*.35))*cos(3y)+.28sin(3z-.4w.time)
                end
            end
            step!(w,(p,),dt)
            s in (steps÷4,steps÷2,3steps÷4,steps) && temporal(s)
            if s%20==0 || s==steps
                all(isfinite,p.phi) || error("nonfinite phase")
                push!(history,(step=s,time=w.time,energy=phase_energy(w,p),mass=sum(Float64,p.phi),range=extrema(p.phi),rms_change=sqrt(mean((p.phi.-phi0).^2)),positive_fraction=mean(p.phi.>0),sign_change_fraction=mean((p.phi.>0).!=(phi0.>0))))
            end
        end
        push!(traces,(label=label,resolution=n,extent=[64,64,64],spacing=64/n,dt=dt,steps=steps,elapsed_seconds=elapsed,history=history,mean_mass_drift=(sum(Float64,p.phi)-m0)/n^3))
        if label=="final"
            record_state(w,(p,),joinpath(out,"state-final");step=steps)
            f=phase_optical_fields(w,p)
            mask=zeros(Float32,n,n,n)
            for I in CartesianIndices(mask)
                xyz=(Tuple(I).-((n+1)/2,))./(n/2)
                mask[I]=sum(abs2,xyz)<.88^2 ? 1 : 0
            end
            fields=merge(f,(activity=.06f0.*f.activity,mask=mask,phi=copy(p.phi),interface=copy(p.interface),chemical_potential=copy(p.chemical_potential),chemical_forcing=copy(p.chemical_forcing),nutrient=copy(w.fields.nutrient),material_activity=copy(w.fields.activity),deposit=copy(w.fields.deposit)))
            export_volume(fields,joinpath(out,"volume"))
        end
        println(slug," ",label," done ",elapsed); flush(stdout)
    end
    hashes=Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in ("src/phase.jl","src/implicit_phase.jl","studies/collection_phase.jl"))
    open(joinpath(out,"run.json"),"w") do io
        JSON3.pretty(io,(id=slug,title=cfg["title"],description=cfg["description"],solver="CPU semi-implicit conserved Cahn-Hilliard; constant mobility; stabilization 3",initialization="Explicit periodic analytic field in studies/collection_phase.jl, initial kind "*cfg["initial"],forcing=cfg["active"] ? "Prescribed time-varying chemical potential, amplitude .42 and .28; conservative Laplacian coupling. No claim of spontaneous coral growth." : "none",pilot="48 cubed and 96 cubed have identical extent64; spacing changes4/3 to2/3; identical physical endtime and dt; pilot checks finiteness before highres bake",phase_parameters=(bulk=1,kappa=cfg["kappa"],mobility=get(cfg,"mobility",1.),stabilization=3,dt=finaldt,t_end=tend),temporal="Quarter-time raw phase fields plus central xy PPM slices in temporal/; full initial and final record_state; initial analytic topology is explicitly prescribed",optical_mapping="structure=(1+tanh(2phi))/2; activity=.06*(.8 interface+.2 material_activity); spherical viewing crop only",traces=traces,source_sha256=hashes))
    end
    job=Dict("id"=>slug,"family"=>"phase","title"=>cfg["title"],"description"=>cfg["description"],"volume"=>"volume","renderer"=>"volume","args"=>["--inclusion-surface","structure","--inclusion-isovalue","0.5","--polished-carrier","--palette",cfg["palette"],"--density-scale","0.18","--emission","0.3","--absorption","0.4","--scattering","0.08","--orbit","0.65","--elevation","0.32"],"size"=>1600,"samples"=>96)
    open(joinpath(out,"render-job.json"),"w") do io; JSON3.pretty(io,job); end
end
for recipe in (isempty(ARGS) ? sort(filter(p->endswith(p,".toml"),readdir(joinpath(ROOT,"experiments","collection","phase");join=true))) : ARGS)
    bake(abspath(recipe))
end
