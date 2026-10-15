# Included by LivingTissueStudy; simulation laws live in living_tissue.jl.
function tissue_backend(kind)
    if kind=="amd"
        Core.eval(@__MODULE__,:(using AMDGPU))
        # Import can introduce a binding newer than this function's world.
        return Core.eval(@__MODULE__,:(AMDGPU.ROCBackend()))
    elseif kind=="metal"
        Core.eval(@__MODULE__,:(using Metal))
        return Core.eval(@__MODULE__,:(Metal.MetalBackend()))
    elseif kind=="cpu"
        return CPU()
    end
    throw(ArgumentError("backend must be cpu, amd or metal"))
end

function write_tissue_frame(directory,state,grid,step,time;native=true)
    mkdir(directory)
    if native
        for (channel,name) in enumerate(CHANNELS)
            open(joinpath(directory,"$name.f32"),"w") do io
                write(io,Array(view(state,:,:,:,channel)))
            end
        end
    end
    open(joinpath(directory,"volume.json"),"w") do io
        JSON3.pretty(io,(shape=grid.size,extent=grid.extent,spacing=grid.spacing,
            dtype="float32",byte_order=ENDIAN_BOM==0x04030201 ? "little" : "big",
            array_order="F",centering="cell",boundary="periodic",step,time,
            fields=native ? string.(CHANNELS) : (),
            optical_mapping="Biomass and matrix are actual simulation fields. Optical properties require a separate declared mapping."))
    end
    # A compact fixed-scale middle section; no per-frame contrast adjustment.
    indices=(round.(Int,range(1,grid.size[1];length=min(160,grid.size[1]))),
        round.(Int,range(1,grid.size[2];length=min(160,grid.size[2]))))
    replay=zeros(UInt8,length(indices[1]),length(indices[2]),3)
    for (c,(channel,low,high)) in enumerate(((2,0f0,1.2f0),(3,0f0,.6f0),(5,-1f0,1f0)))
        for (j,jj) in enumerate(indices[2]),(i,ii) in enumerate(indices[1])
            replay[i,j,c]=round(UInt8,255clamp((state[ii,jj,cld(grid.size[3],2),channel]-low)/(high-low),0,1))
        end
    end
    open(io->write(io,replay),joinpath(directory,"section.u8"),"w")
end

function run_living_tissue(recipe,output;backend_kind="cpu")
    backend=tissue_backend(backend_kind)
    # Keep late-loaded backend constructors outside the compiled runner age.
    Base.invokelatest(run_living_tissue_backend,recipe,output,backend,backend_kind)
end
function run_living_tissue_backend(recipe,output,backend,backend_kind)
    c=TOML.parsefile(recipe);dims=Tuple(Int.(c["dimensions"]));extent=Tuple(Float32.(c["extent"]))
    length(dims)==length(extent)==3 || throw(ArgumentError("three dimensions and extents required"))
    dt=Float32(c["dt"]);steps=Int(c["steps"]);every=Int(c["capture_every"])
    steps>0 && every>0 || throw(ArgumentError("positive steps/capture cadence required"))
    grid=Grid(dims...;extent)
    model=NamedTuple(Symbol(k)=>v for (k,v) in get(c,"model",Dict()))
    q=tissue_parameters(;model...)
    resident=big(prod(BigInt.(dims)))*56
    resident<=get(c,"resident_budget_mb",1024)*1024^2 || throw(ArgumentError("tissue resident memory budget exceeded"))
    frames=union(Set(0:every:steps),Set((steps,)))
    native=get(c,"native_fields",true)
    disk=length(frames)*big(prod(BigInt.(dims)))*(native ? 28 : 0)+3big(prod(BigInt.(dims)))*56+10_000_000
    disk<=get(c,"capture_budget_mb",1024)*1024^2 || throw(ArgumentError("tissue capture budget exceeded"))
    ispath(output) && throw(ArgumentError("use a fresh output directory"))
    # Preflight an empty-buffer process before claiming the output destination.
    pattern=get(c,"seed_pattern","ribbons")
    w=tissue_world(grid;backend,seed=get(c,"seed",911));p=LivingTissue(w;parameters=q,pattern)
    Interior.check_dt(w,p,dt)
    mkpath(dirname(abspath(output)));mkdir(output);mkdir(joinpath(output,"frames"))
    cp(recipe,joinpath(output,"recipe.toml"))
    write(joinpath(output,"study-source.jl"),MODEL_SOURCE)
    runner_source=read(@__FILE__,String);write(joinpath(output,"living_tissue_runner.jl"),runner_source)
    sources=Dict("studies/living_tissue.jl"=>MODEL_SHA,
        "studies/living_tissue_runner.jl"=>bytes2hex(sha256(runner_source)))
    for name in ("core.jl","recurrent_material.jl")
        source=SHARED_SOURCES[name]
        write(joinpath(output,name),source);sources["src/$name"]=SHARED_SHA[name]
    end
    metadata=(id=c["id"],title=c["title"],description=c["description"],dimensions=dims,extent,
        backend=backend_kind,dt,steps,parameters=q,resident_array_bytes=Int(resident),
        source_sha256=sources,channels=string.(CHANNELS),
        section=(shape=(min(160,dims[1]),min(160,dims[2]),3),dtype="uint8",array_order="F",
            channels=("biomass","matrix","excitation"),ranges=((0,1.2),(0,.6),(-1,1))),
        material_convention="Equal-volume concentrations; N+B+M+W conserved. N/B/M potential 1, spent W potential 0. Metabolic N→W funds neural drive. Finite initial supply; no input or output flux.",
        initialization="Explicit $(pattern) seed strands and smooth precursor heterogeneity; anatomy is not claimed to arise from an unstructured seed.",
        simulation="Eulerian continuum biomass; no individual mechanics, inheritance, training or evolutionary selection.")
    open(io->JSON3.pretty(io,metadata),joinpath(output,"configuration.json"),"w")
    initial=Array(p.state);mass=sum(Float64,view(initial,:,:,:,1:4));potential=sum(Float64,view(initial,:,:,:,1:3))
    center=argmax(view(initial,:,:,:,2));probe=Tuple(center)
    initial=nothing # Do not retain a full-size host replica throughout the run.
    history=Any[];probes=Any[]
    function capture(step)
        KernelAbstractions.synchronize(backend);state=Array(p.state)
        statistics=tissue_statistics(w,p;state)
        residual=statistics.material-mass;potential_residual=statistics.potential+statistics.spent-potential
        statistics.finite && statistics.physical_min>=0 && statistics.neural_maxabs<=1.00001 || error("invalid tissue state")
        abs(residual)<=max(1e-5,5e-5mass) || error("tissue material balance failed: $residual")
        push!(history,merge((step=step,material_residual=residual,potential_residual=potential_residual),statistics))
        dir=joinpath("frames",lpad(step,7,'0'));write_tissue_frame(joinpath(output,dir),state,grid,step,w.time;native)
        open(io->println(io,JSON3.write(last(history))),joinpath(output,"metrics.jsonl"),"a")
        println(JSON3.write(last(history)));flush(stdout)
    end
    capture(0);record_state(w,p,joinpath(output,"initial-state");step=0)
    # Warm the real dispatch on a tiny separate world, preserving initial state.
    warm=tissue_world(Grid(3,3,1;extent=(3,3,1));backend);warm_p=LivingTissue(warm;parameters=q)
    step!(warm,warm_p,.001f0);KernelAbstractions.synchronize(backend)
    elapsed=@elapsed for k in 1:steps
        step!(w,p,dt)
        k%64==0 && KernelAbstractions.synchronize(backend)
        if k%max(1,every÷8)==0
            # Probe readback happens only at explicit observation boundaries.
            KernelAbstractions.synchronize(backend)
            probe_state=Array(view(p.state,probe[1]:probe[1],probe[2]:probe[2],probe[3]:probe[3],:))
            push!(probes,(step=k,time=w.time,values=vec(probe_state)))
        end
        k in frames && capture(k)
    end
    KernelAbstractions.synchronize(backend)
    record_state(w,p,joinpath(output,"final-state");step=steps)
    tissue_checkpoint(w,p,joinpath(output,"final.jls"))
    report=merge(metadata,(elapsed_seconds=elapsed,voxel_steps_per_second=prod(dims)*steps/elapsed,
        frame_steps=sort!(collect(frames)),history,probe_site=probe,probes))
    open(io->JSON3.pretty(io,report),joinpath(output,"run.json"),"w")
    println("READY ",output);output
end
