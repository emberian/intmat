# Study-only exact fusion; production source is read, never overwritten.
using Interior
using SHA, Statistics, Random, TOML
using Interior: JSON3, KernelAbstractions
include(joinpath(@__DIR__,"..","scripts","native_constructed_garden.jl"))

function feedback_section(source, firstmarker, nextmarker)
    a=findfirst(firstmarker,source);a===nothing && error("missing $firstmarker")
    b=findnext(nextmarker,source,last(a)+1);b===nothing && error("missing $nextmarker")
    source[first(a):prevind(source,first(b))]
end

function install_feedback_factorization!()
    kernelpath=joinpath(pkgdir(Interior),"src","device_cellular_kernels.jl")
    hostpath=joinpath(pkgdir(Interior),"src","device_cellular.jl")
    sense=feedback_section(read(kernelpath,String),"@kernel function dc_sense!","@kernel function dc_controller!")
    sense=replace(sense,"dc_sense!(s,e,q,g)"=>"ff_fused_sense!(s,e,q,g,feed,channels,scales)")
    sense=replace(sense,"        concentration1=0f0;concentration2=0f0"=>"        raw1=0f0;raw2=0f0;bound1=0f0;bound2=0f0;enzyme1=0f0;enzyme2=0f0\n        concentration1=0f0;concentration2=0f0")
    # Keep receptor and chemical sums separate: raw inventory may be Float64,
    # and averaging before saturation is part of the existing law.
    insertion="""
            raw1+=weight*Float32(feed.raw[site,1]);raw2+=weight*Float32(feed.raw[site,2])
            bound1+=weight*feed.bound[site,1];bound2+=weight*feed.bound[site,2]
            enzyme1+=weight*reinterpret(Float32,feed.enzyme_bits[site,1])
            enzyme2+=weight*reinterpret(Float32,feed.enzyme_bits[site,2])
"""
    sense=replace(sense,"            weight=edgeweight/g.volume\n"=>"            weight=edgeweight/g.volume\n"*insertion)
    receptors="""
        for j in 1:3
            code=channels[j]
            value=code<=4 ? dh_material_channel(code,raw1,raw2) :
                code<=8 ? dh_material_channel(code,bound1,bound2) : dh_material_channel(code,enzyme1,enzyme2)
            s.inputs[i,2Int(g.hidden)+15+j]=value/(scales[j]+abs(value))
            !isfinite(value) && dh_fail!(e.status,8)
        end
"""
    # Insert after all native writes so native zero initialization cannot erase it.
    ending="    end\nend\n"
    endswith(sense,ending)||error("native sense tail changed")
    sense=sense[1:end-length(ending)]*receptors*ending
    Core.eval(Interior,Meta.parse(sense))
    source=read(hostpath,String)
    enqueue=feedback_section(source,"function device_cellular_enqueue!","function device_cellular_mechanics!")
    enqueue=replace(enqueue,"function device_cellular_enqueue!"=>"function ff_enqueue!")
    original="    device_cellular_launch!(p,dc_sense!,s,e,q,g;ndrange=C)\n    device_anabolic_sense!(anabolism,w,p,delta)"
    replacement="""
    feed=(bound=anabolism.arrays.bound,enzyme_bits=anabolism.arrays.enzyme_bits,raw=anabolism.feedstock.arrays.fields)
    device_cellular_launch!(p,ff_fused_sense!,s,e,q,g,feed,anabolism.parameters.material_channels,anabolism.parameters.material_scales;ndrange=C)
"""
    occursin(original,enqueue)||error("native sense schedule changed")
    Core.eval(Interior,Meta.parse(replace(enqueue,original=>replacement)))
    step=feedback_section(source,"function device_cellular_step!","# Internal fixed-capacity launch schedule.")
    step=replace(step,"function device_cellular_step!"=>"function ff_step!","device_cellular_enqueue!("=>"ff_enqueue!(")
    Core.eval(Interior,Meta.parse(step))
    (;kernel_sha256=bytes2hex(sha256(read(kernelpath))),host_sha256=bytes2hex(sha256(read(hostpath))))
end
const FEEDBACK_FACTORIZATION_SOURCE=install_feedback_factorization!()

function feedback_factorized_step!(w,p,dt)
    hasproperty(p.habitat.parameters,:material_channels)||error("factorized study requires explicit material receptors")
    Interior.ff_step!(w,p.cells,dt;anabolism=p.habitat)
end


function feedback_physical_error(a,b)
    error=0.
    for key in (:active,:ids,:bonds,:bond_ids)
        a.cells.state[key]==b.cells.state[key]||return Inf
    end
    for key in (:position,:mass,:hidden,:adaptation,:actions,:metabolites,:energy)
        error=max(error,maximum(abs,Float64.(a.cells.state[key]).-Float64.(b.cells.state[key]);init=0.))
    end
    for (x,y) in ((a.cells.environment.fields,b.cells.environment.fields),(a.feedstock.fields,b.feedstock.fields),
                  (a.arrays.bound,b.arrays.bound),(reinterpret(Float32,a.arrays.enzyme_bits),reinterpret(Float32,b.arrays.enzyme_bits)))
        error=max(error,maximum(abs,Float64.(x).-Float64.(y);init=0.))
    end
    error
end

# Invoke explicitly after inclusion; all measurements cover full production
# schedule and status synchronization, without captures in the timed window.
function benchmark_feedback_factorization(output;steps=20,repeats=5,backend=KernelAbstractions.CPU(),
        profile=joinpath(@__DIR__,"..","experiments","native","constructed-garden-workbench-16.toml"),founders_counts=(4,64),parity_tolerance=0.,candidate=feedback_factorized_step!,optimization="sensing_fusion")
    config=TOML.parsefile(profile)
    get!(config,"material_sensing",Dict("channels"=>["raw_total","bound_total","catalyst_total"],"scales"=>[4.,.02,.004]))
    results=Dict{String,Any}()
    for founders in founders_counts
        cfg=deepcopy(config);cfg["cellular"]["founders"]=founders
        cfg["cellular"]["capacity"]=max(cfg["cellular"]["capacity"],16founders)
        w,a,_=native_constructed_garden(backend,cfg)
        v,b,_=native_constructed_garden(backend,cfg)
        dt=Float32(cfg["native"]["dt"])
        ordinary()=Interior.device_cellular_step!(w,a.cells,dt;anabolism=a.habitat)
        fused()=candidate(v,b,dt)
        ordinary();fused()
        warm_a=Interior.device_habitat_snapshot(a);warm_b=Interior.device_habitat_snapshot(b)
        feedback_physical_error(warm_a,warm_b)<=parity_tolerance||error("warm physical parity failed")
        exact=warm_a==warm_b;errors=Float64[]
        times=[Float64[],Float64[]]
        for repeat in 1:repeats
            for arm in (isodd(repeat) ? (1,2) : (2,1))
                step=arm==1 ? ordinary : fused
                GC.gc();start=time_ns()
                for _ in 1:steps;step();end
                push!(times[arm],(time_ns()-start)/1e6/steps)
            end
            sa=Interior.device_habitat_snapshot(a);sb=Interior.device_habitat_snapshot(b)
            err=feedback_physical_error(sa,sb);push!(errors,err);exact &= sa==sb
            err<=parity_tolerance||error("complete trajectory physical parity failed: $err")
        end
        results[string(founders)]=(;founders,grid=cfg["grid"]["shape"],capacity=cfg["cellular"]["capacity"],steps,repeats,
            ordinary_ms=times[1],fused_ms=times[2],speedup=Statistics.median(times[1])/Statistics.median(times[2]),exact_all_snapshots=exact,physical_max_errors=errors,parity_tolerance)
    end
    evidence=(;source=FEEDBACK_FACTORIZATION_SOURCE,julia=string(VERSION),threads=Threads.nthreads(),backend=string(typeof(backend)),optimization,results)
    open(io->JSON3.pretty(io,evidence),output,"w")
    evidence
end

# Independent follow-up candidate: retain all scratch writes and commit launches,
# but avoid six-neighbor stencils in species whose declared diffusion is zero.
# Explicit installation keeps the completed fusion experiment source unchanged.
const FEEDBACK_ZERO_DIFFUSION_TOKEN_SHA256="3634c451d8c6331adda4173f6f6d034edc4aaf06cd7d3fcc2c2e2231e666f8ca"
function install_zero_diffusion_factorization!()
    source=read(joinpath(pkgdir(Interior),"src","device_cellular_kernels.jl"),String)
    kernel=feedback_section(source,"@kernel function dc_diffuse!","@kernel function dc_geometry!")
    if !occursin("if q.diffusion[k]!=0",kernel)
        marker="        pos=dc_coords(site,g.dims);value=e.fields[site,k];flux=zero(value)"
        occursin(marker,kernel)||error("native diffusion setup changed")
        kernel=replace(kernel,marker=>"        value=e.fields[site,k];flux=zero(value)\n        if q.diffusion[k]!=0\n        pos=dc_coords(site,g.dims)")
        kernel=replace(kernel,"        e.scratch[site,k]=value+dt*flux"=>"        end\n        e.scratch[site,k]=value+dt*flux")
    end
    # Whitespace-independent token hash matches the exact kernel body measured
    # before production integration, including the arithmetic/scratch expression.
    bytes2hex(sha256(replace(kernel,r"\s+"=>"")))==FEEDBACK_ZERO_DIFFUSION_TOKEN_SHA256 || error("zero-D benchmark body changed")
    Core.eval(Interior,Meta.parse(replace(kernel,"dc_diffuse!"=>"ff_zero_diffuse!")))
    # Force the old stencil even for D=0; compiler erases the constant true guard.
    reference=replace(kernel,"dc_diffuse!"=>"ff_reference_diffuse!","if q.diffusion[k]!=0"=>"if true")
    Core.eval(Interior,Meta.parse(reference))
    habitat=read(joinpath(pkgdir(Interior),"src","device_habitat.jl"),String)
    host=read(joinpath(pkgdir(Interior),"src","device_cellular.jl"),String)
    for prefix in ("ff_zero","ff_reference")
        beginstep=feedback_section(habitat,"function device_anabolic_begin!","function device_anabolic_sense!")
        beginstep=replace(beginstep,"device_anabolic_begin!"=>prefix*"_begin!","dc_diffuse!"=>prefix*"_diffuse!")
        Core.eval(Interior,Meta.parse(beginstep))
        Core.eval(Interior,Meta.parse(prefix*"_begin!(::Nothing,w,p,dt)=nothing"))
        enqueue=feedback_section(host,"function device_cellular_enqueue!","function device_cellular_mechanics!")
        enqueue=replace(enqueue,"function device_cellular_enqueue!"=>"function "*prefix*"_enqueue!",
            "dc_diffuse!"=>prefix*"_diffuse!","device_anabolic_begin!"=>prefix*"_begin!")
        Core.eval(Interior,Meta.parse(enqueue))
        step=feedback_section(host,"function device_cellular_step!","# Internal fixed-capacity launch schedule.")
        step=replace(step,"function device_cellular_step!"=>"function "*prefix*"_step!","device_cellular_enqueue!("=>prefix*"_enqueue!(")
        Core.eval(Interior,Meta.parse(step))
    end
    nothing
end
feedback_zero_diffusion_step!(w,p,dt)=Interior.ff_zero_step!(w,p.cells,dt;anabolism=p.habitat)
feedback_reference_diffusion_step!(w,p,dt)=Interior.ff_reference_step!(w,p.cells,dt;anabolism=p.habitat)
