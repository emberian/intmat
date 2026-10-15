using Interior, FFTW, LinearAlgebra, JSON3, TOML, SHA, Statistics
BLAS.set_num_threads(1); FFTW.set_num_threads(Threads.nthreads())
const ROOT=dirname(@__DIR__)
# Original periodic mode counts; extent64 reproduces the old physical wavelengths.
@inline function seedvalue(i,j,k,n)
    a,b,c=2pi*(i-1)/n,2pi*(j-1)/n,2pi*(k-1)/n
    .14*(sin(5a)*cos(4b)+sin(4b)*cos(5c)+sin(5c)*cos(4a))+.035sin(7a+3b-5c)
end
function setup(n;extent=(64,64,64),kappa=.85,mobility=1,stabilization=3,seed=702)
    w=World(Grid(n,n,n;extent=Tuple(extent));seed=seed)
    p=ImplicitPhaseSeparation(w;noise=0,kappa=kappa,mobility=mobility,mineral_pinning=0,nutrient_diffusion=0,interface_activity=.65,activity_relaxation=.3,stabilization=stabilization)
    for k in 1:n,j in 1:n,i in 1:n;p.phi[i,j,k]=seedvalue(i,j,k,n);end
    Interior.phase_observe!(w,p)
    w,p
end
function rss_kib()
    parse(Int,strip(read(`ps -o rss= -p $(getpid())`,String)))
end
function benchmark(sizes=(96,192))
    rows=[]
    for n in sizes
        GC.gc(); start=time();w,p=setup(n); setup_s=time()-start
        step!(w,p,.25f0)
        timings=Float64[];alloc=Int[]
        for _ in 1:3
            t=@timed step!(w,p,.25f0);push!(timings,t.time);push!(alloc,t.bytes)
        end
        row=(n=n,extent=[64,64,64],spacing=64/n,setup_seconds=setup_s,warm_seconds=timings,median_step_seconds=median(timings),allocated_bytes_per_step=alloc,rss_kib=rss_kib(),persistent_array_bytes=sum(sizeof,values(w.fields))+sum(sizeof,values(w.scratch))+sum(sizeof,(p.phi,p.interface,p.chemical_potential,p.scratch,p.chemical_forcing,p.eigenvalues,p.real_buffer,p.spectrum,p.nonlinear_spectrum)),finite=all(isfinite,p.phi))
        push!(rows,row);println(JSON3.write(row));flush(stdout)
        w=nothing;p=nothing;GC.gc()
    end
    out=joinpath(ROOT,"outputs","refinement-20261015");mkpath(out)
    open(io->JSON3.pretty(io,(julia_threads=Threads.nthreads(),fftw_threads=FFTW.get_num_threads(),blas_threads=BLAS.get_num_threads(),solver="Interior.ImplicitPhaseSeparation, CPU real-to-complexF64 FFT; ephemeral native plans",benchmarks=rows)),joinpath(out,"phase-core-realfft-benchmark-"*join(sizes,"-")*".json"),"w")
end
function compare_baseline()
    w,p=setup(96)
    for _ in 1:240;step!(w,p,.25f0);end
    path=joinpath(ROOT,"outputs","collection-20261015","phase-gyroid-cathedral-aged","volume","phi.f32")
    old=reshape(open(io->read!(io,Vector{Float32}(undef,96^3)),path),96,96,96)
    result=(max_error=maximum(abs.(old.-p.phi)),rms_error=sqrt(mean((old.-p.phi).^2)),time=w.time,relative_signed_mass_difference=(sum(Float64,p.phi)-sum(Float64,old))/96^3)
    println(JSON3.write(result))
    open(io->JSON3.pretty(io,result),joinpath(ROOT,"outputs","refinement-20261015","phase-core-realfft-baseline-comparison.json"),"w")
end
function metrics(w,p,m0)
    n=size(p.phi,1);squared=0.;changed=0;positive=0
    for k in 1:n,j in 1:n,i in 1:n
        a=seedvalue(i,j,k,n);b=p.phi[i,j,k];squared+=(b-a)^2;changed+=(b>0)!=(a>0);positive+=b>0
    end
    (time=w.time,energy=phase_energy(w,p),range=extrema(p.phi),rms_change=sqrt(squared/n^3),sign_change_fraction=changed/n^3,positive_fraction=positive/n^3,mean_mass_drift=(sum(Float64,p.phi)-m0)/n^3,finite=all(isfinite,p.phi))
end
function bake(recipe)
    cfg=TOML.parsefile(recipe);n=cfg["resolution"];out=joinpath(ROOT,"outputs","refinement-20261015",cfg["id"])
    ispath(out) && error("Fresh refinement output required; existing observations preserved: $out")
    get(cfg,"real_fft",true) || error("ImplicitPhaseSeparation uses the real-FFT representation")
    dt=Float32(cfg["dt"]); isfinite(dt) && dt>0 || error("Positive finite dt required")
    steps=round(Int,cfg["t_end"]/dt)
    steps>0 && isapprox(steps*Float64(dt),cfg["t_end"];rtol=1e-7) || error("t_end must be a positive multiple of dt")
    mkpath(out);cp(recipe,joinpath(out,"recipe.toml");force=true);cp(@__FILE__,joinpath(out,"study-source.jl");force=true)
    w,p=setup(n;extent=get(cfg,"extent",[64,64,64]),kappa=get(cfg,"kappa",.85),mobility=get(cfg,"mobility",1),stabilization=get(cfg,"stabilization",3),seed=get(cfg,"seed",702))
    m0=sum(Float64,p.phi); observations=Any[metrics(w,p,m0)]
    warm_timings=Float64[]
    elapsed=@elapsed for s in 1:steps
        step_start=time_ns()
        step!(w,p,dt)
        if s<=4
            push!(warm_timings,(time_ns()-step_start)/1e9)
            if s==4
                resource=(resolution=n,threads=Threads.nthreads(),fftw_threads=FFTW.get_num_threads(),rss_kib=rss_kib(),four_step_seconds=warm_timings,median_warm_step_seconds=median(warm_timings[2:4]),estimated_remaining_seconds=(steps-4)*median(warm_timings[2:4]))
                open(io->JSON3.pretty(io,resource),joinpath(out,"resource-measurement.json"),"w")
                println("PHASE PRODUCTION RESOURCE ",JSON3.write(resource));flush(stdout)
                resource.rss_kib<14*1024^2 || error("384 benchmark exceeds coordinated14GiB budget")
            end
        end
        if s==steps÷2
            push!(observations,metrics(w,p,m0))
            open(io->write(io,p.phi),joinpath(out,"phi-midpoint.f32"),"w")
        end
        s%10==0 && (println("$(cfg["id"]) $s/$steps t=$(w.time) rss_kib=$(rss_kib())");flush(stdout))
    end
    push!(observations,metrics(w,p,m0));all(o->o.finite,observations) || error("Nonfinite refinement")
    # Six arrays only; no duplicated full World/FFT checkpoint graph.
    density=fill(.01f0,size(p.phi));activity=zeros(Float32,size(p.phi));structure=similar(p.phi);mask=similar(p.phi)
    for I in CartesianIndices(p.phi)
        structure[I]=(1+tanh(2p.phi[I]))/2
        xyz=(Tuple(I).-((n+1)/2,))./(n/2);mask[I]=sum(abs2,xyz)<.88^2 ? 1 : 0
    end
    export_volume((density=density,activity=activity,structure=structure,mask=mask,phi=p.phi,interface=p.interface),joinpath(out,"volume"))
    hashes=Dict(path=>bytes2hex(sha256(read(joinpath(ROOT,path)))) for path in ("studies/refinement_phase.jl","studies/collection_phase.jl","src/phase.jl","src/implicit_phase.jl"))
    baseline=get(cfg,"baseline","outputs/collection-20261015/phase-gyroid-cathedral-aged")
    baseline_match=Tuple(w.grid.extent)==(64f0,64f0,64f0) && p.kappa==.85f0 && p.mobility==1 && p.stabilization==3 && dt==.25f0 && cfg["t_end"]==60
    metadata=(id=cfg["id"],title=cfg["title"],family="phase",shape=[n,n,n],extent=collect(w.grid.extent),spacing=collect(w.grid.spacing),dt=dt,steps=steps,time=w.time,seed=w.seed,parameters=cfg,solver="Interior.ImplicitPhaseSeparation; CPU real-to-complexF64 FFT; constant mobility",phase_parameters=(bulk=p.bulk,kappa=p.kappa,mobility=p.mobility,mineral_pinning=p.mineral_pinning,nutrient_diffusion=p.nutrient_diffusion,interface_activity=p.interface_activity,activity_relaxation=p.activity_relaxation,stabilization=p.stabilization),initialization="Original weak periodic mode counts evaluated afresh on the native lattice; physical wavelengths follow declared extent; not interpolated",baseline=baseline,baseline_dynamics_match=baseline_match,spatial_refinement_per_axis=n/96,observations=observations,elapsed_seconds=elapsed,source_sha256=hashes,threads=(julia=Threads.nthreads(),fftw=FFTW.get_num_threads(),blas=BLAS.get_num_threads()),solver_equivalence=JSON3.read(read(joinpath(ROOT,"outputs","refinement-20261015","phase-core-realfft-baseline-comparison.json"),String)),retention="Three scalar observations; midpoint raw phi; final six scalar volume channels. Initialization is exactly reconstructible. No full World/FFT buffers retained.",optical_mapping="structure=(1+tanh(2phi))/2; zero emission, spherical viewing crop; constant dilute optical carrier")
    open(io->JSON3.pretty(io,metadata),joinpath(out,"run.json"),"w")
    baseline_job=JSON3.read(read(joinpath(ROOT,baseline,"render-job.json"),String),Dict{String,Any})
    merge!(baseline_job,Dict("id"=>cfg["id"],"title"=>cfg["title"],"description"=>"Fresh $(n)³ spinodal evolution to physical time$(w.time), extent$(Tuple(w.grid.extent)), kappa$(p.kappa), mobility$(p.mobility), stabilization$(p.stabilization); not upsampled. Baseline dynamics match: $baseline_match.","size"=>2400,"samples"=>128,"native_shape"=>[n,n,n]))
    open(io->JSON3.pretty(io,baseline_job),joinpath(out,"render-job.json"),"w")
end
if ARGS==["--compare-baseline"]
    compare_baseline()
elseif ARGS==["--benchmark-384"]
    benchmark((384,))
elseif ARGS==["--benchmark-real"]
    benchmark()
elseif ARGS==["--benchmark"]
    benchmark()
else
    length(ARGS)==1 || error("usage: refinement_phase.jl --benchmark | RECIPE.toml")
    bake(ARGS[1])
end
